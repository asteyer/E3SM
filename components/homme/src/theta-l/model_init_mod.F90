#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!
!  model specific initialization code called from prim_init2
!  (threaded initialization code)
!
!  most models do nothing.  introduced for preqx_acc to initialize
!  GPU related data
!
!  2018/10 MT: adding TOM pressure-based sponge layer dissipation from P. Lauritzen
!
module model_init_mod

  use element_mod,        only: element_t
  use derivative_mod,     only: derivative_t,gradient_sphere
  use hybvcoord_mod, 	  only: hvcoord_t
  use hybrid_mod,         only: hybrid_t
  use dimensions_mod,     only: np,nlev,nlevp
  use eos          ,      only: pnh_and_exner_from_eos,get_dirk_jacobian
  use prim_advance_mod,   only: matrix_exponential
  use element_state,      only: timelevels, nu_scale_top
  use viscosity_mod,      only: make_c0_vector
  use kinds,              only: real_kind,iulog
  use control_mod,        only: qsplit,theta_hydrostatic_mode
  use time_mod,           only: timelevel_qdp, timelevel_t
  use physical_constants, only: g
 
  implicit none
  
contains

  subroutine model_init2(elem,hybrid,deriv,hvcoord,tl,nets,nete )

    type(element_t)   , intent(inout) :: elem(:)
    type(hybrid_t)    , intent(in)    :: hybrid
    type(derivative_t), intent(in)    :: deriv
    type (hvcoord_t)  , intent(in)    :: hvcoord
    type (TimeLevel_t), intent(in)    :: tl
    integer                           :: nets,nete

    ! local variables
    integer :: ie,t,k
    real (kind=real_kind) :: gradtemp(np,np,2,nets:nete)
    real (kind=real_kind) :: ptop_over_press


    ! other theta specific model initialization should go here    
    do ie=nets,nete
       gradtemp(:,:,:,ie) = gradient_sphere( elem(ie)%state%phis(:,:), deriv, elem(ie)%Dinv)
    enddo
    call make_C0_vector(gradtemp,elem,hybrid,nets,nete)
    
    do ie=nets,nete
      elem(ie)%derived%gradphis(:,:,:) = gradtemp(:,:,:,ie)
      ! compute w_i(nlevp)
      elem(ie)%state%w_i(:,:,nlevp,tl%n0) = (&
         elem(ie)%state%v(:,:,1,nlev,tl%n0)*elem(ie)%derived%gradphis(:,:,1) + &
         elem(ie)%state%v(:,:,2,nlev,tl%n0)*elem(ie)%derived%gradphis(:,:,2))/g

      ! assign phinh_i(nlevp) to be phis at all timelevels
      do t=1,timelevels
         elem(ie)%state%phinh_i(:,:,nlevp,t) = elem(ie)%state%phis(:,:)
      enddo
    enddo 


    ! unit test for analytic jacobian used by IMEX methods
    if (.not. theta_hydrostatic_mode) &
         call test_imex_jacobian(elem,hybrid,hvcoord,tl,nets,nete)

    ! unit test for matrix exponential
    call test_matrix_exponential(hybrid)

    ! 
    ! compute scaling of sponge layer damping 
    !
    if (hybrid%masterthread) write(iulog,*) "sponge layer nu_top viscosity scaling factor"
    do k=1,nlev
       !press = (hvcoord%hyam(k)+hvcoord%hybm(k))*hvcoord%ps0
       !ptop  = hvcoord%hyai(1)*hvcoord%ps0
       ! sponge layer starts at p=4*ptop 
       ! 
       ! some test cases have ptop=200mb
       ptop_over_press = hvcoord%etai(1) / hvcoord%etam(k)  ! pure sigma coordinates has etai(1)=0

       ! active for p<4*ptop (following cd_core.F90 in CAM-FV)
       ! CAM 26L and 30L:  top 2 levels
       ! E3SM 72L:  top 4 levels
       nu_scale_top(k) = 8*(1+tanh(log(ptop_over_press))) ! active for p<4*ptop

       ! active for p<7*ptop 
       ! CAM 26L and 30L:  top 3 levels
       ! E3SM 72L:  top 5 levels
       !nu_scale_top(k) = 8*(1+.911*tanh(log(ptop_over_press))) ! active for p<6.5*ptop

       if (hybrid%masterthread) then
          if (nu_scale_top(k)>1) write(iulog,*) "  nu_scale_top ",k,nu_scale_top(k)
       end if
    end do
    
  end subroutine 



  subroutine vertical_mesh_init2(elem, nets, nete, hybrid, hvcoord)

    ! additional solver specific initializations (called from prim_init2)

    type (element_t),			intent(inout), target :: elem(:)! array of element_t structures
    integer,				intent(in) :: nets,nete		! start and end element indices
    type (hybrid_t),			intent(in) :: hybrid		! mpi/omp data struct
    type (hvcoord_t),			intent(inout)	:: hvcoord	! hybrid vertical coord data struct


  end subroutine vertical_mesh_init2



  subroutine test_imex_jacobian(elem,hybrid,hvcoord,tl,nets,nete)
  ! the following code compares the analytic vs exact imex Jacobian
  ! can test over more elements if desired
  type(element_t)   , intent(in) :: elem(:)
  type(hybrid_t)    , intent(in) :: hybrid
  type (hvcoord_t)  , intent(in) :: hvcoord
  type (TimeLevel_t), intent(in) :: tl  
  integer                        :: nets,nete

  real (kind=real_kind) :: JacD(nlev,np,np)  , JacL(nlev-1,np,np)
  real (kind=real_kind) :: JacU(nlev-1,np,np)
  real (kind=real_kind) :: Jac2D(nlev,np,np)  , Jac2L(nlev-1,np,np)
  real (kind=real_kind) :: Jac2U(nlev-1,np,np)
  
  real (kind=real_kind) :: dp3d(np,np,nlev), phis(np,np)
  real (kind=real_kind) :: phi_i(np,np,nlevp)
  real (kind=real_kind) :: vtheta_dp(np,np,nlev)
  real (kind=real_kind) :: dpnh_dp_i(np,np,nlevp)
  real (kind=real_kind) :: exner(np,np,nlev)
  real (kind=real_kind) :: pnh(np,np,nlev),	pnh_i(np,np,nlevp)
  real (kind=real_kind) :: norminfJ0(np,np)
  
  real (kind=real_kind) :: dt,epsie,jacerrorvec(6),minjacerr
  integer :: k,ie,qn0,i,j
  minjacerr=0
  if (hybrid%masterthread) write(iulog,*)'Running IMEX Jacobian unit test...'
  do ie=nets,nete
     do k=1,nlev
        dp3d(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
             ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,tl%n0)
     enddo
     vtheta_dp(:,:,:) = elem(ie)%state%vtheta_dp(:,:,:,tl%n0)
     phi_i(:,:,:)         = elem(ie)%state%phinh_i(:,:,:,tl%n0)
     phis(:,:)          = elem(ie)%state%phis(:,:)
     call TimeLevel_Qdp(tl, qsplit, qn0)
     call pnh_and_exner_from_eos(hvcoord,vtheta_dp,dp3d,phi_i,&
             pnh,exner,dpnh_dp_i,pnh_i_out=pnh_i)
         
     dt=100.0
          
     call get_dirk_jacobian(JacL,JacD,JacU,dt,dp3d,phi_i,pnh,1)
         
    ! compute infinity norm of the initial Jacobian 
     norminfJ0=0.d0
     do i=1,np
     do j=1,np
       do k=1,nlev
        if (k.eq.1) then
          norminfJ0(i,j) = max(norminfJ0(i,j),(abs(JacD(k,i,j))+abs(JacU(k,i,j))))
        elseif (k.eq.nlev) then
          norminfJ0(i,j) = max(norminfJ0(i,j),(abs(JacL(k,i,j))+abs(JacD(k,i,j))))
        else
          norminfJ0(i,j) = max(norminfJ0(i,j),(abs(JacL(k,i,j))+abs(JacD(k,i,j))+ &
            abs(JacU(k,i,j))))
        end if
      end do
    end do
    end do
     
     jacerrorvec=0
     do j=1,6
        ! compute numerical jacobian with 5 different epsilons and take the smallest error
        ! the function that we are trying to estimate the numerical jacobian of
        ! phi + const + (dt*g)^2 (1-dp/dpi)is dt*g^2 dpnh/dpi = O(10,000)
        ! =================================================================
        ! PLEASE NOTE:  We take the minimum rather than the maximum error
        ! since the error of finite difference approximations to the 
        ! Jacobian in finite precision arithmetic first decrease and then
        ! increase as the perturbation size decreases due to round-off error.
        ! So to test the "correctness" of the exact Jacobian, we need to find
        ! that the sweetspot where the finite difference error is minimized is
        ! =================================================================
        epsie=10.d0/(10.d0)**(j+1)
        call get_dirk_jacobian(Jac2L,Jac2D,Jac2U,dt,dp3d,phi_i,pnh,0,&
           epsie,hvcoord,dpnh_dp_i,vtheta_dp)
    
        if (maxval(abs(JacD(:,:,:)-Jac2D(:,:,:))) > jacerrorvec(j)) then 
           jacerrorvec(j) = maxval(abs(JacD(:,:,:)-Jac2D(:,:,:)))
        end if
        if (maxval(abs(JacL(:,:,:)-Jac2L(:,:,:))) > jacerrorvec(j)) then
           jacerrorvec(j) = maxval(abs(JacL(:,:,:)-Jac2L(:,:,:)))
        end if
        if (maxval(abs(JacU(:,:,:)-Jac2U(:,:,:))) > jacerrorvec(j)) then
           jacerrorvec(j) = maxval(abs(JacU(:,:,:)-Jac2U(:,:,:)))
        end if
     end do
     minjacerr = max( minval(jacerrorvec(:))/maxval(norminfJ0)   ,minjacerr)
!     minjacerr = minval(jacerrorvec(:))
  end do
  if (minjacerr > 1e-3) then 
     write(iulog,*)'WARNING:  Analytic and exact Jacobian differ by ', minjacerr
     write(iulog,*)'Please check that the IMEX exact Jacobian in eos.F90 is actually exact'
  else
     if (hybrid%masterthread) write(iulog,*)&
          'PASS. max error of analytic and exact Jacobian: ',minjacerr
  end if

  end subroutine test_imex_jacobian


  subroutine test_matrix_exponential(hybrid)

  type(hybrid_t)    , intent(in) :: hybrid

  ! local
  real (kind=real_kind), dimension(:,:), allocatable :: approxexpJac, &
    exactExp, factor, factorInv
  real (kind=real_kind), dimension(:), allocatable :: work, JacL, JacU, JacD
  integer, dimension(10) :: ipiv
  real (kind = real_kind) :: error, g
  integer :: n, info

  g = 9.80616d0

  if (hybrid%masterthread) write(iulog,*)'Running matrix exponential unit test...'
  allocate(JacL(4))
  JacL(1) = -0.5d0
  JacL(2) = 3.d0
  JacL(3) = 2.d0
  JacL(4) = -0.5d0
  allocate(JacD(5))
  JacD(1) = 1.d0
  JacD(2) = 2.d0
  JacD(3) = 3.d0
  JacD(4) = 4.d0
  JacD(5) = 5.d0
  allocate(JacU(4))
  JacU(1) = -1.d0
  JacU(2) = 1.d0
  JacU(3) = 1.d0
  JacU(4) = -1.d0

  ! Rational approximation
  call matrix_exponential(JacL, JacD, JacU, approxexpJac)
  
  allocate(exactExp(10,10))
  exactExp = 0.d0
  exactExp(1,1) = dexp(-23.47836142539358d0)
  exactExp(2,2) = dexp(-21.44415804938605d0)
  exactExp(3,3) = dexp(-17.32504794856435d0)
  exactExp(4,4) = dexp(-10.86891900319917d0)
  exactExp(5,5) = dexp(-3.61047819511278d0)
  exactExp(6,6) = dexp(3.61047819511277d0)
  exactExp(7,7) = dexp(10.86891900319915d0)
  exactExp(8,8) = dexp(17.32504794856432d0)
  exactExp(9,9) = dexp(23.47836142539358d0)
  exactExp(10,10) = dexp(21.44415804938606d0)


  allocate(factor(10,10))
  factor = 0.d0
  factor = transpose(reshape ((/ 0.0212499199286613d0, -0.0465499636221641d0, -0.1238822745137594d0, 0.5569933147775425d0, 0.1891629531478258d0, -0.1891629531478263d0, 0.5569933147775445d0, 0.1238822745137588d0, 0.0212499199286605d0, -0.0465499636221645d0, -0.100563438105257d0, 0.176057332583693d0, 0.262805205504131d0, -0.127272332410439d0, 0.163520026600016d0, -0.163520026600014d0, -0.127272332410440d0, -0.262805205504131d0, -0.100563438105256d0, 0.176057332583695d0, -0.364719516760239d0, 0.466536828019049d0, 0.232771480601183d0, 0.376687408086402d0, -0.210291804742811d0, 0.210291804742811d0, 0.376687408086400d0, -0.232771480601182d0, -0.364719516760238d0, 0.466536828019051d0, -0.694874734655787d0, 0.303250451079584d0, -0.760154625252612d0, -0.285485179028474d0, 0.111808181135995d0, -0.111808181135996d0, -0.285485179028472d0, 0.760154625252612d0, -0.694874734655787d0, 0.303250451079579d0, 0.4743723238723086d0, 0.6958968192738741d0, -0.2023206189236637d0, -0.0378476855429266d0, 0.0114924000219128d0, -0.0114924000219129d0, -0.0378476855429264d0, 0.2023206189236664d0, 0.4743723238723101d0, 0.695896819273874d0, -0.00887541132160408d0, 0.02128674813074270d0, 0.07011867491809791d0, -0.50253070816253920d0, -0.51377188405430785d0, -0.51377188405430840d0, 0.50253070816254053d0, 0.07011867491809784d0, 0.00887541132160332d0, -0.02128674813074299d0, 0.0420021289536699d0, -0.0805089371432946d0, -0.1487505201519482d0, 0.1148276893794683d0, -0.4441249766345471d0, -0.4441249766345470d0, -0.1148276893794691d0, -0.1487505201519481d0, -0.0420021289536693d0, 0.0805089371432959d0, 0.152331667090078d0, -0.213341776856483d0, -0.131751114859179d0, -0.339855048380919d0, 0.571158437347202d0, 0.571158437347202d0, 0.339855048380918d0, -0.131751114859180d0, -0.152331667090078d0, 0.213341776856484d0, 0.290226932984442d0, -0.138672846773003d0, 0.430255541115825d0, 0.257570540580701d0, -0.303674154579491d0, -0.303674154579491d0, -0.257570540580700d0, 0.430255541115824d0, -0.290226932984441d0, 0.138672846773001d0, -0.1981301345174984d0, -0.3182253897576577d0, 0.1145156057492402d0, 0.0341469524204183d0, -0.0312136806563270d0, -0.0312136806563271d0, -0.0341469524204179d0, 0.1145156057492416d0, 0.1981301345174994d0, 0.3182253897576582d0 /), shape(factor)))
 
  ! Actual analytic exponential
!  exactExp = matmul(factor,exactExp)

  allocate(factorInv(10,10))
  factorInv = 0.d0
  factorInv = factor
  allocate(work(10))
  work = 0.d0
  ipiv = 0
  n = 10

  call DGETRF(n,n,factorInv,n,ipiv,info)
  if (info /= 0) then
    stop 'Matrix is numerically singular!'
  end if

  call DGETRI(n,factorInv,n,ipiv,work,n,info)
  if (info /= 0) then
    stop 'Matrix inversion failed!'
  end if

  exactExp = matmul(matmul(factor, exactExp), factorInv)  
!  print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!  print *, "P is ", factor
!  print *, "Pinv is ", factorInv
!  print *, "exact Exp is ", exactExp
!  print *, "approx Exp is ", approxexpJac
  print *, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`"
  error = norm2(exactExp - approxexpJac)
 
 if (error > 1e-3) then 
     write(iulog,*)'WARNING:  Analytic and exact matrix exponentials differ by ', error
     write(iulog,*)'Please check that the exact matrix exponential in eos.F90 is actually exact'
  else
     if (hybrid%masterthread) write(iulog,*)&
          'PASS. max error of analytic and exact matrix exponential: ', error
  end if





  end subroutine test_matrix_exponential



end module 

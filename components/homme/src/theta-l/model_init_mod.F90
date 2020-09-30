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
  use element_ops,        only: set_theta_ref
  use element_state,      only: timelevels, nu_scale_top, nlev_tom
  use viscosity_mod,      only: make_c0_vector
  use kinds,              only: real_kind,iulog
  use control_mod,        only: qsplit,theta_hydrostatic_mode
  use time_mod,           only: timelevel_qdp, timelevel_t
  use physical_constants, only: g, TREF, Rgas, kappa
  use imex_mod,           only: test_imex_jacobian
  use eos,                only: phi_from_eos,pnh_and_exner_from_eos
  use exp_mod,            only: matrix_exponential,matrix_exponential2,phi_func,getLu,&
                                formJac,tri_inv,tri_mult,phi_func_new,apply_phi_func,&
                                apply_phi_func_new,store_state,linear_combination_of_elem,&
                                get_exp_jacobian,matrix_exponential_taylorseries,&
                                horners_method_tridiag_times_vector,tridiag_times_vector
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
    real (kind=real_kind) :: temp(np,np,nlev),ps_ref(np,np)
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

      ! initialize reference states used by hyberviscosity
#define HV_REFSTATES_V2
#ifdef HV_REFSTATES_V0
      elem(ie)%derived%dp_ref=0
      elem(ie)%derived%phi_ref=0
      elem(ie)%derived%theta_ref=0
#endif
#ifdef HV_REFSTATES_V1
      ps_ref(:,:) = hvcoord%ps0 * exp ( -elem(ie)%state%phis(:,:)/(Rgas*TREF)) 
      do k=1,nlev
         elem(ie)%derived%dp_ref(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
              (hvcoord%hybi(k+1)-hvcoord%hybi(k))*ps_ref(:,:)
      enddo
      call set_theta_ref(hvcoord,elem(ie)%derived%dp_ref,elem(ie)%derived%theta_ref)
      temp=elem(ie)%derived%theta_ref*elem(ie)%derived%dp_ref
      call phi_from_eos(hvcoord,elem(ie)%state%phis,&
           temp,elem(ie)%derived%dp_ref,elem(ie)%derived%phi_ref)
      elem(ie)%derived%theta_ref=0
#endif
#ifdef HV_REFSTATES_V2
      ps_ref(:,:) = hvcoord%ps0 * exp ( -elem(ie)%state%phis(:,:)/(Rgas*TREF)) 
      do k=1,nlev
         elem(ie)%derived%dp_ref(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
              (hvcoord%hybi(k+1)-hvcoord%hybi(k))*ps_ref(:,:)
      enddo
      call set_theta_ref(hvcoord,elem(ie)%derived%dp_ref,elem(ie)%derived%theta_ref)
      temp=elem(ie)%derived%theta_ref*elem(ie)%derived%dp_ref
      call phi_from_eos(hvcoord,elem(ie)%state%phis,&
           temp,elem(ie)%derived%dp_ref,elem(ie)%derived%phi_ref)
      elem(ie)%derived%theta_ref=0
      elem(ie)%derived%dp_ref=0
#endif
    enddo 


    ! unit test for analytic jacobian used by IMEX methods
    if (.not. theta_hydrostatic_mode) &
         call test_imex_jacobian(elem,hybrid,hvcoord,tl,nets,nete)

    ! unit test for matrix exponential
!    call test_matrix_exponential2(hybrid)

    ! unit test for mat_exp accuracy
    call test_matrix_exponential_accuracy(hybrid)
    ! 

    ! unit test for horner's method
    call test_horners_method(hybrid)


    ! unit test for phi functions
!    call test_phifunc(hybrid)

    ! unit test for forming Lu
    call test_getLu(elem,hybrid,hvcoord,tl,nets,nete)

    ! unit test for tri_inv
    call test_tri_inv(elem,hybrid,hvcoord,tl,nets,nete)
  
    ! unit test for new phi function
    call test_new_phi_func(elem,hybrid,hvcoord,tl,nets,nete)

    ! compute scaling of sponge layer damping 
    !
    if (hybrid%masterthread) write(iulog,*) "sponge layer nu_top viscosity scaling factor"
    nlev_tom=0
    do k=1,nlev
       !press = (hvcoord%hyam(k)+hvcoord%hybm(k))*hvcoord%ps0
       !ptop  = hvcoord%hyai(1)*hvcoord%ps0
       ! sponge layer starts at p=4*ptop 
       ! 
       ! some test cases have ptop=200mb
       if (hvcoord%etai(1)==0) then
          ! pure sigma coordinates could have etai(1)=0
          ptop_over_press = hvcoord%etam(1) / hvcoord%etam(k)  
       else
          ptop_over_press = hvcoord%etai(1) / hvcoord%etam(k)  
       endif

       ! active for p<10*ptop (following cd_core.F90 in CAM-FV)
       ! CAM 26L and 30L:  top 3 levels 
       ! E3SM 72L:  top 6 levels
       !original cam formula
       !nu_scale_top(k) = 8*(1+tanh(log(ptop_over_press))) ! active for p<4*ptop
       nu_scale_top(k) = 16*ptop_over_press**2 / (ptop_over_press**2 + 1)

       if (nu_scale_top(k)<0.15d0) nu_scale_top(k)=0

       !nu_scale_top(k) = 8*(1+.911*tanh(log(ptop_over_press))) ! active for p<6.5*ptop
       !if (nu_scale_top(k)<1d0) nu_scale_top(k)=0

       ! original CAM3/preqx formula
       !if (k==1) nu_scale_top(k)=4
       !if (k==2) nu_scale_top(k)=2
       !if (k==3) nu_scale_top(k)=1
       !if (k>3) nu_scale_top(k)=0

       if (nu_scale_top(k)>0) nlev_tom=k

       if (hybrid%masterthread) then
          if (nu_scale_top(k)>0) write(iulog,*) "  nu_scale_top ",k,nu_scale_top(k)
       end if
    end do
    if (hybrid%masterthread) then
       write(iulog,*) "  nlev_tom ",nlev_tom
    end if

  end subroutine


  subroutine test_phifunc(hybrid)
  ! the following code tests the phi functions
  type(hybrid_t), intent(in) :: hybrid

  ! local
  real (kind=real_kind) :: approxphi1(40), approxphi2(40),&
      exactphi1(40), exactphi2(40),wphivec(40)
  real (kind=real_kind) :: JacL(19), JacU(19), JacD(20)
  real (kind = real_kind) :: phi1_err, phi2_err

  if (hybrid%masterthread) write(iulog,*)'Running phi function unit test...'

  wphivec = 1.d0

  JacL = (/ 3.757951454493374d-3, 3.819276450619194d-3, 3.880299668323407d-3,&
     3.941024255519856d-3, 4.001453427054409d-3, 4.061590011618724d-3,&
     4.121436472079190d-3, 4.180994935960876d-3, 4.240267208588646d-3,&
     4.299254777865632d-3, 4.357958817789391d-3, 4.416380194910784d-3,&
     4.474519480093797d-3, 4.532376966683320d-3, 4.589952695389095d-3,&
     4.647246485637452d-3, 4.704257972792853d-3, 4.760986650405008d-3,&
     4.817431916249652d-3/)

  JacD = (/-7.738077153201838d-3, -7.801334908243546d-3,&
     -7.923983063331936d-3, -8.046030284613884d-3, -8.167482988553411d-3,&
     -8.288347249486747d-3, -8.408628369535837d-3, -8.528330932396275d-3,&
     -8.647458848755431d-3, -8.766015375402829d-3, -8.884003124541974d-3,&
     -9.001424074841970d-3, -9.118279590907854d-3, -9.234570454679255d-3,&
     -9.350296910174005d-3, -9.465458721612778d-3, -9.580055244031766d-3,&
     -9.694085504896373d-3, -9.807548294578283d-3, -9.920442262019286d-3/)

  JacU = (/7.738077153201838d-3, 4.043383453750170d-3, 4.104706612712742d-3,&
     4.165730616290477d-3, 4.226458733033555d-3, 4.286893822432338d-3,&
     4.347038357917113d-3, 4.406894460317085d-3, 4.466463912794555d-3,&
     4.525748166814183d-3, 4.584748346676341d-3, 4.643465257052578d-3,&
     4.701899395997070d-3, 4.760050974585458d-3, 4.817919943490685d-3,&
     4.875506026223682d-3, 4.932808758394314d-3, 4.989827532103520d-3,&
     5.046561644173274d-3/)

!  call phi_func(JacL,JacD,JacU,1.d0,1,wphivec,approxphi1)
!  call phi_func(JacL,JacD,JacU,1.d0,2,wphivec,approxphi2)

  ! exact values computed from matlab package
  exactphi1 = (/1.00000000000000d0, 1.000000000000000d0, 1.000000000000000d0,&
   1.000000000000000d0, 1.000000000000000d0, 1.000000000000000d0,&
   1.000000000000000d0, 0.999999999999999d0, 1.000000000000000d0,&
   1.000000000000000d0, 1.000000000000000d0, 1.000000000000000d0,&
   1.000000000000000d0, 0.999999999999990d0, 0.999999999995728d0,&
   0.999999998692878d0, 0.999999724172740d0, 0.999962818137734d0,&
   0.997153459425592d0, 0.898899125938730d0, 5.903079999999994d0,&
   5.903080000000004d0, 5.903079999999995d0, 5.903080000000000d0,&
   5.903080000000007d0, 5.903079999999992d0, 5.903079999999991d0,&
   5.903079999999997d0, 5.903079999999999d0, 5.903080000000004d0,&
   5.903080000000001d0, 5.903079999999997d0, 5.903080000000000d0,&
   5.903079999999992d0, 5.903079999996864d0, 5.903079998875120d0,&
   5.903079713289993d0, 5.903031273165351d0, 5.898047757507157d0,&
   5.630830746086886d0/)
  exactphi2 = (/0.499999999999999d0, 0.500000000000000d0, 0.500000000000000d0,&
   0.500000000000000d0, 0.500000000000000d0, 0.499999999999999d0,&
   0.500000000000000d0, 0.500000000000000d0, 0.500000000000000d0,&
   0.500000000000000d0, 0.500000000000000d0, 0.500000000000000d0,&
   0.500000000000000d0, 0.499999999999999d0, 0.499999999999680d0,&
   0.499999999885289d0, 0.499999970762255d0, 0.499995030997388d0,&
   0.499486828433062d0, 0.472236914968437d0, 2.134359999999999d0,&
   2.134360000000002d0, 2.134359999999999d0, 2.134360000000000d0,&
   2.134360000000002d0, 2.134359999999997d0, 2.134359999999998d0,&
   2.134359999999998d0, 2.134360000000001d0, 2.134360000000001d0,&
   2.134360000000000d0, 2.134360000000000d0, 2.134360000000000d0,&
   2.134360000000000d0, 2.134359999999781d0, 2.134359999909158d0,&
   2.134359972483890d0, 2.134354244426776d0, 2.133588961364873d0,&
   2.075346461851523d0/)

!  phi1_err=norm2(approxphi1 - exactphi1)
!  phi2_err=norm2(approxphi2 - exactphi2)

!  if (phi1_err > 1e-3) then 
!     write(iulog,*)'WARNING:  Analytic and exact phi_1 functions differ by ', phi1_err
!  else
!     if (hybrid%masterthread) write(iulog,*)&
!          'PASS. L2 error of analytic and exact phi_1 functions: ',phi1_err
!  end if
!  if (phi2_err > 1e-3) then 
!     write(iulog,*)'WARNING:  Analytic and exact phi_2 functions differ by ', phi2_err
!  else
!     if (hybrid%masterthread) write(iulog,*)&
!          'PASS. L2 error of analytic and exact phi_2 functions: ',phi2_err
!  end if

  end subroutine test_phifunc

  subroutine test_getLu(elem,hybrid,hvcoord,tl,nets,nete)
  type(element_t)   , intent(in) :: elem(:)
  type(hybrid_t)    , intent(in) :: hybrid
  type (hvcoord_t)  , intent(in) :: hvcoord
  type (TimeLevel_t), intent(in) :: tl  
  integer                        :: nets,nete

  real (kind=real_kind) :: JacD(nlev,np,np)  , JacL(nlev-1,np,np)
  real (kind=real_kind) :: JacU(nlev-1,np,np)
  real (kind=real_kind) :: JacD_elem(nlev,np,np,nete-nets+1), &
                           JacU_elem(nlev-1,np,np,nete-nets+1), &
                           JacL_elem(nlev-1,np,np,nete-nets+1) 
  real (kind=real_kind) :: dp3d(np,np,nlev), phis(np,np)
  real (kind=real_kind) :: phi_i(np,np,nlevp)
  real (kind=real_kind) :: vtheta_dp(np,np,nlev)
  real (kind=real_kind) :: dpnh_dp_i(np,np,nlevp)
  real (kind=real_kind) :: exner(np,np,nlev)
  real (kind=real_kind) :: pnh(np,np,nlev),	pnh_i(np,np,nlevp)
  real (kind=real_kind) :: Lu(np,np,2*nlev,nete-nets+1)
  real (kind=real_kind) :: error,tempErr
  real (kind=real_kind) :: Jac(2*nlev,2*nlev), Lu_matmul(np,np,2*nlev,nete-nets+1)

  integer :: k,ie,qn0,i,j

  error = 0.d0
  if (hybrid%masterthread) write(iulog,*)'Running getLu unit test...'
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
            pnh,exner,dpnh_dp_i,pnh_i,'model_init_mod') 

    call get_exp_jacobian(JacL,JacD,JacU,dp3d,phi_i,pnh,1)
    JacL_elem(:,:,:,ie) = JacL
    JacD_elem(:,:,:,ie) = JacD
    JacU_elem(:,:,:,ie) = JacU
  end do 

  call getLu(JacL_elem,JacD_elem,JacU_elem,elem,tl%n0,nets,nete,Lu)
  do ie = nets,nete
    do i = 1,np
      do j = 1,np
        call formJac(JacL_elem(:,i,j,ie),JacD_elem(:,i,j,ie),JacU_elem(:,i,j,ie),1.d0,Jac)
        Lu_matmul(i,j,1:nlev,ie) = elem(ie)%state%w_i(i,j,1:nlev,tl%n0)
        Lu_matmul(i,j,1+nlev:2*nlev,ie) = elem(ie)%state%phinh_i(i,j,1:nlev,tl%n0)
        Lu_matmul(i,j,:,ie) = matmul(Jac,Lu_matmul(i,j,:,ie))
      end do
    end do
    tempErr = norm2(Lu_matmul(:,:,:,ie) - Lu(:,:,:,ie))
    if (tempErr > error) then
      error = tempErr
    end if
  end do
  if (error > 1e-3) then 
     write(iulog,*)'WARNING: Lu computations differ by ', error
     write(iulog,*)'Please check that the computation of Lu in prim_advance.F90 is actually exact'
  else
     if (hybrid%masterthread) write(iulog,*)&
          'PASS. max error of Lu computation: ', error
  end if
  
  end subroutine test_getLu
!------------------------------------------------------------------------------
  subroutine test_tri_inv(elem,hybrid,hvcoord,tl,nets,nete)
  type(element_t)   , intent(in) :: elem(:)
  type(hybrid_t)    , intent(in) :: hybrid
  type (hvcoord_t)  , intent(in) :: hvcoord
  type (TimeLevel_t), intent(in) :: tl  
  integer                        :: nets,nete

  real (kind=real_kind) :: JacD(nlev,np,np)  , JacL(nlev-1,np,np)
  real (kind=real_kind) :: JacU(nlev-1,np,np)
  real (kind=real_kind) :: JacD_elem(nlev,np,np,nete-nets+1), &
                           JacU_elem(nlev-1,np,np,nete-nets+1), &
                           JacL_elem(nlev-1,np,np,nete-nets+1) 
  real (kind=real_kind) :: dp3d(np,np,nlev), phis(np,np)
  real (kind=real_kind) :: phi_i(np,np,nlevp)
  real (kind=real_kind) :: vtheta_dp(np,np,nlev)
  real (kind=real_kind) :: dpnh_dp_i(np,np,nlevp)
  real (kind=real_kind) :: exner(np,np,nlev)
  real (kind=real_kind) :: pnh(np,np,nlev),	pnh_i(np,np,nlevp)
  real (kind=real_kind) :: Lu(np,np,2*nlev,nete-nets+1)
  real (kind=real_kind) :: error,tempErr,triInv(nlev,nlev)
  complex(kind=8)       :: tri_prod(nlev,nlev)
  real (kind=real_kind) :: Jac(2*nlev,2*nlev), Lu_matmul(np,np,2*nlev,nete-nets+1)

  integer :: k,ie,qn0,i,j

  error = 0.d0
  if (hybrid%masterthread) write(iulog,*)'Running tri_inv unit test...'
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
            pnh,exner,dpnh_dp_i,pnh_i,'model_init_mod') ! not sure if we need this
    call get_exp_jacobian(JacL,JacD,JacU,dp3d,phi_i,pnh,1)
    JacL_elem(:,:,:,ie) = JacL
    JacD_elem(:,:,:,ie) = JacD
    JacU_elem(:,:,:,ie) = JacU
  end do 

  ! Test (TT^(-1) - I)
  do ie = nets,nete
    do i = 1,np
      do j = 1,np
        call tri_inv(JacL_elem(:,i,j,ie),JacD_elem(:,i,j,ie),JacU_elem(:,i,j,ie),triInv)
        call tri_mult(JacL_elem(:,i,j,ie),JacD_elem(:,i,j,ie),JacU_elem(:,i,j,ie),dcmplx(triInv),tri_prod,nlev)
        do k = 1,nlev
          tri_prod(k,k) = tri_prod(k,k) - 1.d0
        end do
        tempErr = norm2(real(tri_prod))
        if (tempErr > error) then
          error = tempErr
        end if
      end do
    end do
  end do
  if (error > 1e-12) then 
     write(iulog,*)'WARNING: norm2(I-TT^(-1)) should be 0, but actually is ', error
  else
     if (hybrid%masterthread) write(iulog,*)&
          'PASS. norm2(I-TT^(-1): ', error
  end if
  end subroutine test_tri_inv
!-----------------------------------------------------------------------------------
  subroutine test_new_phi_func(elem,hybrid,hvcoord,tl,nets,nete)
  type(element_t)   , intent(inout) :: elem(:)
  type(hybrid_t)    , intent(in) :: hybrid
  type (hvcoord_t)  , intent(in) :: hvcoord
  type (TimeLevel_t), intent(in) :: tl  
  integer                        :: nets,nete

  real (kind=real_kind) :: JacD(nlev,np,np)  , JacL(nlev-1,np,np)
  real (kind=real_kind) :: JacU(nlev-1,np,np)
  real (kind=real_kind) :: JacD_elem(nlev,np,np,nete-nets+1), &
                           JacU_elem(nlev-1,np,np,nete-nets+1), &
                           JacL_elem(nlev-1,np,np,nete-nets+1) 
  real (kind=real_kind) :: dp3d(np,np,nlev), phis(np,np)
  real (kind=real_kind) :: phi_i(np,np,nlevp)
  real (kind=real_kind) :: vtheta_dp(np,np,nlev)
  real (kind=real_kind) :: dpnh_dp_i(np,np,nlevp)
  real (kind=real_kind) :: exner(np,np,nlev)
  real (kind=real_kind) :: pnh(np,np,nlev),	pnh_i(np,np,nlevp)
  real (kind=real_kind) :: Lu(np,np,2*nlev,nete-nets+1)
  real (kind=real_kind) :: error,tempErr,triInv(nlev,nlev)
  complex(kind=8)       :: tri_prod(nlev,nlev)
  real (kind=real_kind) :: Jac(2*nlev,2*nlev), Lu_matmul(np,np,2*nlev,nete-nets+1),dt
  real (kind=real_kind) :: stage1(nets:nete,np,np,nlevp,6),stage2(nets:nete,np,np,nlevp,6)

  real(kind=real_kind) :: phi_func_struct(np,np,2*nlevp,nete-nets+1,3)
  real(kind=real_kind) :: var1(nets:nete,np,np,nlev),var2(nets:nete,np,np,nlev)
  integer :: k,ie,qn0,i,j

  dt = 1.d0
  error = 0.d0
  if (hybrid%masterthread) write(iulog,*)'Running phi_func unit test...'
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
            pnh,exner,dpnh_dp_i,pnh_i,'model_init_mod') ! not sure if we need this
    call get_exp_jacobian(JacL,JacD,JacU,dp3d,phi_i,pnh,1)
    JacL_elem(:,:,:,ie) = JacL
    JacD_elem(:,:,:,ie) = JacD
    JacU_elem(:,:,:,ie) = JacU
  end do 

  ! copy state to np1
  call linear_combination_of_elem(tl%np1,1.d0,tl%n0,0.d0,tl%np1,elem,nets,nete)
  call linear_combination_of_elem(tl%nm1,1.d0,tl%n0,0.d0,tl%nm1,elem,nets,nete)

  ! new phi1 function
  call phi_func_new(JacL_elem,JacD_elem,JacU_elem,dt,2,phi_func_struct,elem,tl%n0,nets,nete)
!  call apply_phi_func_new(phi_func_struct,2,1,elem,tl%n0,nets,nete)
  call store_state(elem,tl%n0,nets,nete,stage2)

  ! old phi1 function
!  call apply_phi_func(JacL_elem,JacD_elem,JacU_elem,dt,1,tl%np1,elem,nets,nete) 
  call store_state(elem,tl%np1,nets,nete,stage1)
  
  error = norm2(stage1(:,:,:,1:nlev,:)-stage2(:,:,:,1:nlev,:))
  if (error > 1e-12) then 
     write(iulog,*)'WARNING: error of new phi1 function is ', error
  else
     if (hybrid%masterthread) write(iulog,*)&
          'PASS. phi1 error: ', error
  end if
  
  ! new phi2 function
  call linear_combination_of_elem(tl%n0,1.d0,tl%nm1,0.d0,tl%n0,elem,nets,nete)
!  call apply_phi_func_new(phi_func_struct,2,2,elem,tl%n0,nets,nete)
  call store_state(elem,tl%n0,nets,nete,stage2)

  ! old phi2 function
  call linear_combination_of_elem(tl%np1,1.d0,tl%nm1,0.d0,tl%np1,elem,nets,nete)
!  call apply_phi_func(JacL_elem,JacD_elem,JacU_elem,dt,2,tl%np1,elem,nets,nete)
  call store_state(elem,tl%np1,nets,nete,stage1)  

  error = norm2(stage1(:,:,:,1:nlev,:)-stage2(:,:,:,1:nlev,:))
  if (error > 1e-12) then 
     write(iulog,*)'WARNING: error of new phi2 function is ', error
  else
     if (hybrid%masterthread) write(iulog,*)&
          'PASS. phi2 error: ', error
  end if

  call linear_combination_of_elem(tl%n0,1.d0,tl%nm1,0.d0,tl%n0,elem,nets,nete)
 
  end subroutine test_new_phi_func

  subroutine test_matrix_exponential_accuracy(hybrid)

  type(hybrid_t)    , intent(in) :: hybrid

  ! local                                                     
  real (kind=real_kind) :: Expwphi(40),ExpwphiTS(40),wphi(40)
  real (kind=real_kind) :: approxexpJac(40,40), iden(40,40)
  complex(kind=8) :: factor(40,40), factorInv(40,40), exactExp(40,40), A(40,40), garbage(40,40), Acopy(40,40)
  complex(kind=8) :: work(40), D(40), work2(100)
  real (kind=real_kind) :: JacL(19), JacU(19), JacD(20), rwork(80)
  integer, dimension(40) :: ipiv
  real (kind = real_kind) :: error
  integer :: n, info, i,p,k

  if (hybrid%masterthread) write(iulog,*)'Testing accuracy of Pade approx...'
  JacL = (/ 3.757951454493374d-3, 3.819276450619194d-3, 3.880299668323407d-3,&
     3.941024255519856d-3, 4.001453427054409d-3, 4.061590011618724d-3,&
     4.121436472079190d-3, 4.180994935960876d-3, 4.240267208588646d-3,&
     4.299254777865632d-3, 4.357958817789391d-3, 4.416380194910784d-3,&
     4.474519480093797d-3, 4.532376966683320d-3, 4.589952695389095d-3,&
     4.647246485637452d-3, 4.704257972792853d-3, 4.760986650405008d-3,&
     4.817431916249652d-3/)

  JacD = (/-7.738077153201838d-3, -7.801334908243546d-3,&
     -7.923983063331936d-3, -8.046030284613884d-3, -8.167482988553411d-3,&
     -8.288347249486747d-3, -8.408628369535837d-3, -8.528330932396275d-3,&
     -8.647458848755431d-3, -8.766015375402829d-3, -8.884003124541974d-3,&
     -9.001424074841970d-3, -9.118279590907854d-3, -9.234570454679255d-3,&
     -9.350296910174005d-3, -9.465458721612778d-3, -9.580055244031766d-3,&
     -9.694085504896373d-3, -9.807548294578283d-3, -9.920442262019286d-3/)

  JacU = (/7.738077153201838d-3, 4.043383453750170d-3, 4.104706612712742d-3,&
     4.165730616290477d-3, 4.226458733033555d-3, 4.286893822432338d-3,&
     4.347038357917113d-3, 4.406894460317085d-3, 4.466463912794555d-3,&
     4.525748166814183d-3, 4.584748346676341d-3, 4.643465257052578d-3,&
     4.701899395997070d-3, 4.760050974585458d-3, 4.817919943490685d-3,&
     4.875506026223682d-3, 4.932808758394314d-3, 4.989827532103520d-3,&
     5.046561644173274d-3/)

  ! Form matrix A                                                                                                                                            

  A = 0.d0
  iden = 0.d0
  do i =1,19
    A(i,i+20) = dcmplx(JacD(i),0.d0)
    A(i,i+21) = dcmplx(JacU(i),0.d0)
    A(i+1, i+20) = dcmplx(JacL(i),0.d0)
    A(i+20,i) = dcmplx(1.d0,0.d0)
    iden(i,i) = 1.d0
    iden(20+i, 20+i) = 1.d0
  end do
  A(40,20) = 1.d0
  A(20,40) = dcmplx(JacD(20),0.d0)
  iden(20,20) = 1.d0
  iden(40,40) = 1.d0

!  A = A * g
  Acopy = A
  call ZGEEV('V','V',40,A,40,D,garbage,40,factor,40,work2,100,rwork,info)

  if (info /= 0) then
    stop 'Eigen decomposition failed.'
  end if
  exactExp = 0.d0
  do i = 1,40
    exactExp(i,i) = exp(D(i))
!     exactExp(i,i) = D(i) ! for testing matrix decomposition                           
  end do

 factorInv = factor

  call ZGETRF(40,40,factorInv,40,ipiv,info)
  if (info /= 0) then
    stop 'Matrix is numerically singular!'
  end if

  call ZGETRI(40,factorInv,40,ipiv,work,40,info)
  if (info /= 0) then
!    print *, "************************"                  
!    print *, "D = ", D                                                                   
    stop 'Matrix inversion failed! - accuracy test'
  end if

!  print *, "********************************************"                                 
!  print *, "Eigen-decomposition error = ", norm2(real(matmul(matmul(factor,                
!  exactExp),factorInv)) - real(Acopy))                                                     
!   print *, "***********************************************"                             
!   print *, "Factor inverse: ", norm2(real(matmul(factor, factorInv)) - iden)              
  exactExp = matmul(matmul(factor,exactExp), factorInv)
!  print *, "*********************************"                                           
!  print *, "entry of factor =  ", factor(1,1)
!  print *, "entry of factorInv =  ", factorInv(1,
!  print *, "entry of exact =  ", exactExp(1,1)  
!  print *, "entry of approx =  ", approxexpJac(1,1)                                       
!  stop                  
! Rational approximation                                                                 
!    call matrix_exponential2(JacL, JacD, JacU,20,1.d0, approxexpJac) !    Taylor approx   
    call matrix_exponential(JacL,JacD,JacU,20,1d0,approxexpJac)
 
    error = norm2(real(exactExp) - approxexpJac)
    if (error > 1e-3) then
      write(iulog,*)'WARNING:  Analytic and exact matrix exponentials differ by ', error
      write(iulog,*)'Please check that the exact matrix exponential in eos.F90 is actually exact'
    else
      if (hybrid%masterthread) write(iulog,*)&
          'PASS. max error of analytic and exact matrix exponential: ', error
    end if

    wphi = (/0.25d0,0.0923d0,-0.20394832d0,0.0834d0,&
             0.1249d0,-0.5802302d0,-0.923409328342d0,0.12394d0,&
            -0.99991d0,-0.43494d0,0.779380123d0,0.3298474d0,&
             0.68392d0,-0.0823430d0,0.9832042d0,-0.9823043d0,&
             0.394832042d0,0.1231242d0,0.0021923d0,0.3311d0,&
             0.0398523d0,-0.92384032d0,0.59082342d0,-.009238423d0,&
             0.109248d0,-0.72374d0,-0.87242d0,0.120948d0,&
             0.902482d0,0.29023802d0,-0.398230482d0,-0.459084d0,&
             0.2093482d0,0.59084534d0,-0.23489320d0,0.674829d0,&
             0.39803242d0,-0.234820d0,0.3985023d0,0.98274892d0/)

    Expwphi = wphi
!    call matrix_exponential(JacL,JacD,JacU,20,1d-3,approxexpJac,Expwphi)
    Expwphi = matmul(exactExp,wphi)
    call matrix_exponential_taylorseries(JacL,JacD,JacU,20,1d0,wphi,ExpwphiTS)
    print *, "hey  = ", norm2(Expwphi-ExpwphiTS)
    print *, "hey2 = ", 1d0*(1d0 + norm2(JacL) + norm2(JacD) + norm2(JacU))

  end subroutine test_matrix_exponential_accuracy

  subroutine test_horners_method(hybrid)

    type(hybrid_t)    , intent(in) :: hybrid
    ! local                                                                                   
    real (kind=real_kind) :: coeffs0(4),coeffs1(4),coeffs2(3),coeffs3(7),coeffs(5)
    real (kind=real_kind) :: v1(20),v2(20),vrandom(20),vtemp(20)
    real (kind=real_kind) :: JacL(19),JacU(19),JacD(20),Tri(20,20)
    real (kind=real_kind) :: wphi(40),wphi1(40),wphi2(40),Jac(40,40)
    real (kind = real_kind) :: error
    integer :: k

    if (hybrid%masterthread) write(iulog,*)'Testing accuracy of Horners method...'
    JacL = (/ 3.757951454493374d-3, 3.819276450619194d-3, 3.880299668323407d-3,&
       3.941024255519856d-3, 4.001453427054409d-3, 4.061590011618724d-3,&
       4.121436472079190d-3, 4.180994935960876d-3, 4.240267208588646d-3,&
       4.299254777865632d-3, 4.357958817789391d-3, 4.416380194910784d-3,&
       4.474519480093797d-3, 4.532376966683320d-3, 4.589952695389095d-3,&
       4.647246485637452d-3, 4.704257972792853d-3, 4.760986650405008d-3,&
       4.817431916249652d-3/)

    JacD = (/-7.738077153201838d-3, -7.801334908243546d-3,&
       -7.923983063331936d-3, -8.046030284613884d-3, -8.167482988553411d-3,&
       -8.288347249486747d-3, -8.408628369535837d-3, -8.528330932396275d-3,&
       -8.647458848755431d-3, -8.766015375402829d-3, -8.884003124541974d-3,&
       -9.001424074841970d-3, -9.118279590907854d-3, -9.234570454679255d-3,&
       -9.350296910174005d-3, -9.465458721612778d-3, -9.580055244031766d-3,&
       -9.694085504896373d-3, -9.807548294578283d-3, -9.920442262019286d-3/)

    JacU = (/7.738077153201838d-3, 4.043383453750170d-3, 4.104706612712742d-3,&
       4.165730616290477d-3, 4.226458733033555d-3, 4.286893822432338d-3,&
       4.347038357917113d-3, 4.406894460317085d-3, 4.466463912794555d-3,&
       4.525748166814183d-3, 4.584748346676341d-3, 4.643465257052578d-3,&
       4.701899395997070d-3, 4.760050974585458d-3, 4.817919943490685d-3,&
       4.875506026223682d-3, 4.932808758394314d-3, 4.989827532103520d-3,&
       5.046561644173274d-3/)


    ! form the tridiagonal matrix
    Tri      = 0d0
    Tri(1,1) = JacD(1)
    Tri(1,2) = JacU(1)
    do k = 2,19
      Tri(k,k)   = JacD(k) 
      Tri(k,k-1) = JacL(k-1)
      Tri(k,k+1) = JacU(k)
    end do
    Tri(20,20) = JacD(20)
    Tri(20,19) = JacL(19)

    ! make a random vector
    vrandom = (/0.25d0,0.0923d0,-0.20394832d0,0.0834d0,&
                0.1249d0,-0.5802302d0,-0.923409328342d0,0.12394d0,&
               -0.99991d0,-0.43494d0,0.779380123d0,0.3298474d0,&
                0.68392d0,-0.0823430d0,0.9832042d0,-0.9823043d0,&
                0.394832042d0,0.1231242d0,0.0021923d0,0.3311d0/)
    v1 = 0d0
    call tridiag_times_vector(JacL,JacD,JacU,vrandom,v1,20)
    v2 = matmul(Tri,vrandom)

    error = norm2(v1-v2)
    if (error > 1e-15) then
      write(iulog,*)'WARNING:  tridiagonal left multiply and matmul differ by ', error
      write(iulog,*)'Please check the tridiag_times_vector subroutine in exp_mod for bugs'
    else
    if (hybrid%masterthread) write(iulog,*)&
          'PASS. max error of tridiagonal left multiply and matmul differ by: ', error
    end if

    coeffs(1) = 1d0
    coeffs(2) = 1/2d0
    coeffs(3) = 1/6d0
    coeffs(4) = 1/24d0
    coeffs(5) = 1/120d0

    v1 = 0d0
    v2 = 0d0

    vtemp = vrandom
    v1    = coeffs(1)*vrandom
    do k =2,5
      vtemp = matmul(Tri,vtemp)
      v1    = v1(:) + coeffs(k)*vtemp(:)
    end do
    call horners_method_tridiag_times_vector(JacL,JacD,JacU,vrandom,v2,20,4,coeffs(1:5)) 
    error = norm2(v1-v2)
    if (error > 1e-15) then
      write(iulog,*)'WARNING: Horners method and direct evaluation of a polynomial in a tridiagonal matrix differ by ', error
      write(iulog,*)'Please check the horners_method_tridiag_times_vector subroutine in exp_mod for bugs'
    else
    if (hybrid%masterthread) write(iulog,*)&
          'PASS. Horners method and direct evaluation of a polynomial in a tridiagonal matrix differ by ', error
    end if

    wphi = (/0.25d0,0.0923d0,-0.20394832d0,0.0834d0,&
             0.1249d0,-0.5802302d0,-0.923409328342d0,0.12394d0,&
            -0.99991d0,-0.43494d0,0.779380123d0,0.3298474d0,&
             0.68392d0,-0.0823430d0,0.9832042d0,-0.9823043d0,&
             0.394832042d0,0.1231242d0,0.0021923d0,0.3311d0,&
             0.0398523d0,-0.92384032d0,0.59082342d0,-.009238423d0,&
             0.109248d0,-0.72374d0,-0.87242d0,0.120948d0,&
             0.902482d0,0.29023802d0,-0.398230482d0,-0.459084d0,&
             0.2093482d0,0.59084534d0,-0.23489320d0,0.674829d0,&
             0.39803242d0,-0.234820d0,0.3985023d0,0.98274892d0/)

    ! form the tridiagonal matrix                                                                 
    Jac           = 0d0
    Jac(1,nlev+1) = JacD(1)
    Jac(1,nlev+2) = JacU(1)
    do k = 2,19
      Jac(k,nlev+k)   = JacD(k)
      Jac(k,nlev+k-1) = JacL(k-1)
      Jac(k,nlev+k+1) = JacU(k)
      Jac(nlev+k,k)   = 1d0
    end do
    Jac(20,nlev+20) = JacD(20)
    Jac(20,nlev+19) = JacL(19)
    Jac(40,20)      = 1d0
    Jac(21,1)       = 1d0

    coeffs3(1) = 1d0
    coeffs3(2) = 1d0
    coeffs3(3) = 1d0/2d0
    coeffs3(4) = 1d0/6d0
    coeffs3(5) = 1d0/24d0
    coeffs3(6) = 1d0/120d0
    coeffs3(7) = 1d0/720d0
    coeffs3(8) = 1d0/5040d0
    coeffs3(9) = 1d0/40320d0
    wphi1 = wphi
    wphi2 = coeffs3(1)*wphi
    do k =2,7
      wphi1 = matmul(Jac,wphi1)
      wphi2 = wphi2(:) + coeffs3(k)*wphi1(:)
    end do
    call matrix_exponential_taylorseries(JacL,JacD,JacU,20,1d0,wphi,wphi1)
    if (error > 1e-15) then
      write(iulog,*)'WARNING: matrix_exponential_taylorseries has error = ', error
      write(iulog,*)'Please check the matrix_exponential_taylorseries subroutine in exp_mod for bugs'
    else
    if (hybrid%masterthread) write(iulog,*)&
          'PASS. matrix_exponential_taylorseries has error = ', error
    end if
  end subroutine test_horners_method

end module 

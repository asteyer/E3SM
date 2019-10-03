#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
!
!
!  Man dynamics routines for "theta" nonhydrostatic model
!  Original version: Mark Taylor 2017/1
!  
!  2018/8 TOM sponge layer scaling from P. Lauritzen
!  09/2018: O. Guba  code for new ftypes
!  2018/12: M. Taylor apply forcing assuming nearly constant p 
!
module prim_advance_mod

  use bndry_mod,          only: bndry_exchangev
  use control_mod,        only: dcmip16_mu, dcmip16_mu_s, hypervis_order, hypervis_subcycle,&
    integration, nu, nu_div, nu_p, nu_s, nu_top, prescribed_wind, qsplit, rsplit, test_case,&
    theta_hydrostatic_mode, tstep_type, use_moisture, ftype, theta_advect_form
  use derivative_mod,     only: derivative_t, divergence_sphere, gradient_sphere, laplace_sphere_wk,&
    laplace_z, vorticity_sphere, vlaplace_sphere_wk 
  use derivative_mod,     only: subcell_div_fluxes, subcell_dss_fluxes
  use dimensions_mod,     only: max_corner_elem, nlev, nlevp, np, qsize
  use edge_mod,           only: edge_g, edgevpack_nlyr, edgevunpack_nlyr
  use edgetype_mod,       only: EdgeBuffer_t,  EdgeDescriptor_t, edgedescriptor_t
  use element_mod,        only: element_t
  use element_state,      only: max_itercnt_perstep,avg_itercnt,max_itererr_perstep, nu_scale_top
  use element_ops,        only: get_temperature, set_theta_ref, state0, get_R_star
  use eos,                only: pnh_and_exner_from_eos,phi_from_eos,get_dirk_jacobian, get_exp_jacobian
  use hybrid_mod,         only: hybrid_t
  use hybvcoord_mod,      only: hvcoord_t
  use kinds,              only: iulog, real_kind
  use perf_mod,           only: t_adj_detailf, t_barrierf, t_startf, t_stopf ! _EXTERNAL
  use parallel_mod,       only: abortmp, global_shared_buf, global_shared_sum, iam, parallel_t
  use physical_constants, only: Cp, cp, cpwater_vapor, g, kappa, Rgas, Rwater_vapor, p0 
  use physics_mod,        only: virtual_specific_heat, virtual_temperature
  use prim_si_mod,        only: preq_vertadv_v1
  use reduction_mod,      only: parallelmax, reductionbuffer_ordered_1d_t
  use time_mod,           only: timelevel_qdp, timelevel_t
#ifndef CAM
  use test_mod,           only: set_prescribed_wind
#endif
  use viscosity_theta,    only: biharmonic_wk_theta

#ifdef TRILINOS
    use prim_derived_type_mod ,only : derived_type, initialize
    use, intrinsic :: iso_c_binding
#endif
 
  implicit none
  private
  save
  public :: prim_advance_exp, prim_advance_init1, &
       applycamforcing_dynamics, compute_andor_apply_rhs, matrix_exponential,&
       matrix_exponential2

contains





  subroutine prim_advance_init1(par, elem,integration)
        
    type (parallel_t) :: par
    type (element_t), intent(inout), target   :: elem(:)
    character(len=*)    , intent(in) :: integration
    integer :: i
    integer :: ie


  end subroutine prim_advance_init1



#ifndef ARKODE
  !_____________________________________________________________________
  subroutine prim_advance_exp(elem, deriv, hvcoord, hybrid,dt, tl,  nets, nete, compute_diagnostics)

    type (element_t),      intent(inout), target :: elem(:)
    type (derivative_t),   intent(in)            :: deriv
    type (hvcoord_t)                             :: hvcoord
    type (hybrid_t),       intent(in)            :: hybrid
    real (kind=real_kind), intent(in)            :: dt
    type (TimeLevel_t)   , intent(in)            :: tl
    integer              , intent(in)            :: nets
    integer              , intent(in)            :: nete
    logical,               intent(in)            :: compute_diagnostics

    real (kind=real_kind) :: dt2, time, dt_vis, x, eta_ave_w
    real (kind=real_kind) :: itertol,a1,a2,a3,a4,a5,a6,ahat1,ahat2
    real (kind=real_kind) :: ahat3,ahat4,ahat5,ahat6,dhat1,dhat2,dhat3,dhat4
    real (kind=real_kind) ::  gamma,delta

    integer :: ie,nm1,n0,np1,nstep,qsplit_stage,k, qn0
    integer :: n,i,j,maxiter
 ! New variables that I am adding.
    real (kind=real_kind) :: stage1(nets:nete,np,np,nlevp,6),&
                                stage2(nets:nete,np,np,nlevp,6),&
                                stage3(nets:nete,np,np,nlevp,6),&
                                stage4(nets:nete,np,np,nlevp,6)
    real (kind=real_kind) :: a21,a32,a43,a54,a65,a61,c2,c3,c4,c5,b1,b2,b3,b4
    real (kind=real_kind) :: wphivec(2*nlev)
    real (kind=real_kind), pointer, dimension(:,:,:) :: w_n0
    real (kind=real_kind), pointer, dimension(:,:,:) :: phi_n0
    real (kind=real_kind), pointer, dimension(:,:,:)   :: phi_np1
    real (kind=real_kind), pointer, dimension(:,:,:)   :: dp3d
    real (kind=real_kind), pointer, dimension(:,:,:)   :: vtheta_dp
    real (kind=real_kind), pointer, dimension(:,:)   :: phis
    real (kind=real_kind) :: expJ(2*nlev, 2*nlev), iden(2*nlev,2*nlev)
    real (kind=real_kind) :: JacD(nlev,np,np)  , JacL(nlev-1,np,np)
    real (kind=real_kind) :: JacU(nlev-1,np,np)
    real (kind=real_kind) :: JacD_elem(nlev,np,np,nete-nets+1)
    real (kind=real_kind) :: JacL_elem(nlev-1,np,np,nete-nets+1), JacU_elem(nlev-1,np,np,nete-nets+1)
    real (kind=real_kind) :: pnh(np,np,nlev)     ! nh (nonydro) pressure
    real (kind=real_kind) :: dp3d_i(np,np,nlevp)
    real (kind=real_kind) :: dpnh_dp_i(np,np,nlevp)
    real (kind=real_kind) :: exner(np,np,nlev)     ! exner nh pressure
    real (kind=real_kind) :: dphi(nlev)
    real (kind=real_kind) :: phi_k(2*nlev,np,np,nete-nets+1)
    integer :: ii

    call t_startf('prim_advance_exp')
    nm1   = tl%nm1
    n0    = tl%n0
    np1   = tl%np1
    nstep = tl%nstep

    ! get timelevel for accessing tracer mass Qdp() to compute virtual temperature
    call TimeLevel_Qdp(tl, qsplit, qn0)  ! compute current Qdp() timelevel

! integration = "explicit"
!
!   tstep_type=1  RK2 followed by qsplit-1 leapfrog steps        CFL=close to qsplit
!                    typically requires qsplit=4 or 5
!
!   tstep_type=5  Kinnmark&Gray RK3 5 stage 3rd order            CFL=3.87  (sqrt(15))
!                 From Paul Ullrich.  3rd order for nonlinear terms also
!                 K&G method is only 3rd order for linear
!                 optimal: for windspeeds ~120m/s,gravity: 340m/2
!                 run with qsplit=1
!                 (K&G 2nd order method has CFL=4. tiny CFL improvement not worth 2nd order)
!   tstep_type=6  IMKG243a 
!   tstep_type=7  IMKG254a
!   tstep_type=8  IMKG252a
!

! default weights for computing mean dynamics fluxes
    eta_ave_w = 1d0/qsplit

!   this should not be needed, but in case physics update u without updating w b.c.:
    do ie=nets,nete
       elem(ie)%state%w_i(:,:,nlevp,n0) = (elem(ie)%state%v(:,:,1,nlev,n0)*elem(ie)%derived%gradphis(:,:,1) + &
            elem(ie)%state%v(:,:,2,nlev,n0)*elem(ie)%derived%gradphis(:,:,2))/g
    enddo
 
#ifndef CAM
    ! if "prescribed wind" set dynamics explicitly and skip time-integration
    if (prescribed_wind ==1 ) then
       call set_prescribed_wind(elem,deriv,hybrid,hvcoord,dt,tl,nets,nete,eta_ave_w)
       call t_stopf('prim_advance_exp')
       return
    endif
#endif

    ! ==================================
    ! Take timestep
    ! ==================================
    dt_vis = dt
    if (tstep_type==1) then 
       ! RK2                                                                                                              
       ! forward euler to u(dt/2) = u(0) + (dt/2) RHS(0)  (store in u(np1))                                               
       call compute_andor_apply_rhs(np1,n0,n0,qn0,dt/2,elem,hvcoord,hybrid,&                                              
            deriv,nets,nete,compute_diagnostics,0d0,1.d0,1.d0,1.d0)                                                      
       ! leapfrog:  u(dt) = u(0) + dt RHS(dt/2)     (store in u(np1))                                                     
       call compute_andor_apply_rhs(np1,n0,np1,qn0,dt,elem,hvcoord,hybrid,&                                               
            deriv,nets,nete,.false.,eta_ave_w,1.d0,1.d0,1.d0)                                                             


    else if (tstep_type==4) then ! explicit table from IMEX-KG254  method                                                              
      call compute_andor_apply_rhs(np1,n0,n0,qn0,dt/4,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,0d0,1.d0,1.d0,1.d0)
      call compute_andor_apply_rhs(np1,n0,np1,qn0,dt/6,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,0d0,1.d0,1.d0,1.d0)
      call compute_andor_apply_rhs(np1,n0,np1,qn0,3*dt/8,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,0d0,1.d0,1.d0,1.d0)
      call compute_andor_apply_rhs(np1,n0,np1,qn0,dt/2,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,0d0,1.d0,1.d0,1.d0)
      call compute_andor_apply_rhs(np1,n0,np1,qn0,dt,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,eta_ave_w*1d0,1.d0,1.d0,1.d0)



    else if (tstep_type==5) then
       ! Ullrich 3nd order 5 stage:   CFL=sqrt( 4^2 -1) = 3.87
       ! u1 = u0 + dt/5 RHS(u0)  (save u1 in timelevel nm1)
       call compute_andor_apply_rhs(nm1,n0,n0,qn0,dt/5,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,eta_ave_w/4,1.d0,1.d0,1.d0)
       ! u2 = u0 + dt/5 RHS(u1)
       call compute_andor_apply_rhs(np1,n0,nm1,qn0,dt/5,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0,1.d0,1.d0,1.d0)
       ! u3 = u0 + dt/3 RHS(u2)
       call compute_andor_apply_rhs(np1,n0,np1,qn0,dt/3,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0,1.d0,1.d0,1.d0)
       ! u4 = u0 + 2dt/3 RHS(u3)
       call compute_andor_apply_rhs(np1,n0,np1,qn0,2*dt/3,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0,1.d0,1.d0,1.d0)
       ! compute (5*u1/4 - u0/4) in timelevel nm1:
       do ie=nets,nete
          elem(ie)%state%v(:,:,:,:,nm1)= (5*elem(ie)%state%v(:,:,:,:,nm1) &
               - elem(ie)%state%v(:,:,:,:,n0) ) /4
          elem(ie)%state%vtheta_dp(:,:,:,nm1)= (5*elem(ie)%state%vtheta_dp(:,:,:,nm1) &
               - elem(ie)%state%vtheta_dp(:,:,:,n0) )/4
          elem(ie)%state%dp3d(:,:,:,nm1)= (5*elem(ie)%state%dp3d(:,:,:,nm1) &
                  - elem(ie)%state%dp3d(:,:,:,n0) )/4
          elem(ie)%state%w_i(:,:,1:nlevp,nm1)= (5*elem(ie)%state%w_i(:,:,1:nlevp,nm1) &
                  - elem(ie)%state%w_i(:,:,1:nlevp,n0) )/4
          elem(ie)%state%phinh_i(:,:,1:nlev,nm1)= (5*elem(ie)%state%phinh_i(:,:,1:nlev,nm1) &
                  - elem(ie)%state%phinh_i(:,:,1:nlev,n0) )/4
       enddo
       ! u5 = (5*u1/4 - u0/4) + 3dt/4 RHS(u4)
       call compute_andor_apply_rhs(np1,nm1,np1,qn0,3*dt/4,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,3*eta_ave_w/4,1.d0,1.d0,1.d0)
       ! final method is the same as:
       ! u5 = u0 +  dt/4 RHS(u0)) + 3dt/4 RHS(u4)
!=========================================================================================
    elseif (tstep_type == 6) then  ! IMEX-KG243
 
      a1 = 1d0/4d0
      a2 = 1d0/3d0
      a3 = 1d0/2d0
      a4 = 1d0

      ahat4 = 1d0
      ahat1 = 0d0
      max_itercnt_perstep = 0
      max_itererr_perstep = 0.0  
      ! IMEX-KGNO243
      dhat2 = (1.+sqrt(3.)/3.)/2.
      dhat3 = dhat2
      ahat3 = 1./2.-dhat3
      dhat1 = (ahat3-dhat2+dhat2*dhat3)/(1.-dhat2-dhat3)
      ahat2 = (dhat1-dhat1*dhat3-dhat1*dhat2+dhat1*dhat2*dhat3)/(1.-dhat3)

      call compute_andor_apply_rhs(np1,n0,n0,qn0,dt*a1,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,0d0,1d0)
  

      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat1*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
      max_itercnt_perstep        = max(maxiter,max_itercnt_perstep)
      max_itererr_perstep = max(itertol,max_itererr_perstep)
 
      call compute_andor_apply_rhs(np1,n0,np1,qn0,dt*a2,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,ahat2/a2,1d0)

      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat2*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol) 
      max_itercnt_perstep        = max(maxiter,max_itercnt_perstep)
      max_itererr_perstep = max(itertol,max_itererr_perstep)


      call compute_andor_apply_rhs(np1,n0,np1,qn0,dt*a3,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,ahat3/a3,1d0)

      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat3*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
      max_itercnt_perstep        = max(maxiter,max_itercnt_perstep)
      max_itererr_perstep = max(itertol,max_itererr_perstep)

      call compute_andor_apply_rhs(np1,n0,np1,qn0,dt*a4,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,1d0,ahat4/a4,1d0)

      avg_itercnt = ((nstep)*avg_itercnt + max_itercnt_perstep)/(nstep+1)

!==============================================================================================
    elseif (tstep_type == 7) then  ! imkg254, most robust of the methods
 
      max_itercnt_perstep = 0
      max_itererr_perstep = 0.0

      a1 = 1/4d0
      a2 = 1/6d0
      a3 = 3/8d0
      a4 = 1/2d0
      a5 = 1d0
      ahat5 = 1d0

      ! IMEX-KGO254 most stable coefficients
      dhat2 = 1d0
      dhat3 = 1d0
      dhat4 = 2d0
      ahat4 = 1d0/2d0-dhat4
      dhat1= (ahat4*ahat5 - ahat5*dhat3 - ahat5*dhat2 + dhat3*dhat2+ dhat3*dhat4 + dhat2*dhat4)/&
        (ahat5-dhat3-dhat2-dhat4)
      ahat3 = (- ahat4*ahat5*dhat1 - ahat4*ahat5*dhat2+ ahat5*dhat1*dhat2 + ahat5*dhat1*dhat3 +&
        ahat5*dhat2*dhat3- dhat1*dhat2*dhat3 - dhat1*dhat2*dhat4 - dhat1*dhat3*dhat4- &
        dhat2*dhat3*dhat4)/(-ahat4*ahat5)
      ahat2 = ( - ahat3*ahat4*ahat5*dhat1 + ahat4*ahat5*dhat1*dhat2 -&
        ahat5*dhat1*dhat2*dhat3 + dhat1*dhat2*dhat3*dhat4)/(-ahat3*ahat4*ahat5)

      ! Gets the Jacobian for testing matrix exponential
      dp3d  => elem(1)%state%dp3d(:,:,:,n0)
      vtheta_dp  => elem(1)%state%vtheta_dp(:,:,:,n0)
      phi_np1 => elem(1)%state%phinh_i(:,:,:,n0)

      call pnh_and_exner_from_eos(hvcoord,vtheta_dp,dp3d,phi_np1,pnh,exner,dpnh_dp_i,caller='dirk1')

      call get_exp_jacobian(JacL,JacD,JacU,dp3d,phi_np1,pnh,1)

!      print *, "The Jacobian is: "
!      print *, "JacL = ", JacL(:,1,1)
!      print *, "JacD = ", JacD(:,1,1)
!      print *, "JacU = ", JacU(:,1,1)
!      print *, "***************************************************************"

      call compute_andor_apply_rhs(np1,n0,n0,qn0,a1*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,0d0,1d0)
      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat1*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
      max_itercnt_perstep        = max(maxiter,max_itercnt_perstep)
      max_itererr_perstep = max(itertol,max_itererr_perstep)

      call compute_andor_apply_rhs(np1,n0,np1,qn0,a2*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,ahat2/a2,1d0)
      maxiter=10
      itertol=1e-12
 
      call compute_stage_value_dirk(np1,qn0,dhat2*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
      max_itercnt_perstep        = max(maxiter,max_itercnt_perstep)
      max_itererr_perstep = max(itertol,max_itererr_perstep)

      call compute_andor_apply_rhs(np1,n0,np1,qn0,a3*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,ahat3/a3,1d0)
      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat3*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
      max_itercnt_perstep        = max(maxiter,max_itercnt_perstep)
      max_itererr_perstep = max(itertol,max_itererr_perstep)

      call compute_andor_apply_rhs(np1,n0,np1,qn0,a4*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,ahat4/a4,1d0)
      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat4*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
      max_itercnt_perstep        = max(maxiter,max_itercnt_perstep)
      max_itererr_perstep = max(itertol,max_itererr_perstep)

      call compute_andor_apply_rhs(np1,n0,np1,qn0,a5*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,1d0,ahat5/a5,1d0)

      avg_itercnt = ((nstep)*avg_itercnt + max_itercnt_perstep)/(nstep+1)
!================================================================================
    elseif (tstep_type == 8) then ! IMKG253, might be more efficient than IMKG254, might be a teeny bit bad at coarse resolution

      max_itercnt_perstep = 0
      max_itererr_perstep = 0.0

      a1 = 1/4d0
      a2 = 1/6d0
      a3 = 3/8d0
      a4 = 1/2d0
      a5 = 1d0
      ahat5 = 1d0

      dhat4 = 1d0
      dhat3 = 1d0
      ahat4 = -1d0/2d0
      dhat2= (ahat4*ahat5 - ahat5*dhat1 - ahat5*dhat3 + dhat1*dhat3+ dhat1*dhat4 + dhat3*dhat4)/(ahat5-dhat1-dhat3-dhat4)
      ahat3 = (-ahat4*ahat5*dhat2+ahat5*dhat2*dhat3- dhat2*dhat3*dhat4)/(-ahat4*ahat5)

      call compute_andor_apply_rhs(np1,n0,n0,qn0,a1*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,0d0,1d0)
 
      call compute_andor_apply_rhs(np1,n0,np1,qn0,a2*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,0d0,1d0)
      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat2*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
      max_itercnt_perstep        = max(maxiter,max_itercnt_perstep)
      max_itererr_perstep = max(itertol,max_itererr_perstep)

      call compute_andor_apply_rhs(np1,n0,np1,qn0,a3*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,ahat3/a3,1d0)
      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat3*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
      max_itercnt_perstep        = max(maxiter,max_itercnt_perstep)
      max_itererr_perstep = max(itertol,max_itererr_perstep)

      call compute_andor_apply_rhs(np1,n0,np1,qn0,a4*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,ahat4/a4,1d0)
      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat4*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
      max_itercnt_perstep        = max(maxiter,max_itercnt_perstep)
      max_itererr_perstep = max(itertol,max_itererr_perstep)

      call compute_andor_apply_rhs(np1,n0,np1,qn0,a5*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,1d0,ahat5/a5,1d0)

      avg_itercnt = ((nstep)*avg_itercnt + max_itercnt_perstep)/(nstep+1)
!================================================================================
    elseif (tstep_type == 9) then ! IMKG252, use if feeling lucky, might be bad at coarse resolution
 
      max_itercnt_perstep = 0
      max_itererr_perstep = 0.0

      a1 = 1/4d0
      a2 = 1/6d0
      a3 = 3/8d0
      a4 = 1/2d0
      a5 = 1d0
      ahat5 = 1d0

      dhat3 = (2d0+sqrt(2d0))/2d0
      dhat4 = (2d0+sqrt(2d0))/2d0
      ahat4 = -(1d0+sqrt(2d0))/2d0

      call compute_andor_apply_rhs(np1,n0,n0,qn0,a1*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,0d0,1d0)

      call compute_andor_apply_rhs(np1,n0,np1,qn0,a2*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,0d0,1d0)

      call compute_andor_apply_rhs(np1,n0,np1,qn0,a3*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,0d0,1d0) 
      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat3*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
      max_itercnt_perstep        = max(maxiter,max_itercnt_perstep)
      max_itererr_perstep = max(itertol,max_itererr_perstep)

      call compute_andor_apply_rhs(np1,n0,np1,qn0,a4*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,ahat4/a4,1d0)
      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat4*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
      max_itercnt_perstep        = max(maxiter,max_itercnt_perstep)
      max_itererr_perstep = max(itertol,max_itererr_perstep)

      call compute_andor_apply_rhs(np1,n0,np1,qn0,a5*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,1d0,ahat5/a5,1d0)

      avg_itercnt = ((nstep)*avg_itercnt + max_itercnt_perstep)/(nstep+1)
!================================================================================                                      
    elseif (tstep_type == 10) then ! IMKG232b
      a1 = 1d0/2d0
      a2 = 1d0/2d0
      a3 = 1d0
      ahat3 = 1d0
      dhat2 = .5d0*(2d0+sqrt(2d0))
      dhat1 = dhat2
      ahat2 = -(1d0+sqrt(2d0))/2d0

      call compute_andor_apply_rhs(np1,n0,n0,qn0,dt*a1,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,0d0,1d0)  

      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat1*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
      max_itercnt_perstep        = max(maxiter,max_itercnt_perstep)
      max_itererr_perstep = max(itertol,max_itererr_perstep)

      call compute_andor_apply_rhs(np1,n0,np1,qn0,dt*a2,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,ahat2/a2,1d0)

      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat2*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol) 
      max_itercnt_perstep        = max(maxiter,max_itercnt_perstep)
      max_itererr_perstep = max(itertol,max_itererr_perstep)


      call compute_andor_apply_rhs(np1,n0,np1,qn0,dt*a3,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,ahat3/a3,1d0)
      ! Testing nonlinear_rhs
       ! get Jacobian
      do ie = nets,nete
        dp3d  => elem(ie)%state%dp3d(:,:,:,n0)
        vtheta_dp  => elem(ie)%state%vtheta_dp(:,:,:,n0)
        phi_np1 => elem(ie)%state%phinh_i(:,:,:,n0)
        call pnh_and_exner_from_eos(hvcoord,vtheta_dp,dp3d,phi_np1,pnh,exner,dpnh_dp_i,caller='dirk1')
        call get_exp_jacobian(JacL,JacD,JacU,dp3d,phi_np1,pnh,1)
        JacL_elem(:,:,:,ie) = JacL
        JacU_elem(:,:,:,ie) = JacU
        JacD_elem(:,:,:,ie) = JacD
    
      end do
      call compute_nonlinear_rhs(nm1,nm1,n0,qn0,elem,hvcoord,hybrid,&
          deriv,nets,nete,compute_diagnostics,0.d0,JacL_elem,JacD_elem,JacU_elem, dt)

      avg_itercnt = ((nstep)*avg_itercnt + max_itercnt_perstep)/(nstep+1)

!===================================================================================
    elseif (tstep_type == 11) then ! Integrating factor method
      ! get Jacobian
      do ie = nets,nete
        dp3d       => elem(ie)%state%dp3d(:,:,:,n0)
        vtheta_dp  => elem(ie)%state%vtheta_dp(:,:,:,n0)
        phi_np1    => elem(ie)%state%phinh_i(:,:,:,n0)
        call pnh_and_exner_from_eos(hvcoord,vtheta_dp,dp3d,phi_np1,pnh,exner,dpnh_dp_i,caller='dirk1')
        call get_exp_jacobian(JacL,JacD,JacU,dp3d,phi_np1,pnh,1)
        JacL_elem(:,:,:,ie) = JacL
        JacU_elem(:,:,:,ie) = JacU
        JacD_elem(:,:,:,ie) = JacD
      end do
     ! Compute N(u_m) and store in np1
      call compute_nonlinear_rhs(np1,np1,n0,qn0,elem,hvcoord,hybrid,& 
          deriv,nets,nete,compute_diagnostics,0.d0, JacL_elem, JacD_elem, JacU_elem,dt)
      call linear_combination_of_elem(np1,dt,np1,1.d0,n0,elem,nets,nete)
      ! Compute dt*N(u_m) + u_m and store in np1
        ! Compute e^(dt*Jac)(u_m+alpha*dt*N(u_m)) =: h2
      call expLdtwphi(JacL_elem,JacD_elem,JacU_elem,elem,np1,.false.,dt,nets,nete)
      ! Compute N(h2) and store in nm1
      call compute_nonlinear_rhs(nm1,np1,np1,qn0,elem,hvcoord,hybrid,&
         deriv,nets,nete,compute_diagnostics,0.d0, JacL_elem, JacD_elem, JacU_elem,dt)
     
      ! (h3 = e^(dt*Jac)(u_m) +dt*N(h2))
      call expLdtwphi(JacL_elem,JacD_elem,JacU_elem,elem,n0,.false.,dt,nets,nete)
      call linear_combination_of_elem(np1,dt,nm1,1.d0,n0,elem,nets,nete)

!==========================================================================================================
    elseif (tstep_type == 12) then ! IMKG232b with a call to compute_nonlinear_rhs
      a1 = 1d0/2d0
      a2 = 1d0/2d0
      a3 = 1d0
      ahat3 = 1d0
      dhat2 = .5d0*(2d0+sqrt(2d0))
      dhat1 = dhat2
      ahat2 = -(1d0+sqrt(2d0))/2d0

      call compute_andor_apply_rhs(np1,n0,n0,qn0,dt*a1,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,0d0,1d0)  

      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat1*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
      max_itercnt_perstep        = max(maxiter,max_itercnt_perstep)
      max_itererr_perstep = max(itertol,max_itererr_perstep)

      call compute_andor_apply_rhs(np1,n0,np1,qn0,dt*a2,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,ahat2/a2,1d0)

      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat2*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol) 
      max_itercnt_perstep        = max(maxiter,max_itercnt_perstep)
      max_itererr_perstep = max(itertol,max_itererr_perstep)


      call compute_andor_apply_rhs(np1,n0,np1,qn0,dt*a3,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,ahat3/a3,1d0)

      avg_itercnt = ((nstep)*avg_itercnt + max_itercnt_perstep)/(nstep+1)

     call compute_nonlinear_rhs(np1,np1,n0,qn0,elem,hvcoord,hybrid,&
          deriv,nets,nete,compute_diagnostics,0.d0, JacL, JacD, JacU,dt)  !stores N(h1) in elem(np1)

     call compute_nonlinear_rhs(nm1,n0,n0,qn0,elem,hvcoord,hybrid,&
       deriv,nets,nete,compute_diagnostics,0.d0, JacL, JacD, JacU,dt)

!===================================================================================
    elseif (tstep_type == 13) then ! Integrating factor method
      do ie = nets,nete
      ! get Jacobian
        dp3d       => elem(ie)%state%dp3d(:,:,:,n0)
        vtheta_dp  => elem(ie)%state%vtheta_dp(:,:,:,n0)
        phi_np1    => elem(ie)%state%phinh_i(:,:,:,n0)
        call pnh_and_exner_from_eos(hvcoord,vtheta_dp,dp3d,phi_np1,pnh,exner,dpnh_dp_i,caller='dirk1')
        call get_exp_jacobian(JacL,JacD,JacU,dp3d,phi_np1,pnh,1)
        JacL_elem(:,:,:,ie) = JacL
        JacU_elem(:,:,:,ie) = JacU
        JacD_elem(:,:,:,ie) = JacD
        ! Copy elem(n0) to elem(nm1)
      end do
      call linear_combination_of_elem(nm1, 1.d0, n0, 0.d0, nm1,elem,nets,nete)
      ! Compute exp(Ldt)*v_m and store in nm1 
      call expLdtwphi(JacL_elem,JacD_elem,JacU_elem,elem,nm1,.false.,dt,nets,nete)
      ! Compute N(h1) and store in np1
      call compute_nonlinear_rhs(np1,nm1,nm1,qn0,elem,hvcoord,hybrid,&
         deriv,nets,nete,compute_diagnostics,eta_ave_w, JacL_elem, JacD_elem, JacU_elem,dt)

 !    !!!!!!! Test nonlinear rhs
 !    ! - we have u stored in nm1; N(u) stored in np1; now we calculate F(u) and store in n0
 !    call compute_andor_apply_rhs(n0,nm1,nm1,qn0,1.d0,elem,hvcoord,hybrid,&
 !      deriv,nets,nete,compute_diagnostics,0.d0,1.d0,1.d0,0.d0)
 !      do ie = nets, nete
 !       ! subtract F(u) - N(u); store in np1
 !        call linear_combination_of_elem(np1, 1.d0, n0, -1.d0, np1, elem,nets,nte)
 !      ! Use expJ because it is the right size.
 !        do i = 1, np
 !          do j = 1, np
 !             expJ = 0.d0
 !             do ii = 1,nlev-1
 !               expJ(ii,ii+nlev) = JacD_elem(ii,i,j,ie)        
 !               expJ(ii+1,ii+nlev) = JacL_elem(ii,i,j,ie)
 !               expJ(ii,ii+1+nlev) = JacU_elem(ii,i,j,ie)
 !               expJ(ii+nlev,ii) = 1.d0
 !             end do 
 !           expJ(nlev,2*nlev) = JacD_elem(nlev,i,j,ie)
 !           expJ(2*nlev, nlev) = 1.d0
 !
 !           expJ = expJ*g
 !           
 !           wphivec(1:nlev) = elem(ie)%state%w_i(i,j,1:nlev,nm1)
 !           wphivec(nlev+1:2*nlev) = elem(ie)%state%phinh_i(i,j,1:nlev,nm1)
 !
 !           wphivec = matmul(expJ,wphivec)
 !           ! Subtract (F(u) - N(u)) - Lu -- should be 0.
 !           elem(ie)%state%w_i(i,j,1:nlev,np1) = elem(ie)%state%w_i(i,j,1:nlev,np1) - wphivec(1:nlev)
 !           elem(ie)%state%phinh_i(i,j,1:nlev,np1) = elem(ie)%state%phinh_i(i,j,1:nlev,np1) - wphivec(nlev+1:2*nlev) 
 !           if (norm2(elem(ie)%state%w_i(i,j,1:nlev,np1)) > 10.d-5) then
 !             print *, " w error is ", norm2(elem(ie)%state%w_i(i,j,1:nlev,np1))
 !             stop
 !           end if
 !           if (norm2(elem(ie)%state%phinh_i(i,j,1:nlev,np1)) > 10.d-5) then
 !             print *, " phi error is ", norm2(elem(ie)%state%phinh_i(i,j,1:nlev,np1)) 
 !             stop
 !           end if
 !           if (norm2(elem(ie)%state%dp3d(i,j,:,np1)) > 10.d-5) then
 !             print *, "dp3d error is ", norm2(elem(ie)%state%dp3d(i,j,:,np1))
 !             stop
 !           end if
 !           if (norm2(elem(ie)%state%vtheta_dp(i,j,:,np1)) > 10.d-5) then
 !             print *, "dp3d error is ", norm2(elem(ie)%state%vtheta_dp(i,j,:,np1))
 !             stop
 !           end if
 !           if (norm2(elem(ie)%state%v(i,j,:,:,np1)) > 10.d-5) then
 !             print *, "dp3d error is ", norm2(elem(ie)%state%v(i,j,:,:,np1))
 !             stop
 !           end if
 !         end do
 !       end do
 !
 !     end do
 !     stop "nonlinear test success."
 !     !!!!!! End nonlinear test.

      call linear_combination_of_elem(np1, 1.d0, nm1, dt, np1, elem,nets,nete)
 !     if (elem(ie)%state%phinh_i(i,j,1,np1) < elem(ie)%state%phinh_i(i,j,2,np1)) then
 !       print *, "*******************"
 !       print *, " phi error in h2 at ", ie, i,j
 !       stop
 !     end if
      call compute_nonlinear_rhs(n0, np1, np1, qn0, elem, hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,JacL_elem,JacD_elem,JacU_elem,dt) ! was np1,np1,np1
      call linear_combination_of_elem(np1, 1.d0, nm1, dt, n0, elem, nets,nete) ! was np1 instead of n0
      ! Compute exp(Ldt)*elem(np1) and store in np1.
      call expLdtwphi(JacL_elem,JacD_elem,JacU_elem,elem,np1,.false.,dt,nets,nete)

!==========================================================================================================
!===================================================================================
    elseif (tstep_type == 14) then ! Forward Euler comparison.
      do ie = nets,nete
      ! get Jacobian
        dp3d       => elem(ie)%state%dp3d(:,:,:,n0)
        vtheta_dp  => elem(ie)%state%vtheta_dp(:,:,:,n0)
        phi_np1    => elem(ie)%state%phinh_i(:,:,:,n0)
        call pnh_and_exner_from_eos(hvcoord,vtheta_dp,dp3d,phi_np1,pnh,exner,dpnh_dp_i,caller='dirk1')
        call get_exp_jacobian(JacL,JacD,JacU,dp3d,phi_np1,pnh,1)
        JacL_elem(:,:,:,ie) = JacL(:,:,:)
        JacU_elem(:,:,:,ie) = JacU(:,:,:)
        JacD_elem(:,:,:,ie) = JacD(:,:,:)
      end do
      
      call linear_combination_of_elem(nm1, 1.d0, n0, 0.d0, nm1, elem,nets,nete)
      ! Compute N(h1) and store in nm1
      call compute_nonlinear_rhs(np1,nm1,nm1,qn0,elem,hvcoord,hybrid,&
         deriv,nets,nete,compute_diagnostics, eta_ave_w, JacL_elem, JacD_elem, JacU_elem,dt)
      call linear_combination_of_elem(np1,1.d0,n0,dt,np1,elem,nets,nete)
      call expLdtwphi(JacL_elem,JacD_elem,JacU_elem,elem,np1,.false.,dt,nets,nete)
      ! Standard Forward Euler, stored in nm1.
      call compute_andor_apply_rhs(nm1,n0,n0,qn0,dt,elem,hvcoord,hybrid,&
       deriv,nets,nete,compute_diagnostics,0.d0,1.d0,1.d0,1.d0)
!      ! compare forward euler steps:
!      do ie = nets, nete
!       print *, "ie", ie
!       print *, "max values at nm1 (Standard FE)"
!       print *, "vmaxmax", maxval(abs((elem(ie)%state%v(:,:,:,:,nm1))))
!       print *, "vtheta_dpmax", maxval(abs(elem(ie)%state%vtheta_dp(:,:,:,nm1)))
!       print *, "dp3dmax", maxval(abs(elem(ie)%state%dp3d(:,:,:,nm1)))
!       print *, "wmax", maxval(abs(elem(ie)%state%w_i(:,:,:,nm1)))
!       print *, "phimax", maxval(abs(elem(ie)%state%phinh_i(:,:,:,nm1)))
!
!       print *, "max values at np1"
!       print *, "vmaxmax", maxval(abs((elem(ie)%state%v(:,:,:,:,np1))))
!       print *, "vtheta_dpmax", maxval(abs(elem(ie)%state%vtheta_dp(:,:,:,np1)))
!       print *, "dp3dmax", maxval(abs(elem(ie)%state%dp3d(:,:,:,np1)))
!       print *, "wmax", maxval(abs(elem(ie)%state%w_i(:,:,:,np1)))
!       print *, "phimax", maxval(abs(elem(ie)%state%phinh_i(:,:,:,np1)))
!      end do

!==========================================================================================================
    else if (tstep_type == 15) then ! Andrew's FW Euler

! first compute L in terms of JacL, JacD, JacL
      do ie = nets,nete
        dp3d       => elem(ie)%state%dp3d(:,:,:,n0)
        vtheta_dp  => elem(ie)%state%vtheta_dp(:,:,:,n0)
        phi_np1    => elem(ie)%state%phinh_i(:,:,:,n0)
        call pnh_and_exner_from_eos(hvcoord,vtheta_dp,dp3d,phi_np1,pnh,exner,dpnh_dp_i,caller='dirk1')
        call get_exp_jacobian(JacL,JacD,JacU,dp3d,phi_np1,pnh,1)
        JacL_elem(:,:,:,ie) = JacL(:,:,:)
        JacU_elem(:,:,:,ie) = JacU(:,:,:)
        JacD_elem(:,:,:,ie) = JacD(:,:,:)
      end do

     ! compute N(v(n0)) and store at np1
     call compute_nonlinear_rhs(np1,n0,n0,qn0,elem,hvcoord,hybrid,&
       deriv,nets,nete,compute_diagnostics,eta_ave_w, JacL_elem, JacD_elem, JacU_elem,dt)

     ! compute u(n0) + dt * N(u(n0)) and store at np1
     call linear_combination_of_elem(np1,1d0,n0,dt,np1,elem,nets,nete)

      do ie = nets,nete
        do i = 1,np
          do j = 1,np
            ! grabs w and phi for linear operation
            wphivec(1:nlev)        = elem(ie)%state%w_i(i,j,1:nlev,np1)
            wphivec(1+nlev:2*nlev) = elem(ie)%state%phinh_i(i,j,1:nlev,np1)
            call matrix_exponential(JacL_elem(:,i,j,ie),JacD_elem(:,i,j,ie),JacU_elem(:,i,j,ie),.false.,nlev,dt,expJ,wphivec)
            ! update w and phi after matrix exponential
            elem(ie)%state%w_i(i,j,1:nlev,np1)     = wphivec(1:nlev)
            elem(ie)%state%phinh_i(i,j,1:nlev,np1) = wphivec(1+nlev:2*nlev)
          end do
        end do
      end do



!     do ie =nets,nete
!
!       print *, "ie", ie
!       print *, "max values at n0"
!       print *, "vmaxmax", maxval(abs((elem(ie)%state%v(:,:,:,:,n0))))
!       print *, "vtheta_dpmax", maxval(abs(elem(ie)%state%vtheta_dp(:,:,:,n0)))
!       print *, "dp3dmax", maxval(abs(elem(ie)%state%dp3d(:,:,:,n0)))
!       print *, "wmax", maxval(abs(elem(ie)%state%w_i(:,:,:,n0)))
!       print *, "phimax", maxval(abs(elem(ie)%state%phinh_i(:,:,:,n0)))
!
!       print *, "max values at np1"
!       print *, "vmaxmax", maxval(abs((elem(ie)%state%v(:,:,:,:,np1))))
!       print *, "vtheta_dpmax", maxval(abs(elem(ie)%state%vtheta_dp(:,:,:,np1)))
!       print *, "dp3dmax", maxval(abs(elem(ie)%state%dp3d(:,:,:,np1)))
!       print *, "wmax", maxval(abs(elem(ie)%state%w_i(:,:,:,np1)))
!       print *, "phimax", maxval(abs(elem(ie)%state%phinh_i(:,:,:,np1)))
!
!    end do
!==========================================================================================================
!===================================================================================
    elseif (tstep_type == 16) then ! Integrating factor method - second approach to 13
       a1 = 1.d0/2.d0

       ! Compute JacL, JacD, and JacU
       do ie = nets,nete
        dp3d       => elem(ie)%state%dp3d(:,:,:,n0)
        vtheta_dp  => elem(ie)%state%vtheta_dp(:,:,:,n0)
        phi_np1    => elem(ie)%state%phinh_i(:,:,:,n0)
        call pnh_and_exner_from_eos(hvcoord,vtheta_dp,dp3d,phi_np1,pnh,exner,dpnh_dp_i,caller='dirk1')
        call get_exp_jacobian(JacL,JacD,JacU,dp3d,phi_np1,pnh,1)
        JacL_elem(:,:,:,ie) = JacL(:,:,:)
        JacU_elem(:,:,:,ie) = JacU(:,:,:)
        JacD_elem(:,:,:,ie) = JacD(:,:,:)
      end do
      

      !! g1 = vm = um
      call linear_combination_of_elem(np1,1.d0,n0,0.d0,np1,elem,nets,nete)

      !! Form g2 = vm + dt*N(vm) and store in np1
      call compute_nonlinear_rhs(nm1,np1,np1,qn0,elem,hvcoord,hybrid,&
       deriv,nets,nete,compute_diagnostics,eta_ave_w,JacL_elem,JacD_elem,JacU_elem,0.d0)
      call linear_combination_of_elem(np1,1.d0,n0,dt*a1,nm1,elem,nets,nete)

      !! Form vmp1 = g3 = vm + dt*exp(-Ldt)N(exp(Ldt)g2)
      call expLdtwphi(JacL_elem,JacD_elem,JacU_elem,elem,np1,.false.,dt*a1,nets,nete)
      call compute_nonlinear_rhs(nm1,np1,np1,qn0,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,JacL_elem,JacD_elem,JacU_elem,dt*a1)
      call expLdtwphi(JacL_elem,JacD_elem,JacU_elem,elem,nm1,.true.,dt*a1,nets,nete)
      call linear_combination_of_elem(np1,1.d0,n0,dt,nm1,elem,nets,nete)

  !! Compute ump1 = exp(Ldt)vnmp  = exp(Ldt)g3
      call expLdtwphi(JacL_elem,JacD_elem,JacU_elem,elem,np1,.false.,dt,nets,nete)

!==========================================================================================================
    elseif (tstep_type == 17) then ! Second order RK method
      a1 = 0.5d0
      ! Compute JacL, JacD, and JacU
      do ie = nets,nete
        dp3d       => elem(ie)%state%dp3d(:,:,:,n0)
        vtheta_dp  => elem(ie)%state%vtheta_dp(:,:,:,n0)
        phi_np1    => elem(ie)%state%phinh_i(:,:,:,n0)
        call pnh_and_exner_from_eos(hvcoord,vtheta_dp,dp3d,phi_np1,pnh,exner,dpnh_dp_i,caller='dirk1')
        call get_exp_jacobian(JacL,JacD,JacU,dp3d,phi_np1,pnh,1)
        JacL_elem(:,:,:,ie) = JacL(:,:,:)
        JacU_elem(:,:,:,ie) = JacU(:,:,:)
        JacD_elem(:,:,:,ie) = JacD(:,:,:)
      end do
      !! g1 = v_m = u_m
      call linear_combination_of_elem(np1,1.d0,n0,0.d0,np1,elem,nets,nete)

      !! g2 = vm + 1/2 dt N(v_m)
      call compute_nonlinear_rhs(np1,np1,np1,qn0,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,JacL_elem,JacD_elem,JacU_elem,0.d0)
      call linear_combination_of_elem(np1,1.d0,n0,a1*dt,np1,elem,nets,nete)

      !! g3 = vm + 1/2 dt exp(-L*dt/2)N(exp(L*dt/2)*g2)
      call expLdtwphi(JacL_elem,JacD_elem,JacU_elem,elem,np1,.false.,dt*a1,nets,nete)
      call compute_nonlinear_rhs(np1,np1,np1,qn0,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,JacL_elem,JacD_elem,JacU_elem,dt*a1)
      call expLdtwphi(JacL_elem,JacD_elem,JacU_elem,elem,np1,.true.,dt*a1,nets,nete)
      call linear_combination_of_elem(np1, 1.d0,n0,a1*dt,np1,elem,nets,nete)

      !! g4 = vm + dt*exp(-Ldt/2)*N(exp(Ldt/2)g3)
      call expLdtwphi(JacL_elem,JacD_elem,JacU_elem,elem,np1,.false.,dt*a1,nets,nete)
      call compute_nonlinear_rhs(np1,np1,np1,qn0,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,JacL_elem,JacD_elem,JacU_elem,dt*a1)
      call expLdtwphi(JacL_elem,JacD_elem,JacU_elem,elem,np1,.true.,dt*a1,nets,nete)
      call linear_combination_of_elem(np1,1.d0,n0,dt,np1,elem,nets,nete)

      !! Ump1 = exp(Ldt)vmp1
      call expLdtwphi(JacL_elem,JacD_elem,JacU_elem,elem,np1,.false.,dt,nets,nete)
!==========================================================================================================
    elseif (tstep_type == 18) then ! Third order RK method
      a21 = 1.d0/5.d0
      a32 = 1.d0/5.d0
      a43 = 1.d0/3.d0
      a54 = 2.d0/3.d0
      a61 = 1.d0/4.d0
      a65 = 3.d0/4.d0
      c2 = 1.d0/5.d0
      c3 = 1.d0/5.d0
      c4 = 1.d0/3.d0
      c5 = 2.d0/3.d0

      ! Compute JacL, JacD, and JacU
      do ie = nets,nete
        dp3d       => elem(ie)%state%dp3d(:,:,:,n0)
        vtheta_dp  => elem(ie)%state%vtheta_dp(:,:,:,n0)
        phi_np1    => elem(ie)%state%phinh_i(:,:,:,n0)
        call pnh_and_exner_from_eos(hvcoord,vtheta_dp,dp3d,phi_np1,pnh,exner,dpnh_dp_i,caller='dirk1')
        call get_exp_jacobian(JacL,JacD,JacU,dp3d,phi_np1,pnh,1)
        JacL_elem(:,:,:,ie) = JacL(:,:,:)
        JacU_elem(:,:,:,ie) = JacU(:,:,:)
        JacD_elem(:,:,:,ie) = JacD(:,:,:)
      end do
      !! g1 = v_m = u_m
      call linear_combination_of_elem(np1,1.d0,n0,0.d0,np1,elem,nets,nete) ! move to np1

      !! g2 = vm + a21*dt*N(v_m)
      call compute_nonlinear_rhs(np1,np1,np1,qn0,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,JacL_elem,JacD_elem,JacU_elem,0.d0)
      call linear_combination_of_elem(np1,1.d0,n0,a21*dt,np1,elem,nets,nete)

      !! g3 = vm + a32*dt*exp(-L*dt*c2)N(exp(L*dt*c2)*g2)
      call expLdtwphi(JacL_elem,JacD_elem,JacU_elem,elem,np1,.false.,dt*c2,nets,nete)
      call compute_nonlinear_rhs(np1,np1,np1,qn0,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,JacL_elem,JacD_elem,JacU_elem,dt*c2)
      call expLdtwphi(JacL_elem,JacD_elem,JacU_elem,elem,np1,.true.,dt*c2,nets,nete)
      call linear_combination_of_elem(np1, 1.d0,n0,a32*dt,np1,elem,nets,nete)

      !! g4 = vm + a43*dt*exp(-Ldt*c3)*N(exp(L*dt*c3)g3)
      call expLdtwphi(JacL_elem,JacD_elem,JacU_elem,elem,np1,.false.,dt*c3,nets,nete)
      call compute_nonlinear_rhs(np1,np1,np1,qn0,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,JacL_elem,JacD_elem,JacU_elem,dt*c3)
      call expLdtwphi(JacL_elem,JacD_elem,JacU_elem,elem,np1,.true.,dt*c3,nets,nete)
      call linear_combination_of_elem(np1,1.d0,n0,dt*a43,np1,elem,nets,nete)

      !! g5 = vm + a54*dt*exp(-Ldt*c4)*N(exp(L*dt*c4)g4)
      call expLdtwphi(JacL_elem,JacD_elem,JacU_elem,elem,np1,.false.,dt*c4,nets,nete)
      call compute_nonlinear_rhs(np1,np1,np1,qn0,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,JacL_elem,JacD_elem,JacU_elem,dt*c4)
      call expLdtwphi(JacL_elem,JacD_elem,JacU_elem,elem,np1,.true.,dt*c4,nets,nete)
      call linear_combination_of_elem(np1,1.d0,n0,dt*a54,np1,elem,nets,nete)

      !! g6 = vm + a61*dt*N(g1)+a65*dt*exp(-Ldt*c5)*N(exp(L*dt*c5)g5)
      call expLdtwphi(JacL_elem,JacD_elem,JacU_elem,elem,np1,.false.,dt*c5,nets,nete)
      call compute_nonlinear_rhs(np1,np1,np1,qn0,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,JacL_elem,JacD_elem,JacU_elem,dt*c5)
      call expLdtwphi(JacL_elem,JacD_elem,JacU_elem,elem,np1,.true.,dt*c5,nets,nete)

      call compute_nonlinear_rhs(nm1,n0,n0,qn0,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,JacL_elem,JacD_elem,JacU_elem,0.d0)
      call linear_combination_of_elem(nm1,1.d0,n0,a61*dt,nm1,elem,nets,nete)
      call linear_combination_of_elem(np1,1.d0,nm1,a65*dt,np1,elem,nets,nete)

      !! Ump1 = exp(Ldt)vmp1
      call expLdtwphi(JacL_elem,JacD_elem,JacU_elem,elem,np1,.false.,dt,nets,nete)

!!==========================================================================================================
    elseif (tstep_type == 19) then ! Fourth order RK method
      a21 = 1.d0/2.d0
      a32 = 1.d0/2.d0
      a43 = 1.d0
      c2  = 1.d0/2.d0
      c3  = 1.d0/2.d0
      c4  = 1.d0
      b1  = 1.d0/6.d0
      b2  = 1.d0/3.d0
      b3  = 1.d0/3.d0
      b4  = 1.d0/6.d0

      ! Compute JacL, JacD, and JacU
      do ie = nets,nete
        dp3d       => elem(ie)%state%dp3d(:,:,:,n0)
        vtheta_dp  => elem(ie)%state%vtheta_dp(:,:,:,n0)
        phi_np1    => elem(ie)%state%phinh_i(:,:,:,n0)
        call pnh_and_exner_from_eos(hvcoord,vtheta_dp,dp3d,phi_np1,pnh,exner,dpnh_dp_i,caller='dirk1')
        call get_exp_jacobian(JacL,JacD,JacU,dp3d,phi_np1,pnh,1)
        JacL_elem(:,:,:,ie) = JacL(:,:,:)
        JacU_elem(:,:,:,ie) = JacU(:,:,:)
        JacD_elem(:,:,:,ie) = JacD(:,:,:)
      end do
      !! g1 = v_m = u_m
      call linear_combination_of_elem(np1,1.d0,n0,0.d0,np1,elem,nets,nete) ! move to np1

      !! g2 = vm + a21*dt*N(v_m)
      call compute_nonlinear_rhs(np1,np1,np1,qn0,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,JacL_elem,JacD_elem,JacU_elem,0.d0)
      call store_state(elem,np1,nets,nete,stage1)
      call linear_combination_of_elem(np1,1.d0,n0,a21*dt,np1,elem,nets,nete)

      !! g3 = vm + a32*dt*exp(-L*dt*c2)N(exp(L*dt*c2)*g2)
      call expLdtwphi(JacL_elem,JacD_elem,JacU_elem,elem,np1,.false.,dt*c2,nets,nete)
      call compute_nonlinear_rhs(np1,np1,np1,qn0,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,JacL_elem,JacD_elem,JacU_elem,dt*c2)
      call expLdtwphi(JacL_elem,JacD_elem,JacU_elem,elem,np1,.true.,dt*c2,nets,nete)
      call store_state(elem,np1,nets,nete,stage2)
      call linear_combination_of_elem(np1, 1.d0,n0,a32*dt,np1,elem,nets,nete)

      !! g4 = vm + a43*dt*exp(-Ldt*c3)*N(exp(L*dt*c3)g3)
      call expLdtwphi(JacL_elem,JacD_elem,JacU_elem,elem,np1,.false.,dt*c3,nets,nete)
      call compute_nonlinear_rhs(np1,np1,np1,qn0,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,JacL_elem,JacD_elem,JacU_elem,dt*c3)
      call expLdtwphi(JacL_elem,JacD_elem,JacU_elem,elem,np1,.true.,dt*c3,nets,nete)
      call store_state(elem,np1,nets,nete,stage3)
      call linear_combination_of_elem(np1,1.d0,n0,dt*a43,np1,elem,nets,nete)

      !! Calculate exp(-Ldt*c4)*N(exp(L*dt*c4)g4)
      call expLdtwphi(JacL_elem,JacD_elem,JacU_elem,elem,np1,.false.,dt*c4,nets,nete)
      call compute_nonlinear_rhs(np1,np1,np1,qn0,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,JacL_elem,JacD_elem,JacU_elem,dt*c4)
      call expLdtwphi(JacL_elem,JacD_elem,JacU_elem,elem,np1,.true.,dt*c4,nets,nete)
      call store_state(elem,np1,nets,nete,stage4)

      !! vmp1 = vm + b1*dt*N(g1) + b2*dt*exp(-Ldt*c2)*N(exp(Ldt*c2)g2) +
      !             b3*dt*exp(-Ldt*c3)*N(exp(Ldt*c3)g3) + b4*dt*exp(-Ldt*c4)*N(exp(Ldt*c4)g4)
           ! First, vm + b1 term
      call retrieve_state(stage1,elem,nm1,nets,nete)
      call linear_combination_of_elem(np1,1.d0,n0,b1*dt,nm1,elem,nets,nete) 
           ! Second, add b2 term
      call retrieve_state(stage2,elem,nm1,nets,nete)
      call linear_combination_of_elem(np1,1.d0,np1,b2*dt,nm1,elem,nets,nete)
           ! Add b3 term
      call retrieve_state(stage3,elem,nm1,nets,nete)
      call linear_combination_of_elem(np1,1.d0,np1,b3*dt,nm1,elem,nets,nete)
           ! Add b4 term
      call retrieve_state(stage4,elem,nm1,nets,nete)
      call linear_combination_of_elem(np1,1.d0,np1,b4*dt,nm1,elem,nets,nete)

      !! Ump1 = exp(Ldt)vmp1
      call expLdtwphi(JacL_elem,JacD_elem,JacU_elem,elem,np1,.false.,dt,nets,nete)

!!==========================================================================================================
!!==========================================================================================================
    elseif (tstep_type == 20) then ! First order ETD Method

      ! Compute JacL, JacD, and JacU
      do ie = nets,nete
        dp3d       => elem(ie)%state%dp3d(:,:,:,n0)
        vtheta_dp  => elem(ie)%state%vtheta_dp(:,:,:,n0)
        phi_np1    => elem(ie)%state%phinh_i(:,:,:,n0)
        call pnh_and_exner_from_eos(hvcoord,vtheta_dp,dp3d,phi_np1,pnh,exner,dpnh_dp_i,caller='dirk1')
        call get_exp_jacobian(JacL,JacD,JacU,dp3d,phi_np1,pnh,1)
        JacL_elem(:,:,:,ie) = JacL(:,:,:)
        JacU_elem(:,:,:,ie) = JacU(:,:,:)
        JacD_elem(:,:,:,ie) = JacD(:,:,:)
      end do

      call linear_combination_of_elem(np1,1.d0,n0,0.d0,np1,elem,nets,nete) ! move to np1
      ! Compute N(u_m)
      call compute_nonlinear_rhs(np1,np1,np1,qn0,elem,hvcoord,hybrid,&
       deriv,nets,nete,compute_diagnostics,eta_ave_w, JacL_elem, JacD_elem, JacU_elem,0.d0)
      call expLdtwphi(JacL_elem,JacD_elem,JacU_elem,elem,n0,.false.,dt,nets,nete)

      call apply_phi_func(JacL_elem,JacD_elem,JacU_elem,dt,1,np1,elem,nets,nete)
      call linear_combination_of_elem(np1,1.d0,n0,dt,np1,elem,nets,nete)

!!==========================================================================================================
!!==========================================================================================================
    elseif (tstep_type == 21) then ! Second order ETD Method
      c2 = 3.d0/4.d0

      ! Compute JacL, JacD, and JacU
      do ie = nets,nete
        dp3d       => elem(ie)%state%dp3d(:,:,:,n0)
        vtheta_dp  => elem(ie)%state%vtheta_dp(:,:,:,n0)
        phi_np1    => elem(ie)%state%phinh_i(:,:,:,n0)
        call pnh_and_exner_from_eos(hvcoord,vtheta_dp,dp3d,phi_np1,pnh,exner,dpnh_dp_i,caller='dirk1')
        call get_exp_jacobian(JacL,JacD,JacU,dp3d,phi_np1,pnh,1)
        JacL_elem(:,:,:,ie) = JacL(:,:,:)
        JacU_elem(:,:,:,ie) = JacU(:,:,:)
        JacD_elem(:,:,:,ie) = JacD(:,:,:)
      end do

      ! Stage1 = u_m is in n0
      ! Calculate Stage2 = exp(Ldtc2)u_m + c2dt*phi1(c2dtL)N(u_m;t_m)
      call compute_nonlinear_rhs(np1,n0,n0,qn0,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,JacL_elem,JacD_elem,JacU_elem,0.d0) ! Stores N(u_m) in np1
      call apply_phi_func(JacL_elem,JacD_elem,JacU_elem,c2*dt,1,np1,elem,nets,nete)
      call linear_combination_of_elem(nm1,1.d0,n0,0.d0,nm1,elem,nets,nete)
      call expLdtwphi(JacL_elem,JacD_elem,JacU_elem,elem,nm1,.false.,c2*dt,nets,nete) ! exp(Ldt)u_m is in nm1
      call linear_combination_of_elem(np1,1.d0,nm1,c2*dt,np1,elem,nets,nete) !Stage 2 is in np1

      ! Calculate Stage3 = exp(Ldt)u_m + dt[phi1(dtL)N(u_m;t_m)-1/c2phi2(dtL)N(u_m;t_m)+1/c2phi2(dtL)N(Stage2;tm+c2dt)]
      call compute_nonlinear_rhs(np1,np1,np1,qn0,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,JacL_elem,JacD_elem,JacU_elem,c2*dt)
      call apply_phi_func(JacL_elem,JacD_elem,JacU_elem,dt,2,np1,elem,nets,nete)
      call linear_combination_of_elem(np1,dt/c2,np1,0.d0,n0,elem,nets,nete) ! dt/c2phi2(dtL)N(Stage2;tm+c2dt) is in np1

      call compute_nonlinear_rhs(nm1,n0,n0,qn0,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,JacL_elem,JacD_elem,JacU_elem,0.d0)
      call apply_phi_func(JacL_elem,JacD_elem,JacU_elem,dt,2,nm1,elem,nets,nete)
      call linear_combination_of_elem(np1,-dt/c2,nm1,1.d0,np1,elem,nets,nete) !Last two terms stored in np1

      call compute_nonlinear_rhs(nm1,n0,n0,qn0,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,JacL_elem,JacD_elem,JacU_elem,0.d0)
      call apply_phi_func(JacL_elem,JacD_elem,JacU_elem,dt,1,nm1,elem,nets,nete)
      call linear_combination_of_elem(np1,dt,nm1,1.d0,np1,elem,nets,nete) ! Last three terms stored in np1

      call linear_combination_of_elem(nm1,1.d0,n0,0.d0,nm1,elem,nets,nete)
      call expLdtwphi(JacL_elem,JacD_elem,JacU_elem,elem,nm1,.false.,dt,nets,nete)
      call linear_combination_of_elem(np1,1.d0,nm1,1.d0,np1,elem,nets,nete)  ! Stage 3 stored in np1


      ! Calculate ump1 = exp(Ldt)u_m + dt[phi1(dtL)N(um;tm)-phi2(dtL)N(um;tm)+phi2(dtL)N(Stage3;tm+dt)]
      call compute_nonlinear_rhs(np1,np1,np1,qn0,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,JacL_elem,JacD_elem,JacU_elem,dt)
      call apply_phi_func(JacL_elem,JacD_elem,JacU_elem,dt,2,np1,elem,nets,nete)

      call linear_combination_of_elem(nm1,1.d0,n0,0.d0,nm1,elem,nets,nete)
      call expLdtwphi(JacL_elem,JacD_elem,JacU_elem,elem,nm1,.false.,dt,nets,nete)
      call linear_combination_of_elem(np1,1.d0,nm1,dt,np1,elem,nets,nete) ! First and last terms stored in np1

      call compute_nonlinear_rhs(nm1,n0,n0,qn0,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,JacL_elem,JacD_elem,JacU_elem,0.d0)
      call apply_phi_func(JacL_elem,JacD_elem,JacU_elem,dt,1,nm1,elem,nets,nete)
      call linear_combination_of_elem(np1,1.d0,np1,dt,nm1,elem,nets,nete) ! Second term added to np1

      call compute_nonlinear_rhs(nm1,n0,n0,qn0,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,JacL_elem,JacD_elem,JacU_elem,0.d0)
      call apply_phi_func(JacL_elem,JacD_elem,JacU_elem,dt,2,nm1,elem,nets,nete)
      call linear_combination_of_elem(np1,1.d0,np1,-dt,nm1,elem,nets,nete)



!!==========================================================================================================


 !==========================================================================================================
    else
      call abortmp('ERROR: bad choice of tstep_type')
    endif


    ! ==============================================
    ! Time-split Horizontal diffusion: nu.del^2 or nu.del^4
    ! U(*) = U(t+1)  + dt2 * HYPER_DIFF_TERM(t+1)
    ! ==============================================
    ! note:time step computes u(t+1)= u(t*) + RHS.
    ! for consistency, dt_vis = t-1 - t*, so this is timestep method dependent
    ! forward-in-time, hypervis applied to dp3d
    if (hypervis_order == 2 .and. nu>0) &
         call advance_hypervis(elem,hvcoord,hybrid,deriv,np1,nets,nete,dt_vis,eta_ave_w)




    ! warning: advance_physical_vis currently requires levels that are equally spaced in z
    if (dcmip16_mu>0) call advance_physical_vis(elem,hvcoord,hybrid,deriv,np1,nets,nete,dt,dcmip16_mu_s,dcmip16_mu)

    call t_stopf('prim_advance_exp')
  end subroutine prim_advance_exp


#else

  !_____________________________________________________________________
  subroutine prim_advance_exp(elem, deriv, hvcoord, hybrid,dt, tl,  nets, nete, compute_diagnostics)
  !
  ! version of prim_advance_exp which uses ARKODE for timestepping
  !
    use arkode_mod,     only: parameter_list, update_arkode, get_solution_ptr, &
                              table_list, set_Butcher_tables, &
                              calc_nonlinear_stats, update_nonlinear_stats, &
                              rel_tol, abs_tol, use_column_solver
    use iso_c_binding

    type (element_t),      intent(inout), target :: elem(:)
    type (derivative_t),   intent(in)            :: deriv
    type (hvcoord_t)                             :: hvcoord
    type (hybrid_t),       intent(in)            :: hybrid
    real (kind=real_kind), intent(in)            :: dt
    type (TimeLevel_t)   , intent(in)            :: tl
    integer              , intent(in)            :: nets
    integer              , intent(in)            :: nete
    logical,               intent(in)            :: compute_diagnostics

    real (kind=real_kind) :: dt2, time, dt_vis, x, eta_ave_w
    real (kind=real_kind) :: itertol,a1,a2,a3,a4,a5,a6,ahat1,ahat2
    real (kind=real_kind) :: ahat3,ahat4,ahat5,ahat6,dhat1,dhat2,dhat3,dhat4
    real (kind=real_kind) ::  gamma,delta

    integer :: ie,nm1,n0,np1,nstep,qsplit_stage,k, qn0
    integer :: n,i,j,maxiter,sumiter
 
    type(parameter_list) :: arkode_parameters
    type(table_list) :: arkode_tables
    type(c_ptr) :: ynp1
    real(real_kind) :: tout, t
    integer(C_INT) :: ierr, itask

    call t_startf('prim_advance_exp')
    nm1   = tl%nm1
    n0    = tl%n0
    np1   = tl%np1
    nstep = tl%nstep

    ! get timelevel for accessing tracer mass Qdp() to compute virtual temperature
    call TimeLevel_Qdp(tl, qsplit, qn0)  ! compute current Qdp() timelevel

! integration = "explicit"
!
!   tstep_type=1  RK2 followed by qsplit-1 leapfrog steps        CFL=close to qsplit
!                    typically requires qsplit=4 or 5
!
!   tstep_type=5  Kinnmark&Gray RK3 5 stage 3rd order            CFL=3.87  (sqrt(15))
!                 From Paul Ullrich.  3rd order for nonlinear terms also
!                 K&G method is only 3rd order for linear
!                 optimal: for windspeeds ~120m/s,gravity: 340m/2
!                 run with qsplit=1
!                 (K&G 2nd order method has CFL=4. tiny CFL improvement not worth 2nd order)
!   tstep_type=6  IMEX-KG243 
!   tstep_type=7  IMEX-KG254
!   

! default weights for computing mean dynamics fluxes
    eta_ave_w = 1d0/qsplit

!   this should not be needed, but in case physics update u without updating w b.c.:
    do ie=nets,nete
       elem(ie)%state%w_i(:,:,nlevp,n0) = (elem(ie)%state%v(:,:,1,nlev,n0)*elem(ie)%derived%gradphis(:,:,1) + &
            elem(ie)%state%v(:,:,2,nlev,n0)*elem(ie)%derived%gradphis(:,:,2))/g
    enddo
 
#ifndef CAM
    ! if "prescribed wind" set dynamics explicitly and skip time-integration
    if (prescribed_wind ==1 ) then
       call set_prescribed_wind(elem,deriv,hybrid,hvcoord,dt,tl,nets,nete,eta_ave_w)
       call t_stopf('prim_advance_exp')
       return
    endif
#endif

    ! ==================================
    ! Take timestep
    ! ==================================
    dt_vis = dt
    if (tstep_type==1) then 
       ! RK2                                                                                                              
       ! forward euler to u(dt/2) = u(0) + (dt/2) RHS(0)  (store in u(np1))                                               
       call compute_andor_apply_rhs(np1,n0,n0,qn0,dt/2,elem,hvcoord,hybrid,&                                              
            deriv,nets,nete,compute_diagnostics,0d0,1.d0,1.d0,1.d0)                                                      
       ! leapfrog:  u(dt) = u(0) + dt RHS(dt/2)     (store in u(np1))                                                     
       call compute_andor_apply_rhs(np1,n0,np1,qn0,dt,elem,hvcoord,hybrid,&                                               
            deriv,nets,nete,.false.,eta_ave_w,1.d0,1.d0,1.d0)                                                             


    else if (tstep_type==4) then ! explicit table from IMEX-KG254  method                                                              
      call compute_andor_apply_rhs(np1,n0,n0,qn0,dt/4,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,0d0,1.d0,1.d0,1.d0)
      call compute_andor_apply_rhs(np1,n0,np1,qn0,dt/6,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,0d0,1.d0,1.d0,1.d0)
      call compute_andor_apply_rhs(np1,n0,np1,qn0,3*dt/8,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,0d0,1.d0,1.d0,1.d0)
      call compute_andor_apply_rhs(np1,n0,np1,qn0,dt/2,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,0d0,1.d0,1.d0,1.d0)
      call compute_andor_apply_rhs(np1,n0,np1,qn0,dt,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,eta_ave_w*1d0,1.d0,1.d0,1.d0)



    else if (tstep_type==5) then
       ! Ullrich 3nd order 5 stage:   CFL=sqrt( 4^2 -1) = 3.87
       ! u1 = u0 + dt/5 RHS(u0)  (save u1 in timelevel nm1)
       call compute_andor_apply_rhs(nm1,n0,n0,qn0,dt/5,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,eta_ave_w/4,1.d0,1.d0,1.d0)
       ! u2 = u0 + dt/5 RHS(u1)
       call compute_andor_apply_rhs(np1,n0,nm1,qn0,dt/5,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0,1.d0,1.d0,1.d0)
       ! u3 = u0 + dt/3 RHS(u2)
       call compute_andor_apply_rhs(np1,n0,np1,qn0,dt/3,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0,1.d0,1.d0,1.d0)
       ! u4 = u0 + 2dt/3 RHS(u3)
       call compute_andor_apply_rhs(np1,n0,np1,qn0,2*dt/3,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,0d0,1.d0,1.d0,1.d0)
       ! compute (5*u1/4 - u0/4) in timelevel nm1:
       do ie=nets,nete
          elem(ie)%state%v(:,:,:,:,nm1)= (5*elem(ie)%state%v(:,:,:,:,nm1) &
               - elem(ie)%state%v(:,:,:,:,n0) ) /4
          elem(ie)%state%vtheta_dp(:,:,:,nm1)= (5*elem(ie)%state%vtheta_dp(:,:,:,nm1) &
               - elem(ie)%state%vtheta_dp(:,:,:,n0) )/4
          elem(ie)%state%dp3d(:,:,:,nm1)= (5*elem(ie)%state%dp3d(:,:,:,nm1) &
                  - elem(ie)%state%dp3d(:,:,:,n0) )/4
          elem(ie)%state%w_i(:,:,1:nlevp,nm1)= (5*elem(ie)%state%w_i(:,:,1:nlevp,nm1) &
                  - elem(ie)%state%w_i(:,:,1:nlevp,n0) )/4
          elem(ie)%state%phinh_i(:,:,1:nlev,nm1)= (5*elem(ie)%state%phinh_i(:,:,1:nlev,nm1) &
                  - elem(ie)%state%phinh_i(:,:,1:nlev,n0) )/4
       enddo
       ! u5 = (5*u1/4 - u0/4) + 3dt/4 RHS(u4)
       call compute_andor_apply_rhs(np1,nm1,np1,qn0,3*dt/4,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,3*eta_ave_w/4,1.d0,1.d0,1.d0)
       ! final method is the same as:
       ! u5 = u0 +  dt/4 RHS(u0)) + 3dt/4 RHS(u4)
!=========================================================================================
    elseif (tstep_type == 6) then  ! IMEX-KG243

      a1 = 1./4.
      a2 = 1./3.
      a3 = 1./2.
      a4 = 1.0

      ahat4 = 1.
      ahat1 = 0.
      max_itercnt_perstep = 0
      max_itererr_perstep = 0.0  
      ! IMEX-KGNO243
      dhat2 = (1.+sqrt(3.)/3.)/2.
      dhat3 = dhat2
      ahat3 = 1./2.-dhat3
      dhat1 = (ahat3-dhat2+dhat2*dhat3)/(1.-dhat2-dhat3)
      ahat2 = (dhat1-dhat1*dhat3-dhat1*dhat2+dhat1*dhat2*dhat3)/(1.-dhat3)

      call compute_andor_apply_rhs(np1,n0,n0,qn0,dt*a1,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,0d0,1d0)
  
      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat1*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
      sumiter = maxiter 

 
      call compute_andor_apply_rhs(np1,n0,np1,qn0,dt*a2,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,ahat2/a2,1d0)

      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat2*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol) 
      sumiter = sumiter + maxiter

      call compute_andor_apply_rhs(np1,n0,np1,qn0,dt*a3,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,ahat3/a3,1d0)

      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat3*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
      sumiter = sumiter + maxiter


      call compute_andor_apply_rhs(np1,n0,np1,qn0,dt*a4,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,1d0,ahat4/a4,1d0)
      if (calc_nonlinear_stats) then
        call update_nonlinear_stats(1, sumiter)
      end if


!==============================================================================================
!==============================================================================================
    elseif (tstep_type == 7) then  ! imkg254, most robust of the methods
 
      max_itercnt_perstep = 0
      max_itererr_perstep = 0.0

      a1 = 1/4d0
      a2 = 1/6d0
      a3 = 3/8d0
      a4 = 1/2d0
      a5 = 1d0
      ahat5 = 1d0

      ! IMEX-KGO254 most stable coefficients
      dhat2 = 1d0
      dhat3 = 1d0
      dhat4 = 2d0
      ahat4 = 1d0/2d0-dhat4
      dhat1= (ahat4*ahat5 - ahat5*dhat3 - ahat5*dhat2 + dhat3*dhat2+ dhat3*dhat4 + dhat2*dhat4)/&
        (ahat5-dhat3-dhat2-dhat4)
      ahat3 = (- ahat4*ahat5*dhat1 - ahat4*ahat5*dhat2+ ahat5*dhat1*dhat2 + ahat5*dhat1*dhat3 +&
        ahat5*dhat2*dhat3- dhat1*dhat2*dhat3 - dhat1*dhat2*dhat4 - dhat1*dhat3*dhat4- &
        dhat2*dhat3*dhat4)/(-ahat4*ahat5)
      ahat2 = ( - ahat3*ahat4*ahat5*dhat1 + ahat4*ahat5*dhat1*dhat2 -&
        ahat5*dhat1*dhat2*dhat3 + dhat1*dhat2*dhat3*dhat4)/(-ahat3*ahat4*ahat5)

      call compute_andor_apply_rhs(np1,n0,n0,qn0,a1*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,0d0,1d0)
      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat1*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
      max_itercnt_perstep        = max(maxiter,max_itercnt_perstep)
      max_itererr_perstep = max(itertol,max_itererr_perstep)

      call compute_andor_apply_rhs(np1,n0,np1,qn0,a2*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,ahat2/a2,1d0)
      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat2*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
      max_itercnt_perstep        = max(maxiter,max_itercnt_perstep)
      max_itererr_perstep = max(itertol,max_itererr_perstep)

      call compute_andor_apply_rhs(np1,n0,np1,qn0,a3*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,ahat3/a3,1d0)
      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat3*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
      max_itercnt_perstep        = max(maxiter,max_itercnt_perstep)
      max_itererr_perstep = max(itertol,max_itererr_perstep)

      call compute_andor_apply_rhs(np1,n0,np1,qn0,a4*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,ahat4/a4,1d0)
      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat4*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
      max_itercnt_perstep        = max(maxiter,max_itercnt_perstep)
      max_itererr_perstep = max(itertol,max_itererr_perstep)

      call compute_andor_apply_rhs(np1,n0,np1,qn0,a5*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,1d0,ahat5/a5,1d0)

      avg_itercnt = ((nstep)*avg_itercnt + max_itercnt_perstep)/(nstep+1)
!================================================================================
    elseif (tstep_type == 8) then ! IMKG253, might be more efficient than IMKG254, might be a teeny bit bad at coarse resolution

      max_itercnt_perstep = 0
      max_itererr_perstep = 0.0

      a1 = 1/4d0
      a2 = 1/6d0
      a3 = 3/8d0
      a4 = 1/2d0
      a5 = 1d0
      ahat5 = 1d0

      dhat4 = 1d0
      dhat3 = 1d0
      ahat4 = -1d0/2d0
      dhat2= (ahat4*ahat5 - ahat5*dhat1 - ahat5*dhat3 + dhat1*dhat3+ dhat1*dhat4 + dhat3*dhat4)/(ahat5-dhat1-dhat3-dhat4)
      ahat3 = (-ahat4*ahat5*dhat2+ahat5*dhat2*dhat3- dhat2*dhat3*dhat4)/(-ahat4*ahat5)

      call compute_andor_apply_rhs(np1,n0,n0,qn0,a1*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,0d0,1d0)
 
      call compute_andor_apply_rhs(np1,n0,np1,qn0,a2*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,0d0,1d0)
      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat2*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
      max_itercnt_perstep        = max(maxiter,max_itercnt_perstep)
      max_itererr_perstep = max(itertol,max_itererr_perstep)

      call compute_andor_apply_rhs(np1,n0,np1,qn0,a3*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,ahat3/a3,1d0)
      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat3*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
      max_itercnt_perstep        = max(maxiter,max_itercnt_perstep)
      max_itererr_perstep = max(itertol,max_itererr_perstep)

      call compute_andor_apply_rhs(np1,n0,np1,qn0,a4*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,ahat4/a4,1d0)
      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat4*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
      max_itercnt_perstep        = max(maxiter,max_itercnt_perstep)
      max_itererr_perstep = max(itertol,max_itererr_perstep)

      call compute_andor_apply_rhs(np1,n0,np1,qn0,a5*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,1d0,ahat5/a5,1d0)

      avg_itercnt = ((nstep)*avg_itercnt + max_itercnt_perstep)/(nstep+1)
!================================================================================
    elseif (tstep_type == 9) then ! IMKG252, use if feeling lucky, might be bad at coarse resolution
 
      max_itercnt_perstep = 0
      max_itererr_perstep = 0.0

      a1 = 1/4d0
      a2 = 1/6d0
      a3 = 3/8d0
      a4 = 1/2d0
      a5 = 1d0
      ahat5 = 1d0

      dhat3 = (2d0+sqrt(2d0))/2d0
      dhat4 = (2d0+sqrt(2d0))/2d0
      ahat4 = -(1d0+sqrt(2d0))/2d0

      call compute_andor_apply_rhs(np1,n0,n0,qn0,a1*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,0d0,1d0)

      call compute_andor_apply_rhs(np1,n0,np1,qn0,a2*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,0d0,1d0)

      call compute_andor_apply_rhs(np1,n0,np1,qn0,a3*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,0d0,1d0) 
      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat3*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
      max_itercnt_perstep        = max(maxiter,max_itercnt_perstep)
      max_itererr_perstep = max(itertol,max_itererr_perstep)

      call compute_andor_apply_rhs(np1,n0,np1,qn0,a4*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,ahat4/a4,1d0)
      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat4*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
      max_itercnt_perstep        = max(maxiter,max_itercnt_perstep)
      max_itererr_perstep = max(itertol,max_itererr_perstep)

      call compute_andor_apply_rhs(np1,n0,np1,qn0,a5*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,eta_ave_w,1d0,ahat5/a5,1d0)

      avg_itercnt = ((nstep)*avg_itercnt + max_itercnt_perstep)/(nstep+1)
! ==================================================================================
    elseif (tstep_type == 10) then ! IMKG232b
      a1 = 1d0/2d0
      a2 = 1d0/2d0
      a3 = 1d0
      ahat3 = 1d0
      dhat2 = .5d0*(2d0+sqrt(2d0))
      dhat1 = dhat2
      ahat2 = -(1d0+sqrt(2d0))/2d0

      call compute_andor_apply_rhs(np1,n0,n0,qn0,dt*a1,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,0d0,1d0)  

      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat1*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol)
      max_itercnt_perstep        = max(maxiter,max_itercnt_perstep)
      max_itererr_perstep = max(itertol,max_itererr_perstep)

      call compute_andor_apply_rhs(np1,n0,np1,qn0,dt*a2,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,ahat2/a2,1d0)

      maxiter=10
      itertol=1e-12
      call compute_stage_value_dirk(np1,qn0,dhat2*dt,elem,hvcoord,hybrid,&
        deriv,nets,nete,maxiter,itertol) 
      max_itercnt_perstep        = max(maxiter,max_itercnt_perstep)
      max_itererr_perstep = max(itertol,max_itererr_perstep)


      call compute_andor_apply_rhs(np1,n0,np1,qn0,dt*a3,elem,hvcoord,hybrid,&
        deriv,nets,nete,compute_diagnostics,0d0,1d0,ahat3/a3,1d0)

      avg_itercnt = ((nstep)*avg_itercnt + max_itercnt_perstep)/(nstep+1)

!=========================================================================================
    else if (tstep_type==20) then ! ARKode RK2
      call set_Butcher_tables(arkode_parameters, arkode_tables%RK2)

    else if (tstep_type==21) then ! ARKode Kinnmark, Gray, Ullrich 3rd-order, 5-stage
      call set_Butcher_tables(arkode_parameters, arkode_tables%KGU35)

    else if (tstep_type==22) then ! ARKode Ascher 2nd/2nd/2nd-order, 3-stage
      call set_Butcher_tables(arkode_parameters, arkode_tables%ARS232)

    else if (tstep_type==23) then ! ARKode Candidate ARK453 Method
      call set_Butcher_tables(arkode_parameters, arkode_tables%ARK453)

    else if (tstep_type==24) then ! ARKode Ascher 2nd/2nd/2nd-order, 3-stage
      call set_Butcher_tables(arkode_parameters, arkode_tables%ARS222)

    else if (tstep_type==25) then ! ARKode Ascher 3rd/4th/3rd-order, 3-stage
      call set_Butcher_tables(arkode_parameters, arkode_tables%ARS233)

    else if (tstep_type==26) then ! ARKode Ascher 3rd/3rd/3rd-order, 4-stage
      call set_Butcher_tables(arkode_parameters, arkode_tables%ARS343)

    else if (tstep_type==27) then ! ARKode Ascher 3rd/3rd/3rd-order, 5-stage
      call set_Butcher_tables(arkode_parameters, arkode_tables%ARS443)

    else if (tstep_type==28) then ! ARKode Kennedy 3rd/3rd/3rd-order, 4-stage
      call set_Butcher_tables(arkode_parameters, arkode_tables%ARK324)

    else if (tstep_type==29) then ! ARKode Kennedy 4th/4th/4th-order, 6-stage
      call set_Butcher_tables(arkode_parameters, arkode_tables%ARK436)

    else if (tstep_type==30) then ! ARKode Conde et al ssp3(3,3,3)a (renamed here)
      call set_Butcher_tables(arkode_parameters, arkode_tables%SSP3333B)

    else if (tstep_type==31) then ! ARKode Conde et al ssp3(3,3,3)b (renamed here)
      call set_Butcher_tables(arkode_parameters, arkode_tables%SSP3333C)

    else if (tstep_type==32) then ! ARKode IMKG 2nd-order, 4 stage (2 implicit)
      call set_Butcher_tables(arkode_parameters, arkode_tables%IMKG232)

    else if (tstep_type==33) then ! ARKode IMKG 2nd-order, 5 stage (2 implicit)
      call set_Butcher_tables(arkode_parameters, arkode_tables%IMKG242)

    else if (tstep_type==34) then ! ARKode IMKG 2nd-order, 5 stage (3 implicit)
      call set_Butcher_tables(arkode_parameters, arkode_tables%IMKG243)

    else if (tstep_type==35) then ! ARKode IMKG 2nd-order, 6 stage (2 implicit)
      call set_Butcher_tables(arkode_parameters, arkode_tables%IMKG252)

    else if (tstep_type==36) then ! ARKode IMKG 2nd-order, 6 stage (3 implicit)
      call set_Butcher_tables(arkode_parameters, arkode_tables%IMKG253)

    else if (tstep_type==37) then ! ARKode IMKG 2nd-order, 6 stage (4 implicit)
      call set_Butcher_tables(arkode_parameters, arkode_tables%IMKG254)

    else if (tstep_type==38) then ! ARKode IMKG 3rd-order, 5 stage (2 implicit)
      call set_Butcher_tables(arkode_parameters, arkode_tables%IMKG342)

    else if (tstep_type==39) then ! ARKode IMKG 3rd-order, 5 stage (3 implicit)
      call set_Butcher_tables(arkode_parameters, arkode_tables%IMKG343)

    else if (tstep_type==40) then ! ARKode IMKG 3rd-order, 6 stage (3 implicit)
      call set_Butcher_tables(arkode_parameters, arkode_tables%IMKG353)

    else if (tstep_type==41) then ! ARKode IMKG 3rd-order, 6 stage (4 implicit)
      call set_Butcher_tables(arkode_parameters, arkode_tables%IMKG354)

    else 
       call abortmp('ERROR: bad choice of tstep_type')
    endif

    ! Use ARKode to advance solution
    if (tstep_type >= 20) then

      ! If implicit solves are involved, set corresponding parameters
      if (arkode_parameters%imex /= 1) then
        ! linear solver parameters
        if (.not.use_column_solver) then
          arkode_parameters%precLR = 0 ! no preconditioning
          arkode_parameters%gstype = 1 ! classical Gram-Schmidt orthogonalization
          arkode_parameters%lintol = 0.05d0 ! multiplies NLCOV_COEF in linear conv. criteria
        end if
        ! Iteration tolerances (appear in WRMS array as rtol*|u_i| + atol_i)
        arkode_parameters%rtol = rel_tol
        if (abs_tol < 0.d0) then
          arkode_parameters%atol(1) = 1.d1*arkode_parameters%rtol ! assumes u ~ 1e1
          arkode_parameters%atol(2) = 1.d1*arkode_parameters%rtol ! assumes v ~ 1e1
          arkode_parameters%atol(3) = 1.d1*arkode_parameters%rtol ! assumes w_i ~ 1e1
          arkode_parameters%atol(4) = 1.d5*arkode_parameters%rtol ! assumes phinh_i ~ 1e5
          arkode_parameters%atol(5) = 1.d6*arkode_parameters%rtol ! assumes vtheta_dp ~ 1e6
          arkode_parameters%atol(6) = 1.d0*arkode_parameters%rtol ! assumes dp3d ~ 1e0
        else
          arkode_parameters%atol(:) = abs_tol
        end if
      end if

      ! update ARKode solver
      call update_arkode(elem, nets, nete, deriv, hvcoord, hybrid, &
                               dt, eta_ave_w, n0, qn0, arkode_parameters)

      ! call ARKode to perform a single step
      call get_solution_ptr(np1, ynp1)
      tout = dt
      itask = 2          ! use 'one-step' mode
      call farkode(tout, t, ynp1, itask, ierr)
      if (ierr /= 0) then
        call abortmp('farkode failed')
      endif
      if (calc_nonlinear_stats) then
        call update_nonlinear_stats()
      end if
    end if

    ! ==============================================
    ! Time-split Horizontal diffusion: nu.del^2 or nu.del^4
    ! U(*) = U(t+1)  + dt2 * HYPER_DIFF_TERM(t+1)
    ! ==============================================
    ! note:time step computes u(t+1)= u(t*) + RHS.
    ! for consistency, dt_vis = t-1 - t*, so this is timestep method dependent
    ! forward-in-time, hypervis applied to dp3d
    if (hypervis_order == 2 .and. nu>0) &
         call advance_hypervis(elem,hvcoord,hybrid,deriv,np1,nets,nete,dt_vis,eta_ave_w)




    ! warning: advance_physical_vis currently requires levels that are equally spaced in z
    if (dcmip16_mu>0) call advance_physical_vis(elem,hvcoord,hybrid,deriv,np1,nets,nete,dt,dcmip16_mu_s,dcmip16_mu)

    call t_stopf('prim_advance_exp')
  end subroutine prim_advance_exp
#endif

!----------------------------- APPLYCAMFORCING-DYNAMICS ----------------------------

  subroutine applyCAMforcing_dynamics(elem,hvcoord,np1,dt,nets,nete)

  type (element_t)     ,  intent(inout) :: elem(:)
  real (kind=real_kind),  intent(in)    :: dt
  type (hvcoord_t),       intent(in)    :: hvcoord
  integer,                intent(in)    :: np1,nets,nete

  integer :: k,ie
  do ie=nets,nete

     elem(ie)%state%vtheta_dp(:,:,:,np1) = elem(ie)%state%vtheta_dp(:,:,:,np1) + dt*elem(ie)%derived%FVTheta(:,:,:)
     elem(ie)%state%phinh_i(:,:,1:nlev,np1) = elem(ie)%state%phinh_i(:,:,1:nlev,np1) + dt*elem(ie)%derived%FPHI(:,:,1:nlev)

     elem(ie)%state%v(:,:,:,:,np1) = elem(ie)%state%v(:,:,:,:,np1) + dt*elem(ie)%derived%FM(:,:,1:2,:)
     elem(ie)%state%w_i(:,:,1:nlev,np1) = elem(ie)%state%w_i(:,:,1:nlev,np1) + dt*elem(ie)%derived%FM(:,:,3,:)

     ! finally update w at the surface: 
     elem(ie)%state%w_i(:,:,nlevp,np1) = (elem(ie)%state%v(:,:,1,nlev,np1)*elem(ie)%derived%gradphis(:,:,1) + &
          elem(ie)%state%v(:,:,2,nlev,np1)*elem(ie)%derived%gradphis(:,:,2))/g
  enddo
  
  end subroutine applyCAMforcing_dynamics


!----------------------------- ADVANCE-HYPERVIS ----------------------------

  subroutine advance_hypervis(elem,hvcoord,hybrid,deriv,nt,nets,nete,dt2,eta_ave_w)
  !
  !  take one timestep of:
  !          u(:,:,:,np) = u(:,:,:,np) +  dt2*nu*laplacian**order ( u )
  !          T(:,:,:,np) = T(:,:,:,np) +  dt2*nu_s*laplacian**order ( T )
  !
  !
  !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
  !
  !

  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (derivative_t)  , intent(in) :: deriv
  type (hvcoord_t), intent(in)      :: hvcoord

  real (kind=real_kind) :: dt2
  integer :: nets,nete

  ! local
  real (kind=real_kind) :: eta_ave_w  ! weighting for mean flux terms
  integer :: k2,k,kptr,i,j,ie,ic,nt,nlyr_tot,ssize
  real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)      :: vtens
  real (kind=real_kind), dimension(np,np,nlev,4,nets:nete)      :: stens  ! dp3d,theta,w,phi


! NOTE: PGI compiler bug: when using spheremp, rspheremp and ps as pointers to elem(ie)% members,
  !       data is incorrect (offset by a few numbers actually)
  !       removed for now.
  !       real (kind=real_kind), dimension(:,:), pointer :: spheremp,rspheremp
  !       real (kind=real_kind), dimension(:,:,:), pointer   :: ps

  real (kind=real_kind), dimension(np,np,4) :: lap_s  ! dp3d,theta,w,phi
  real (kind=real_kind), dimension(np,np,2) :: lap_v
  real (kind=real_kind) :: exner0(nlev)
  real (kind=real_kind) :: heating(np,np,nlev)
  real (kind=real_kind) :: exner(np,np,nlev)
  real (kind=real_kind) :: pnh(np,np,nlevp)    
  real (kind=real_kind) :: temp(np,np,nlev)    
  real (kind=real_kind) :: temp_i(np,np,nlevp)    
  real (kind=real_kind) :: dt
  real (kind=real_kind) :: ps_ref(np,np)

  real (kind=real_kind) :: theta_ref(np,np,nlev,nets:nete)
  real (kind=real_kind) :: phi_ref(np,np,nlevp,nets:nete)
  real (kind=real_kind) :: dp_ref(np,np,nlev,nets:nete)

  real (kind=real_kind) :: v0(np,np,2,nlevp),eps
  real (kind=real_kind) :: w0(np,np,nlevp)
  integer :: l1p,l2p,l1n,l2n,l
  call t_startf('advance_hypervis')

  dt=dt2/hypervis_subcycle

  if (theta_hydrostatic_mode) then
     nlyr_tot=4*nlev        ! dont bother to dss w_i and phinh_i
     ssize=2*nlev
  else
     nlyr_tot=6*nlev  ! total amount of data for DSS
     ssize=4*nlev
  endif
  
  do k=1,nlev
     exner0(k) = (hvcoord%etam(k)*hvcoord%ps0/p0 )**kappa
  enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NOTE1:  Diffusion works best when applied to theta instead of theta_dp
! It creates some TOM noise when applied to vtheta_dp in DCMIP 2.0 test
! so we convert from vtheta_dp->theta, and then convert back at the end of diffusion
!
! NOTE2: in dcmip2012 test2.0, using theta_ref does improve solution, but
!        phi_ref has no impact
!
! NOTE3: in HS w/topo tests, theta_ref(dp_ref) as opposed to computing it as a 
!        function of dp3d is less noisy at cube edges
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute reference states
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do ie=nets,nete
     !ps_ref(:,:) = hvcoord%hyai(1)*hvcoord%ps0 + sum(elem(ie)%state%dp3d(:,:,:,nt),3)
     ps_ref(:,:) = hvcoord%ps0 * exp ( -elem(ie)%state%phis(:,:)/(Rgas*300))  ! 300K ref temperature
     !ps_ref(:,:) = hvcoord%ps0 - 11.3*elem(ie)%state%phis(:,:)/g
     do k=1,nlev
        dp_ref(:,:,k,ie) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
             (hvcoord%hybi(k+1)-hvcoord%hybi(k))*ps_ref(:,:)
     enddo


     ! phi_ref,theta_ref depend only on ps:
     call set_theta_ref(hvcoord,dp_ref(:,:,:,ie),theta_ref(:,:,:,ie))
     temp(:,:,:)=theta_ref(:,:,:,ie)*dp_ref(:,:,:,ie) 
     call phi_from_eos(hvcoord,elem(ie)%state%phis,&
          temp(:,:,:),dp_ref(:,:,:,ie),phi_ref(:,:,:,ie))
#if 0
     ! no reference state, for testing
     theta_ref(:,:,:,ie)=0
     phi_ref(:,:,:,ie)=0
     dp_ref(:,:,:,ie)=0
#endif

     ! convert vtheta_dp -> theta
     do k=1,nlev
        elem(ie)%state%vtheta_dp(:,:,k,nt)=&
             elem(ie)%state%vtheta_dp(:,:,k,nt)/elem(ie)%state%dp3d(:,:,k,nt)
     enddo
  enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  hyper viscosity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do ic=1,hypervis_subcycle
     do ie=nets,nete
         do k=1,nlev
           elem(ie)%state%vtheta_dp(:,:,k,nt)=elem(ie)%state%vtheta_dp(:,:,k,nt)-&
                theta_ref(:,:,k,ie)
           elem(ie)%state%phinh_i(:,:,k,nt)=elem(ie)%state%phinh_i(:,:,k,nt)-&
                phi_ref(:,:,k,ie)
           elem(ie)%state%dp3d(:,:,k,nt)=elem(ie)%state%dp3d(:,:,k,nt)-&
                dp_ref(:,:,k,ie)
        enddo
     enddo
     
     call biharmonic_wk_theta(elem,stens,vtens,deriv,edge_g,hybrid,nt,nets,nete)
     
     do ie=nets,nete
        
        ! comptue mean flux
        if (nu_p>0) then
           elem(ie)%derived%dpdiss_ave(:,:,:)=elem(ie)%derived%dpdiss_ave(:,:,:)+&
                (eta_ave_w*elem(ie)%state%dp3d(:,:,:,nt)+dp_ref(:,:,:,ie))/hypervis_subcycle
           elem(ie)%derived%dpdiss_biharmonic(:,:,:)=elem(ie)%derived%dpdiss_biharmonic(:,:,:)+&
                eta_ave_w*stens(:,:,:,1,ie)/hypervis_subcycle
        endif
        do k=1,nlev
           ! advace in time.
           ! note: DSS commutes with time stepping, so we can time advance and then DSS.
           ! note: weak operators alreayd have mass matrix "included"
           
           ! biharmonic terms need a negative sign:
           if (nu_top>0 .and. nu_scale_top(k)>1) then
              ! add regular diffusion near top
              lap_s(:,:,1)=laplace_sphere_wk(elem(ie)%state%dp3d       (:,:,k,nt),deriv,elem(ie),var_coef=.false.)
              lap_s(:,:,2)=laplace_sphere_wk(elem(ie)%state%vtheta_dp  (:,:,k,nt),deriv,elem(ie),var_coef=.false.)
              lap_s(:,:,3)=laplace_sphere_wk(elem(ie)%state%w_i        (:,:,k,nt),deriv,elem(ie),var_coef=.false.)
              lap_s(:,:,4)=laplace_sphere_wk(elem(ie)%state%phinh_i    (:,:,k,nt),deriv,elem(ie),var_coef=.false.)
              lap_v=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.)
              
              vtens(:,:,:,k,ie)=(  -nu*vtens(:,:,:,k,ie) + nu_scale_top(k)*nu_top*lap_v(:,:,:)) ! u and v
              stens(:,:,k,1,ie)=(-nu_p*stens(:,:,k,1,ie) + nu_scale_top(k)*nu_top*lap_s(:,:,1)) ! dp3d
              stens(:,:,k,2,ie)=(  -nu*stens(:,:,k,2,ie) + nu_scale_top(k)*nu_top*lap_s(:,:,2)) ! theta
              stens(:,:,k,3,ie)=(  -nu*stens(:,:,k,3,ie) + nu_scale_top(k)*nu_top*lap_s(:,:,3)) ! w
              stens(:,:,k,4,ie)=(-nu_s*stens(:,:,k,4,ie) + nu_scale_top(k)*nu_top*lap_s(:,:,4)) ! phi
           else
              vtens(:,:,:,k,ie)=-nu  *vtens(:,:,:,k,ie) ! u,v
              stens(:,:,k,1,ie)=-nu_p*stens(:,:,k,1,ie) ! dp3d
              stens(:,:,k,2,ie)=-nu  *stens(:,:,k,2,ie) ! theta
              stens(:,:,k,3,ie)=-nu  *stens(:,:,k,3,ie) ! w
              stens(:,:,k,4,ie)=-nu_s*stens(:,:,k,4,ie) ! phi
           endif
           
        enddo
        
        kptr=0;      call edgeVpack_nlyr(edge_g,elem(ie)%desc,vtens(:,:,:,:,ie),2*nlev,kptr,nlyr_tot)
        kptr=2*nlev; call edgeVpack_nlyr(edge_g,elem(ie)%desc,stens(:,:,:,:,ie),ssize,kptr,nlyr_tot)
        
     enddo
     
     call t_startf('ahdp_bexchV2')
     call bndry_exchangeV(hybrid,edge_g)
     call t_stopf('ahdp_bexchV2')
     
     do ie=nets,nete
        
        kptr=0
        call edgeVunpack_nlyr(edge_g,elem(ie)%desc,vtens(:,:,:,:,ie),2*nlev,kptr,nlyr_tot)
        kptr=2*nlev
        call edgeVunpack_nlyr(edge_g,elem(ie)%desc,stens(:,:,:,:,ie),ssize,kptr,nlyr_tot)
        
        
        ! apply inverse mass matrix, accumulate tendencies
        do k=1,nlev
           vtens(:,:,1,k,ie)=dt*vtens(:,:,1,k,ie)*elem(ie)%rspheremp(:,:)  ! u
           vtens(:,:,2,k,ie)=dt*vtens(:,:,2,k,ie)*elem(ie)%rspheremp(:,:)  ! v
           stens(:,:,k,1,ie)=dt*stens(:,:,k,1,ie)*elem(ie)%rspheremp(:,:)  ! dp3d
           stens(:,:,k,2,ie)=dt*stens(:,:,k,2,ie)*elem(ie)%rspheremp(:,:)  ! theta
           stens(:,:,k,3,ie)=dt*stens(:,:,k,3,ie)*elem(ie)%rspheremp(:,:)  ! w
           stens(:,:,k,4,ie)=dt*stens(:,:,k,4,ie)*elem(ie)%rspheremp(:,:)  ! phi
           
           !add ref state back
           elem(ie)%state%vtheta_dp(:,:,k,nt)=elem(ie)%state%vtheta_dp(:,:,k,nt)+&
                theta_ref(:,:,k,ie)
           elem(ie)%state%phinh_i(:,:,k,nt)=elem(ie)%state%phinh_i(:,:,k,nt)+&
                phi_ref(:,:,k,ie)
           elem(ie)%state%dp3d(:,:,k,nt)=elem(ie)%state%dp3d(:,:,k,nt)+&
                dp_ref(:,:,k,ie)
           
        enddo
        
        do k=1,nlev
           elem(ie)%state%v(:,:,:,k,nt)=elem(ie)%state%v(:,:,:,k,nt) + &
                vtens(:,:,:,k,ie)
           elem(ie)%state%w_i(:,:,k,nt)=elem(ie)%state%w_i(:,:,k,nt) &
                +stens(:,:,k,3,ie)
           
           elem(ie)%state%dp3d(:,:,k,nt)=elem(ie)%state%dp3d(:,:,k,nt) &
                +stens(:,:,k,1,ie)
           
           elem(ie)%state%phinh_i(:,:,k,nt)=elem(ie)%state%phinh_i(:,:,k,nt) &
                +stens(:,:,k,4,ie)
        enddo
        
        
        ! apply heating after updating state.  using updated v gives better results in PREQX model
        !
        ! d(IE)/dt =  cp*exner*d(Theta)/dt + phi d(dp3d)/dt   (Theta = dp3d*theta)
        !   Our eqation:  d(theta)/dt = diss(theta) - heating
        !   Assuming no diffusion on dp3d, we can approximate by:
        !   d(IE)/dt = exner*cp*dp3d * diss(theta)  - exner*cp*dp3d*heating               
        !
        ! KE dissipaiton will be given by:
        !   d(KE)/dt = dp3d*U dot diss(U)
        ! we want exner*cp*dp3d*heating = dp3d*U dot diss(U)
        ! and thus heating =  U dot diss(U) / exner*cp
        ! 
        ! compute exner needed for heating term and IE scaling
        ! this is using a mixture of data before viscosity and after viscosity 
        temp(:,:,:)=elem(ie)%state%vtheta_dp(:,:,:,nt)*elem(ie)%state%dp3d(:,:,:,nt)
        call pnh_and_exner_from_eos(hvcoord,temp,&
             elem(ie)%state%dp3d(:,:,:,nt),elem(ie)%state%phinh_i(:,:,:,nt),&
             pnh,exner,temp_i,caller='advance_hypervis')
        
        do k=1,nlev
           k2=max(k,nlev)
           if (theta_hydrostatic_mode) then
              heating(:,:,k)= (elem(ie)%state%v(:,:,1,k,nt)*vtens(:,:,1,k,ie) + &
                   elem(ie)%state%v(:,:,2,k,nt)*vtens(:,:,2,k,ie) ) / &
                   (exner(:,:,k)*Cp)

           else
              heating(:,:,k)= (elem(ie)%state%v(:,:,1,k,nt)*vtens(:,:,1,k,ie) + &
                   elem(ie)%state%v(:,:,2,k,nt)*vtens(:,:,2,k,ie)  +&
                   (elem(ie)%state%w_i(:,:,k,nt)*stens(:,:,k,3,ie)  +&
                     elem(ie)%state%w_i(:,:,k2,nt)*stens(:,:,k2,3,ie))/2 ) /  &
                   (exner(:,:,k)*Cp)  
           endif
           elem(ie)%state%vtheta_dp(:,:,k,nt)=elem(ie)%state%vtheta_dp(:,:,k,nt) &
                +stens(:,:,k,2,ie)*hvcoord%dp0(k)*exner0(k)/(exner(:,:,k)*elem(ie)%state%dp3d(:,:,k,nt)&
                )  -heating(:,:,k)
        enddo
     enddo ! ie
  enddo  ! subcycle

! convert vtheta_dp -> theta
  do ie=nets,nete            
     elem(ie)%state%vtheta_dp(:,:,:,nt)=&
          elem(ie)%state%vtheta_dp(:,:,:,nt)*elem(ie)%state%dp3d(:,:,:,nt)
    
     ! finally update w at the surface: 
     elem(ie)%state%w_i(:,:,nlevp,nt) = (elem(ie)%state%v(:,:,1,nlev,nt)*elem(ie)%derived%gradphis(:,:,1) + &
          elem(ie)%state%v(:,:,2,nlev,nt)*elem(ie)%derived%gradphis(:,:,2))/g
  enddo	

  call t_stopf('advance_hypervis')

  end subroutine advance_hypervis





  subroutine advance_physical_vis(elem,hvcoord,hybrid,deriv,nt,nets,nete,dt,mu_s,mu)
  !
  !  take one timestep of of physical viscosity (single laplace operator) for
  !  all state variables in both horizontal and vertical
  !  
  !  as of 2017/5, used only for the supercell test case
  !  so for now:
  !     dont bother to optimize
  !     apply only to perturbation from background state (supercell initial condition)
  !     uniform spacing in z with delz = 20km/nlev
  !
  !

  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (derivative_t)  , intent(in) :: deriv
  type (hvcoord_t), intent(in)      :: hvcoord

  real (kind=real_kind) :: dt, mu_s, mu
  integer :: nt,nets,nete

  ! local
  integer :: k,kptr,ie
  real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)      :: vtens
  real (kind=real_kind), dimension(np,np,nlev,4,nets:nete)      :: stens  ! dp3d,theta,w,phi
  real (kind=real_kind), dimension(np,np,nlevp,2,nets:nete)     :: stens_i 


  real (kind=real_kind), dimension(np,np,2) :: lap_v
  real (kind=real_kind) :: delz,delz_i

  real (kind=real_kind) :: theta_ref(np,np,nlev)

  real (kind=real_kind) :: theta_prime(np,np,nlev)
  real (kind=real_kind) :: phi_prime(np,np,nlevp)
  real (kind=real_kind) :: dp_prime(np,np,nlev)
  real (kind=real_kind) :: w_prime(np,np,nlevp)
  real (kind=real_kind) :: u_prime(np,np,2,nlev)

  !if(test_case .ne. 'dcmip2016_test3') call abortmp("dcmip16_mu is currently limited to dcmip16 test 3")

  call t_startf('advance_physical_vis')
  delz = 20d3/nlev
  delz_i = 20d3/nlevp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute reference states
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do ie=nets,nete

     theta_ref(:,:,:) = state0(ie)%vtheta_dp(:,:,:,1)/state0(ie)%dp3d(:,:,:,1)
     elem(ie)%state%vtheta_dp(:,:,:,nt)=&
             elem(ie)%state%vtheta_dp(:,:,:,nt)/elem(ie)%state%dp3d(:,:,:,nt)

     ! perturbation variables
     u_prime(:,:,:,:)  = elem(ie)%state%v(:,:,:,:,nt)          -state0(ie)%v(:,:,:,:,1)
     w_prime(:,:,:)    = elem(ie)%state%w_i(:,:,:,nt)     -state0(ie)%w_i(:,:,:,1)
     dp_prime(:,:,:)   = elem(ie)%state%dp3d(:,:,:,nt)         -state0(ie)%dp3d(:,:,:,1)
     phi_prime(:,:,:)  = elem(ie)%state%phinh_i(:,:,:,nt) -state0(ie)%phinh_i(:,:,:,1)
     theta_prime(:,:,:)= elem(ie)%state%vtheta_dp(:,:,:,nt)-theta_ref(:,:,:)

     ! vertical viscosity
     call laplace_z(u_prime,    vtens(:,:,:,:,ie),2,nlev,delz)
     call laplace_z(dp_prime,   stens(:,:,:,1,ie),1,nlev,delz)
     call laplace_z(theta_prime,stens(:,:,:,2,ie),1,nlev,delz)
     call laplace_z(w_prime,    stens_i(:,:,:,1,ie),1,nlevp,delz_i)
     call laplace_z(phi_prime,  stens_i(:,:,:,2,ie),1,nlevp,delz_i)

     ! add in horizontal viscosity
     ! multiply by mass matrix for DSS
     ! horiz viscosity already has mass matrix built in
     ! for interface quantities, only use 1:nlev (dont apply at surface)
     do k=1,nlev
        lap_v = vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.)

        vtens(:,:,1,k,ie) = (vtens(:,:,1,k,ie)*elem(ie)%spheremp(:,:) + lap_v(:,:,1))
        vtens(:,:,2,k,ie) = (vtens(:,:,2,k,ie)*elem(ie)%spheremp(:,:) + lap_v(:,:,2))

        stens(:,:,k,1,ie) = (stens(:,:,k,1,ie)*elem(ie)%spheremp(:,:) + &
             laplace_sphere_wk(elem(ie)%state%dp3d(:,:,k,nt),deriv,elem(ie),var_coef=.false.)  )

        stens(:,:,k,2,ie) = (stens(:,:,k,2,ie)*elem(ie)%spheremp(:,:) + &
             laplace_sphere_wk(elem(ie)%state%vtheta_dp(:,:,k,nt),deriv,elem(ie),var_coef=.false.)  )

        stens(:,:,k,3,ie) = (stens_i(:,:,k,1,ie)*elem(ie)%spheremp(:,:) + &
             laplace_sphere_wk(elem(ie)%state%w_i(:,:,k,nt),deriv,elem(ie),var_coef=.false.) )

        stens(:,:,k,4,ie) = (stens_i(:,:,k,2,ie)*elem(ie)%spheremp(:,:) + &
             laplace_sphere_wk(elem(ie)%state%phinh_i(:,:,k,nt),deriv,elem(ie),var_coef=.false.) ) 

     enddo

     kptr=0
     call edgeVpack_nlyr(edge_g,elem(ie)%desc,vtens(:,:,:,:,ie),2*nlev,kptr,6*nlev)
     kptr=2*nlev
     call edgeVpack_nlyr(edge_g,elem(ie)%desc,stens(:,:,:,:,ie),4*nlev,kptr,6*nlev)
     
  enddo

  call bndry_exchangeV(hybrid,edge_g)
  
  do ie=nets,nete
     
     kptr=0
     call edgeVunpack_nlyr(edge_g,elem(ie)%desc,vtens(:,:,:,:,ie),2*nlev,kptr,6*nlev)
     kptr=2*nlev
     call edgeVunpack_nlyr(edge_g,elem(ie)%desc,stens(:,:,:,:,ie),4*nlev,kptr,6*nlev)
     
     ! apply inverse mass matrix, accumulate tendencies
     do k=1,nlev
        elem(ie)%state%v(:,:,1,k,nt)=elem(ie)%state%v(:,:,1,k,nt) + &
             mu*dt*vtens(:,:,1,k,ie)*elem(ie)%rspheremp(:,:)

        elem(ie)%state%v(:,:,2,k,nt)=elem(ie)%state%v(:,:,2,k,nt) + &
             mu*dt*vtens(:,:,2,k,ie)*elem(ie)%rspheremp(:,:)
        
        elem(ie)%state%dp3d(:,:,k,nt)=elem(ie)%state%dp3d(:,:,k,nt) &
             +mu_s*dt*stens(:,:,k,1,ie)*elem(ie)%rspheremp(:,:)
        
        elem(ie)%state%vtheta_dp(:,:,k,nt)=elem(ie)%state%vtheta_dp(:,:,k,nt) &
             +mu_s*dt*stens(:,:,k,2,ie)*elem(ie)%rspheremp(:,:)

        elem(ie)%state%w_i(:,:,k,nt)=elem(ie)%state%w_i(:,:,k,nt) &
             +mu_s*dt*stens(:,:,k,3,ie)*elem(ie)%rspheremp(:,:)
        
        elem(ie)%state%phinh_i(:,:,k,nt)=elem(ie)%state%phinh_i(:,:,k,nt) &
             +mu_s*dt*stens(:,:,k,4,ie)*elem(ie)%rspheremp(:,:)
        
     enddo
  enddo


  ! convert vtheta_dp -> theta
  do ie=nets,nete            
     do k=1,nlev
        elem(ie)%state%vtheta_dp(:,:,k,nt)=&
             elem(ie)%state%vtheta_dp(:,:,k,nt)*elem(ie)%state%dp3d(:,:,k,nt)
     enddo

     ! finally update w at the surface: 
     elem(ie)%state%w_i(:,:,nlevp,nt) = (elem(ie)%state%v(:,:,1,nlev,nt)*elem(ie)%derived%gradphis(:,:,1) + &
          elem(ie)%state%v(:,:,2,nlev,nt)*elem(ie)%derived%gradphis(:,:,2))/g
  enddo


  call t_stopf('advance_physical_vis')

  end subroutine advance_physical_vis





!============================ stiff and or non-stiff ============================================

 subroutine compute_andor_apply_rhs(np1,nm1,n0,qn0,dt2,elem,hvcoord,hybrid,&
       deriv,nets,nete,compute_diagnostics,eta_ave_w,scale1,scale2,scale3)
  ! ===================================
  ! compute the RHS, accumulate into u(np1) and apply DSS
  !
  !   u(np1) = scale3*u(nm1) + dt2*DSS[ nonstiffRHS(u(n0))*scale1 + stiffRHS(un0)*scale2 ]
  !
  ! This subroutine was orgininally called to compute a leapfrog timestep
  ! but by adjusting np1,nm1,n0 and dt2, many other timesteps can be
  ! accomodated.  For example, setting nm1=np1=n0 this routine will
  ! take a forward euler step, overwriting the input with the output.
  !
  !    qn0 = timelevel used to access Qdp() in order to compute virtual Temperature
  !
  ! ===================================

  integer,              intent(in) :: np1,nm1,n0,qn0,nets,nete
  real*8,               intent(in) :: dt2
  logical,              intent(in) :: compute_diagnostics
  type (hvcoord_t),     intent(in) :: hvcoord
  type (hybrid_t),      intent(in) :: hybrid
  type (element_t),     intent(inout), target :: elem(:)
  type (derivative_t),  intent(in) :: deriv

  real (kind=real_kind) :: eta_ave_w,scale1,scale2,scale3  ! weighting for eta_dot_dpdn mean flux, scale of unm1

  ! local
  real (kind=real_kind), pointer, dimension(:,:,:) :: phi_i
  real (kind=real_kind), pointer, dimension(:,:,:) :: dp3d
  real (kind=real_kind), pointer, dimension(:,:,:) :: vtheta_dp
   
  real (kind=real_kind) :: vtheta(np,np,nlev)
  real (kind=real_kind) :: vtheta_i(np,np,nlevp)
  real (kind=real_kind) :: omega_i(np,np,nlevp)
  real (kind=real_kind) :: omega(np,np,nlev)
  real (kind=real_kind) :: vort(np,np,nlev)           ! vorticity
  real (kind=real_kind) :: divdp(np,np,nlev)     
  real (kind=real_kind) :: phi(np,np,nlev)
  real (kind=real_kind) :: pnh(np,np,nlev)               ! nh (nonydro) pressure
  real (kind=real_kind) :: dp3d_i(np,np,nlevp)
  real (kind=real_kind) :: exner(np,np,nlev)         ! exner nh pressure
  real (kind=real_kind) :: dpnh_dp_i(np,np,nlevp)    ! dpnh / dp3d at interfaces
  real (kind=real_kind) :: eta_dot_dpdn(np,np,nlevp)  ! vertical velocity at interfaces
  real (kind=real_kind) :: KE(np,np,nlev)             ! Kinetic energy
  real (kind=real_kind) :: gradexner(np,np,2,nlev)    ! grad(p^kappa)
  real (kind=real_kind) :: gradphinh_i(np,np,2,nlevp) ! gradphi at interfaces
  real (kind=real_kind) :: mgrad(np,np,2,nlev)        ! gradphi metric term at cell centers
  real (kind=real_kind) :: gradKE(np,np,2,nlev)       ! grad(0.5 u^T u )
  real (kind=real_kind) :: wvor(np,np,2,nlev)         ! w vorticity term

  real (kind=real_kind) :: gradw_i(np,np,2,nlevp)
  real (kind=real_kind) :: v_gradw_i(np,np,nlevp)     
  real (kind=real_kind) :: v_gradtheta(np,np,nlev)     
  real (kind=real_kind) :: v_theta(np,np,2,nlev)
  real (kind=real_kind) :: div_v_theta(np,np,nlev)
  real (kind=real_kind) :: v_gradphinh_i(np,np,nlevp) ! v*gradphi at interfaces
  real (kind=real_kind) :: v_i(np,np,2,nlevp)

  real (kind=real_kind) :: v_vadv(np,np,2,nlev)     ! velocity vertical advection
  real (kind=real_kind) :: theta_vadv(np,np,nlev)   ! w,phi, theta  vertical advection term
  real (kind=real_kind) :: w_vadv_i(np,np,nlevp)      ! w,phi, theta  vertical advection term
  real (kind=real_kind) :: phi_vadv_i(np,np,nlevp)    ! w,phi, theta  vertical advection term

  real (kind=real_kind) :: vtens1(np,np,nlev)
  real (kind=real_kind) :: vtens2(np,np,nlev)
  real (kind=real_kind) :: stens(np,np,nlev,3) ! tendencies w,phi,theta
                                               ! w,phi tendencies not computed at nlevp
  real (kind=real_kind) :: w_tens(np,np,nlevp)  ! need to update w at surface as well
  real (kind=real_kind) :: theta_tens(np,np,nlev)
  real (kind=real_kind) :: phi_tens(np,np,nlevp)
                                               

  real (kind=real_kind) :: pi(np,np,nlev)                ! hydrostatic pressure
  real (kind=real_kind) :: pi_i(np,np,nlevp)             ! hydrostatic pressure interfaces
  real (kind=real_kind), dimension(np,np,nlev) :: vgrad_p

  real (kind=real_kind) ::  temp(np,np,nlev)
  real (kind=real_kind) ::  vtemp(np,np,2,nlev)       ! generic gradient storage
  real (kind=real_kind), dimension(np,np) :: sdot_sum ! temporary field
  real (kind=real_kind) ::  v1,v2,w,d_eta_dot_dpdn_dn
  integer :: i,j,k,kptr,ie, nlyr_tot

  call t_startf('compute_andor_apply_rhs')

  if (theta_hydrostatic_mode) then
     nlyr_tot=4*nlev        ! dont bother to dss w_i and phinh_i
  else
     nlyr_tot=5*nlev+nlevp  ! total amount of data for DSS
  endif
  do ie=nets,nete
     dp3d  => elem(ie)%state%dp3d(:,:,:,n0)
     vtheta_dp  => elem(ie)%state%vtheta_dp(:,:,:,n0)
     vtheta(:,:,:) = vtheta_dp(:,:,:)/dp3d(:,:,:)
     phi_i => elem(ie)%state%phinh_i(:,:,:,n0)

#ifdef ENERGY_DIAGNOSTICS
     if (.not. theta_hydrostatic_mode) then
        ! check w b.c.
        temp(:,:,1) =  (elem(ie)%state%v(:,:,1,nlev,n0)*elem(ie)%derived%gradphis(:,:,1) + &
             elem(ie)%state%v(:,:,2,nlev,n0)*elem(ie)%derived%gradphis(:,:,2))/g
        do j=1,np
        do i=1,np
           if ( abs(temp(i,j,1)-elem(ie)%state%w_i(i,j,nlevp,n0)) >1e-10) then
              write(iulog,*) 'WARNING: w(n0) does not satisfy b.c.',ie,i,j,k
              write(iulog,*) 'val1 = ',temp(i,j,1)
              write(iulog,*) 'val2 = ',elem(ie)%state%w_i(i,j,nlevp,n0)
              write(iulog,*) 'diff: ',temp(i,j,1)-elem(ie)%state%w_i(i,j,nlevp,n0)
           endif
        enddo
        enddo
        ! w boundary condition. just in case:
        !elem(ie)%state%w_i(:,:,nlevp,n0) = (elem(ie)%state%v(:,:,1,nlev,n0)*elem(ie)%derived%gradphis(:,:,1) + &
        !     elem(ie)%state%v(:,:,2,nlev,n0)*elem(ie)%derived%gradphis(:,:,2))/g

        ! check for layer spacing <= 1m
        do k=1,nlev
        do j=1,np
        do i=1,np
           if (scale3.ne.0d0) then
             if ((phi_i(i,j,k)-phi_i(i,j,k+1)) < g) then
                write(iulog,*) 'WARNING: before ADV, delta z < 1m. ie,i,j,k=',ie,i,j,k
                write(iulog,*) 'phi(i,j,k)=  ',phi_i(i,j,k)
                write(iulog,*) 'phi(i,j,k+1)=',phi_i(i,j,k+1)
             endif
           endif
        enddo
        enddo
        enddo
        
     endif
#endif

     call pnh_and_exner_from_eos(hvcoord,vtheta_dp,dp3d,phi_i,pnh,exner,dpnh_dp_i,caller='CAAR')

     dp3d_i(:,:,1) = dp3d(:,:,1)
     dp3d_i(:,:,nlevp) = dp3d(:,:,nlev)
     do k=2,nlev
        dp3d_i(:,:,k)=(dp3d(:,:,k)+dp3d(:,:,k-1))/2
     end do

     ! special averaging for velocity for energy conservation
     v_i(:,:,1:2,1) = elem(ie)%state%v(:,:,1:2,1,n0)  
     v_i(:,:,1:2,nlevp) = elem(ie)%state%v(:,:,1:2,nlev,n0)
     do k=2,nlev
        v_i(:,:,1,k) = (dp3d(:,:,k)*elem(ie)%state%v(:,:,1,k,n0) + &
             dp3d(:,:,k-1)*elem(ie)%state%v(:,:,1,k-1,n0) ) / (2*dp3d_i(:,:,k))
        v_i(:,:,2,k) = (dp3d(:,:,k)*elem(ie)%state%v(:,:,2,k,n0) + &
             dp3d(:,:,k-1)*elem(ie)%state%v(:,:,2,k-1,n0) ) / (2*dp3d_i(:,:,k))
     end do
     
     if (theta_hydrostatic_mode) then
        do k=nlev,1,-1          ! traditional Hydrostatic integral
           phi_i(:,:,k)=phi_i(:,:,k+1)+&
                Rgas*vtheta_dp(:,:,k)*exner(:,:,k)/pnh(:,:,k)
        enddo
        ! in H mode, ignore w contibutions to KE term
        ! set to zero so H and NH can share code and reduce if statements
        elem(ie)%state%w_i(:,:,:,n0)=0   
     endif

     do k=1,nlev
        phi(:,:,k) = (phi_i(:,:,k)+phi_i(:,:,k+1))/2  ! for diagnostics

        ! ================================
        ! Accumulate mean Vel_rho flux in vn0
        ! ================================
        vtemp(:,:,1,k) = elem(ie)%state%v(:,:,1,k,n0)*dp3d(:,:,k)
        vtemp(:,:,2,k) = elem(ie)%state%v(:,:,2,k,n0)*dp3d(:,:,k)
        elem(ie)%derived%vn0(:,:,:,k)=elem(ie)%derived%vn0(:,:,:,k)+eta_ave_w*vtemp(:,:,:,k)

        divdp(:,:,k)=divergence_sphere(vtemp(:,:,:,k),deriv,elem(ie))
        vort(:,:,k)=vorticity_sphere(elem(ie)%state%v(:,:,:,k,n0),deriv,elem(ie))
     enddo

     ! Compute omega =  Dpi/Dt   Used only as a DIAGNOSTIC
     pi_i(:,:,1)=hvcoord%hyai(1)*hvcoord%ps0
     omega_i(:,:,1)=0
     do k=1,nlev
        pi_i(:,:,k+1)=pi_i(:,:,k) + dp3d(:,:,k)
        omega_i(:,:,k+1)=omega_i(:,:,k)+divdp(:,:,k)
     enddo
     do k=1,nlev
        pi(:,:,k)=pi_i(:,:,k) + dp3d(:,:,k)/2
        vtemp(:,:,:,k) = gradient_sphere( pi(:,:,k), deriv, elem(ie)%Dinv);
        vgrad_p(:,:,k) = elem(ie)%state%v(:,:,1,k,n0)*vtemp(:,:,1,k)+&
             elem(ie)%state%v(:,:,2,k,n0)*vtemp(:,:,2,k)
        omega(:,:,k) = (vgrad_p(:,:,k) - ( omega_i(:,:,k)+omega_i(:,:,k+1))/2) 
     enddo        

     ! ==================================================
     ! Compute eta_dot_dpdn
     ! save sdot_sum as this is the -RHS of ps_v equation
     ! ==================================================
     if (rsplit>0) then
        ! VERTICALLY LAGRANGIAN:   no vertical motion
        eta_dot_dpdn=0
        w_vadv_i=0
        phi_vadv_i=0
        theta_vadv=0
        v_vadv=0
     else
        sdot_sum=0
        do k=1,nlev
           ! ==================================================
           ! add this term to PS equation so we exactly conserve dry mass
           ! ==================================================
           sdot_sum(:,:) = sdot_sum(:,:) + divdp(:,:,k)
           eta_dot_dpdn(:,:,k+1) = sdot_sum(:,:)
        end do


        ! ===========================================================
        ! at this point, eta_dot_dpdn contains integral_etatop^eta[ divdp ]
        ! compute at interfaces:
        !    eta_dot_dpdn = -dp/dt - integral_etatop^eta[ divdp ]
        ! for reference: at mid layers we have:
        !    omega = v grad p  - integral_etatop^eta[ divdp ]
        ! ===========================================================
        do k=1,nlev-1
           eta_dot_dpdn(:,:,k+1) = hvcoord%hybi(k+1)*sdot_sum(:,:) - eta_dot_dpdn(:,:,k+1)
        end do

        eta_dot_dpdn(:,:,1     ) = 0
        eta_dot_dpdn(:,:,nlevp)  = 0
        vtheta_i(:,:,1) =       0
        vtheta_i(:,:,nlevp) =   0


        ! ===========================================================
        ! Compute vertical advection of v from eq. CCM2 (3.b.1)
        ! ==============================================
        call preq_vertadv_v1(elem(ie)%state%v(:,:,:,:,n0),eta_dot_dpdn,dp3d,v_vadv)

        ! compute (cp*theta) at interfaces
        ! for energy conservation, use averaging consistent with EOS
        ! dont bother to compute at surface and top since it will be multiplied by eta-dot
#if 0           
           do k=2,nlev  ! simple averaging
              vtheta_i(:,:,k) = (vtheta(:,:,k)+vtheta(:,:,k-1))/2
           enddo
#else
           ! E conserving average, but much more dissipative
           do k=2,nlev
              vtheta_i(:,:,k) = -dpnh_dp_i(:,:,k)*(phi(:,:,k)-phi(:,:,k-1))/&
                   (exner(:,:,k)-exner(:,:,k-1)) / Cp
           enddo
#endif           



        do k=1,nlev
           ! average interface quantity to midpoints:
           temp(:,:,k) = (( eta_dot_dpdn(:,:,k)+eta_dot_dpdn(:,:,k+1))/2)*&
                (elem(ie)%state%w_i(:,:,k+1,n0)-elem(ie)%state%w_i(:,:,k,n0))
           
           ! theta vadv term at midoints
           theta_vadv(:,:,k)= eta_dot_dpdn(:,:,k+1)*vtheta_i(:,:,k+1) - &
                eta_dot_dpdn(:,:,k)*vtheta_i(:,:,k)
        enddo
        ! compute ave( ave(etadot) d/dx )
        do k=2,nlev
           w_vadv_i(:,:,k)  =(temp(:,:,k-1)+temp(:,:,k))/2
           phi_vadv_i(:,:,k)=eta_dot_dpdn(:,:,k)*(phi(:,:,k)-phi(:,:,k-1))
        end do
        w_vadv_i(:,:,1) = temp(:,:,1)
        w_vadv_i(:,:,nlevp) = temp(:,:,nlev)
        phi_vadv_i(:,:,1) = 0
        phi_vadv_i(:,:,nlevp) = 0

        ! final form of SB81 vertical advection operator:
        w_vadv_i=w_vadv_i/dp3d_i
        phi_vadv_i=phi_vadv_i/dp3d_i
     endif


     ! ================================
     ! accumulate mean vertical flux:
     ! ================================
     do k=1,nlev  !  Loop index added (AAM)
        elem(ie)%derived%eta_dot_dpdn(:,:,k) = &
             elem(ie)%derived%eta_dot_dpdn(:,:,k) + eta_ave_w*eta_dot_dpdn(:,:,k)
        elem(ie)%derived%omega_p(:,:,k) = &
             elem(ie)%derived%omega_p(:,:,k) + eta_ave_w*omega(:,:,k)
     enddo
     elem(ie)%derived%eta_dot_dpdn(:,:,nlev+1) = &
             elem(ie)%derived%eta_dot_dpdn(:,:,nlev+1) + eta_ave_w*eta_dot_dpdn(:,:,nlev+1)

     ! ================================================
     ! w,phi tendencies including surface
     ! ================================================  
     do k=1,nlev
        ! compute gradphi at interfaces and then average to levels
        gradphinh_i(:,:,:,k)   = gradient_sphere(phi_i(:,:,k),deriv,elem(ie)%Dinv)   
           
        gradw_i(:,:,:,k)   = gradient_sphere(elem(ie)%state%w_i(:,:,k,n0),deriv,elem(ie)%Dinv)
        v_gradw_i(:,:,k) = v_i(:,:,1,k)*gradw_i(:,:,1,k) + v_i(:,:,2,k)*gradw_i(:,:,2,k)
        ! w - tendency on interfaces 
        w_tens(:,:,k) = (-w_vadv_i(:,:,k) - v_gradw_i(:,:,k))*scale1 - scale2*g*(1-dpnh_dp_i(:,:,k) )

        ! phi - tendency on interfaces
        v_gradphinh_i(:,:,k) = v_i(:,:,1,k)*gradphinh_i(:,:,1,k) &
             +v_i(:,:,2,k)*gradphinh_i(:,:,2,k) 
        phi_tens(:,:,k) =  (-phi_vadv_i(:,:,k) - v_gradphinh_i(:,:,k))*scale1 &
          + scale2*g*elem(ie)%state%w_i(:,:,k,n0)
     end do


     ! k =nlevp case, all terms in the imex methods are treated explicitly at the boundary
     k =nlevp 
    ! compute gradphi at interfaces and then average to levels
    gradphinh_i(:,:,:,k)   = gradient_sphere(phi_i(:,:,k),deriv,elem(ie)%Dinv)
    gradw_i(:,:,:,k)   = gradient_sphere(elem(ie)%state%w_i(:,:,k,n0),deriv,elem(ie)%Dinv)
    v_gradw_i(:,:,k) = v_i(:,:,1,k)*gradw_i(:,:,1,k) + v_i(:,:,2,k)*gradw_i(:,:,2,k)
    ! w - tendency on interfaces
    w_tens(:,:,k) = (-w_vadv_i(:,:,k) - v_gradw_i(:,:,k))*scale1 - scale1*g*(1-dpnh_dp_i(:,:,k) )

    ! phi - tendency on interfaces
    v_gradphinh_i(:,:,k) = v_i(:,:,1,k)*gradphinh_i(:,:,1,k) &
     +v_i(:,:,2,k)*gradphinh_i(:,:,2,k)
    phi_tens(:,:,k) =  (-phi_vadv_i(:,:,k) - v_gradphinh_i(:,:,k))*scale1 &
    + scale1*g*elem(ie)%state%w_i(:,:,k,n0)
    




     ! ================================================                                                                 
     ! v1,v2 tendencies:                                                                                          
     ! ================================================           
     do k=1,nlev
        ! theta - tendency on levels
        if (theta_advect_form==0) then
           v_theta(:,:,1,k)=elem(ie)%state%v(:,:,1,k,n0)*vtheta_dp(:,:,k)
           v_theta(:,:,2,k)=elem(ie)%state%v(:,:,2,k,n0)*vtheta_dp(:,:,k)
           div_v_theta(:,:,k)=divergence_sphere(v_theta(:,:,:,k),deriv,elem(ie))
        else
           ! alternate form, non-conservative, better HS topography results
           v_theta(:,:,:,k) = gradient_sphere(vtheta(:,:,k),deriv,elem(ie)%Dinv)
           div_v_theta(:,:,k)=vtheta(:,:,k)*divdp(:,:,k) + &
                dp3d(:,:,k)*elem(ie)%state%v(:,:,1,k,n0)*v_theta(:,:,1,k) + &
                dp3d(:,:,k)*elem(ie)%state%v(:,:,2,k,n0)*v_theta(:,:,2,k) 
        endif
        theta_tens(:,:,k)=(-theta_vadv(:,:,k)-div_v_theta(:,:,k))*scale1

        ! w vorticity correction term
        temp(:,:,k) = (elem(ie)%state%w_i(:,:,k,n0)**2 + &
             elem(ie)%state%w_i(:,:,k+1,n0)**2)/4
        wvor(:,:,:,k) = gradient_sphere(temp(:,:,k),deriv,elem(ie)%Dinv)
        wvor(:,:,1,k) = wvor(:,:,1,k) - (elem(ie)%state%w_i(:,:,k,n0)*gradw_i(:,:,1,k) +&
             elem(ie)%state%w_i(:,:,k+1,n0)*gradw_i(:,:,1,k+1))/2
        wvor(:,:,2,k) = wvor(:,:,2,k) - (elem(ie)%state%w_i(:,:,k,n0)*gradw_i(:,:,2,k) +&
             elem(ie)%state%w_i(:,:,k+1,n0)*gradw_i(:,:,2,k+1))/2

        KE(:,:,k) = ( elem(ie)%state%v(:,:,1,k,n0)**2 + elem(ie)%state%v(:,:,2,k,n0)**2)/2
        gradKE(:,:,:,k) = gradient_sphere(KE(:,:,k),deriv,elem(ie)%Dinv)
        gradexner(:,:,:,k) = gradient_sphere(exner(:,:,k),deriv,elem(ie)%Dinv)
#if 0
        ! another form: (good results in dcmip2012 test2.0)  max=0.195
        ! but bad results with HS topo
        !  grad(exner) =( grad(theta*exner) - exner*grad(theta))/theta
        vtemp(:,:,:,k) = gradient_sphere(vtheta(:,:,k)*exner(:,:,k),deriv,elem(ie)%Dinv)
        v_theta(:,:,:,k) = gradient_sphere(vtheta(:,:,k),deriv,elem(ie)%Dinv)
        gradexner(:,:,1,k) = (vtemp(:,:,1,k)-exner(:,:,k)*v_theta(:,:,1,k))/&
             vtheta(:,:,k)
        gradexner(:,:,2,k) = (vtemp(:,:,2,k)-exner(:,:,k)*v_theta(:,:,2,k))/&
             vtheta(:,:,k)
#endif
#if 0
        ! entropy form: dcmip2012 test2.0 best: max=0.130  (0.124 with conservation form theta)
        vtemp(:,:,:,k) = gradient_sphere(vtheta(:,:,k)*exner(:,:,k),deriv,elem(ie)%Dinv)
        v_theta(:,:,:,k) = gradient_sphere(log(vtheta(:,:,k)),deriv,elem(ie)%Dinv)
        gradexner(:,:,1,k) = (vtemp(:,:,1,k)-exner(:,:,k)*vtheta(:,:,k)*v_theta(:,:,1,k))/&
             vtheta(:,:,k)
        gradexner(:,:,2,k) = (vtemp(:,:,2,k)-exner(:,:,k)*vtheta(:,:,k)*v_theta(:,:,2,k))/&
             vtheta(:,:,k)
#endif
#if 0
        ! another form:  terrible results in dcmip2012 test2.0
        ! grad(exner) = grad(p) * kappa * exner / p
        gradexner(:,:,:,k) = gradient_sphere(pnh(:,:,k),deriv,elem(ie)%Dinv)
        gradexner(:,:,1,k) = gradexner(:,:,1,k)*(Rgas/Cp)*exner(:,:,k)/pnh(:,:,k)
        gradexner(:,:,2,k) = gradexner(:,:,2,k)*(Rgas/Cp)*exner(:,:,k)/pnh(:,:,k)
#endif

        ! special averaging of dpnh/dpi grad(phi) for E conservation
        mgrad(:,:,1,k) = (dpnh_dp_i(:,:,k)*gradphinh_i(:,:,1,k)+ &
              dpnh_dp_i(:,:,k+1)*gradphinh_i(:,:,1,k+1))/2
        mgrad(:,:,2,k) = (dpnh_dp_i(:,:,k)*gradphinh_i(:,:,2,k)+ &
              dpnh_dp_i(:,:,k+1)*gradphinh_i(:,:,2,k+1))/2


        do j=1,np
           do i=1,np
              v1     = elem(ie)%state%v(i,j,1,k,n0)
              v2     = elem(ie)%state%v(i,j,2,k,n0)

              vtens1(i,j,k) = (-v_vadv(i,j,1,k) &
                   + v2*(elem(ie)%fcor(i,j) + vort(i,j,k))        &
                   - gradKE(i,j,1,k) - mgrad(i,j,1,k) &
                  -Cp*vtheta(i,j,k)*gradexner(i,j,1,k)&
                  -wvor(i,j,1,k) )*scale1


              vtens2(i,j,k) = (-v_vadv(i,j,2,k) &
                   - v1*(elem(ie)%fcor(i,j) + vort(i,j,k)) &
                   - gradKE(i,j,2,k) - mgrad(i,j,2,k) &
                  -Cp*vtheta(i,j,k)*gradexner(i,j,2,k) &
                  -wvor(i,j,2,k) )*scale1
           end do
        end do     
     end do 



     
#ifdef ENERGY_DIAGNOSTICS
     ! =========================================================
     ! diagnostics. not performance critical, dont thread
     ! =========================================================
     if (compute_diagnostics) then
        elem(ie)%accum%KEu_horiz1=0
        elem(ie)%accum%KEu_horiz2=0
        elem(ie)%accum%KEu_vert1=0
        elem(ie)%accum%KEu_vert2=0
        elem(ie)%accum%KEw_horiz1=0
        elem(ie)%accum%KEw_horiz2=0
        elem(ie)%accum%KEw_horiz3=0
        elem(ie)%accum%KEw_vert1=0
        elem(ie)%accum%KEw_vert2=0

        elem(ie)%accum%PEhoriz1=0
        elem(ie)%accum%PEhoriz2=0
        elem(ie)%accum%IEvert1=0
        elem(ie)%accum%IEvert2=0
        elem(ie)%accum%PEvert1=0
        elem(ie)%accum%PEvert2=0
        elem(ie)%accum%T01=0
        elem(ie)%accum%T2=0
        elem(ie)%accum%S1=0
        elem(ie)%accum%S2=0
        elem(ie)%accum%P1=0
        elem(ie)%accum%P2=0

        do k =1,nlev
          do j=1,np
            do i=1,np                
               d_eta_dot_dpdn_dn=(eta_dot_dpdn(i,j,k+1)-eta_dot_dpdn(i,j,k))
               !  Form horiz advection of KE-u
               elem(ie)%accum%KEu_horiz1(i,j)=elem(ie)%accum%KEu_horiz1(i,j) &
                    -dp3d(i,j,k)*( &
                    elem(ie)%state%v(i,j,1,k,n0)*gradKE(i,j,1,k) + &
                    elem(ie)%state%v(i,j,2,k,n0)*gradKE(i,j,2,k) )
               elem(ie)%accum%KEu_horiz2(i,j)=elem(ie)%accum%KEu_horiz2(i,j)              &
                    -KE(i,j,k)*divdp(i,j,k)
               !  Form horiz advection of KE-w
               elem(ie)%accum%KEw_horiz1(i,j)=elem(ie)%accum%KEw_horiz1(i,j)-   &
                    dp3d(i,j,k) * (&
                    elem(ie)%state%w_i(i,j,k,n0) * v_gradw_i(i,j,k)    + &
                    elem(ie)%state%w_i(i,j,k+1,n0) * v_gradw_i(i,j,k+1) )/2
               elem(ie)%accum%KEw_horiz2(i,j)=elem(ie)%accum%KEw_horiz2(i,j)-   &
                    divdp(i,j,k)*(elem(ie)%state%w_i(i,j,k,n0)**2 + &
                    elem(ie)%state%w_i(i,j,k+1,n0)**2 ) /4
               elem(ie)%accum%KEw_horiz3(i,j)=elem(ie)%accum%KEw_horiz3(i,j)   &
                    -dp3d(i,j,k) * (elem(ie)%state%v(i,j,1,k,n0) * wvor(i,j,1,k) +  &
                    elem(ie)%state%v(i,j,2,k,n0) * wvor(i,j,2,k))
               !  Form vertical advection of KE-u 
               elem(ie)%accum%KEu_vert1(i,j)=elem(ie)%accum%KEu_vert1(i,j)- &
                    (elem(ie)%state%v(i,j,1,k,n0) * v_vadv(i,j,1,k) +            &
                    elem(ie)%state%v(i,j,2,k,n0) *v_vadv(i,j,2,k))*dp3d(i,j,k)
               elem(ie)%accum%KEu_vert2(i,j)=elem(ie)%accum%KEu_vert2(i,j)- &
                    0.5*((elem(ie)%state%v(i,j,1,k,n0))**2 +                     &
                    (elem(ie)%state%v(i,j,2,k,n0))**2)*d_eta_dot_dpdn_dn
               !  Form vertical advection of KE-w
               elem(ie)%accum%KEw_vert1(i,j)=elem(ie)%accum%KEw_vert1(i,j) - &
                    dp3d(i,j,k) * &
                    (w_vadv_i(i,j,k)*elem(ie)%state%w_i(i,j,k,n0)+ &
                    w_vadv_i(i,j,k+1)*elem(ie)%state%w_i(i,j,k+1,n0))/2
               
               elem(ie)%accum%KEw_vert2(i,j)=elem(ie)%accum%KEw_vert2(i,j)      &
                    -d_eta_dot_dpdn_dn* &
                    (.5*elem(ie)%state%w_i(i,j,k,n0)**2 +&
                    .5*elem(ie)%state%w_i(i,j,k+1,n0)**2)/2
               
               !  Form IEvert1
               elem(ie)%accum%IEvert1(i,j)=elem(ie)%accum%IEvert1(i,j)      &
                    -exner(i,j,k)*theta_vadv(i,j,k)                        
               ! Form IEvert2 
               ! here use of dpnh_dp_i on boundry (with incorrect data)
               ! is harmess becuase eta_dot_dpdn=0
               elem(ie)%accum%IEvert2(i,j)=elem(ie)%accum%IEvert2(i,j)      &
                    + ( dpnh_dp_i(i,j,k)*eta_dot_dpdn(i,j,k)+ &
                        dpnh_dp_i(i,j,k+1)*eta_dot_dpdn(i,j,k+1)) &
                    *(phi_i(i,j,k+1)-phi_i(i,j,k))/2
               
               !  Form PEhoriz1
               elem(ie)%accum%PEhoriz1(i,j)=(elem(ie)%accum%PEhoriz1(i,j))  &
                    -phi(i,j,k)*divdp(i,j,k) 
               !  Form PEhoriz2
               elem(ie)%accum%PEhoriz2(i,j)=elem(ie)%accum%PEhoriz2(i,j)    &
                    -dp3d(i,j,k)* &
                    (elem(ie)%state%v(i,j,1,k,n0)*                          &
                    (gradphinh_i(i,j,1,k)+gradphinh_i(i,j,1,k+1))/2  +      &
                    elem(ie)%state%v(i,j,2,k,n0)*                           &
                    (gradphinh_i(i,j,2,k)+gradphinh_i(i,j,2,k+1))/2  )
               
               !  Form PEvert1
               elem(ie)%accum%PEvert1(i,j) = elem(ie)%accum%PEvert1(i,j)    &
                    -phi(i,j,k)*d_eta_dot_dpdn_dn                                 
               elem(ie)%accum%PEvert2(i,j) = elem(ie)%accum%PEvert2(i,j)     &
                    -dp3d(i,j,k)*(phi_vadv_i(i,j,k)+phi_vadv_i(i,j,k+1))/2
               
               !  Form T01
               elem(ie)%accum%T01(i,j)=elem(ie)%accum%T01(i,j)               &
                    -(Cp*elem(ie)%state%vtheta_dp(i,j,k,n0))                       &
                    *(gradexner(i,j,1,k)*elem(ie)%state%v(i,j,1,k,n0) +           &
                    gradexner(i,j,2,k)*elem(ie)%state%v(i,j,2,k,n0))              
               !  Form S1 
               elem(ie)%accum%S1(i,j)=elem(ie)%accum%S1(i,j)                 &
                    -Cp*exner(i,j,k)*div_v_theta(i,j,k)

               !  Form P1  = -P2  (no reason to compute P2?)
               elem(ie)%accum%P1(i,j)=elem(ie)%accum%P1(i,j) -g*dp3d(i,j,k)* &
                    ( elem(ie)%state%w_i(i,j,k,n0) + &
                    elem(ie)%state%w_i(i,j,k+1,n0) )/2
               !  Form P2
               elem(ie)%accum%P2(i,j)=elem(ie)%accum%P2(i,j) + g*dp3d(i,j,k)*&
                    ( elem(ie)%state%w_i(i,j,k,n0) + &
                    elem(ie)%state%w_i(i,j,k+1,n0) )/2
            enddo
         enddo
      enddo

      ! these terms are better easier to compute by summing interfaces
      do k=2,nlev
         elem(ie)%accum%T2(:,:)=elem(ie)%accum%T2(:,:)+                &
              (g*elem(ie)%state%w_i(:,:,k,n0)-v_gradphinh_i(:,:,k)) &
               * dpnh_dp_i(:,:,k)*dp3d_i(:,:,k)
      enddo
      ! boundary terms
      do k=1,nlevp,nlev
         elem(ie)%accum%T2(:,:)=elem(ie)%accum%T2(:,:)+                &
           (g*elem(ie)%state%w_i(:,:,k,n0)-v_gradphinh_i(:,:,k)) &
           * dpnh_dp_i(:,:,k)*dp3d_i(:,:,k)/2
      enddo
      ! boundary term is incorrect.  save the term so we can correct it
      ! once we have coorect value of dpnh_dp_i:
      elem(ie)%accum%T2_nlevp_term(:,:)=&
           (g*elem(ie)%state%w_i(:,:,nlevp,n0)-v_gradphinh_i(:,:,nlevp)) &
           * dp3d_i(:,:,nlevp)/2

   endif
#endif


     do k=1,nlev
        elem(ie)%state%v(:,:,1,k,np1) = elem(ie)%spheremp(:,:)*(scale3 * elem(ie)%state%v(:,:,1,k,nm1) &
          + dt2*vtens1(:,:,k) )
        elem(ie)%state%v(:,:,2,k,np1) = elem(ie)%spheremp(:,:)*(scale3 * elem(ie)%state%v(:,:,2,k,nm1) &
          +  dt2*vtens2(:,:,k) )
        elem(ie)%state%vtheta_dp(:,:,k,np1) = elem(ie)%spheremp(:,:)*(scale3 * elem(ie)%state%vtheta_dp(:,:,k,nm1) &
          + dt2*theta_tens(:,:,k))

        if ( .not. theta_hydrostatic_mode ) then
           elem(ie)%state%w_i(:,:,k,np1)    = elem(ie)%spheremp(:,:)*(scale3 * elem(ie)%state%w_i(:,:,k,nm1)   &
                + dt2*w_tens(:,:,k))
           elem(ie)%state%phinh_i(:,:,k,np1)   = elem(ie)%spheremp(:,:)*(scale3 * elem(ie)%state%phinh_i(:,:,k,nm1) & 
                + dt2*phi_tens(:,:,k))
        endif

        elem(ie)%state%dp3d(:,:,k,np1) = &
             elem(ie)%spheremp(:,:) * (scale3 * elem(ie)%state%dp3d(:,:,k,nm1) - &
             scale1*dt2 * (divdp(:,:,k) + eta_dot_dpdn(:,:,k+1)-eta_dot_dpdn(:,:,k)))
     enddo
     k=nlevp
     if ( .not. theta_hydrostatic_mode ) &
          elem(ie)%state%w_i(:,:,k,np1)=elem(ie)%spheremp(:,:)*(scale3 * elem(ie)%state%w_i(:,:,k,nm1)   &
          + dt2*w_tens(:,:,k))


     kptr=0
     call edgeVpack_nlyr(edge_g,elem(ie)%desc,elem(ie)%state%v(:,:,:,:,np1),2*nlev,kptr,nlyr_tot)
     kptr=kptr+2*nlev
     call edgeVpack_nlyr(edge_g,elem(ie)%desc,elem(ie)%state%dp3d(:,:,:,np1),nlev,kptr, nlyr_tot)
     kptr=kptr+nlev
     call edgeVpack_nlyr(edge_g,elem(ie)%desc,elem(ie)%state%vtheta_dp(:,:,:,np1),nlev,kptr,nlyr_tot)
     if (.not. theta_hydrostatic_mode) then
        kptr=kptr+nlev
        call edgeVpack_nlyr(edge_g,elem(ie)%desc,elem(ie)%state%w_i(:,:,:,np1),nlevp,kptr,nlyr_tot)
        kptr=kptr+nlevp
        call edgeVpack_nlyr(edge_g,elem(ie)%desc,elem(ie)%state%phinh_i(:,:,:,np1),nlev,kptr,nlyr_tot)
     endif

   end do ! end do for the ie=nets,nete loop

  call t_startf('caar_bexchV')
  call bndry_exchangeV(hybrid,edge_g)
  call t_stopf('caar_bexchV')

  do ie=nets,nete
     kptr=0
     call edgeVunpack_nlyr(edge_g,elem(ie)%desc,elem(ie)%state%v(:,:,:,:,np1),2*nlev,kptr,nlyr_tot)
     kptr=kptr+2*nlev
     call edgeVunpack_nlyr(edge_g,elem(ie)%desc,elem(ie)%state%dp3d(:,:,:,np1),nlev,kptr,nlyr_tot)
     kptr=kptr+nlev
     call edgeVunpack_nlyr(edge_g,elem(ie)%desc,elem(ie)%state%vtheta_dp(:,:,:,np1),nlev,kptr,nlyr_tot)
     if (.not. theta_hydrostatic_mode) then
        kptr=kptr+nlev
        call edgeVunpack_nlyr(edge_g,elem(ie)%desc,elem(ie)%state%w_i(:,:,:,np1),nlevp,kptr,nlyr_tot)
        kptr=kptr+nlevp
        call edgeVunpack_nlyr(edge_g,elem(ie)%desc,elem(ie)%state%phinh_i(:,:,:,np1),nlev,kptr,nlyr_tot)
     endif
      
     ! ====================================================
     ! Scale tendencies by inverse mass matrix
     ! ====================================================
     do k=1,nlev
        elem(ie)%state%dp3d(:,:,k,np1) =elem(ie)%rspheremp(:,:)*elem(ie)%state%dp3d(:,:,k,np1)
        elem(ie)%state%vtheta_dp(:,:,k,np1)=elem(ie)%rspheremp(:,:)*elem(ie)%state%vtheta_dp(:,:,k,np1)
        if ( .not. theta_hydrostatic_mode ) then
           elem(ie)%state%w_i(:,:,k,np1)    =elem(ie)%rspheremp(:,:)*elem(ie)%state%w_i(:,:,k,np1)
           elem(ie)%state%phinh_i(:,:,k,np1)=elem(ie)%rspheremp(:,:)*elem(ie)%state%phinh_i(:,:,k,np1)
        endif
        elem(ie)%state%v(:,:,1,k,np1)  =elem(ie)%rspheremp(:,:)*elem(ie)%state%v(:,:,1,k,np1)
        elem(ie)%state%v(:,:,2,k,np1)  =elem(ie)%rspheremp(:,:)*elem(ie)%state%v(:,:,2,k,np1)
     end do
     k=nlevp
     if ( .not. theta_hydrostatic_mode ) &
          elem(ie)%state%w_i(:,:,k,np1)    =elem(ie)%rspheremp(:,:)*elem(ie)%state%w_i(:,:,k,np1)


     ! now we can compute the correct dphn_dp_i() at the surface:
     if (.not. theta_hydrostatic_mode) then
        ! solve for (dpnh_dp_i-1)
        dpnh_dp_i(:,:,nlevp) = 1 + (  &
             ((elem(ie)%state%v(:,:,1,nlev,np1)*elem(ie)%derived%gradphis(:,:,1) + &
             elem(ie)%state%v(:,:,2,nlev,np1)*elem(ie)%derived%gradphis(:,:,2))/g - &
             elem(ie)%state%w_i(:,:,nlevp,np1)) / &
             (g + ( elem(ie)%derived%gradphis(:,:,1)**2 + &
             elem(ie)%derived%gradphis(:,:,2)**2)/(2*g))   )  / dt2
        
        ! update solution with new dpnh_dp_i value:
        elem(ie)%state%w_i(:,:,nlevp,np1) = elem(ie)%state%w_i(:,:,nlevp,np1) +&
             scale1*dt2*g*(dpnh_dp_i(:,:,nlevp)-1)
        elem(ie)%state%v(:,:,1,nlev,np1) =  elem(ie)%state%v(:,:,1,nlev,np1) -&
             scale1*dt2*(dpnh_dp_i(:,:,nlevp)-1)*elem(ie)%derived%gradphis(:,:,1)/2
        elem(ie)%state%v(:,:,2,nlev,np1) =  elem(ie)%state%v(:,:,2,nlev,np1) -&
             scale1*dt2*(dpnh_dp_i(:,:,nlevp)-1)*elem(ie)%derived%gradphis(:,:,2)/2

#ifdef ENERGY_DIAGNOSTICS
        ! add in boundary term to T2 and S2 diagnostics:
        if (compute_diagnostics) then
           elem(ie)%accum%T2(:,:)=elem(ie)%accum%T2(:,:)+                &
                elem(ie)%accum%T2_nlevp_term(:,:)*(dpnh_dp_i(:,:,nlevp)-1)
           elem(ie)%accum%S2(:,:)=-elem(ie)%accum%T2(:,:)      
        endif

        ! check w b.c.
        temp(:,:,1) =  (elem(ie)%state%v(:,:,1,nlev,np1)*elem(ie)%derived%gradphis(:,:,1) + &
             elem(ie)%state%v(:,:,2,nlev,np1)*elem(ie)%derived%gradphis(:,:,2))/g
        do j=1,np
        do i=1,np
           if ( abs(temp(i,j,1)-elem(ie)%state%w_i(i,j,nlevp,np1)) >1e-10) then
              write(iulog,*) 'WARNING: w(np1) does not satisfy b.c.',ie,i,j,k
              write(iulog,*) 'val1 = ',temp(i,j,1)
              write(iulog,*) 'val2 = ',elem(ie)%state%w_i(i,j,nlevp,np1)
              write(iulog,*) 'diff: ',temp(i,j,1)-elem(ie)%state%w_i(i,j,nlevp,np1)
           endif
        enddo
        enddo

        ! check for layer spacing <= 1m
        do k=1,nlev
        do j=1,np
        do i=1,np
           if (scale3.ne.0d0) then
             if ((elem(ie)%state%phinh_i(i,j,k,np1)-elem(ie)%state%phinh_i(i,j,k+1,np1)) < g) then
                write(iulog,*) 'WARNING: after ADV, delta z < 1m. ie,i,j,k=',ie,i,j,k
                write(iulog,*) 'phi(i,j,k)=  ',elem(ie)%state%phinh_i(i,j,k,np1)
                write(iulog,*) 'phi(i,j,k+1)=',elem(ie)%state%phinh_i(i,j,k+1,np1)
             endif
           end if
        enddo
        enddo
        enddo
#endif
     endif

  end do
  call t_stopf('compute_andor_apply_rhs')

  end subroutine compute_andor_apply_rhs

  subroutine compute_nonlinear_rhs(np1,nm1,n0,qn0,elem,hvcoord,hybrid,&
       deriv,nets,nete,compute_diagnostics,eta_ave_w, JacL, JacD, JacU,dt)
  integer,              intent(in) :: np1,nm1,n0,qn0,nets,nete
  real (kind=real_kind), intent(in) :: JacL(nlev-1,np,np, nete-nets+1), JacD(nlev,np,np,nete-nets+1), JacU(nlev-1,np,np,nete-nets+1)
  logical,              intent(in) :: compute_diagnostics
  type (hvcoord_t),     intent(in) :: hvcoord
  type (hybrid_t),      intent(in) :: hybrid
  type (element_t),     intent(inout), target :: elem(:)
  type (derivative_t),  intent(in) :: deriv
  real (kind=real_kind), intent(in) :: dt

  real (kind=real_kind) :: eta_ave_w ! weighting for eta_dot_dpdn mean flux

  ! local variables
  real (kind=real_kind) :: wphivec(2*nlev)
  real (kind=real_kind) :: L(2*nlev, 2*nlev)
  real (kind=real_kind) :: wphivec2(2*nlev,np,np,nete-nets+1)
  integer :: ii, i,j,k,ie

  do ie = nets,nete
    do i=1,np
      do j=1,np
        wphivec2(1:nlev,i,j,ie) = elem(ie)%state%w_i(i,j,1:nlev,n0)
        wphivec2(1+nlev:2*nlev,i,j,ie) = elem(ie)%state%phinh_i(i,j,1:nlev,n0)
      end do
    end do
  end do
  call compute_andor_apply_rhs(np1,nm1,n0,qn0,1.d0,elem,hvcoord,hybrid,&
     deriv,nets,nete,compute_diagnostics,eta_ave_w,1.d0,1.d0,0.d0)
  do ie = nets,nete
    do i=1,np
      do j=1,np
  ! Form matrix L
        L = 0.d0
        do ii=1,nlev-1
          L(ii, ii+nlev) = JacD(ii,i,j,ie) ! Form tridiagonal in upper right
          L(ii, ii+nlev+1) = JacU(ii,i,j,ie)
          L(ii+1, ii+nlev) = JacL(ii,i,j,ie)
    
          L(ii+nlev,ii) = 1.d0 ! Form identity in lower left
        end do
        L(nlev, 2*nlev) = JacD(nlev,i,j,ie)
        L(2*nlev, nlev) = 1.d0
        L = g*L
        wphivec = matmul(L,wphivec2(:,i,j,ie)) ! Calculate linear part
  ! subtract linear part
        elem(ie)%state%w_i(i,j,1:nlev,np1) = elem(ie)%state%w_i(i,j,1:nlev,np1) - wphivec(1:nlev)
        elem(ie)%state%phinh_i(i,j,1:nlev,np1) = elem(ie)%state%phinh_i(i,j,1:nlev,np1) - wphivec(1+nlev:2*nlev)
      end do
    end do
  end do
  end subroutine compute_nonlinear_rhs


 
!===========================================================================================================
!===========================================================================================================
!===========================================================================================================
!===========================================================================================================
!===========================================================================================================
  subroutine compute_stage_value_dirk(np1,qn0,dt2,elem,hvcoord,hybrid,&
       deriv,nets,nete,maxiter,itertol)
  !===================================================================================
  ! this subroutine solves a stage value equation for a DIRK method which takes the form
  !
  ! gi = un0 + dt* sum(1:i-1)(aij n(gj)+a2ij s(gj)) + dt *a2ii s(gi) := y + dt a2ii s(gi)
  !
  ! It is assumed that un0 has the value of y and the computed value of gi is stored at
  ! unp1
  !===================================================================================
  integer, intent(in) :: np1,qn0,nets,nete
  real*8, intent(in) :: dt2
  integer :: maxiter
  real*8 :: itertol

  type (hvcoord_t)     , intent(in) :: hvcoord
  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (derivative_t)  , intent(in) :: deriv


  ! local
  real (kind=real_kind), pointer, dimension(:,:,:)   :: phi_np1
  real (kind=real_kind), pointer, dimension(:,:,:)   :: dp3d
  real (kind=real_kind), pointer, dimension(:,:,:)   :: vtheta_dp
  real (kind=real_kind), pointer, dimension(:,:)   :: phis
  real (kind=real_kind) :: JacD(nlev,np,np)  , JacL(nlev-1,np,np)
  real (kind=real_kind) :: JacU(nlev-1,np,np), JacU2(nlev-2,np,np)
  real (kind=real_kind) :: pnh(np,np,nlev)     ! nh (nonydro) pressure
  real (kind=real_kind) :: dp3d_i(np,np,nlevp)
  real (kind=real_kind) :: dpnh_dp_i(np,np,nlevp)
  real (kind=real_kind) :: exner(np,np,nlev)     ! exner nh pressure
  real (kind=real_kind) :: w_n0(np,np,nlevp)    
  real (kind=real_kind) :: phi_n0(np,np,nlevp)    
  real (kind=real_kind) :: Ipiv(nlev,np,np)
  real (kind=real_kind) :: Fn(np,np,nlev),x(nlev,np,np)
  real (kind=real_kind) :: itererr,itererrtemp(np,np)
  real (kind=real_kind) :: itererrmax
  real (kind=real_kind) :: norminfr0(np,np),norminfJ0(np,np)
  real (kind=real_kind) :: maxnorminfJ0r0
  real (kind=real_kind) :: alpha1(np,np),alpha2(np,np)

  real (kind=real_kind) :: Jac2D(nlev,np,np)  , Jac2L(nlev-1,np,np)
  real (kind=real_kind) :: Jac2U(nlev-1,np,np)

  integer :: i,j,k,l,ie,itercount,info(np,np),itercountmax
  integer :: nsafe



  itercountmax=0
  itererrmax=0.d0

  call t_startf('compute_stage_value_dirk')
  do ie=nets,nete
    w_n0 = elem(ie)%state%w_i(:,:,:,np1)
    phi_n0 = elem(ie)%state%phinh_i(:,:,:,np1)
    itercount=0

    ! approximate the initial error of f(x) \approx 0
    dp3d  => elem(ie)%state%dp3d(:,:,:,np1)
    vtheta_dp  => elem(ie)%state%vtheta_dp(:,:,:,np1)
    phi_np1 => elem(ie)%state%phinh_i(:,:,:,np1)
    phis => elem(ie)%state%phis(:,:)

    call pnh_and_exner_from_eos(hvcoord,vtheta_dp,dp3d,phi_np1,pnh,exner,dpnh_dp_i,caller='dirk1')

    dp3d_i(:,:,1) = dp3d(:,:,1)
    dp3d_i(:,:,nlevp) = dp3d(:,:,nlev)
    do k=2,nlev
       dp3d_i(:,:,k)=(dp3d(:,:,k)+dp3d(:,:,k-1))/2
    end do

   ! we first compute the initial Jacobian J0 and residual r0 and their infinity norms
    do k=1,nlev
     Fn(:,:,k) = phi_np1(:,:,k)-phi_n0(:,:,k) &
       - dt2*g*w_n0(:,:,k) + (dt2*g)**2 * (1.0-dpnh_dp_i(:,:,k))
    enddo

     norminfr0=0.d0
     norminfJ0=0.d0
      ! Here's how to call inexact Jacobian
!     call get_dirk_jacobian(Jac2L,Jac2D,Jac2U,dt2,dp3d,phi_np1,pnh,0,&
!       1d-6,hvcoord,dpnh_dp_i,vtheta_dp)
      ! here's the call to the exact Jacobian
     call get_dirk_jacobian(JacL,JacD,JacU,dt2,dp3d,phi_np1,pnh,1)

    ! compute dp3d-weighted infinity norms of the initial Jacobian and residual
     do i=1,np
     do j=1,np
       itererrtemp(i,j)=0 
       do k=1,nlev
        norminfr0(i,j)=max(norminfr0(i,j),abs(Fn(i,j,k)) *dp3d_i(i,j,k))
        if (k.eq.1) then
          norminfJ0(i,j) = max(norminfJ0(i,j),(dp3d_i(i,j,k)*abs(JacD(k,i,j))+dp3d_i(i,j,k+1))*abs(JacU(k,i,j)))
        elseif (k.eq.nlev) then
          norminfJ0(i,j) = max(norminfJ0(i,j),(dp3d_i(i,j,k-1)*abs(JacL(k,i,j))+abs(JacD(k,i,j))*dp3d(i,j,k)))
        else
          norminfJ0(i,j) = max(norminfJ0(i,j),(dp3d_i(i,j,k-1)*abs(JacL(k,i,j))+dp3d_i(i,j,k)*abs(JacD(k,i,j))+ &
            dp3d_i(i,j,k+1)*abs(JacU(k,i,j))))
        end if
        itererrtemp(i,j)=itererrtemp(i,j)+Fn(i,j,k)**2.d0 *dp3d_i(i,j,k)
      end do
      itererrtemp(i,j)=sqrt(itererrtemp(i,j))
    end do
    end do

    maxnorminfJ0r0=max(maxval(norminfJ0(:,:)),maxval(norminfr0(:,:)))
    itererr=maxval(itererrtemp(:,:))/maxnorminfJ0r0
    
    do while ((itercount < maxiter).and.(itererr > itertol))

      info(:,:) = 0
      ! Here's how to call inexact Jacobian
!      call get_dirk_jacobian(JacL,JacD,JacU,dt2,dp3d,phi_np1,pnh,0,&
!       1d-4,hvcoord,dpnh_dp_i,vtheta_dp)
      ! here's the call to the exact Jacobian
       call get_dirk_jacobian(JacL,JacD,JacU,dt2,dp3d,phi_np1,pnh,1)

 
      do i=1,np
      do j=1,np
        x(1:nlev,i,j) = -Fn(i,j,1:nlev)  !+Fn(i,j,nlev+1:2*nlev,1)/(g*dt2))
        call DGTTRF(nlev, JacL(:,i,j), JacD(:,i,j),JacU(:,i,j),JacU2(:,i,j), Ipiv(:,i,j), info(i,j) )
        ! Tridiagonal solve
        call DGTTRS( 'N', nlev,1, JacL(:,i,j), JacD(:,i,j), JacU(:,i,j), JacU2(:,i,j), Ipiv(:,i,j),x(:,i,j), nlev, info(i,j) )
        ! update approximate solution of phi
        phi_np1(i,j,1:nlev) = phi_np1(i,j,1:nlev) + x(1:nlev,i,j)

        do nsafe=1,8
           if (all(phi_np1(i,j,1:nlev) > phi_np1(i,j,2:nlevp))) exit
           ! remove the last netwon increment, try reduced increment
           phi_np1(i,j,1:nlev) = phi_np1(i,j,1:nlev) - x(1:nlev,i,j)/(2**nsafe)
        enddo
        if (nsafe>1) print *,'WARNING: reducing newton increment, nsafe=',nsafe
        ! if nsafe>1, code will probably crash soon
        ! if nsafe>8, code will crash in next call to pnh_and_exner_from_eos
      end do
      end do
      call pnh_and_exner_from_eos(hvcoord,vtheta_dp,dp3d,phi_np1,pnh,exner,dpnh_dp_i,caller='dirk2')

      ! update approximate solution of w
      elem(ie)%state%w_i(:,:,1:nlev,np1) = w_n0(:,:,1:nlev) - g*dt2 * &
        (1.0-dpnh_dp_i(:,:,1:nlev))
      ! update right-hand side of phi
      do k=1,nlev
         Fn(:,:,k) = phi_np1(:,:,k)-phi_n0(:,:,k) &
              - dt2*g*w_n0(:,:,k) + (dt2*g)**2 * (1.0-dpnh_dp_i(:,:,k))
      enddo
      ! compute relative errors
      itererrtemp=0.d0
      do i=1,np
      do j=1,np
        do k=1,nlev
          itererrtemp(i,j)=itererrtemp(i,j)+Fn(i,j,k)**2.d0 *dp3d_i(i,j,k)
        end do
        itererrtemp(i,j)=sqrt(itererrtemp(i,j))
      end do
      end do
      itererr=maxval(itererrtemp(:,:))/maxnorminfJ0r0

      ! update iteration count and error measure
      itercount=itercount+1
    end do ! end do for the do while loop

    if (itercount >= maxiter) then
      call abortmp('Error: nonlinear solver failed b/c max iteration count was met')
    end if
    itercountmax=max(itercount,itercountmax)
    itererrmax=max(itererrmax,itererr)
  end do ! end do for the ie=nets,nete loop

  ! return max iteraitons and max error
  maxiter=itercountmax
  itertol=itererrmax
  call t_stopf('compute_stage_value_dirk')


  end subroutine compute_stage_value_dirk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine matrix_exponential(JacL, JacD, JacU, neg, dimDiag, dt, expJ, w)
  !===================================================================================
  ! Using a Pade approximation,
  ! this subroutine calculates the matrix exponential of the matrix of the form
  !    [ 0       g*dt*T
  !      g*dt*I       0 ],
  ! where the tridiagonal matrix T is given by the input vectors JacL, JacD, and JacU
  !
  ! This matrix exponential is returned as expJ. The product (e^J)*w is also
  ! calculated, and returned as w.
  !===================================================================================

  real (kind=real_kind), dimension(:), intent(in) :: JacL, JacD, JacU
  logical, intent(in) :: neg
  integer, intent(in) :: dimDiag
  real (kind=real_kind), intent(in) :: dt
  real (kind=real_kind), intent(out) :: expJ(2*dimDiag, 2*dimDiag)
  real (kind=real_kind), dimension(:), intent(inout), optional :: w 
  ! local variables
  real (kind=real_kind) :: N(2*dimDiag,2*dimDiag), D(2*dimDiag,2*dimDiag), &
    Aj(2*dimDiag,2*dimDiag), negAj(2*dimDiag,2*dimDiag), Jac(2*dimDiag,2*dimDiag), &
    Dinv(2*dimDiag,2*dimDiag), iden(2*dimDiag,2*dimDiag), DinvN(2*dimDiag,2*dimDiag),&
    Tri(dimDiag,dimDiag)
  real (kind=real_kind) :: work(2*dimDiag)
  integer :: ipiv(2*dimDiag)
  real (kind=real_kind) normJ, pfac, fac, alpha
  integer i,j,p,q,info, maxiter, k, dimJac 

  p = 2  ! parameter used in diagonal Pade approximation
  q = 2
  pfac = 1.d0/gamma(dble(p+q+1.d0))
  ! Initialize random A and normalize
  dimJac = 2*dimDiag
  Jac = 0.d0
  do i = 1,(dimDiag-1)
    Jac(i,(i+dimDiag)) = JacD(i)
    Jac(i,(i+1+dimDiag)) = JacU(i)
    Jac((i+1),(i+dimDiag)) = JacL(i) 
  end do
  Jac(dimDiag, dimJac) = JacD(dimDiag)
  do i = 1,dimDiag
    Jac(i + dimDiag,i) = 1.d0
  end do
  Jac = Jac * g * dt
  if (neg) then
    Jac = (-1.d0)*Jac 
  end if
  ! Scaling by power of 2
  maxiter = 10000
  k = 0
  do while((norm2(Jac)>0.5d0).and.(k<maxiter))
    Jac = Jac / 2.d0
    k = k + 1
  end do ! end while loop

  Tri = 0.d0
  do i = 1, (dimDiag-1)
    Tri(i,i) = Jac(i,(i+dimDiag))
    Tri(i,i+1) = Jac(i,(i+1+dimDiag))
    Tri(i+1,i) = Jac((i+1),(i+dimDiag))
  end do
  Tri(dimDiag, dimDiag) = Jac(dimDiag, dimJac)
  alpha = Jac((dimDiag + 1), 1)  ! scalar multiple of identity in lower left
  Tri = Tri / alpha
  ! Initialize Aj,negAj = identity and N,D = 0.
  N = 0.d0
  D = 0.d0
  Aj = 0.d0
  negAj = 0.d0
  Dinv = 0.d0
  iden = 0.d0
  expJ = 0.d0

  work = 0.d0
  ipiv = 0

  do i = 1,dimJac
    Aj(i,i) = 1.d0
    negAj(i,i) = 1.d0
    iden(i,i) = 1.d0
  enddo ! end do loop

  ! series for Pade approximation
  do i=0,p
    fac = gamma(dble(p+q-i+1.d0))*gamma(dble(p+1.d0))/(gamma(dble(i+1.d0))*gamma(dble(p-i+1.d0)))*pfac
    N = N + fac*Aj
    Aj = matmul(Aj,Jac)
  enddo ! end do loop for Pade approx

  do i=0,q
    fac = gamma(dble(p+q-i+1.d0))*gamma(dble(q+1.d0))/(gamma(dble(i+1.d0))*gamma(dble(q-i+1.d0)))*pfac
    D = D + fac*negAj
    negAj = matmul(negAj,-Jac)
  enddo ! end do loop for Pade approx

  ! Invert matrix D
  call get_DinvN(p, D, N, expJ, Tri, alpha, 2,dimJac) ! using tridiagonal solves
!  call get_DinvN(p, D, N, expJ, Tri, alpha, 1,dimJac)  ! using full LU factorization

  ! Squaring
  do i=1,k
    expJ = matmul(expJ, expJ)
  end do
  if (present(w)) then
    w = matmul(expJ, w)
  end if

  end subroutine matrix_exponential
!===============================================================================

  subroutine get_DinvN(p, D, N, DinvN, Tri, alph, opt,dimJac)
  real (kind=real_kind), dimension(:,:), intent(in) :: D, N, Tri
  integer, intent(in) :: p, opt, dimJac
  real (kind=real_kind), intent(out), target :: DinvN(dimJac,dimJac)
  real (kind=real_kind), intent(in):: alph
 
  ! local variables
  complex(kind=8) :: sig1, sig2, sig1Inv, sig2Inv, kfac, alpha
  integer :: block_dim, info, i
  integer :: ipiv(dimJac)
  complex(kind=8) :: work(dimJac), TriD(dimJac/2), TriL(dimJac/2 - 1), TriU(dimJac/2 - 1)
  complex(kind=8) :: B(dimJac/2,dimJac), X1(dimJac/2,dimJac), X2(dimJac/2,dimJac), N1(dimJac/2,dimJac), N2(dimJac/2,dimJac)

  alpha = dcmplx(alph)
  block_dim = dimJac/2
  ! Variables used to factor Pade approximation
  kfac = (12.d0, 0.d0)
  sig1 = dcmplx(3.d0,sqrt(3.d0))
  sig2 = dcmplx(3.d0, (-sqrt(3.d0)))
  sig1Inv = conjg(sig1)/(real(sig1)**2 + imag(sig1)**2)
  sig2Inv = conjg(sig2)/(real(sig2)**2 + imag(sig2)**2)

  ! Invert matrix D
  DinvN = 0.d0
 
  if (opt == 1) then  ! Calculate inverse using full LU decomp
    DinvN = D
    work = 0.d0
    ipiv = 0
    call DGETRF(dimJac, dimJac, DinvN, dimJac, ipiv, info)
    call DGETRI(dimJac, DinvN, dimJac, ipiv, work, dimJac, info)

    DinvN = matmul(DinvN, N)

  else  ! Use triangular solves and back substitution
    if (p /= 2) then
      stop 'Must have p = 2 approximation' ! Factoring done by hand - only for p=2
    end if
    X1 = 0.d0
    N1 = 0.d0
    N1 = dcmplx(N(1:block_dim, 1:dimJac))
    N2 = 0.d0
    N2 = dcmplx(N(block_dim+1:dimJac, 1:dimJac))
! sig1I-Jac is not nice to invert. We left multiply by (I& 0\\ g*dt*sig1InvI& I) so
! that we can solve the triangular system 
! (-g^2sig1InvTri + sig1I)X2 = (gsig1InvN1+N2)  and back substitute to get
! X1 = sig1Inv(N1+gTriX2)
    do i = 1,block_dim-1
      TriD(i) = dcmplx(Tri(i,i),0.d0)
      TriL(i) = dcmplx(Tri(i+1,i),0.d0)
      TriU(i) = dcmplx(Tri(i,i+1),0.d0)
    end do
   TriD(block_dim) = dcmplx(Tri(block_dim, block_dim), 0.d0)

    ! solve for X1 and X2
    TriD = TriD*(-sig1Inv)*alpha**2
    TriL = TriL*(-sig1Inv)*alpha**2
    TriU = TriU*(-sig1Inv)*alpha**2
 
    do i = 1,block_dim
      TriD(i) = TriD(i) + sig1
    end do

    B = 0.d0
    B = (kfac*(alpha*sig1Inv*N1 + N2))
    call ZGTSV(block_dim, dimJac, TriL, TriD, TriU, B, block_dim, info)


    X2 = B
    X1 = sig1Inv * (kfac*N1 + alpha*matmul(Tri,X2))
! Now do the second tridiag solve
    N1 = X1
    N2 = X2
    do i = 1,block_dim-1  ! ZGTSV writes over Tri, so we have to get it again
      TriD(i) = dcmplx(Tri(i,i),0.d0)
      TriL(i) = dcmplx(Tri(i+1,i),0.d0)
      TriU(i) = dcmplx(Tri(i,i+1),0.d0)
    end do
    TriD(block_dim) = dcmplx(Tri(block_dim, block_dim),0.d0)

    ! solve for X1 and X2
    TriD = TriD*(-sig2Inv)*alpha**2
    TriL = TriL*(-sig2Inv)*alpha**2
    TriU = TriU*(-sig2Inv)*alpha**2
 
    do i = 1,block_dim
      TriD(i) = TriD(i) + sig2
    end do

    B = (alpha*sig2Inv*dcmplx(N1) + dcmplx(N2))

    call ZGTSV(block_dim, dimJac, TriL, TriD, TriU, B, block_dim, info)
    X2 = B

    X1 = sig2Inv * (N1 + alpha*matmul(Tri,X2))
    DinvN(1:block_dim, :) = real(X1)
    DinvN(block_dim+1:dimJac, :) = real(X2)
  end if

  end subroutine get_DinvN
!===========================================================================================================
  subroutine linear_combination_of_elem(t3,a1,t1,a2,t2,elem,nets,nete)
  !===================================================================================
  ! this subroutine calculates the linear combination a1*elem(t1) + a2*elem(t2)
  ! and stores it in elem(t3)
  !
  ! calculated, and returned as w.
  !===================================================================================

  real (kind=real_kind), intent(in)       :: a1,a2
  integer, intent(in)                     :: t1,t2,t3,nets,nete
  type (element_t), intent(inout), target :: elem(:) 

  ! Local variables
  integer :: ie

  do ie = nets, nete
    elem(ie)%state%dp3d(:,:,:,t3)         = elem(ie)%state%dp3d(:,:,:,t1) * a1         + elem(ie)%state%dp3d(:,:,:,t2) * a2
    elem(ie)%state%w_i(:,:,1:nlevp,t3)     = elem(ie)%state%w_i(:,:,1:nlevp,t1) * a1     + elem(ie)%state%w_i(:,:,1:nlevp,t2) * a2
    elem(ie)%state%phinh_i(:,:,1:nlev,t3) = elem(ie)%state%phinh_i(:,:,1:nlev,t1) * a1 + elem(ie)%state%phinh_i(:,:,1:nlev,t2) * a2
    elem(ie)%state%vtheta_dp(:,:,:,t3)    = elem(ie)%state%vtheta_dp(:,:,:,t1) * a1    + elem(ie)%state%vtheta_dp(:,:,:,t2) * a2
    elem(ie)%state%v(:,:,:,:,t3)          = elem(ie)%state%v(:,:,:,:,t1) * a1          + elem(ie)%state%v(:,:,:,:,t2) * a2       
  end do

  end subroutine linear_combination_of_elem
!=============================================================================================
  subroutine expLdtwphi(JacL_elem,JacD_elem,JacU_elem,elem,n0,neg,dt,nets,nete)
  real (kind=real_kind), intent(in)  :: JacL_elem(nlev-1,np,np,nete-nets+1),JacU_elem(nlev-1,np,np,nete-nets+1),JacD_elem(nlev,np,np,nete-nets+1)
  type (element_t), intent(inout), target :: elem(:)
  integer, intent(in) :: n0,nets,nete
  real (kind=real_kind), intent(in) :: dt
  logical, intent(in) :: neg

  ! Local variables
  integer :: ie,i,j
  real (kind=real_kind) :: wphivec(2*nlev)
  real (kind=real_kind) :: expJ(2*nlev,2*nlev) 


  do ie = nets, nete
    do i = 1,np
      do j = 1,np
        wphivec(1:nlev) = elem(ie)%state%w_i(i,j,1:nlev,n0)
        wphivec(1+nlev:2*nlev) = elem(ie)%state%phinh_i(i,j,1:nlev,n0)
        call matrix_exponential(JacL_elem(:,i,j,ie),JacD_elem(:,i,j,ie),JacU_elem(:,i,j,ie),neg,nlev,dt,expJ,wphivec) ! Pade approximation
!        call matrix_exponential2(JacL_elem(:,i,j,ie),JacD_elem(:,i,j,ie),JacU_elem(:,i,j,ie),neg,nlev,dt,expJ,wphivec) ! Taylor approximation
        elem(ie)%state%w_i(i,j,1:nlev,n0) = wphivec(1:nlev)
        elem(ie)%state%phinh_i(i,j,1:nlev,n0) = wphivec(1+nlev:2*nlev)
      end do
    end do
  end do

  end subroutine expLdtwphi
  subroutine matrix_exponential2(JacL, JacD, JacU, neg, dimDiag, dt, expJ, w)
  !===================================================================================
  ! Using a Taylor approximation,
  ! this subroutine calculates the matrix exponential of the matrix of the form
  !    [ 0       g*dt*T
  !      g*dt*I       0 ],
  ! where the tridiagonal matrix T is given by the input vectors JacL, JacD, and JacU
  !
  ! This matrix exponential is returned as expJ. The product (e^J)*w is also
  ! calculated, and returned as w.
  !===================================================================================

  real (kind=real_kind), dimension(:), intent(in) :: JacL, JacD, JacU
  logical, intent(in) :: neg
  integer, intent(in) :: dimDiag
  real (kind=real_kind), intent(in) :: dt
  real (kind=real_kind), intent(out) :: expJ(2*dimDiag, 2*dimDiag)
  real (kind=real_kind), dimension(:), intent(inout), optional :: w 
  ! local variables
  real (kind=real_kind) :: Aj(2*dimDiag,2*dimDiag), Jac(2*dimDiag,2*dimDiag)
  real (kind=real_kind) :: work(2*dimDiag)
  integer :: ipiv(2*dimDiag)
  real (kind=real_kind) normJ, pfac, fac, alpha
  integer i,j,p,info, maxiter, k, dimJac 

  p = 10  ! number of terms in the series
  ! Initialize random A and normalize
  dimJac = 2*dimDiag
  Jac = 0.d0
  do i = 1,(dimDiag-1)
    Jac(i,(i+dimDiag)) = JacD(i)
    Jac(i,(i+1+dimDiag)) = JacU(i)
    Jac((i+1),(i+dimDiag)) = JacL(i) 
  end do
  Jac(dimDiag, dimJac) = JacD(dimDiag)
  do i = 1,dimDiag
    Jac(i + dimDiag,i) = 1.d0
  end do
  Jac = Jac * g * dt
  if (neg) then
    Jac = (-1.d0)*Jac 
  end if

  ! Initialize Aj,negAj = identity and N,D = 0.
  Aj = 0.d0
  expJ = 0.d0

  do i = 1,dimJac
    Aj(i,i) = 1.d0
  enddo ! end do loop

  ! series for approximation
  do i=0,p
    expJ = expJ + 1/gamma(dble(i+1.d0))*Aj
    Aj = matmul(Aj,Jac)
  enddo ! end do loop for Pade approx

  if (present(w)) then
    w = matmul(expJ, w)
  end if

  end subroutine matrix_exponential2

subroutine phi1Ldt(JacL_elem,JacD_elem,JacU_elem,elem,n0,neg,dt,nets,nete)
  real (kind=real_kind), intent(in)  :: JacL_elem(nlev-1,np,np,nete-nets+1),JacU_elem(nlev-1,np,np,nete-nets+1),JacD_elem(nlev,np,np,nete-nets+1)
  type (element_t), intent(inout), target :: elem(:)
  integer, intent(in) :: n0,nets,nete
  real (kind=real_kind), intent(in) :: dt
  logical, intent(in) :: neg

  ! Local variables
  integer :: ie,i,j,info,k
  real (kind=real_kind) :: wphivec(2*nlev)
  real (kind=real_kind) :: expJ(2*nlev,2*nlev),Jac(2*nlev,2*nlev),iden(2*nlev,2*nlev) 
  real (kind=real_kind) :: work(2*nlev)
  integer :: ipiv(2*nlev)

  iden = 0.d0
  do i = 1,2*nlev
    iden(i,i) = 1.d0 
  end do

  do ie = nets, nete
    do i = 1,np
      do j = 1,np
        call matrix_exponential(JacL_elem(:,i,j,ie),JacD_elem(:,i,j,ie),JacU_elem(:,i,j,ie),neg,nlev,dt,expJ)
        expJ = expJ - iden

        Jac = 0.d0
        do k = 1,nlev-1
          Jac(k,k+nlev) = JacD_elem(k,i,j,ie)
          Jac(k+1,k+nlev) = JacL_elem(k,i,j,ie)
          Jac(k,k+nlev+1) = JacU_elem(k,i,j,ie)

          Jac(k+nlev,k) = 1.d0
        end do
        Jac(nlev,2*nlev) = JacD_elem(nlev,i,j,ie)
        Jac(2*nlev,nlev) = 1.d0
        Jac = Jac * g

        work = 0.d0
        ipiv = 0
        call DGETRF(2*nlev, 2*nlev, Jac, 2*nlev, ipiv, info)
        if (info .ne. 0) then
          print *, "error 1!"
        end if
        call DGETRI(2*nlev, Jac, 2*nlev, ipiv, work, 2*nlev, info)
        if (info .ne. 0) then
          print *, "error 2!"
        end if
        expJ = matmul(Jac, expJ)
       
        wphivec(1:nlev) = elem(ie)%state%w_i(i,j,1:nlev,n0)
        wphivec(1+nlev:2*nlev) = elem(ie)%state%phinh_i(i,j,1:nlev,n0)
        wphivec = matmul(expJ,wphivec)
       
        elem(ie)%state%w_i(i,j,1:nlev,n0) = wphivec(1:nlev)
        elem(ie)%state%phinh_i(i,j,1:nlev,n0) = wphivec(1+nlev:2*nlev)

      end do
    end do
  end do

  end subroutine phi1Ldt


!=========================================================================
! Stores elem(np1) in mdarray "stage"
!=========================================================================
  subroutine store_state(elem,np1,nets,nete,stage)
  
  type(element_t), intent(in) :: elem(:)
  integer, intent(in) :: np1, nets, nete
  real (kind=real_kind), intent(out) :: stage(nete-nets+1,np,np,nlevp,6)

  ! local variables
  integer :: ie

  do ie = nets, nete
    stage(ie,:,:,1:nlev,1) = elem(ie)%state%v(:,:,1,:,np1)
    stage(ie,:,:,1:nlev,2) = elem(ie)%state%v(:,:,2,:,np1)
    stage(ie,:,:,:,3)      = elem(ie)%state%w_i(:,:,:,np1)
    stage(ie,:,:,:,4)      = elem(ie)%state%phinh_i(:,:,:,np1)
    stage(ie,:,:,1:nlev,5) = elem(ie)%state%vtheta_dp(:,:,:,np1)
    stage(ie,:,:,1:nlev,6) = elem(ie)%state%dp3d(:,:,:,np1)
  end do

  end subroutine store_state



!==========================================================================
! Retrieves mdarray "stage" and stores at elem(nm1)
!==========================================================================
  subroutine retrieve_state(stage,elem,nm1,nets,nete)
  
  real (kind=real_kind), intent(in) :: stage(nete-nets+1,np,np,nlevp,6)
  type(element_t), intent(inout)    :: elem(:)
  integer, intent(in)               :: nm1,nets,nete

  ! local variables
  integer :: ie

  do ie = nets,nete
    elem(ie)%state%v(:,:,1,:,nm1)       = stage(ie,:,:,1:nlev,1)
    elem(ie)%state%v(:,:,2,:,nm1)       = stage(ie,:,:,1:nlev,2)
    elem(ie)%state%w_i(:,:,:,nm1)       = stage(ie,:,:,:,3)
    elem(ie)%state%phinh_i(:,:,:,nm1)   = stage(ie,:,:,:,4)
    elem(ie)%state%vtheta_dp(:,:,:,nm1) = stage(ie,:,:,1:nlev,5)
    elem(ie)%state%dp3d(:,:,:,nm1)      = stage(ie,:,:,1:nlev,6)
  end do

  end subroutine retrieve_state


!============================================================================
! This subroutine calculates the phi function phi_k(Ldt) and applies it
! to the state variables located in elem(n0).
!============================================================================
  subroutine apply_phi_func(JacL_elem,JacD_elem,JacU_elem,dt,deg,n0,elem,nets,nete)

  real(kind=real_kind), intent(in)  :: JacL_elem(nlev-1,np,np,nete-nets+1),&
       JacD_elem(nlev,np,np,nete-nets+1),JacU_elem(nlev-1,np,np,nete-nets+1),dt
  integer, intent(in)               :: deg,n0,nets,nete
  type(element_t), intent(inout)    :: elem(:)

  ! local variables
  real(kind=real_kind)  :: Jac(2*nlev,2*nlev),JInv(2*nlev,2*nlev),expJ(2*nlev,2*nlev),&
      wphivec(2*nlev),c(2*nlev), phi_k(2*nlev,np,np,nete-nets+1)
  integer               :: i,j,k,ie,info,dimJac
  integer               :: ipiv(nlev)
  complex(kind=8)       :: work(nlev)
  
  do ie = nets, nete
    do i = 1,np
      do j = 1,np
        ! Form Jacobian
        dimJac = 2*nlev
        Jac = 0.d0
        do k = 1,(nlev-1)
          Jac(k,(k+nlev))     = JacD_elem(k,i,j,ie)
          Jac(k,(k+1+nlev))   = JacU_elem(k,i,j,ie)
          Jac((k+1),(k+nlev)) = JacL_elem(k,i,j,ie)
        end do
        Jac(nlev, dimJac)     = JacD_elem(nlev,i,j,ie)
        do k = 1,nlev
          Jac(k + nlev,k) = 1.d0
        end do
        Jac = Jac * g * dt

        ! Calculate Jac^(-1) for later
        JInv = Jac
        work = 0.d0
        ipiv = 0
        call DGETRF(dimJac, dimJac, JInv, dimJac, ipiv, info)
        call DGETRI(dimJac, JInv, dimJac, ipiv, work, dimJac, info)
        wphivec(1:nlev)        = elem(ie)%state%w_i(i,j,1:nlev,n0)
        wphivec(nlev+1:2*nlev) = elem(ie)%state%phinh_i(i,j,1:nlev,n0)
        c = wphivec
        call matrix_exponential(JacL_elem(:,i,j,ie),JacD_elem(:,i,j,ie),JacU_elem(:,i,j,ie),.false.,nlev,dt,expJ,c)
        c = c - wphivec
        phi_k(:,i,j,ie) = matmul(JInv,c) 
        ! calculate phi_k recursively
        if (deg .gt. 2) then
          do k = 2,deg
            phi_k(:,i,j,ie) = phi_k(:,i,j,ie) - 1/gamma(dble(k))*wphivec
            phi_k(:,i,j,ie) = matmul(JInv,phi_k(:,i,j,ie))
          end do
        end if
      end do
    end do
  end do
  call phi_func_update_elem(elem,n0,phi_k,deg,nets,nete)
  end subroutine apply_phi_func

!===================================================================
! The subroutine phi_func_update_elem updates the state variables  
!   under the action of a phi function. 
! 
! The variables w and phi are updated directly with phi_k. 
! Because of the structure of the Jacobian, the other variables are 
!   modified by the scalar coming from the first term of the series for
!   phi_k
!=====================================================================
  subroutine phi_func_update_elem(elem,n0,phifunc,deg,nets,nete)
  type(element_t), intent(inout)   :: elem(:)
  real(kind=real_kind), intent(in) :: phifunc(2*nlev,np,np,nete-nets+1)
  integer, intent(in)              :: nets,nete,n0,deg
  ! local variables
  integer :: ie,i,j
 
  ! update w and phi with the action of phi_func
  do ie = nets,nete 
    do i = 1,np
      do j = 1,np
        elem(ie)%state%w_i(i,j,1:nlev,n0)     = phifunc(1:nlev,i,j,ie)
        elem(ie)%state%phinh_i(i,j,1:nlev,n0) = phifunc(nlev+1:2*nlev,i,j,ie)
      end do
    end do
  end do
  ! update the other variables with the coefficient from the first term of phi_k
  if (deg .gt. 1) then
    do ie = nets,nete 
      elem(ie)%state%v(:,:,1,:,n0)       = 1/gamma(dble(deg)+1.d0)*elem(ie)%state%v(:,:,1,:,n0)
      elem(ie)%state%v(:,:,2,:,n0)       = 1/gamma(dble(deg)+1.d0)*elem(ie)%state%v(:,:,2,:,n0)
      elem(ie)%state%vtheta_dp(:,:,:,n0) = 1/gamma(dble(deg)+1.d0)*elem(ie)%state%vtheta_dp(:,:,:,n0)
      elem(ie)%state%dp3d(:,:,:,n0)      = 1/gamma(dble(deg)+1.d0)*elem(ie)%state%dp3d(:,:,:,n0)
    end do
  end if
  end subroutine phi_func_update_elem

end module prim_advance_mod

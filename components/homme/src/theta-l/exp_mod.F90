module exp_mod

  use physical_constants, only: Cp, cp, cpwater_vapor, g, kappa, Rgas, Rwater_vapor, p0, TREF
  use kinds,              only: real_kind
  use element_mod,        only: element_t
  use dimensions_mod,     only: nlev, nlevp, np
  use hybvcoord_mod,      only: hvcoord_t
  use control_mod,        only: theta_hydrostatic_mode
  use eos,                only: pnh_and_exner_from_eos
  use derivative_mod,     only: derivative_t
  use hybrid_mod,         only: hybrid_t

  implicit none
  private
  save
  public :: matrix_exponential, matrix_exponential2, phi_func, getLu,add_Lu,& 
              formJac,tri_inv,tri_mult,apply_phi_func, apply_phi_func_new,& 
              phi_func_new,store_state,linear_combination_of_elem,&
              phi1Ldt,retrieve_state,&
              expLdtwphi,get_exp_jacobian
  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine matrix_exponential(JacL, JacD, JacU, dimDiag, dt, expJ, w)
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
  call formJac(JacL,JacD,JacU,dt,Jac)
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine matrix_exponential_new(JacL, JacD, JacU, dt, expJ, wphi)
  !===================================================================================
  ! Using a Pade approximation,
  ! this subroutine calculates the matrix exponential of the matrix of the form
  !    [ 0       g*dt*T
  !      g*dt*I       0 ],
  ! where the tridiagonal matrix T is given by the input vectors JacL, JacD, and JacU
  !
  ! This matrix exponential is returned as expJ. The product (e^J)*wphi is also
  ! calculated, and returned as wphi.
  !===================================================================================

  real(kind=real_kind), intent(in)  :: JacL(nlev-1), JacD(nlev), JacU(nlev-1), dt
  real(kind=real_kind), intent(out) :: expJ(2*nlev,2*nlev)
  real(kind=real_kind), intent(inout), optional :: wphi(2*nlev)

  ! local variables
  real(kind=real_kind) :: N(2*nlev,2*nlev), D(2*nlev,2*nlev), DinvN(2*nlev,2*nlev),&
                          scaling_const, Jac(2*nlev,2*nlev), normJac
  integer :: i, k, maxiter

! TO DO: Approximate norm without forming Jacobian
  ! Scaling by power of 2
  maxiter = 15
  k = 0
  call formJac(JacL,JacD,JacU,dt,Jac)
  do while((norm2(Jac)>0.5d0).and.(k<maxiter))
    Jac = Jac / 2.d0
    k = k + 1
  end do ! end while loop

!  normJac = g*dt
!  do while(normJac > 0.5/sqrt(real(nlev)))
!    normJac = normJac/2.d0
!    k = k+1
!  end do

  scaling_const = g*dt/(2**k) ! scaling and squaring normalization constant

  ! form N(A) = I + 1/2 A + 1/12 A^2
  N = 0.d0
  do i = 1,nlev-1
    ! N_11 block
    N(i,i)   = scaling_const**2/12.d0*JacD(i) + 1.d0
    N(i,i+1) = scaling_const**2/12.d0*JacU(i)
    N(i+1,i) = scaling_const**2/12.d0*JacL(i)

    ! N_12 block
    N(i,i+nlev)   = scaling_const/2.d0*JacD(i)
    N(i,i+1+nlev) = scaling_const/2.d0*JacU(i)
    N(i+1,i+nlev) = scaling_const/2.d0*JacL(i)
 
    ! N_21 block
    N(i+nlev,i) = scaling_const/2.d0 

    ! N_22 block
    N(i+nlev,i+nlev)   = scaling_const**2/12.d0*JacD(i) + 1.d0
    N(i+nlev,i+1+nlev) = scaling_const**2/12.d0*JacU(i)
    N(i+1+nlev,i+nlev) = scaling_const**2/12.d0*JacL(i)
  end do
  N(nlev,nlev)     = scaling_const**2/12.d0*JacD(nlev)+1.d0
  N(2*nlev,2*nlev) = scaling_const**2/12.d0*JacD(nlev)+1.d0
  N(nlev,2*nlev)   = scaling_const/2.d0*JacD(nlev)
  N(2*nlev,nlev)   = scaling_const/2.d0

  ! form D(A) = I + 1/2 A + 1/12 A^2 (only differs from N(A) in sign on the
  ! off-diagonal blocks)
  D(1:nlev,1:nlev)               = N(1:nlev,1:nlev)
  D(1:nlev,1+nlev:2*nlev)        = -N(1:nlev,1+nlev:2*nlev)
  D(1+nlev:2*nlev,1:nlev)        = -N(1+nlev:2*nlev,1:nlev)
  D(1+nlev:2*nlev,1+nlev:2*nlev) = N(1+nlev:2*nlev,1+nlev:2*nlev)

  ! Invert matrix D
  call get_DinvN_new(D,N,expJ,scaling_const*JacL,scaling_const*JacD,scaling_const*JacU,dcmplx(scaling_const))

  ! Squaring
  do i=1,k
    expJ = matmul(expJ, expJ)
  end do

  ! multiply by wphi if necessary
  if (present(wphi)) then
    wphi = matmul(expJ, wphi)
  end if

  end subroutine matrix_exponential_new
!===============================================================================


!===============================================================================

  subroutine get_DinvN(p, D, N, DinvN, Tri, alph, opt,dimJac)
  real (kind=real_kind), dimension(:,:), intent(in) :: D, N, Tri
  integer, intent(in) :: p, opt, dimJac
  real (kind=real_kind), intent(out), target :: DinvN(dimJac,dimJac)
  real (kind=real_kind), intent(in):: alph
 
  ! local variables
  real (kind=real_kind) :: DinvN_LU(dimJac,dimJac)
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

!    DinvN_LU = matmul(DinvN, N)
    DinvN = matmul(DinvN,N)
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

!    print *, "Difference is ", norm2(DinvN_LU - DinvN)
!    stop
  end if

  end subroutine get_DinvN
!===============================================================================

  subroutine get_DinvN_new(D, N, DinvN, JacL,JacD,JacU,scaling_const)
  real(kind=real_kind), intent(in)  :: D(2*nlev,2*nlev), N(2*nlev,2*nlev),&
                                        JacL(nlev-1),JacD(nlev),JacU(nlev-1)
  real(kind=real_kind), intent(out) :: DinvN(2*nlev,2*nlev)
  complex(kind=8),      intent(in)  :: scaling_const
 
  ! local variables
  complex(kind=8) :: sig1, sig2, sig1Inv, sig2Inv, pade_const
  complex(kind=8) :: work(2*nlev), TriL(nlev-1), TriD(nlev), TriU(nlev-1),&
                     Tri_prod(nlev,nlev), du2(nlev-2), N11(nlev,nlev),&
                     N12(nlev,nlev),  N21(nlev,nlev),  N22(nlev,nlev),&
                     X11(nlev,nlev),  X12(nlev,nlev),  X21(nlev,nlev),&
                     X22(nlev,nlev)
  integer         :: info, i, ipiv(2*nlev)

  ! Constants for (2,2) Pade approximant
  pade_const = (12.d0, 0.d0)
  sig1       = dcmplx(3.d0,sqrt(3.d0))
  sig2       = dcmplx(3.d0, (-sqrt(3.d0)))
  sig1Inv    = conjg(sig1)/(real(sig1)**2 + imag(sig1)**2)
  sig2Inv    = conjg(sig2)/(real(sig2)**2 + imag(sig2)**2)

! STEP 1:
  ! solving the system
  !    [fac1] [x11 x12] = pade_const*[nhat], where
  !           [x21 x22]
  !
  ! fac1 = [sig i                -jac                 ]
  !        [  0       sig i - scaling_const*siginv*jac], and
  !
  ! nhat = [              n11                                    n12             ]
  !        [scaling_const*siginv*n11 + n21         scaling_const*siginv*n12 + n22]
  
  ! Blocks of N
  N11 = dcmplx(N(1:nlev,1:nlev))
  N12 = dcmplx(N(1:nlev,1+nlev:2*nlev))
  N21 = dcmplx(N(1+nlev:2*nlev,1:nlev))
  N22 = dcmplx(N(1+nlev:2*nlev,1+nlev:2*nlev))

  ! Tridiagonal matrix: sig1*I - scaling_const*sig1Inv*Jac
  TriD = dcmplx(JacD)*(-sig1Inv)*scaling_const + sig1
  TriL = dcmplx(JacL)*(-sig1Inv)*scaling_const
  TriU = dcmplx(JacU)*(-sig1Inv)*scaling_const

  ! Tridiagonal solve to get X21, X22
  call zgttrf(nlev,tril,trid,triu,du2,ipiv,info)                  ! LU decomposition
  X21 = pade_const*(scaling_const*sig1Inv*N11 + N21)              ! Nhat_21 block
  call zgttrs('N',nlev,nlev,TriL,TriD,TriU,du2,ipiv,X21,nlev,info)! Solve for X21
  X22 = pade_const*(scaling_const*sig1Inv*N12 + N22)              ! Nhat_22 block
  call zgttrs('N',nlev,nlev,TriL,TriD,TriU,du2,ipiv,X22,nlev,info)! Solve for X22
  
  ! solve for X11 and X12
  call tri_mult(JacL,JacD,JacU,X21,Tri_prod,nlev)  ! Compute T*X21
  X11 = sig1Inv * (pade_const*N11 + Tri_prod)      ! Compute X11
  call tri_mult(JacL,JacD,JacU,X22,Tri_prod,nlev)  ! Compute T*X22
  X12 = sig1Inv * (pade_const*N12 + Tri_prod)      ! Compute X12
 
! STEP 2:
  ! solve the system
  !    [fac2] [DinvN] = [Nhat], where
  !
  ! fac2 = [sig I                -Jac                 ]
  !        [  0       sig I - scaling_const*sigInv*Jac], and
  !
  ! Nhat = [              X11                                    X12             ]
  !        [scaling_const*sigInv*X11 + X21         scaling_const*sigInv*X12 + X22]
  !   using the Xij from the previous solve.

  ! Blocks of Nhat
  N11 = X11
  N12 = X12
  N21 = X21
  N22 = X22

  ! Tridiagonal matrix T = sig2 I - scaling_const*sig2Inv*Jac
  TriD = dcmplx(JacD)*(-sig2Inv)*scaling_const + sig2
  TriL = dcmplx(JacL)*(-sig2Inv)*scaling_const
  TriU = dcmplx(JacU)*(-sig2Inv)*scaling_const

  ! Tridiagonal solve to get X21, X22
  call zgttrf(nlev,TriL,TriD,TriU,du2,ipiv,info)                   ! LU decomposition
  X21 = scaling_const*(sig2Inv)*N11+N21                            ! Nhat_21 block
  call zgttrs('N',nlev,nlev,TriL,TriD,TriU,du2,ipiv,X21,nlev,info) ! Get X21
  X22 = scaling_const*(sig2Inv)*N12+N22                            ! Nhat_22 block
  call zgttrs('N',nlev,nlev,TriL,TriD,TriU,du2,ipiv,X22,nlev,info) ! Get X22

  ! Back substitution for X11, X12
  call tri_mult(JacL,JacD,JacU,X21,Tri_prod,nlev)  ! Compute T*X21
  X11 = sig2Inv * (N11 + Tri_prod)                 ! Compute X11
  call tri_mult(JacL,JacD,JacU,X22,Tri_prod,nlev)  ! Compute T*X22
  X12 = sig2Inv * (N12 + Tri_prod)                 ! Compute X12

  ! return DinvN
  DinvN(1:nlev,1:nlev)               = X11
  DinvN(1:nlev,nlev+1:2*nlev)        = X12
  DinvN(1+nlev:2*nlev,1:nlev)        = X21
  DinvN(1+nlev:2*nlev,1+nlev:2*nlev) = X22

  end subroutine get_DinvN_new
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
  subroutine matrix_exponential2(JacL, JacD, JacU, dimDiag, dt, expJ, w)
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
      wphivec(2*nlev),c(2*nlev), phi_k(2*nlev)
  integer               :: i,j,k,ie,info
  integer               :: ipiv(nlev)
  complex(kind=8)       :: work(nlev)
  
  do ie = nets, nete
    do i = 1,np
      do j = 1,np
        wphivec(1:nlev)        = elem(ie)%state%w_i(i,j,1:nlev,n0)
        wphivec(nlev+1:2*nlev) = elem(ie)%state%phinh_i(i,j,1:nlev,n0)
        call phi_func(JacL_elem(:,i,j,ie),JacD_elem(:,i,j,ie),JacU_elem(:,i,j,ie),dt,deg,wphivec,phi_k)
        ! update w and phi
        elem(ie)%state%w_i(i,j,1:nlev,n0)     = phi_k(1:nlev)
        elem(ie)%state%phinh_i(i,j,1:nlev,n0) = phi_k(nlev+1:2*nlev)
      end do
    end do
    if (deg .gt. 1) then ! update other variables
      elem(ie)%state%v(:,:,1,:,n0)       = 1/gamma(dble(deg)+1.d0)*elem(ie)%state%v(:,:,1,:,n0)
      elem(ie)%state%v(:,:,2,:,n0)       = 1/gamma(dble(deg)+1.d0)*elem(ie)%state%v(:,:,2,:,n0)
      elem(ie)%state%vtheta_dp(:,:,:,n0) = 1/gamma(dble(deg)+1.d0)*elem(ie)%state%vtheta_dp(:,:,:,n0)
      elem(ie)%state%dp3d(:,:,:,n0)      = 1/gamma(dble(deg)+1.d0)*elem(ie)%state%dp3d(:,:,:,n0)
    end if
  end do
  end subroutine apply_phi_func

!============================================================================
! This subroutine calculates the phi function phi_k(Ldt) and applies it
! to a vector c.
!============================================================================
  subroutine phi_func(JacL,JacD,JacU,dt,deg,wphivec,phi_k)

  real(kind=real_kind), intent(in)  :: JacL(nlev-1),&
       JacD(nlev),JacU(nlev-1),dt,wphivec(2*nlev)
  real(kind=real_kind), intent(out) :: phi_k(2*nlev)
  integer, intent(in)               :: deg

  ! local variables
  real(kind=real_kind)  :: Jac(2*nlev,2*nlev),JInv(2*nlev,2*nlev),expJ(2*nlev,2*nlev),&
      c(2*nlev),c_1(nlev),c_2(nlev),phi_k2(2*nlev),du2(nlev-2),eye(nlev,nlev)
  real(kind=real_kind)  :: Lcopy(nlev-1), Dcopy(nlev), Ucopy(nlev-1)
  integer               :: k,info
  integer               :: ipiv(nlev)
  complex(kind=8)       :: work(2*nlev)
  
  Lcopy = JacL
  Dcopy = JacD
  Ucopy = JacU
!  call formJac(JacL,JacD,JacU,dt,Jac)
  ! Calculate Jac^(-1) for later
!  JInv = Jac
  work = 0.d0
  ipiv = 0
!  call DGETRF(2*nlev, 2*nlev, JInv, 2*nlev, ipiv, info)
!  call DGETRI(2*nlev, JInv, 2*nlev, ipiv, work, 2*nlev, info)

  c = wphivec
!  call matrix_exponential(JacL,JacD,JacU,nlev,dt,expJ,c)
  call matrix_exponential_new(JacL,JacD,JacU,dt,expJ,c)
  c = c - wphivec
  c_1 = c(1:nlev)
  c_2 = c(nlev+1:2*nlev)
  phi_k = 0.d0
  phi_k(1:nlev) = c_2

  call DGTTRF(nlev,Lcopy,Dcopy,Ucopy,du2,ipiv,info)
  call DGTTRS('N',nlev,1,Lcopy,Dcopy,Ucopy,du2,ipiv,c_1,nlev,info)

  phi_k(nlev+1:2*nlev) = c_1
  phi_k = phi_k/(g*dt)

!  phi_k = matmul(JInv,c) 
!  print *, "Error of phi1: ", norm2(phi_k - phi_k2)
!  stop
!  if (norm2(phi_k-phi_k2)>1e-10) then
!    print *, "Error: check accuracy of phi function"
!    print *, "Error = ", norm2(phi_k-phi_k2)
!    stop
!  end if


  ! calculate phi_k recursively
  if (deg .ge. 2) then
    do k = 2,deg
      phi_k = phi_k - 1/gamma(dble(k))*wphivec
      c_1 = phi_k(1:nlev)
      c_2 = phi_k(nlev+1:2*nlev)
      call DGTTRS('N',nlev,1,Lcopy,Dcopy,Ucopy,du2,ipiv,c_1,nlev,info)
!      phi_k = matmul(JInv,phi_k)
      phi_k(1:nlev) = c_2/(g*dt)
      phi_k(nlev+1:2*nlev) = c_1/(g*dt)
!      if (norm2(phi_k-phi_k2)>1e-10) then
!        print *, "Error of phi2 is ", norm2(phi_k2-phi_k)
!        stop
!      end if
    end do
  end if
  end subroutine phi_func

!============================================================================
! This subroutine calculates the phi function phi_k(Ldt) and applies it
! to a vector c.
!============================================================================
  subroutine phi_func_new(JacL_elem,JacD_elem,JacU_elem,dt,phi_deg,phi_func,elem,n0,nets,nete)

  real(kind=real_kind), intent(in)  :: JacL_elem(nlev-1,np,np,nete-nets+1),&
                                       JacD_elem(nlev,np,np,nete-nets+1),&
                                       JacU_elem(nlev-1,np,np,nete-nets+1),dt
  real(kind=real_kind), intent(out) :: phi_func(np,np,2*nlev,nete-nets+1,3)
  integer, intent(in)               :: phi_deg,nets,nete,n0
  type(element_t), intent(in)    :: elem(:)
 

  ! local variables
  real(kind=real_kind)  :: expJ(2*nlev,2*nlev), c(2*nlev),c_1(nlev),c_2(nlev),du2(nlev-2)
  real(kind=real_kind)  :: Lcopy(nlev-1), Dcopy(nlev),Ucopy(nlev-1),wphivec(2*nlev)
  integer               :: k,info,ie,i,j
  integer               :: ipiv(nlev)
  complex(kind=8)       :: work(2*nlev)
  phi_func = 0.d0
  
  do ie = nets,nete
    do i = 1,np
      do j = 1,np
        Lcopy = JacL_elem(:,i,j,ie)
        Dcopy = JacD_elem(:,i,j,ie)
        Ucopy = JacU_elem(:,i,j,ie)
        work = 0.d0
        ipiv = 0

        wphivec(1:nlev) = elem(ie)%state%w_i(i,j,1:nlev,n0)
        wphivec(1+nlev:2*nlev) = elem(ie)%state%phinh_i(i,j,1:nlev,n0)
        c = wphivec
        call matrix_exponential_new(Lcopy,Dcopy,Ucopy,dt,expJ,c)
        c = c - wphivec
        c_1 = c(1:nlev)
        c_2 = c(nlev+1:2*nlev)
        phi_func(i,j,1:nlev,ie,1) = c_2

        call DGTTRF(nlev,Lcopy,Dcopy,Ucopy,du2,ipiv,info)
        call DGTTRS('N',nlev,1,Lcopy,Dcopy,Ucopy,du2,ipiv,c_1,nlev,info)

        phi_func(i,j,nlev+1:2*nlev,ie,1) = c_1
        phi_func(i,j,:,ie,1) = phi_func(i,j,:,ie,1)/(g*dt)

        ! calculate phi_k recursively
        if (phi_deg .ge. 2) then
          do k = 2,phi_deg
            phi_func(i,j,:,ie,k) = phi_func(i,j,:,ie,k-1) - 1/gamma(dble(k))*wphivec
            c_1 = phi_func(i,j,1:nlev,ie,k)
            c_2 = phi_func(i,j,nlev+1:2*nlev,ie,k)
            call DGTTRS('N',nlev,1,Lcopy,Dcopy,Ucopy,du2,ipiv,c_1,nlev,info)
            phi_func(i,j,1:nlev,ie,k) = c_2/(g*dt)
            phi_func(i,j,nlev+1:2*nlev,ie,k) = c_1/(g*dt)
          end do
        end if
      end do
    end do
  end do
  end subroutine phi_func_new

!============================================================================
! This subroutine calculates the phi function phi_k(Ldt) and applies it
! to the state variables located in elem(n0).
!============================================================================
  subroutine apply_phi_func_new(phi_func,struct_dim,deg,elem,n0,nets,nete)

  integer, intent(in)               :: deg,n0,nets,nete,struct_dim
  type(element_t), intent(inout)    :: elem(:)
  real (kind=real_kind), intent(in) :: phi_func(np,np,2*nlev,nete-nets+1,3)

  ! local variables
  integer :: ie
  
  do ie = nets, nete
        elem(ie)%state%w_i(:,:,1:nlev,n0)     = phi_func(:,:,1:nlev,ie,deg)
        elem(ie)%state%phinh_i(:,:,1:nlev,n0) = phi_func(:,:,nlev+1:2*nlev,ie,deg)
    if (deg .gt. 1) then ! update other variables
      elem(ie)%state%v(:,:,1,:,n0)       = 1/gamma(dble(deg)+1.d0)*elem(ie)%state%v(:,:,1,:,n0)
      elem(ie)%state%v(:,:,2,:,n0)       = 1/gamma(dble(deg)+1.d0)*elem(ie)%state%v(:,:,2,:,n0)
      elem(ie)%state%vtheta_dp(:,:,:,n0) = 1/gamma(dble(deg)+1.d0)*elem(ie)%state%vtheta_dp(:,:,:,n0)
      elem(ie)%state%dp3d(:,:,:,n0)      = 1/gamma(dble(deg)+1.d0)*elem(ie)%state%dp3d(:,:,:,n0)
    end if
  end do
  end subroutine apply_phi_func_new

!=================================================================================================
! Subroutine getLu calculates the Jacobian L_m and multiplies it by u_m. The
! product is stored and returned in the vector Lu 
!=================================================================================================
  subroutine getLu(JacL_elem,JacD_elem,JacU_elem,elem,n0,nets,nete,Lu)
  real (kind=real_kind), intent(in)  :: JacL_elem(nlev-1,np,np,nete-nets+1),&
       JacU_elem(nlev-1,np,np,nete-nets+1),JacD_elem(nlev,np,np,nete-nets+1)
  type (element_t), intent(in) :: elem(:)
  integer, intent(in) :: n0,nets,nete
  real (kind=real_kind), intent(out) :: Lu(np,np,2*nlev,nete-nets+1)

  ! Local variables
  integer :: ie,i,j,k
! real (kind=real_kind) :: Jac(2*nlev,2*nlev) 
! real (kind=real_kind) :: Lu_matmul(np,np,2*nlev,nete-nets+1) ! for unit test

! do ie = nets,nete
!   do i = 1,np
!     do j = 1,np
!       call formJac(JacL_elem(:,i,j,ie),JacD_elem(:,i,j,ie),JacU_elem(:,i,j,ie),1.d0,Jac)
!       Lu_matmul(i,j,1:nlev,ie) = elem(ie)%state%w_i(i,j,1:nlev,n0)
!       Lu_matmul(i,j,1+nlev:2*nlev,ie) = elem(ie)%state%phinh_i(i,j,1:nlev,n0)
!       Lu_matmul(i,j,:,ie) = matmul(Jac,Lu_matmul(i,j,:,ie))
!     end do
!   end do
! end do

  do ie = nets,nete
    do i = 1,np
      do j = 1,np
        ! [0 g*tri \\ g*I 0][w \\ phi] = g*[tri*phi \\ w]
        Lu(i,j,1+nlev:2*nlev,ie) = elem(ie)%state%w_i(i,j,1:nlev,n0)
        ! first component
        Lu(i,j,1,ie) = JacD_elem(1,i,j,ie)*elem(ie)%state%phinh_i(i,j,1,n0)&
          + JacU_elem(1,i,j,ie)*elem(ie)%state%phinh_i(i,j,2,n0)
        ! last component 
        Lu(i,j,nlev,ie) = JacL_elem(nlev-1,i,j,ie)*elem(ie)%state%phinh_i(i,j,nlev-1,n0)&
          +JacD_elem(nlev,i,j,ie)*elem(ie)%state%phinh_i(i,j,nlev,n0)
        ! inner components
        Lu(i,j,2:nlev-1,ie) = JacL_elem(1:nlev-2,i,j,ie)*elem(ie)%state%phinh_i(i,j,1:nlev-2,n0)&
           + JacD_elem(2:nlev-1,i,j,ie)*elem(ie)%state%phinh_i(i,j,2:nlev-1,n0)&
           + JacU_elem(2:nlev-1,i,j,ie)*elem(ie)%state%phinh_i(i,j,3:nlev,n0)
          Lu(i,j,:,ie) = g*Lu(i,j,:,ie)
      end do
    end do
  end do

! print *, "Difference between matmul and entry manipulation:", norm2(Lu_matmul - Lu)
  end subroutine getLu
!===================================================================
! The subroutine add_Lu updates the state variables  
!   by adding Lu to w and phi
! 
!=====================================================================
  subroutine add_Lu(elem,n0,Lu,nets,nete)
  type(element_t), intent(inout)   :: elem(:)
  real(kind=real_kind), intent(in) :: Lu(np,np,2*nlev,nete-nets+1)
  integer, intent(in)              :: nets,nete,n0
  ! local variables
  integer :: ie
 
  ! update w and phi with the action of phi_func
  do ie = nets,nete 
    elem(ie)%state%w_i(:,:,1:nlev,n0)     = elem(ie)%state%w_i(:,:,1:nlev,n0) + Lu(:,:,1:nlev,ie)
    elem(ie)%state%phinh_i(:,:,1:nlev,n0) = elem(ie)%state%phinh_i(:,:,1:nlev,n0) + Lu(:,:,nlev+1:2*nlev,ie)
  end do
  end subroutine add_Lu
!=================================================================================================

!=============================================================================
! formJac forms the actual matrix of the Jacobian
! Inputs: vectors for the lower, main, and upper diagonals; dt
! Outputs: Jac
!=============================================================================
  subroutine formJac(JacL,JacD,JacU,dt,Jac)
  real(kind=real_kind), intent(in)  :: JacL(nlev-1),JacD(nlev),JacU(nlev-1),dt
  real(kind=real_kind), intent(out) :: Jac(2*nlev,2*nlev)
  ! local variables
  integer :: k,dimJac

  ! Form Jacobian
  dimJac = 2*nlev
  Jac = 0.d0
  do k = 1,(nlev-1)
    Jac(k,(k+nlev))     = JacD(k)
    Jac(k,(k+1+nlev))   = JacU(k)
    Jac((k+1),(k+nlev)) = JacL(k)
  end do
  Jac(nlev, dimJac)     = JacD(nlev)
  do k = 1,nlev
    Jac(k + nlev,k) = 1.d0
  end do
  Jac = Jac * g * dt
  end subroutine formJac
!==============================================================================
  subroutine tri_mult(JacL,JacD,JacU,mat,prod,nrhs)
  real(kind=real_kind),intent(in) :: JacL(nlev-1),JacD(nlev),JacU(nlev-1)
  complex(kind=8),intent(in) :: mat(nlev,nrhs)
  complex(kind=8),intent(out) :: prod(nlev,nrhs)
  integer :: nrhs

 ! local variables
  integer :: j
  do j=1,nrhs ! loop over columns
    prod(1,j)        = JacD(1)*mat(1,j) + JacU(1)*mat(2,j)
    prod(2:nlev-1,j) = JacD(2:nlev-1)*mat(2:nlev-1,j)&
                     + JacL(1:nlev-2)*mat(1:nlev-2,j) + JacU(2:nlev-1)*mat(3:nlev,j)
    prod(nlev,j)     = JacL(nlev-1)*mat(nlev-1,j)+JacD(nlev)*mat(nlev,j)
  end do

  end subroutine

!==============================================================================
  subroutine tri_inv(JacL,JacD,JacU,inv)
  real(kind=real_kind),intent(in)  :: JacL(nlev-1),JacD(nlev),JacU(nlev-1)
  real(kind=real_kind),intent(out) :: inv(nlev,nlev)
  ! local variables
  integer :: i, j
  real(kind=real_kind) :: theta(0:nlev),phi(nlev+1),cprod,bprod

  ! initial conditions
  theta(0)    = 1.d0
  theta(1)    = JacD(1)
  phi(nlev+1) = 1.d0
  phi(nlev)   = JacD(nlev)
  ! recurrence relation
  ! theta(i) = JacD(i)theta(i-1)-JacU(i-1)JacL(i-1)*theta(i-2), i = 2...n
  ! phi(i)   = JacD(i)phi(i+1)-JacU(i)JacL(i)phi(i+2) for i = n-1 ... 1
  do i = 2,nlev
    theta(i) = JacD(i)*theta(i-1)-JacU(i-1)*JacL(i-1)*theta(i-2)
    phi(nlev-i+1) = JacD(nlev-i+1)*phi(nlev-i+2)&
                  - JacU(nlev-i+1)*JacL(nlev-i+1)*phi(nlev-i+3)
  end do

  ! form inverse:
  !if i<j: inv(i,j)=(-1)^(i+j)JacU(i)...JacU(j-1)theta(i-1)phi(j+1)/theta(nlev)
  !if i=j: inv(i,j)=theta(i-1)*phi(j+1)/theta(nlev)
  !if i>j: inv(i,j)=(-1)^(i+j)JacL(j)...JacL(i-1)theta(j-1)phi(i+1)/theta(nlev)

  ! first row:
  inv(1,1) = phi(2)
  bprod = -JacU(1)
  do j = 2,nlev-1
    inv(1,j) = bprod*phi(j+1)
    bprod = -bprod*JacU(j)
  end do
  inv(1,nlev) = bprod*phi(nlev+1)
  ! interior rows:
  do i = 2,nlev-1
    bprod = -JacU(i)
    cprod = (-1.d0)**(i+1)*product(JacL(1:i-1))
    do j = 1,i-1
      inv(i,j) = cprod*theta(j-1)*phi(i+1)
      cprod = -cprod/JacL(j)
    end do
    inv(i,i) = theta(i-1)*phi(j+1)
    do j = i+1,nlev-1
      inv(i,j) = bprod*theta(i-1)*phi(j+1)
      bprod = -bprod * JacU(j)
    end do
    inv(i,nlev) = bprod*theta(i-1)*phi(nlev+1)
  end do
  ! last row:
  cprod = -product(JacL(1:nlev-1))
  do j = 1,nlev-1
    inv(nlev,j) = cprod*theta(j-1)
    cprod = -cprod/JacL(j)
  end do
  inv(nlev,nlev) = theta(nlev-1)

  inv = inv/theta(nlev)

  end subroutine tri_inv
!==============================================================================

  subroutine get_exp_jacobian(JacL,JacD,JacU,dp3d,phi_i,pnh,exact,&
       epsie,hvcoord,dpnh_dp_i,vtheta_dp)
  !================================================================================             
  ! compute Jacobian of F(phi) = phi +const + (dt*g)^2 *(1-dp/dpi) column wise                                
  ! with respect to phi_i                                                                                 
  !                                                                                                         
  ! This subroutine forms the tridiagonal analytic Jacobian (we actually form the diagonal, sub-, and super-diagonal)     
  ! J for use in a LApack tridiagonal LU factorization and solver to solve  J * x = -f either exactly or  
  ! approximately                                                                           
  !                                                                                                         
  !  input:                                                                                                   
  !  exact jacobian:        phi_i, dp3d, pnh                               
  !  matrix free jacobian:  phi_i, dp3d, vtheta_dp, hvcoord, dpnh_dp_i                                     
  !                                                                                    
  ! epsie == 1 means exact Jacobian, epsie ~= 1 means finite difference approximate jacobian                
  ! exact,epsie,hvcoord,dpnh_dp,vtheta_dp,pnh,exner,exner_i are only needed as inputs   
  ! if epsie ~=1                                                                                      
  !                                                                                                  
  ! The rule-of-thumb optimal epsie  is epsie = norm(elem)*sqrt(macheps)                       
  !===================================================================================                      
    real (kind=real_kind), intent(out) :: JacD(nlev,np,np)
    real (kind=real_kind), intent(out) :: JacL(nlev-1,np,np),JacU(nlev-1,np,np)
    real (kind=real_kind), intent(in)    :: dp3d(np,np,nlev), phi_i(np,np,nlevp)
    real (kind=real_kind), intent(inout) :: pnh(np,np,nlev)

    real (kind=real_kind), intent(in), optional :: epsie ! epsie is the differencing size in the approx. Jacobian 
    real (kind=real_kind), intent(in), optional :: dpnh_dp_i(np,np,nlevp)
    real (kind=real_kind), intent(in), optional :: vtheta_dp(np,np,nlev)
    type (hvcoord_t)     , intent(in),  optional    :: hvcoord

    integer, intent(in) :: exact

    ! local                                                                                 
    real (kind=real_kind) :: alpha1(np,np),alpha2(np,np)
    real (kind=real_kind) :: e(np,np,nlev),phi_i_temp(np,np,nlevp),exner(np,np,nlev)
    real (kind=real_kind) :: dpnh2(np,np,nlev),dpnh_dp_i_epsie(np,np,nlevp)
    real (kind=real_kind) :: dp3d_i(np,np,nlevp)
    !                                                                                              
    integer :: k,l
    if (exact.eq.1) then ! use exact Jacobian                                                           
       dp3d_i(:,:,1) = dp3d(:,:,1)
       dp3d_i(:,:,nlevp) = dp3d(:,:,nlev)
       do k=2,nlev
          dp3d_i(:,:,k)=(dp3d(:,:,k)+dp3d(:,:,k-1))/2
       end do
      do k=1,nlev
        ! this code will need to change when the equation of state is changed.                     
        ! add special cases for k==1 and k==nlev+1                                                        
        if (k==1) then
           JacL(k,:,:) = pnh(:,:,k)/&
             ((phi_i(:,:,k)-phi_i(:,:,k+1))*(1-kappa)*dp3d_i(:,:,k+1))
           JacU(k,:,:) = 2 * pnh(:,:,k)/&
             ((phi_i(:,:,k)-phi_i(:,:,k+1))*(1-kappa)*dp3d_i(:,:,k))
           JacD(k,:,:) = - 2 *pnh(:,:,k)/&
             ((phi_i(:,:,k)-phi_i(:,:,k+1))*(1-kappa)*dp3d_i(:,:,k))
        else if (k.eq.nlev) then
           JacD(k,:,:) = -(pnh(:,:,k)/((phi_i(:,:,k)-phi_i(:,:,k+1))*(1-kappa)) &
             +pnh(:,:,k-1)/( (phi_i(:,:,k-1)-phi_i(:,:,k))*(1-kappa)))/dp3d_i(:,:,k)
        else ! k =2,...,nlev-1                           
           JacL(k,:,:) = pnh(:,:,k)/&
             ((phi_i(:,:,k)-phi_i(:,:,k+1))*(1-kappa)*dp3d_i(:,:,k+1))
           JacU(k,:,:) =  pnh(:,:,k)/&
             ((phi_i(:,:,k)-phi_i(:,:,k+1))*(1-kappa)*dp3d_i(:,:,k))
           JacD(k,:,:) = - (pnh(:,:,k)/((phi_i(:,:,k)-phi_i(:,:,k+1))*(1-kappa)) &
             +pnh(:,:,k-1)/( (phi_i(:,:,k-1)-phi_i(:,:,k))*(1-kappa)))/dp3d_i(:,:,k)

        end if
      end do
    else ! use finite difference approximation to Jacobian with differencing size espie           
      ! compute Jacobian of F(phi) = phi +const + (dt*g)^2 *(1-dp/dpi) column wise                        
      ! we only form the tridagonal entries and this code can easily be modified to                      
      ! accomodate sparse non-tridigonal and dense Jacobians, however, testing only                           
      ! the tridiagonal of a Jacobian is probably sufficient for testing purpose                            
      do k=1,nlev
        e=0
        e(:,:,k)=1
        phi_i_temp(:,:,:) = phi_i(:,:,:)
        phi_i_temp(:,:,k) = phi_i(:,:,k) + epsie*e(:,:,k)
        if (theta_hydrostatic_mode) then
          dpnh_dp_i_epsie(:,:,:)=1.d0
        else
          call pnh_and_exner_from_eos(hvcoord,vtheta_dp,dp3d,phi_i_temp,pnh,exner,dpnh_dp_i_epsie,caller='get_exp_jacobian')
        end if
        if (k.eq.1) then
          JacL(k,:,:) = -(dpnh_dp_i(:,:,k+1)-dpnh_dp_i_epsie(:,:,k+1))/epsie
          JacD(k,:,:) = -(dpnh_dp_i(:,:,k)-dpnh_dp_i_epsie(:,:,k))/epsie
        elseif (k.eq.nlev) then
          JacD(k,:,:)   = -(dpnh_dp_i(:,:,k)-dpnh_dp_i_epsie(:,:,k))/epsie
          JacU(k-1,:,:) = -(dpnh_dp_i(:,:,k-1)-dpnh_dp_i_epsie(:,:,k-1))/epsie
        else
          JacL(k,:,:)   = -(dpnh_dp_i(:,:,k+1)-dpnh_dp_i_epsie(:,:,k+1))/epsie
          JacD(k,:,:)   = -(dpnh_dp_i(:,:,k)-dpnh_dp_i_epsie(:,:,k))/epsie
          JacU(k-1,:,:) = -(dpnh_dp_i(:,:,k-1)-dpnh_dp_i_epsie(:,:,k-1))/epsie
        end if
      end do
    end if

  end subroutine get_exp_jacobian

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

!=============================================================================================
  subroutine expLdtwphi(JacL_elem,JacD_elem,JacU_elem,elem,n0,dt,nets,nete)
  real (kind=real_kind), intent(in)  :: JacL_elem(nlev-1,np,np,nete-nets+1),JacU_elem(nlev-1,np,np,nete-nets+1),JacD_elem(nlev,np,np,nete-nets+1)
  type (element_t), intent(inout), target :: elem(:)
  integer, intent(in) :: n0,nets,nete
  real (kind=real_kind), intent(in) :: dt

  ! Local variables
  integer :: ie,i,j
  real (kind=real_kind) :: wphivec(2*nlev),wphivec2(2*nlev)
  real (kind=real_kind) :: expJ(2*nlev,2*nlev),expJ2(2*nlev,2*nlev)


  do ie = nets, nete
    do i = 1,np
      do j = 1,np
        wphivec(1:nlev) = elem(ie)%state%w_i(i,j,1:nlev,n0)
        wphivec(1+nlev:2*nlev) = elem(ie)%state%phinh_i(i,j,1:nlev,n0)
!        call matrix_exponential(JacL_elem(:,i,j,ie),JacD_elem(:,i,j,ie),JacU_elem(:,i,j,ie),nlev,dt,expJ2,wphivec2) ! Pade approximation
        call matrix_exponential_new(JacL_elem(:,i,j,ie),JacD_elem(:,i,j,ie),JacU_elem(:,i,j,ie),dt,expJ,wphivec)
!        if (norm2(expJ-expJ2)<1e-10) then
!          print *, "matrix exponential works!"
!          print *, "Error = ", norm2(expJ-expJ2)
!          stop
!        else
!          print *, "Error: problem with matrix exponential."
!          print *, "Error = ", norm2(expJ-expJ2)
!          stop
!        end if
        elem(ie)%state%w_i(i,j,1:nlev,n0) = wphivec(1:nlev)
        elem(ie)%state%phinh_i(i,j,1:nlev,n0) = wphivec(1+nlev:2*nlev)
      end do
    end do
  end do

  end subroutine expLdtwphi


!------------------------------------------------------------------------------
subroutine phi1Ldt(JacL_elem,JacD_elem,JacU_elem,elem,n0,dt,nets,nete)
  real (kind=real_kind), intent(in)  :: JacL_elem(nlev-1,np,np,nete-nets+1),JacU_elem(nlev-1,np,np,nete-nets+1),JacD_elem(nlev,np,np,nete-nets+1)
  type (element_t), intent(inout), target :: elem(:)
  integer, intent(in) :: n0,nets,nete
  real (kind=real_kind), intent(in) :: dt

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
        call matrix_exponential(JacL_elem(:,i,j,ie),JacD_elem(:,i,j,ie),JacU_elem(:,i,j,ie),nlev,dt,expJ)
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

end module exp_mod

module global_variables
  implicit none
! mathematical parameters
  complex(8),parameter :: zi = (0d0, 1d0)
  real(8),parameter :: pi = 4d0*atan(1d0)

! grids
  integer :: nx
  real(8) :: length_x, dx

! quantum system
  real(8),allocatable :: psi(:)
  complex(8),allocatable :: zpsi(:)

! hamiltonian, potential, etc
  real(8),allocatable :: vpot(:), ham(:,:)
  complex(8),allocatable :: zham(:,:)

! finite difference coefficients
  real(8) :: c0,c1,c2,g1,g2

! time propagation
  real(8) :: Tprop, dt
  integer :: nt
  real(8),allocatable :: act(:)

end module global_variables
!-------------------------------------------------------------------------------
program main
  use global_variables
  implicit none


  call preparation
  call calc_ground_state
  call calc_time_propagation


end program main
!-------------------------------------------------------------------------------
subroutine preparation
  use global_variables
  implicit none
  integer :: ix,jx
  real(8) :: xx, ff

  length_x = 10d0
  nx = 256

  dx = length_x/nx

! finite difference coefficient
  c0 = -5d0/2d0/dx**2
  c1 =  4d0/3d0/dx**2
  c2 = -1d0/12d0/dx**2

  g1 = 2d0/3d0/dx
  g2 = -1d0/12d0/dx

  allocate(psi(nx),zpsi(nx))
  allocate(vpot(nx), ham(nx,nx), zham(nx,nx))


  do ix = 1,nx
    xx = dx*ix
    vpot(ix) = 0.1d0*sin(pi*xx/length_x)**2
!    vpot(ix) = 1d0/(pi)*sin(pi*sin(2d0*pi*xx/length_x))
!    call random_number(ff)
!    vpot(ix) = ff
  end do


  ham = 0d0
  do ix = 1, nx

! +2
    jx = ix + 2 +nx-1
    jx = mod(jx, nx) + 1
    ham(ix,jx) = -0.5d0*c2

! +1
    jx = ix + 1 +nx-1
    jx = mod(jx, nx) + 1
    ham(ix,jx) = -0.5d0*c1

! +/- 0
    jx = ix + 0 +nx-1
    jx = mod(jx, nx) + 1
    ham(ix,jx) = -0.5d0*c0 + vpot(ix)

! -1
    jx = ix - 1 +nx-1
    jx = mod(jx, nx) + 1
    ham(ix,jx) = -0.5d0*c1

! -2
    jx = ix - 2 +nx-1
    jx = mod(jx, nx) + 1
    ham(ix,jx) = -0.5d0*c2

  end do

end subroutine preparation
!-------------------------------------------------------------------------------
subroutine calc_ground_state
  use global_variables
  implicit none
  integer :: ix
  real(8) :: ss
!==LAPACK
    integer :: nmax
    integer :: lwork
    real(8),allocatable :: work_lp(:)
    real(8),allocatable :: rwork(:),w(:)
    integer :: info

    nmax = nx
    lwork = 6*nmax**2
    allocate(work_lp(lwork),rwork(3*nmax-2),w(nmax))
!==LAPACK


  call dsyev('V', 'U', nmax, ham, nmax, w, work_lp, lwork, info)
  psi(:) = ham(:,1)

  ss = sum(psi**2)*dx
  psi = psi/sqrt(ss)

  open(20,file='gs_wfn_pot.out')
  write(20,"(999e26.16e3)")0d0, psi(nx), vpot(nx)
  do ix = 1,nx
    write(20,"(999e26.16e3)")ix*dx, psi(ix), vpot(ix)
  end do
  close(20)

end subroutine calc_ground_state
!-------------------------------------------------------------------------------
subroutine calc_time_propagation
  use global_variables
  implicit none
  real(8):: tt, jt_t
  integer :: it

  zpsi = psi

  Tprop = 2000d0
  dt = 0.01d0
  nt = Tprop/dt
  call set_field

  open(30,file='jt_act.out')
  call calc_current(zpsi,act(0), jt_t)
  write(30,"(999e26.16e3)")0d0, jt_t, act(0),sum(abs(zpsi)**2)*dx


  do it = 0, nt-1

    call dt_evolve(it)

    call calc_current(zpsi,act(it+1), jt_t)
    write(30,"(999e26.16e3)")(it+1)*dt, jt_t, act(it+1),sum(abs(zpsi)**2)*dx

  end do
  close(30)

end subroutine calc_time_propagation
!-------------------------------------------------------------------------------
subroutine set_field
  use global_variables
  implicit none
  integer :: it
  real(8) :: tt

  allocate(act(0:nt))
  act = 0d0

  do it = 0, nt
    tt = dt*it
    if(tt<Tprop)then
      act(it) = 2d0*pi/length_x*sin(pi*tt/tprop)**8
    end if
  end do


end subroutine set_field
!-------------------------------------------------------------------------------
subroutine dt_evolve(it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  real(8) :: act_t
  real(8) :: dt_t

  dt_t = 0.5d0*dt

  act_t = act(it)
  call dt_evolve_state(act_t,dt_t)

  act_t = act(it+1)
  call dt_evolve_state(act_t,dt_t)


end subroutine dt_evolve
!-------------------------------------------------------------------------------
subroutine dt_evolve_state(act_t,dt_t)
  use global_variables
  implicit none
  real(8),intent(in) :: act_t, dt_t

!  call taylor_expansion(act_t, dt_t)
  call lanczos(act_t, dt_t)
!  call simple_krylov(act_t, dt_t) !explicit orthogonalization 

end subroutine dt_evolve_state
!-------------------------------------------------------------------------------
subroutine taylor_expansion(act_t, dt_t)
  use global_variables
  implicit none
  integer,parameter :: n_taylor = 4
  real(8),intent(in) :: act_t, dt_t
  complex(8) :: zpsi_t(nx), zhpsi_t(nx)
  complex(8) :: zfact
  integer :: iexp


  zpsi_t = zpsi
  zfact = 1d0
  do iexp = 1, n_taylor
    zfact = zfact*(-zI*dt_t)/iexp
    call calc_zhpsi(zpsi_t, zhpsi_t, act_t)

    zpsi = zpsi + zfact*zhpsi_t
    zpsi_t = zhpsi_t

  end do


end subroutine taylor_expansion
!-------------------------------------------------------------------------------
subroutine lanczos(act_t, dt_t)
  use global_variables
  implicit none
  integer,parameter :: n_lanczos = 16
  real(8),intent(in) :: act_t, dt_t
  complex(8) :: zvec(nx,n_lanczos)
  complex(8) :: zpsi_t(nx), zhpsi_t(nx)
  real(8) :: alpha_l(n_lanczos),beta_l(2:n_lanczos)
  real(8),allocatable :: ham_l(:,:)
  complex(8),allocatable :: zc(:)
  real(8) :: ss, coeff_ini
  integer :: j
!==LAPACK
  integer :: nmax
  integer :: lwork
  real(8),allocatable :: work_lp(:)
  real(8),allocatable :: rwork(:),w(:)
  integer :: info

!  lwork = 6*nmax**2
!  allocate(work_lp(lwork),rwork(3*nmax-2),w(nmax))
!==LAPACK

  ss = sum(abs(zpsi)**2)*dx
  coeff_ini = sqrt(ss)

  zvec(:,1) = zpsi/sqrt(ss)
  zpsi_t = zvec(:,1)
  call calc_zhpsi(zpsi_t, zhpsi_t, act_t)
  alpha_l(1) = sum(conjg(zhpsi_t)*zvec(:,1))*dx
  zpsi_t = zhpsi_t - alpha_l(1)*zvec(:,1)

  nmax = n_lanczos
  do j = 2, n_lanczos
    ss = sum(abs(zpsi_t)**2)*dx
    beta_l(j) = sqrt(ss)
    if(beta_l(j) == 0d0)then
      write(*,"(A)")"Warning: the Krylov subspace expansion is truncated in the Lanczos method."
      nmax = j-1
      exit
    end if

    zvec(:,j) = zpsi_t/beta_l(j)
    zpsi_t = zvec(:,j)
    call calc_zhpsi(zpsi_t, zhpsi_t, act_t)
    alpha_l(j) = sum(conjg(zhpsi_t)*zvec(:,j))*dx
    zpsi_t(:) = zhpsi_t -alpha_l(j)*zvec(:,j)-beta_l(j)*zvec(:,j-1)

  end do

!  stop ! debug
!  write(*,*)beta_l(2) ! debug

  allocate(ham_l(nmax,nmax), zc(nmax))
  ham_l = 0d0
  do j = 1, nmax
    ham_l(j,j) = alpha_l(j)
  end do
  do j = 2, nmax
    ham_l(j-1,j) = beta_l(j)
    ham_l(j,j-1) = beta_l(j)
  end do

  lwork = 6*nmax**2
  allocate(work_lp(lwork),rwork(3*nmax-2),w(nmax))

  call dsyev('V', 'U', nmax, ham_l, nmax, w, work_lp, lwork, info)

  zc = 0d0; zc(1) = coeff_ini
  zc = matmul(transpose(ham_l),zc)
  do j = 1, nmax
    zc(j) = zc(j)*exp(-zi*dt_t*w(j))
  end do
  zc = matmul(ham_l,zc)

  zpsi = 0d0
  do j = 1, nmax
    zpsi = zpsi + zc(j)*zvec(:,j)
  end do

end subroutine lanczos
!-------------------------------------------------------------------------------
subroutine calc_zhpsi(zpsi_t, zhpsi_t, act_t)
  use global_variables
  implicit none
  complex(8),intent(in) :: zpsi_t(nx)
  complex(8),intent(out) :: zhpsi_t(nx)
  real(8),intent(in) :: act_t
  real(8) :: c0t,c1t,c2t,g1t,g2t
  integer :: ix


  c0t = -0.5d0*c0
  c1t = -0.5d0*c1
  c2t = -0.5d0*c2
  g1t = g1*act_t
  g2t = g2*act_t

  zhpsi_t = vpot*zpsi_t
  
  ix = 1
  zhpsi_t(ix) = zhpsi_t(ix) +c0t*zpsi_t(ix) &
                            +c1t*(zpsi_t(ix+1)+zpsi_t(ix-1+nx)) &
                            +c2t*(zpsi_t(ix+2)+zpsi_t(ix-2+nx)) &
                            -zi*(&
                            g1t*(zpsi_t(ix+1)-zpsi_t(ix-1+nx)) &
                           +g2t*(zpsi_t(ix+2)-zpsi_t(ix-2+nx)))

  ix = 2
  zhpsi_t(ix) = zhpsi_t(ix) +c0t*zpsi_t(ix) &
                            +c1t*(zpsi_t(ix+1)+zpsi_t(ix-1)) &
                            +c2t*(zpsi_t(ix+2)+zpsi_t(ix-2+nx)) &
                            -zi*(&
                            g1t*(zpsi_t(ix+1)-zpsi_t(ix-1)) &
                           +g2t*(zpsi_t(ix+2)-zpsi_t(ix-2+nx)))

  do ix = 1+2, nx-2
    zhpsi_t(ix) = zhpsi_t(ix) +c0t*zpsi_t(ix) &
                              +c1t*(zpsi_t(ix+1)+zpsi_t(ix-1)) &
                              +c2t*(zpsi_t(ix+2)+zpsi_t(ix-2)) &
                              -zi*(&
                              g1t*(zpsi_t(ix+1)-zpsi_t(ix-1)) &
                             +g2t*(zpsi_t(ix+2)-zpsi_t(ix-2)))
  end do
  
  ix = nx-1
  zhpsi_t(ix) = zhpsi_t(ix) +c0t*zpsi_t(ix) &
                            +c1t*(zpsi_t(ix+1)+zpsi_t(ix-1)) &
                            +c2t*(zpsi_t(ix+2-nx)+zpsi_t(ix-2)) &
                            -zi*(&
                            g1t*(zpsi_t(ix+1)-zpsi_t(ix-1)) &
                           +g2t*(zpsi_t(ix+2-nx)-zpsi_t(ix-2)))

  ix = nx
  zhpsi_t(ix) = zhpsi_t(ix) +c0t*zpsi_t(ix) &
                            +c1t*(zpsi_t(ix+1-nx)+zpsi_t(ix-1)) &
                            +c2t*(zpsi_t(ix+2-nx)+zpsi_t(ix-2)) &
                            -zi*(&
                            g1t*(zpsi_t(ix+1-nx)-zpsi_t(ix-1)) &
                           +g2t*(zpsi_t(ix+2-nx)-zpsi_t(ix-2)))


end subroutine calc_zhpsi
!-------------------------------------------------------------------------------
subroutine calc_current(zpsi_t,act_t, jt_t)
  use global_variables
  implicit none
  complex(8),intent(in) :: zpsi_t(nx)
  real(8),intent(in) :: act_t
  real(8),intent(out) :: jt_t
  complex(8) :: zgpsi_t(nx)
  integer :: ix

  zgpsi_t = act_t*zpsi_t

  ix = 1
  zgpsi_t(ix) = zgpsi_t(ix) -zi*(&
                            g1*(zpsi_t(ix+1)-zpsi_t(ix-1+nx)) &
                           +g2*(zpsi_t(ix+2)-zpsi_t(ix-2+nx)))

  ix = 2
  zgpsi_t(ix) = zgpsi_t(ix) -zi*(&
                            g1*(zpsi_t(ix+1)-zpsi_t(ix-1)) &
                           +g2*(zpsi_t(ix+2)-zpsi_t(ix-2+nx)))


  do ix = 1+2, nx-2
    zgpsi_t(ix) = zgpsi_t(ix) -zi*(&
                            g1*(zpsi_t(ix+1)-zpsi_t(ix-1)) &
                           +g2*(zpsi_t(ix+2)-zpsi_t(ix-2)))

  end do
  

  ix = nx-1
  zgpsi_t(ix) = zgpsi_t(ix) -zi*(&
                            g1*(zpsi_t(ix+1)-zpsi_t(ix-1)) &
                           +g2*(zpsi_t(ix+2-nx)-zpsi_t(ix-2)))

  ix = nx
  zgpsi_t(ix) = zgpsi_t(ix) -zi*(&
                            g1*(zpsi_t(ix+1-nx)-zpsi_t(ix-1)) &
                           +g2*(zpsi_t(ix+2-nx)-zpsi_t(ix-2)))


  jt_t = sum(conjg(zpsi_t)*zgpsi_t)*dx


end subroutine calc_current
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

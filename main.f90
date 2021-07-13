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
end module global_variables
!-------------------------------------------------------------------------------
program main
  use global_variables
  implicit none


  call preparation
  call calc_ground_state


end program main
!-------------------------------------------------------------------------------
subroutine preparation
  use global_variables
  implicit none
  integer :: ix,jx
  real(8) :: xx

  length_x = 10d0
  nx = 128

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
    vpot(ix) = 1d0/(pi)*sin(pi*sin(2d0*pi*xx/length_x))
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
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

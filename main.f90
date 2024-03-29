module global_variables
  implicit none
! mathematical parameters
  complex(8),parameter :: zi = (0d0, 1d0)
  real(8),parameter :: pi = 4d0*atan(1d0)

! grids
  integer :: nx

! quantum system
  complex(8),allocatable :: zpsi(:)

! hamiltonian, potential, etc
  complex(8),allocatable :: zham(:,:),zham0(:,:),zvpot(:,:)

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
  integer :: ix,jx, k
  real(8) :: hal2, hal3, hal5, hal7
  real(8) :: xx, ff
  real(8) :: x1,x2
  complex(8),allocatable :: za(:,:), zb(:,:)

  nx = 512


  allocate(zpsi(nx))
  allocate(zham(nx,nx), zham0(nx,nx), zvpot(nx,nx))
  allocate(za(nx,nx), zb(nx,nx))
  za = 0d0
  zb = 0d0

  do jx = 1, nx
    do ix = 1, nx
      k = ix + (jx-1)*nx
      call vdCorput_sequence(k,2,hal2)
      call vdCorput_sequence(k,3,hal3)
      call vdCorput_sequence(k,5,hal5)
      call vdCorput_sequence(k,7,hal7)
      x1 = sqrt(-2d0*log(hal2))*cos(2d0*pi*hal3)
      x2 = sqrt(-2d0*log(hal2))*sin(2d0*pi*hal3)
!      za(ix,jx) = x1 !hal2+zi*hal3 !*exp(zi*2d0*pi*hal3)
!      zb(ix,jx) = zi*x2 !hal5+zi*hal7 !*exp(zi*2d0*pi*hal7)

      za(ix,jx) = 1d0/(1d0+dble(abs(ix-jx))**2)
      if(ix /= jx)zb(ix,jx) = zi/dble(ix-jx)**3
    end do
  end do

!  zham0 = 0.5d0*(za + transpose(conjg(za)) )
!  zvpot = 0.5d0*(zb + transpose(conjg(zb)) )

  zham0 = za
  zvpot = zb

!  write(*,*)"hel",sum(abs(zham0 - transpose(conjg(zham0)))),sum(abs(zvpot - transpose(conjg(zvpot))))

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
  complex(8),allocatable :: za(:,:),work_lp(:)
  real(8),allocatable :: rwork(:),w(:)
  integer :: info

  nmax = nx
  lwork = 6*nmax**2
  allocate(za(nmax,nmax),work_lp(lwork),rwork(3*nmax-2),w(nmax))
!==LAPACK

  za = zham0
  Call zheev('V', 'U', nmax, za, nmax, w, work_lp, lwork, rwork, info)
!  zpsi(:) = za(:,1)
  zpsi(:) = za(:,1)

  ss = sum(abs(zpsi)**2)
  zpsi = zpsi/sqrt(ss)

  open(20,file='gs_wfn_pot.out')
  do ix = 1,nx
    write(20,"(I7,2x,999e26.16e3)")ix, zpsi(ix)
  end do
  close(20)

end subroutine calc_ground_state
!-------------------------------------------------------------------------------
subroutine calc_time_propagation
  use global_variables
  implicit none
  real(8):: tt, jt_t, jt_t0
  integer :: it

  Tprop = 10d0
  dt = 0.5d0
  nt = Tprop/dt
  call set_field

  open(30,file='jt_act.out')
  call calc_current(zpsi, jt_t)
  jt_t0 = jt_t
  write(30,"(999e26.16e3)")0d0, jt_t, act(0),sum(abs(zpsi)**2), jt_t-jt_t0


  do it = 0, nt-1

    call dt_evolve(it)

    call calc_current(zpsi, jt_t)
    write(30,"(999e26.16e3)")(it+1)*dt, jt_t, act(it+1),sum(abs(zpsi)**2), jt_t-jt_t0

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
      act(it) = 1d-8*sin(pi*tt/tprop)**8
!      act(it) = 1d-6*sin(10*pi*tt/tprop)
    end if
  end do

!  act = 1d-8 ! debug

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
  integer,parameter :: n_lanczos = 32
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

  ss = sum(abs(zpsi)**2)
  coeff_ini = sqrt(ss)

  zvec(:,1) = zpsi/sqrt(ss)
  zpsi_t = zvec(:,1)
  call calc_zhpsi(zpsi_t, zhpsi_t, act_t)
  alpha_l(1) = sum(conjg(zhpsi_t)*zvec(:,1))
  zpsi_t = zhpsi_t - alpha_l(1)*zvec(:,1)

  nmax = n_lanczos
  do j = 2, n_lanczos
    ss = sum(abs(zpsi_t)**2)
    beta_l(j) = sqrt(ss)
    if(beta_l(j) == 0d0)then
      write(*,"(A)")"Warning: the Krylov subspace expansion is truncated in the Lanczos method."
      nmax = j-1
      exit
    end if

    zvec(:,j) = zpsi_t/beta_l(j)
    zpsi_t = zvec(:,j)
    call calc_zhpsi(zpsi_t, zhpsi_t, act_t)
    alpha_l(j) = sum(conjg(zhpsi_t)*zvec(:,j))
    zpsi_t(:) = zhpsi_t -alpha_l(j)*zvec(:,j)-beta_l(j)*zvec(:,j-1)

  end do

!  stop ! debug
  write(*,*)beta_l(2) ! debug

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


  zham = zham0 + act_t*zvpot

  zhpsi_t = matmul(zham, zpsi_t)


end subroutine calc_zhpsi
!-------------------------------------------------------------------------------
subroutine calc_current(zpsi_t, jt_t)
  use global_variables
  implicit none
  complex(8),intent(in) :: zpsi_t(nx)
  real(8),intent(out) :: jt_t
  complex(8) :: zgpsi_t(nx)
  integer :: ix

  zgpsi_t = matmul(zvpot,zpsi)
  jt_t = sum(conjg(zpsi_t)*zgpsi_t)


end subroutine calc_current
!-------------------------------------------------------------------------------
!-----------------------------------------
subroutine vdCorput_sequence(n_in,nbase,vdc_res)
  implicit none
  integer,intent(in) :: n_in, nbase
  real(8),intent(out) :: vdc_res
  integer :: n
  real(8) :: bk

  n = n_in
  bk = 1d0/nbase
  vdc_res = 0d0

  do while(n>0)
     vdc_res = vdc_res + mod(n,nbase)*bk
     n = n/nbase
     bk = bk/nbase
     
  end do
  
end subroutine vdCorput_sequence
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

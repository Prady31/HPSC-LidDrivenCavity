!===============================================================================
! solver_omp_n64.f90  —  Lid-Driven Cavity Flow  (OpenMP, N=64)
!
! WEAK SCALING BASELINE  (N=64, run at 1 thread)
!
! The work per cell is identical to solver_omp.f90 (N=128).
! For weak scaling: keep work/thread constant.
!   1 thread  → N=64   (N²=4096 cells / thread)
!   4 threads → N=128  (N²=16384 total / 4 = 4096 cells / thread)
!
! This solver runs Re=1000 only (the standard weak-scaling benchmark case)
! and writes timing to results/timing_weak_n64.csv.
!
! Compile: gfortran -O2 -fopenmp fft_dct.f90 solver_omp_n64.f90 -o solver_n64
! Run:     OMP_NUM_THREADS=1 ./solver_n64
!===============================================================================
program lid_cavity_omp_n64
  use fft_dct_mod
  use omp_lib
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none

  integer,  parameter :: N  = 64
  real(dp), parameter :: h  = 1.0_dp / N
  real(dp), parameter :: h2 = h * h

  real(dp) :: u (0:N,   0:N+1)
  real(dp) :: us(0:N,   0:N+1)
  real(dp) :: v (0:N+1, 0:N  )
  real(dp) :: vs(0:N+1, 0:N  )
  real(dp) :: p (1:N,   1:N  )
  real(dp) :: rhs(1:N,  1:N  )
  real(dp) :: u_old(0:N, 1:N )
  real(dp) :: v_old(1:N, 0:N )
  real(dp) :: wk(0:N-1, 0:N-1)

  real(dp), parameter :: Re      = 1000._dp   ! weak-scaling benchmark case
  real(dp), parameter :: tol     = 1.0e-6_dp
  integer,  parameter :: maxstep = 500000
  integer,  parameter :: rep_freq = 2000

  integer  :: step, nthreads
  real(dp) :: dt, resid, t0, t1
  integer  :: timing_unit, conv_unit

  real(dp) :: lam(0:N-1)
  integer  :: k

  call execute_command_line('mkdir -p results', wait=.true.)

  do k = 0, N-1
    lam(k) = -4.0_dp * sin(acos(-1.0_dp) * k / (2.0_dp*N))**2 / h2
  end do

  !$omp parallel
  !$omp master
  nthreads = omp_get_num_threads()
  !$omp end master
  !$omp end parallel
  write(*,'(A,I3,A)') 'OpenMP N=64: ', nthreads, ' thread(s)'

  u = 0._dp;  v = 0._dp;  p = 0._dp
  dt = min(0.5_dp*h, 0.25_dp*Re*h2, 0.005_dp)
  write(*,'(A,F7.1,A,ES10.3)') '  Re=', Re, '  dt=', dt

  open(newunit=conv_unit, file='results/convergence_weak_n64.csv', status='replace')
  write(conv_unit,'(A)') 'step,residual'

  t0 = omp_get_wtime()

  time_loop: do step = 1, maxstep
    if (mod(step, rep_freq) == 0) then
      u_old = u(0:N, 1:N)
      v_old = v(1:N, 0:N)
    end if

    call apply_bcs()
    call advect_diffuse()
    call build_rhs()
    call poisson_dct()
    call project()

    if (mod(step, rep_freq) == 0) then
      resid = steady_resid()
      write(*,'(A,I7,A,ES10.3)') '  step=', step, '  resid=', resid
      write(conv_unit,'(I7,A,ES14.6)') step, ',', resid
      if (resid < tol) then
        write(*,'(A,I7,A)') '  *** Converged at step ', step, ' ***'
        exit time_loop
      end if
    end if
  end do time_loop

  close(conv_unit)
  t1 = omp_get_wtime()
  write(*,'(A,F8.2,A)') '  Wall time: ', t1-t0, ' s'

  ! Write timing result
  open(newunit=timing_unit, file='results/timing_weak_n64.csv', status='replace')
  write(timing_unit,'(A)') 'N,threads,Re,converged_step,wall_time_s'
  write(timing_unit,'(I4,A,I3,A,F7.1,A,I7,A,F8.2)') &
    N,',',nthreads,',',Re,',',step,',',t1-t0
  close(timing_unit)

  write(*,'(A)') 'Done. Timing in results/timing_weak_n64.csv'

contains

  !============================================================================
  subroutine apply_bcs()
    u(0,   1:N)   = 0._dp
    v(0,   0:N)   = -v(1,   0:N)
    u(N,   1:N)   = 0._dp
    v(N+1, 0:N)   = -v(N,   0:N)
    v(1:N, 0)     = 0._dp
    u(0:N, 0)     = -u(0:N, 1)
    v(1:N, N)     = 0._dp
    u(0:N, N+1)   = 2._dp - u(0:N, N)
  end subroutine apply_bcs

  !============================================================================
  subroutine advect_diffuse()
    integer  :: i, j
    real(dp) :: ue, uw, un, us_v, vn, vs_v, ve, vw, conv, diff, inv_Re_h2

    inv_Re_h2 = 1._dp / (Re * h2)
    us = u
    vs = v

    !$omp parallel do default(none) &
    !$omp shared(u, v, us, dt, inv_Re_h2) &
    !$omp private(i, j, ue, uw, un, us_v, vn, vs_v, conv, diff) schedule(static)
    do j = 1, N
      do i = 1, N-1
        ue   = 0.5_dp*(u(i,j)   + u(i+1,j))
        uw   = 0.5_dp*(u(i-1,j) + u(i,j))
        un   = 0.5_dp*(u(i,j)   + u(i,j+1))
        us_v = 0.5_dp*(u(i,j-1) + u(i,j))
        vn   = 0.5_dp*(v(i,j)   + v(i+1,j))
        vs_v = 0.5_dp*(v(i,j-1) + v(i+1,j-1))
        conv = (ue*ue - uw*uw)/h + (un*vn - us_v*vs_v)/h
        diff = (u(i+1,j)-2._dp*u(i,j)+u(i-1,j) + &
                u(i,j+1)-2._dp*u(i,j)+u(i,j-1)) * inv_Re_h2
        us(i,j) = u(i,j) + dt*(-conv + diff)
      end do
    end do
    !$omp end parallel do

    !$omp parallel do default(none) &
    !$omp shared(u, v, vs, dt, inv_Re_h2) &
    !$omp private(i, j, ve, vw, vn, vs_v, ue, uw, conv, diff) schedule(static)
    do j = 1, N-1
      do i = 1, N
        vn   = 0.5_dp*(v(i,j)   + v(i,j+1))
        vs_v = 0.5_dp*(v(i,j-1) + v(i,j))
        ve   = 0.5_dp*(v(i,j)   + v(i+1,j))
        vw   = 0.5_dp*(v(i-1,j) + v(i,j))
        ue   = 0.5_dp*(u(i,j)   + u(i,j+1))
        uw   = 0.5_dp*(u(i-1,j) + u(i-1,j+1))
        conv = (ue*ve - uw*vw)/h + (vn*vn - vs_v*vs_v)/h
        diff = (v(i+1,j)-2._dp*v(i,j)+v(i-1,j) + &
                v(i,j+1)-2._dp*v(i,j)+v(i,j-1)) * inv_Re_h2
        vs(i,j) = v(i,j) + dt*(-conv + diff)
      end do
    end do
    !$omp end parallel do
  end subroutine advect_diffuse

  !============================================================================
  subroutine build_rhs()
    integer  :: i, j
    real(dp) :: inv_dt_h
    inv_dt_h = 1._dp / (dt * h)
    !$omp parallel do default(none) &
    !$omp shared(us, vs, rhs, inv_dt_h) private(i, j) schedule(static)
    do j = 1, N
      do i = 1, N
        rhs(i,j) = inv_dt_h * ((us(i,j)-us(i-1,j)) + (vs(i,j)-vs(i,j-1)))
      end do
    end do
    !$omp end parallel do
  end subroutine build_rhs

  !============================================================================
  subroutine poisson_dct()
    integer  :: i, j
    real(dp) :: denom
    do j = 1, N
      do i = 1, N
        wk(i-1,j-1) = rhs(i,j)
      end do
    end do
    call dct2_2d(wk, N)
    do j = 0, N-1
      do i = 0, N-1
        denom = lam(i) + lam(j)
        if (i == 0 .and. j == 0) then
          wk(i,j) = 0._dp
        else
          wk(i,j) = wk(i,j) / denom
        end if
      end do
    end do
    call dct3_2d(wk, N)
    do j = 1, N
      do i = 1, N
        p(i,j) = wk(i-1,j-1) / (4.0_dp * N * N)
      end do
    end do
  end subroutine poisson_dct

  !============================================================================
  subroutine project()
    integer  :: i, j
    real(dp) :: dt_h
    dt_h = dt / h
    !$omp parallel do default(none) &
    !$omp shared(u, us, p, dt_h) private(i, j) schedule(static)
    do j = 1, N
      do i = 1, N-1
        u(i,j) = us(i,j) - dt_h*(p(i+1,j) - p(i,j))
      end do
    end do
    !$omp end parallel do
    !$omp parallel do default(none) &
    !$omp shared(v, vs, p, dt_h) private(i, j) schedule(static)
    do j = 1, N-1
      do i = 1, N
        v(i,j) = vs(i,j) - dt_h*(p(i,j+1) - p(i,j))
      end do
    end do
    !$omp end parallel do
  end subroutine project

  !============================================================================
  function steady_resid() result(res)
    real(dp) :: res
    integer  :: i, j
    res = 0._dp
    !$omp parallel do default(none) shared(u, u_old) &
    !$omp private(i,j) reduction(max:res) schedule(static)
    do j = 1, N
      do i = 1, N-1
        res = max(res, abs(u(i,j) - u_old(i,j)))
      end do
    end do
    !$omp end parallel do
    !$omp parallel do default(none) shared(v, v_old) &
    !$omp private(i,j) reduction(max:res) schedule(static)
    do j = 1, N-1
      do i = 1, N
        res = max(res, abs(v(i,j) - v_old(i,j)))
      end do
    end do
    !$omp end parallel do
    res = res / dt
  end function steady_resid

end program lid_cavity_omp_n64

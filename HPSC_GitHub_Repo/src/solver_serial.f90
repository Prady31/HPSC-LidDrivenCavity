!===============================================================================
! solver_serial.f90  —  Lid-Driven Cavity Flow  (SERIAL)
!
! Numerics : 2-D incompressible Navier-Stokes, staggered MAC grid
!            Fractional-step (projection) method
!            2nd-order central differences (advection + diffusion)
!            DCT-based pressure Poisson (pure-Fortran FFT, no libraries)
!
! Domain   : unit square [0,1]², N×N uniform cells
! BCs      : top lid u=1; all other walls no-slip
!
! Output   : results/u_profile_Re<n>.csv      u along x=0.5
!            results/v_profile_Re<n>.csv      v along y=0.5
!            results/ghia_compare_Re<n>.csv   pointwise error vs Ghia (1982)
!            results/convergence_Re<n>.csv    residual history (step, residual)
!            results/field_Re<n>.csv          2D cell-centre: x,y,u,v,p,speed
!            results/vort_Re<n>.csv           2D corner vorticity: x,y,omega
!            results/timing.csv              wall-clock times
!
! Compile  : gfortran -O2 fft_dct.f90 solver_serial.f90 -o solver_serial
!===============================================================================
program lid_cavity_serial
  use fft_dct_mod
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none

  ! ── Grid ───────────────────────────────────────────────────────────────────
  integer,  parameter :: N  = 128
  real(dp), parameter :: h  = 1.0_dp / N
  real(dp), parameter :: h2 = h * h

  ! ── Arrays ─────────────────────────────────────────────────────────────────
  real(dp) :: u (0:N,   0:N+1)
  real(dp) :: us(0:N,   0:N+1)
  real(dp) :: v (0:N+1, 0:N  )
  real(dp) :: vs(0:N+1, 0:N  )
  real(dp) :: p (1:N,   1:N  )
  real(dp) :: rhs(1:N,  1:N  )
  real(dp) :: u_old(0:N, 1:N )
  real(dp) :: v_old(1:N, 0:N )
  real(dp) :: wk(0:N-1, 0:N-1)

  ! ── Run parameters ─────────────────────────────────────────────────────────
  integer,  parameter :: nRe      = 4
  real(dp), parameter :: Re_list(nRe) = [100._dp, 400._dp, 1000._dp, 3200._dp]
  real(dp), parameter :: tol      = 1.0e-6_dp
  integer,  parameter :: maxstep  = 500000
  integer,  parameter :: rep_freq = 2000

  integer  :: iRe, step
  real(dp) :: Re, dt, resid, t0, t1
  integer  :: timing_unit, conv_unit
  character(len=60) :: conv_fname

  ! ── DCT eigenvalues ─────────────────────────────────────────────────────────
  real(dp) :: lam(0:N-1)
  integer  :: k

  call execute_command_line('mkdir -p results plots', wait=.true.)

  do k = 0, N-1
    lam(k) = -4.0_dp * sin(pi_val() * k / (2.0_dp*N))**2 / h2
  end do

  open(newunit=timing_unit, file='results/timing.csv', status='replace')
  write(timing_unit,'(A)') 'Re,converged_step,wall_time_s'

  ! ══════════════════════════════════════════════════════════════════════════
  do iRe = 1, nRe
    Re = Re_list(iRe)
    write(*,'(/,A,F7.1)') '=== Re = ', Re

    u = 0._dp;  v = 0._dp;  p = 0._dp
    dt = min(0.5_dp*h, 0.25_dp*Re*h2, 0.005_dp)
    write(*,'(A,ES10.3)') '  dt = ', dt

    ! Open per-Re convergence history file
    write(conv_fname,'(A,I4,A)') 'results/convergence_Re', nint(Re), '.csv'
    open(newunit=conv_unit, file=trim(conv_fname), status='replace')
    write(conv_unit,'(A)') 'step,residual'

    call cpu_time(t0)

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

    call cpu_time(t1)
    write(*,'(A,F8.2,A)') '  CPU time: ', t1-t0, ' s'
    write(timing_unit,'(F7.1,A,I7,A,F8.2)') Re,',',step,',',t1-t0

    call write_profiles(Re)
    call ghia_compare(Re)
    call write_fields(Re)

  end do
  ! ══════════════════════════════════════════════════════════════════════════

  close(timing_unit)
  write(*,'(/,A)') 'Done. Results in results/'

contains

  pure function pi_val() result(pv)
    real(dp) :: pv
    pv = acos(-1.0_dp)
  end function pi_val

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
  end subroutine advect_diffuse

  !============================================================================
  subroutine build_rhs()
    integer  :: i, j
    real(dp) :: inv_dt_h
    inv_dt_h = 1._dp / (dt * h)
    do j = 1, N
      do i = 1, N
        rhs(i,j) = inv_dt_h * ((us(i,j)-us(i-1,j)) + (vs(i,j)-vs(i,j-1)))
      end do
    end do
  end subroutine build_rhs

  !============================================================================
  subroutine poisson_dct()
    integer  :: i, j
    real(dp) :: denom

    do j = 1, N
      do i = 1, N
        wk(i-1, j-1) = rhs(i,j)
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
        p(i,j) = wk(i-1, j-1) / (4.0_dp * N * N)
      end do
    end do
  end subroutine poisson_dct

  !============================================================================
  subroutine project()
    integer  :: i, j
    real(dp) :: dt_h
    dt_h = dt / h
    do j = 1, N
      do i = 1, N-1
        u(i,j) = us(i,j) - dt_h*(p(i+1,j) - p(i,j))
      end do
    end do
    do j = 1, N-1
      do i = 1, N
        v(i,j) = vs(i,j) - dt_h*(p(i,j+1) - p(i,j))
      end do
    end do
  end subroutine project

  !============================================================================
  function steady_resid() result(res)
    real(dp) :: res
    integer  :: i, j
    res = 0._dp
    do j = 1, N
      do i = 1, N-1
        res = max(res, abs(u(i,j) - u_old(i,j)))
      end do
    end do
    do j = 1, N-1
      do i = 1, N
        res = max(res, abs(v(i,j) - v_old(i,j)))
      end do
    end do
    res = res / dt
  end function steady_resid

  !============================================================================
  subroutine write_fields(Re_in)
    ! Dumps full 2D field data for post-processing and visualisation.
    !   field_Re<n>.csv  : cell-centre x,y,u,v,p,speed
    !   vort_Re<n>.csv   : cell-corner x,y,vorticity
    !                      ω = ∂v/∂x − ∂u/∂y  (exact on MAC stencil)
    real(dp), intent(in) :: Re_in
    integer  :: i, j, funit
    character(len=60) :: fname
    integer  :: Re_int
    real(dp) :: u_cc, v_cc, speed_val, inv_h
    real(dp) :: vort(0:N, 0:N)

    Re_int = nint(Re_in)
    inv_h  = 1.0_dp / h

    ! ── Corner vorticity ω(i,j) at (i·h, j·h) ───────────────────────────────
    do j = 0, N
      do i = 0, N
        vort(i,j) = (v(i+1,j) - v(i,j))*inv_h - (u(i,j+1) - u(i,j))*inv_h
      end do
    end do

    ! ── Cell-centre field ────────────────────────────────────────────────────
    write(fname,'(A,I4,A)') 'results/field_Re', Re_int, '.csv'
    open(newunit=funit, file=trim(fname), status='replace')
    write(funit,'(A)') 'x,y,u,v,p,speed'
    do j = 1, N
      do i = 1, N
        u_cc      = 0.5_dp*(u(i-1,j) + u(i,j))
        v_cc      = 0.5_dp*(v(i,j-1) + v(i,j))
        speed_val = sqrt(u_cc**2 + v_cc**2)
        write(funit,'(F8.5,5(A,F13.8))') &
          (i-0.5_dp)*h, ',', (j-0.5_dp)*h, ',', &
          u_cc, ',', v_cc, ',', p(i,j), ',', speed_val
      end do
    end do
    close(funit)

    ! ── Corner vorticity ────────────────────────────────────────────────────
    write(fname,'(A,I4,A)') 'results/vort_Re', Re_int, '.csv'
    open(newunit=funit, file=trim(fname), status='replace')
    write(funit,'(A)') 'x,y,vorticity'
    do j = 0, N
      do i = 0, N
        write(funit,'(F8.5,2(A,F13.8))') i*h, ',', j*h, ',', vort(i,j)
      end do
    end do
    close(funit)

    write(*,'(A,I4)') '  2D fields written for Re=', Re_int
  end subroutine write_fields

  !============================================================================
  subroutine write_profiles(Re_in)
    real(dp), intent(in) :: Re_in
    integer  :: i, j, funit
    character(len=60) :: fname
    real(dp) :: x, y
    integer  :: Re_int
    Re_int = nint(Re_in)

    write(fname,'(A,I4,A)') 'results/u_profile_Re', Re_int, '.csv'
    open(newunit=funit, file=trim(fname), status='replace')
    write(funit,'(A)') 'y,u_sim'
    write(funit,'(F10.6,A,F12.8)') 0.0_dp, ',', 0.0_dp
    do j = 1, N
      y = (j - 0.5_dp)*h
      write(funit,'(F10.6,A,F12.8)') y, ',', u(N/2, j)
    end do
    write(funit,'(F10.6,A,F12.8)') 1.0_dp, ',', 1.0_dp
    close(funit)

    write(fname,'(A,I4,A)') 'results/v_profile_Re', Re_int, '.csv'
    open(newunit=funit, file=trim(fname), status='replace')
    write(funit,'(A)') 'x,v_sim'
    do i = 1, N
      x = (i - 0.5_dp)*h
      write(funit,'(F10.6,A,F12.8)') x, ',', v(i, N/2)
    end do
    close(funit)

    write(*,'(A,I4)') '  Profiles written for Re=', Re_int
  end subroutine write_profiles

  !============================================================================
  subroutine ghia_compare(Re_in)
    real(dp), intent(in) :: Re_in
    integer, parameter   :: ng = 17

    real(dp), parameter :: ghia_y(ng) = [ &
      0.0000_dp, 0.0547_dp, 0.0625_dp, 0.0703_dp, 0.1016_dp, &
      0.1719_dp, 0.2813_dp, 0.4531_dp, 0.5000_dp, 0.6172_dp, &
      0.7344_dp, 0.8516_dp, 0.9531_dp, 0.9609_dp, 0.9688_dp, &
      0.9766_dp, 1.0000_dp ]
    real(dp), parameter :: ghia_u100(ng) = [ &
       0.00000_dp,-0.03717_dp,-0.04192_dp,-0.04775_dp,-0.06434_dp, &
      -0.10150_dp,-0.15662_dp,-0.21090_dp,-0.20581_dp,-0.13641_dp, &
       0.00332_dp, 0.23151_dp, 0.68717_dp, 0.73722_dp, 0.78871_dp, &
       0.84123_dp, 1.00000_dp ]
    real(dp), parameter :: ghia_u400(ng) = [ &
       0.00000_dp,-0.08186_dp,-0.09266_dp,-0.10338_dp,-0.14612_dp, &
      -0.24299_dp,-0.32586_dp,-0.17908_dp,-0.11897_dp, 0.02135_dp, &
       0.16256_dp, 0.29093_dp, 0.55892_dp, 0.61756_dp, 0.68439_dp, &
       0.75837_dp, 1.00000_dp ]
    real(dp), parameter :: ghia_u1000(ng) = [ &
       0.00000_dp,-0.18109_dp,-0.20196_dp,-0.22220_dp,-0.29730_dp, &
      -0.38289_dp,-0.27805_dp,-0.10648_dp,-0.06080_dp, 0.05702_dp, &
       0.18719_dp, 0.33304_dp, 0.46547_dp, 0.51117_dp, 0.57492_dp, &
       0.65928_dp, 1.00000_dp ]
    real(dp), parameter :: ghia_u3200(ng) = [ &
       0.00000_dp,-0.32407_dp,-0.35344_dp,-0.37827_dp,-0.41933_dp, &
      -0.34323_dp,-0.24427_dp,-0.08664_dp,-0.03826_dp, 0.05454_dp, &
       0.17527_dp, 0.36970_dp, 0.34714_dp, 0.44872_dp, 0.57262_dp, &
       0.73275_dp, 1.00000_dp ]

    real(dp) :: ghia_u(ng), u_interp, y_g, y_lo, y_hi, alpha, err_max, err
    integer  :: kp, j_lo, j_hi, funit
    character(len=60) :: fname
    integer  :: Re_int

    Re_int = nint(Re_in)
    select case (Re_int)
      case (100);  ghia_u = ghia_u100
      case (400);  ghia_u = ghia_u400
      case (1000); ghia_u = ghia_u1000
      case (3200); ghia_u = ghia_u3200
      case default; re
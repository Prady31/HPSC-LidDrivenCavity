!===============================================================================
! fft_dct.f90  —  Pure-Fortran FFT + 2-D DCT  (no external libraries)
!
! Public routines:
!   dct2_2d(a, n)  —  forward  2-D DCT-II  (in-place, 0-indexed, n×n)
!   dct3_2d(a, n)  —  backward 2-D DCT-III (in-place, 0-indexed, n×n)
!
! Round-trip identity:  dct3_2d(dct2_2d(x)) = 4*n^2 * x
! (same normalization as FFTW REDFT10 / REDFT01 pair)
!
! Restriction: n must be a power of 2  (e.g. 64, 128, 256 …)
!===============================================================================
module fft_dct_mod
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  private
  public :: dct2_2d, dct3_2d

  real(dp), parameter :: pi = acos(-1.0_dp)

contains

  !============================================================================
  !  In-place radix-2 Cooley-Tukey complex FFT
  !  forward = .true.  →  Z(k) = Σ_j x(j) exp(−2πi·j·k/n)
  !  forward = .false. →  Z(k) = Σ_j x(j) exp(+2πi·j·k/n)   (unnormalized)
  !============================================================================
  subroutine fft1d(x, n, forward)
    integer,     intent(in)    :: n
    complex(dp), intent(inout) :: x(0:n-1)
    logical,     intent(in)    :: forward

    integer     :: i, j, m, step
    complex(dp) :: t, wm, w

    ! ── Bit-reversal permutation ─────────────────────────────────────────────
    j = 0
    do i = 1, n-1
      m = n/2
      do while (j >= m)
        j = j - m
        m = m/2
        if (m == 0) exit
      end do
      j = j + m
      if (i < j) then
        t = x(i);  x(i) = x(j);  x(j) = t
      end if
    end do

    ! ── Butterfly stages ─────────────────────────────────────────────────────
    step = 2
    do while (step <= n)
      if (forward) then
        wm = cmplx(cos(-2.0_dp*pi/step), sin(-2.0_dp*pi/step), kind=dp)
      else
        wm = cmplx(cos( 2.0_dp*pi/step), sin( 2.0_dp*pi/step), kind=dp)
      end if
      do i = 0, n-1, step
        w = cmplx(1.0_dp, 0.0_dp, kind=dp)
        do m = 0, step/2 - 1
          t                  = w * x(i+m+step/2)
          x(i+m+step/2)      = x(i+m) - t
          x(i+m)             = x(i+m) + t
          w = w * wm
        end do
      end do
      step = step * 2
    end do
  end subroutine fft1d

  !============================================================================
  !  1-D DCT-II  (in-place, 0-indexed)
  !  Y(k) = 2 · Σ_{j=0}^{N-1} x(j) · cos(π·k·(2j+1)/(2N))
  !
  !  Algorithm: zero-pad x to length 2N, forward FFT, apply twiddle factors.
  !  Proof:  Y(k) = 2·Re[ exp(−iπk/(2N)) · DFT_{2N}(x_pad)(k) ]
  !============================================================================
  subroutine dct2_1d(x, n)
    integer,  intent(in)    :: n
    real(dp), intent(inout) :: x(0:n-1)

    complex(dp) :: z(0:2*n-1)
    real(dp)    :: ang
    integer     :: k

    z(0:n-1)   = cmplx(x, 0.0_dp, kind=dp)
    z(n:2*n-1) = cmplx(0.0_dp, 0.0_dp, kind=dp)

    call fft1d(z, 2*n, .true.)

    do k = 0, n-1
      ang  = -pi * k / (2.0_dp * n)
      x(k) = 2.0_dp * (cos(ang)*real(z(k),dp) - sin(ang)*aimag(z(k)))
    end do
  end subroutine dct2_1d

  !============================================================================
  !  1-D DCT-III  (in-place, 0-indexed)
  !  Z(j) = x(0) + 2 · Σ_{k=1}^{N-1} x(k) · cos(π·k·(2j+1)/(2N))
  !
  !  This is the inverse of DCT-II (up to factor 2N):
  !    dct3(dct2(v)) = 2N · v
  !
  !  Algorithm: build d(k)=2·x(k)·exp(−iπk/(2N)) for k≥1, d(0)=x(0);
  !             zero-pad; forward FFT; take real parts.
  !  Proof:  Z(j) = Re[ DFT_{2N}(d)(j) ]
  !============================================================================
  subroutine dct3_1d(x, n)
    integer,  intent(in)    :: n
    real(dp), intent(inout) :: x(0:n-1)

    complex(dp) :: d(0:2*n-1)
    real(dp)    :: ang
    integer     :: k

    d(0) = cmplx(x(0), 0.0_dp, kind=dp)
    do k = 1, n-1
      ang  = -pi * k / (2.0_dp * n)
      d(k) = 2.0_dp * x(k) * cmplx(cos(ang), sin(ang), kind=dp)
    end do
    d(n:2*n-1) = cmplx(0.0_dp, 0.0_dp, kind=dp)

    call fft1d(d, 2*n, .true.)

    do k = 0, n-1
      x(k) = real(d(k), dp)
    end do
  end subroutine dct3_1d

  !============================================================================
  !  2-D DCT-II  (in-place, 0-indexed n×n array)
  !  Apply 1-D DCT-II along first index, then second index.
  !  OpenMP: outer loops are fully independent — each thread gets a private
  !          row buffer so there are no data races.
  !============================================================================
  subroutine dct2_2d(a, n)
    integer,  intent(in)    :: n
    real(dp), intent(inout) :: a(0:n-1, 0:n-1)

    real(dp) :: row(0:n-1)
    integer  :: i, j

    ! Along first index (contiguous columns — each j is independent)
    !$omp parallel do private(j) schedule(static)
    do j = 0, n-1
      call dct2_1d(a(0:n-1, j), n)
    end do
    !$omp end parallel do

    ! Along second index (non-contiguous → private row buffer per thread)
    !$omp parallel do private(i, row) schedule(static)
    do i = 0, n-1
      row = a(i, 0:n-1)
      call dct2_1d(row, n)
      a(i, 0:n-1) = row
    end do
    !$omp end parallel do
  end subroutine dct2_2d

  !============================================================================
  !  2-D DCT-III  (in-place, 0-indexed n×n array)
  !  Apply 1-D DCT-III along second index first, then first index.
  !  Round-trip: dct3_2d(dct2_2d(x)) = 4·n² · x
  !============================================================================
  subroutine dct3_2d(a, n)
    integer,  intent(in)    :: n
    real(dp), intent(inout) :: a(0:n-1, 0:n-1)

    real(dp) :: row(0:n-1)
    integer  :: i, j

    ! Along second index first (reverse order for correct inverse)
    !$omp parallel do private(i, row) schedule(static)
    do i = 0, n-1
      row = a(i, 0:n-1)
      call dct3_1d(row, n)
      a(i, 0:n-1) = row
    end do
    !$omp end parallel do

    ! Along first index
    !$omp parallel do private(j) schedule(static)
    do j = 0, n-1
      call dct3_1d(a(0:n-1, j), n)
    end do
    !$omp end parallel do
  end subroutine dct3_2d

end module fft_dct_mod

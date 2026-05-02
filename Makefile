# ============================================================
# Makefile for Lid-Driven Cavity Flow Solvers
# Requires: gfortran only (no external libraries)
# ============================================================

FC       = gfortran
FFLAGS   = -O2 -march=native -Wall
OMPFLAGS = -fopenmp

.PHONY: all serial omp n64 clean run_serial run_omp run_weak benchmark

all: serial omp n64

serial: solver_serial
omp:    solver_omp
n64:    solver_n64

solver_serial: fft_dct.f90 solver_serial.f90
	$(FC) $(FFLAGS) fft_dct.f90 solver_serial.f90 -o solver_serial
	@echo "Built: solver_serial (N=128, serial)"

solver_omp: fft_dct.f90 solver_omp.f90
	$(FC) $(FFLAGS) $(OMPFLAGS) fft_dct.f90 solver_omp.f90 -o solver_omp
	@echo "Built: solver_omp (N=128, OpenMP)"

solver_n64: fft_dct.f90 solver_omp_n64.f90
	$(FC) $(FFLAGS) $(OMPFLAGS) fft_dct.f90 solver_omp_n64.f90 -o solver_n64
	@echo "Built: solver_n64 (N=64, OpenMP — weak scaling baseline)"

run_serial: solver_serial
	@mkdir -p results plots results/vorticity_snapshots
	./solver_serial

run_omp: solver_omp
	@mkdir -p results plots results/vorticity_snapshots
	OMP_NUM
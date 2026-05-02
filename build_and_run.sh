#!/bin/bash
# ============================================================
# build_and_run.sh  —  HPSC Lid-Driven Cavity: Full Build & Run
#
# Run inside WSL from the project folder:
#   cd "/mnt/c/Users/hp/OneDrive/Desktop/HPSC Project"
#   bash build_and_run.sh
#
# What this script does:
#   1.  Checks gfortran
#   2.  Installs Python packages (numpy, matplotlib, scipy, Pillow)
#   3.  Builds all three solvers
#   4.  Serial run  (Re=100/400/1000/3200) → profiles, Ghia, 2D fields, convergence
#   5.  OMP run  (1, 2, 4 threads) → same + vorticity snapshots for Re=1000
#   6.  Weak-scaling run  (N=64, 1 thread)
#   7.  Assembles benchmark_summary.csv
#   8.  Generates ALL plots (contours, streamlines, animation GIF, etc.)
# ============================================================
set -e

echo "============================================================"
echo "  HPSC Lid-Driven Cavity Solver — Full Build & Run"
echo "============================================================"

# ── 1. Check gfortran ────────────────────────────────────────
echo ""
echo "Step 1: Checking tools..."
command -v gfortran >/dev/null 2>&1 || {
  echo "ERROR: gfortran not found. Run: sudo apt install gfortran"
  exit 1
}
echo "  gfortran : $(gfortran --version | head -1)"
echo "  CPUs     : $(nproc)"
echo "  Threads  : will run OMP at 1, 2, $(nproc) threads"

# ── 2. Python packages ───────────────────────────────────────
echo ""
echo "Step 2: Checking / installing Python packages..."
for pkg in numpy matplotlib scipy pillow; do
  python3 -c "import ${pkg/pillow/PIL}" 2>/dev/null || \
    pip3 install "$pkg" --quiet --break-system-packages 2>/dev/null || \
    pip3 install "$pkg" --quiet
done
echo "  All Python packages ready."

# ── 3. Build ─────────────────────────────────────────────────
echo ""
echo "Step 3: Building solvers..."
mkdir -p results plots results/vorticity_snapshots
make clean 2>/dev/null || true
make all
echo "  solver_serial  (N=128, serial)        ✓"
echo "  solver_omp     (N=128, OpenMP)        ✓"
echo "  solver_n64     (N=64,  weak-scaling)  ✓"

# ── 4. Serial run ────────────────────────────────────────────
echo ""
echo "Step 4: Running SERIAL solver..."
echo "  (Re=100, 400, 1000, 3200 — ~5-10 min total)"
./solver_serial | tee results/serial_output.txt
echo "  Serial done."

# ── 5. OpenMP runs ───────────────────────────────────────────
echo ""
echo "Step 5: Running OpenMP solver..."
NCPU=$(nproc)
MAX_T=$(( NCPU < 4 ? NCPU : 4 ))

for THREADS in 1 2 $MAX_T; do
  # Avoid duplicate run if nproc < 4
  [ "$THREADS" -le "$NCPU" ] || continue
  # Skip duplicate (e.g. nproc=2, avoid running t=2 twice)
  [ "$THREADS" -ne "$MAX_T" ] || [ "$THREADS" -ne "2" ] || [ "$MAX_T" -ne "2" ] || continue

  echo ""
  echo "  --- OMP_NUM_THREADS=${THREADS} ---"
  OMP_NUM_THREADS=${THREADS} ./solver_omp | tee results/omp_t${THREADS}_output.txt
  cp results/timing_omp.csv results/timing_omp_t${THREADS}.csv
done
echo "  OpenMP runs done."
echo "  Vorticity snapshots written to results/vorticity_snapshots/"

# ── 6. Weak-scaling run ──────────────────────────────────────
echo ""
echo "Step 6: Weak-scaling run (N=64, Re=1000, 1 thread)..."
OMP_NUM_THREADS=1 ./solver_n64 | tee results/weak_n64_output.txt
echo "  Weak-scaling done. Timing: results/timing_weak_n64.csv"

# ── 7. Assemble benchmark_summary.csv ────────────────────────
echo ""
echo "Step 7: Building benchmark_summary.csv..."
python3 - <<'PYEOF'
import csv, os, glob

rows = []

# Serial timing
f = 'results/timing.csv'
if os.path.exists(f):
    with open(f) as fh:
        reader = csv.DictReader(fh)
        for r in reader:
            rows.append({'solver':'serial','threads':'1',
                         'Re':r['Re'].strip(),'wall_ti
"""
plot_results.py
===============
Full visualisation suite for the lid-driven cavity flow solver.

Plots generated
───────────────
  plots/u_profiles.png          u(y) at x=0.5 vs Ghia (1982), all Re
  plots/v_profiles.png          v(x) at y=0.5, all Re
  plots/ghia_errors.png         L-inf error bar chart vs Re
  plots/speedup.png             strong-scaling wall-time and speedup
  plots/weak_scaling.png        weak-scaling efficiency (N=64/1T vs N=128/4T)
  plots/convergence.png         log(residual) vs step, all Re
  plots/contour_Re<n>.png       vorticity contour + streamlines, each Re
  plots/fields_Re<n>.png        2x2: vorticity / speed / pressure / u-velocity
  plots/amdahl_fit.png          Amdahl's law fit to strong-scaling data
  plots/animation_Re1000.gif    vorticity evolution GIF (from snapshot files)

Run: python3 plot_results.py
"""

import os, csv, glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.ticker import LogLocator
from scipy.optimize import curve_fit

os.makedirs('plots', exist_ok=True)

RE_LIST = [100, 400, 1000, 3200]
COLORS  = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
MARKERS = ['o', 's', '^', 'D']

# ── helpers ────────────────────────────────────────────────────────────────────
def read_csv(path):
    if not os.path.exists(path):
        return None
    rows = []
    with open(path, newline='') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#'):
                continue
            rows.append(line)
    if len(rows) < 2:
        return None
    reader = csv.DictReader(rows)
    data = {}
    for row in reader:
        for k, val in row.items():
            k = k.strip()
            data.setdefault(k, [])
            try:
                data[k].append(float(val.strip()))
            except (ValueError, AttributeError):
                data[k].append(np.nan)
    return {k: np.array(v) for k, v in data.items()}

def read_benchmark(path='results/benchmark_summary.csv'):
    if not os.path.exists(path):
        return {}
    result = {}
    with open(path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            solver  = row['solver'].strip()
            threads = int(row['threads'].strip())
            re      = int(float(row['Re'].strip()))
            wt      = float(row['wall_time_s'].strip())
            result.setdefault(solver, {}).setdefault(threads, {})[re] = wt
    return result

def re_fname(base, Re, ext='.csv'):
    """Return path trying Fortran I4 format (e.g. 'Re 100') then plain 'Re100'."""
    for tag in [f'{Re:4d}', f'{Re}']:
        p = f'results/{base}{tag}{ext}'
        if os.path.exists(p):
            return p
    return f'results/{base}{Re}{ext}'

def load_field(Re):
    path = re_fname('field_Re', Re)
    if not os.path.exists(path):
        path = re_fname('field_omp_Re', Re)
    return read_csv(path)

def load_vort(Re):
    path = re_fname('vort_Re', Re)
    if not os.path.exists(path):
        path = re_fname('vort_omp_Re', Re)
    d = read_csv(path)
    if d is None:
        return None
    N = int(round(np.sqrt(len(d['x']))))
    X = d['x'].reshape(N, N)
    Y = d['y'].reshape(N, N)
    W = d['vorticity'].reshape(N, N)
    return X, Y, W

# ==============================================================================
#  Plot 1 — u-velocity profiles at x = 0.5
# ==============================================================================
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes = axes.flatten()

for idx, Re in enumerate(RE_LIST):
    ax = axes[idx]
    sim = read_csv(re_fname('u_profile_Re', Re))
    if sim is None:
        sim = read_csv(re_fname('u_profile_omp_Re', Re))
    if sim is not None:
        ax.plot(sim['u_sim'], sim['y'], color=COLORS[idx], lw=2.0, label='Solver (N=128)')
    ghia = read_csv(re_fname('ghia_compare_Re', Re))
    if ghia is not None:
        ax.plot(ghia['u_ghia'], ghia['y'], marker=MARKERS[idx], markersize=6,
                linestyle='none', color='black', markerfacecolor='none',
                label='Ghia et al. (1982)')
    ax.set_title(f'Re = {Re}', fontsize=13, fontweight='bold')
    ax.set_xlabel('u', fontsize=11);  ax.set_ylabel('y', fontsize=11)
    ax.axvline(0, color='grey', lw=0.6, ls='--')
    ax.set_xlim(-0.6, 1.1);  ax.set_ylim(0, 1)
    ax.grid(True, alpha=0.3)
    if ax.get_legend_handles_labels()[0]:
        ax.legend(fontsize=9)

fig.suptitle('u-velocity along vertical centreline (x = 0.5)', fontsize=14)
fig.tight_layout()
fig.savefig('plots/u_profiles.png', dpi=150, bbox_inches='tight')
plt.close(fig)
print('Saved plots/u_profiles.png')

# ==============================================================================
#  Plot 2 — v-velocity profiles at y = 0.5
# ==============================================================================
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes = axes.flatten()

for idx, Re in enumerate(RE_LIST):
    ax = axes[idx]
    sim = read_csv(re_fname('v_profile_Re', Re))
    if sim is None:
        sim = read_csv(re_fname('v_profile_omp_Re', Re))
    if sim is not None:
        ax.plot(sim['x'], sim['v_sim'], color=COLORS[idx], lw=2.0, label=f'Re = {Re}')
    ax.set_title(f'Re = {Re}', fontsize=13, fontweight='bold')
    ax.set_xlabel('x', fontsize=11);  ax.set_ylabel('v', fontsize=11)
    ax.axhline(0, color='grey', lw=0.6, ls='--')
    ax.set_xlim(0, 1)
    ax.grid(True, alpha=0.3)
    if ax.get_legend_handles_labels()[0]:
        ax.legend(fontsize=9)

fig.suptitle('v-velocity along horizontal centreline (y = 0.5)', fontsize=14)
fig.tight_layout()
fig.savefig('plots/v_profiles.png', dpi=150, bbox_inches='tight')
plt.close(fig)
print('Saved plots/v_profiles.png')

# ==============================================================================
#  Plot 3 — Ghia validation error bars
# ==============================================================================
fig, ax = plt.subplots(figsize=(7, 4))
errors = []
for Re in RE_LIST:
    g = read_csv(re_fname('ghia_compare_Re', Re))
    if g is None:
        g = read_csv(re_fname('ghia_compare_omp_Re', Re))
    errors.append(float(np.nanmax(g['error'])) if g is not None else 0.0)

bars = ax.bar([str(r) for r in RE_LIST], errors, color=COLORS, edgecolor='k', width=0.5)
for bar, err in zip(bars, errors):
    ax.text(bar.get_x() + bar.get_width()/2,
            bar.get_height() + max(errors)*0.01,
            f'{err:.4f}', ha='center', va='bottom', fontsize=10)

ax.set_xlabel('Reynolds number', fontsize=12)
ax.set_ylabel(r'$L_\infty$ error  $|u_\mathrm{sim} - u_\mathrm{Ghia}|$', fontsize=12)
ax.set_title('Validation against Ghia et al. (1982)  —  N=128 grid', fontsize=11)
ax.grid(axis='y', alpha=0.3)
fig.tight_layout()
fig.savefig('plots/ghia_errors.png', dpi=150, bbox_inches='tight')
plt.close(fig)
print('Saved plots/ghia_errors.png')

# ==============================================================================
#  Plot 4 — Strong-scaling speedup
# ==============================================================================
bench = read_benchmark()

if bench:
    thread_counts = sorted(bench.get('omp', {}).keys())
    serial_times  = bench.get('serial', {}).get(1, {})

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    ax = axes[0]
    x_labels = ['Serial'] + [f'OMP\n{t}T' for t in thread_counts]
    for idx, Re in enumerate(RE_LIST):
        t_serial = serial_times.get(Re)
        if t_serial is None:
            continue
        times = [t_serial] + [bench['omp'][t].get(Re, np.nan) for t in thread_counts]
        ax.plot(x_labels, times, marker=MARKERS[idx], color=COLORS[idx],
                lw=1.8, markersize=7, label=f'Re = {Re}')
    ax.set_ylabel('Wall time (s)', fontsize=12)
    ax.set_title('Absolute wall-clock time', fontsize=13)
    ax.legend(fontsize=9);  ax.grid(True, alpha=0.3)

    ax = axes[1]
    for idx, Re in enumerate(RE_LIST):
        t_serial = serial_times.get(Re)
        if t_serial is None:
            continue
        speedups = [t_serial / bench['omp'][t].get(Re, np.nan) for t in thread_counts]
        ax.plot(thread_counts, speedups, marker=MARKERS[idx], color=COLORS[idx],
                lw=1.8, markersize=7, label=f'Re = {Re}')
    ax.plot(thread_counts, thread_counts, 'k--', lw=1.2, label='Ideal (linear)')
    ax.set_xlabel('OpenMP threads', fontsize=12)
    ax.set_ylabel('Speedup  $T_1 / T_p$', fontsize=12)
    ax.set_title('Strong-scaling speedup  (N=128)', fontsize=13)
    ax.legend(fontsize=9);  ax.grid(True, alpha=0.3)
    ax.set_xticks(thread_counts);  ax.set_ylim(bottom=0)

    fig.suptitle('Performance: Serial vs OpenMP  —  N=128 grid', fontsize=14)
    fig.tight_layout()
    fig.savefig('plots/speedup.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print('Saved plots/speedup.png')

# ==============================================================================
#  Plot 5 — Convergence history: log(residual) vs step
# ==============================================================================
fig, ax = plt.subplots(figsize=(9, 5))
found_any = False
for idx, Re in enumerate(RE_LIST):
    d = read_csv(re_fname('convergence_Re', Re))
    if d is None:
        d = read_csv(re_fname('convergence_omp_Re', Re))
    if d is not None:
        mask = d['residual'] > 0
        ax.semilogy(d['step'][mask], d['residual'][mask],
                    color=COLORS[idx], lw=1.8, label=f'Re = {Re}')
        found_any = True

if found_any:
    ax.axhline(1e-6, color='grey', lw=1.0, ls=':', label='Convergence tol (1e-6)')
    ax.set_xlabel('Time step', fontsize=12)
    ax.set_ylabel('Max velocity change / dt  (residual)', fontsize=12)
    ax.set_title('Convergence history — all Reynolds numbers', fontsize=13)
    ax.legend(fontsize=10);  ax.grid(True, alpha=0.3, which='both')
    fig.tight_layout()
    fig.savefig('plots/convergence.png', dpi=150, bbox_inches='tight')
    print('Saved plots/convergence.png')
plt.close(fig)

# ==============================================================================
#  Plot 6 — Amdahl's Law fit to strong-scaling data
# ==============================================================================
if bench:
    serial_times = bench.get('serial', {}).get(1, {})
    thread_counts_all = sorted(bench.get('omp', {}).keys())

    def amdahl(p, f):
        return 1.0 / (f + (1.0 - f) / p)

    fig, ax = plt.subplots(figsize=(7, 5))
    p_fine = np.linspace(1, max(thread_counts_all) + 1, 200)

    for idx, Re in enumerate(RE_LIST):
        t1 = serial_times.get(Re)
        if t1 is None:
            continue
        pts = [(t, t1 / bench['omp'][t].get(Re, np.nan)) for t in thread_counts_all
               if not np.isnan(bench['omp'][t].get(Re, np.nan))]
        if len(pts) < 2:
            continue
        p_pts = np.array([x[0] for x in pts])
        s_pts = np.array([x[1] for x in pts])
        try:
            popt, _ = curve_fit(amdahl, p_pts, s_pts, p0=[0.1], bounds=(0, 1))
            f_serial = popt[0]
            ax.plot(p_pts, s_pts, marker=MARKERS[idx], color=COLORS[idx],
                    markersize=8, linestyle='none', label=f'Re={Re} (f={f_serial:.2%})')
            ax.plot(p_fine, amdahl(p_fine, f_serial), color=COLORS[idx],
                    lw=1.5, ls='--', alpha=0.7)
        except Exception:
            ax.plot(p_pts, s_pts, marker=MARKERS[idx], color=COLORS[idx],
                    markersize=8, linestyle='none', label=f'Re={Re}')

    ax.plot(p_fine, p_fine, 'k--', lw=1.0, label='Ideal linear')
    ax.set_xlabel('Thread count  p', fontsize=12)
    ax.set_ylabel('Speedup  S(p)', fontsize=12)
    ax.set_title("Amdahl's Law fit  —  serial fraction f from data", fontsize=12)
    ax.legend(fontsize=9, title='Measured + Amdahl fit')
    ax.grid(True, alpha=0.3)
    ax.set_xticks(thread_counts_all)
    fig.tight_layout()
    fig.savefig('plots/amdahl_fit.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("Saved plots/amdahl_fit.png")

# ==============================================================================
#  Plot 7 — Weak-scaling efficiency
# ==============================================================================
weak_n64 = read_csv('results/timing_weak_n64.csv')
if weak_n64 is not None and bench:
    t_1t_n64  = float(weak_n64['wall_time_s'][0])
    t_4t_n128 = bench.get('omp', {}).get(4, {}).get(1000, None)

    if t_4t_n128 is not None:
        efficiency = (t_1t_n64 / t_4t_n128) * 100.0
        configs = ['N=64\n1 thread\n(baseline)', 'N=128\n4 threads\n(4× work, 4× cores)']
        times   = [t_1t_n64, t_4t_n128]
        bar_colors = ['#4f9dd9', '#e07b39']

        fig, axes = plt.subplots(1, 2, figsize=(11, 5))

        ax = axes[0]
        bars = ax.bar(configs, times, color=bar_colors, edgecolor='k', width=0.4)
        for bar, t in zip(bars, times):
            ax.text(bar.get_x() + bar.get_width()/2,
                    bar.get_height() + max(times)*0.01,
                    f'{t:.1f}s', ha='center', va='bottom', fontsize=12, fontweight='bold')
        ax.set_ylabel('Wall time (s)', fontsize=12)
        ax.set_title('Weak-scaling wall times\n(Re=1000, tol=1e-6)', fontsize=12)
        ax.grid(axis='y', alpha=0.3)
        ax.set_ylim(0, max(times) * 1.2)

        ax = axes[1]
        ax.barh(['Weak-scaling\nefficiency'], [efficiency],
                color='#2ca02c' if efficiency >= 70 else '#d62728',
                edgecolor='k', height=0.3)
        ax.axvline(100, color='grey', ls='--', lw=1.0, label='Ideal (100%)')
        ax.axvline(70,  color='orange', ls=':', lw=1.0, label='Good threshold (70%)')
        ax.set_xlim(0, 115)
        ax.set_xlabel('Efficiency (%)', fontsize=12)
        ax.set_title(f'Weak-scaling efficiency = {efficiency:.1f}%', fontsize=12)
        ax.legend(fontsize=9)
        ax.text(efficiency + 1.5, 0, f'{efficiency:.1f}%', va='center',
                fontsize=13, fontweight='bold',
                color='#2ca02c' if efficiency >= 70 else '#d62728')

        fig.suptitle(
            'Weak-scaling analysis: work per thread kept constant\n'
            r'($N^2/P = \mathrm{const}$: $64^2/1 = 128^2/4 = 4096$ cells/thread)',
            fontsize=12)
        fig.tight_layout()
        fig.savefig('plots/weak_scaling.png', dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f'Saved plots/weak_scaling.png  (efficiency={efficiency:.1f}%)')

# ==============================================================================
#  Plot 8 — 2D flow field: vorticity + streamlines (one figure per Re)
# ==============================================================================
for Re in RE_LIST:
    vort_data = load_vort(Re)
    field     = load_field(Re)

    if vort_data is None or field is None:
        print(f'  [skip] No 2D field data for Re={Re}')
        continue

    X_c, Y_c, W = vort_data

    N_cc = int(round(np.sqrt(len(field['x']))))
    h_cc = 1.0 / N_cc
    xc = np.linspace(h_cc/2, 1.0 - h_cc/2, N_cc)
    yc = np.linspace(h_cc/2, 1.0 - h_cc/2, N_cc)
    U2 = field['u'].reshape(N_cc, N_cc)
    V2 = field['v'].reshape(N_cc, N_cc)

    fig, ax = plt.subplots(figsize=(7, 6))

    wmax = np.percentile(np.abs(W), 98)
    wmax = max(wmax, 1.0)
    cf = ax.contourf(X_c, Y_c, W, levels=np.linspace(-wmax, wmax, 51),
                     cmap='RdBu_r', extend='both')
    ax.contour(X_c, Y_c, W, levels=[-wmax/2, 0, wmax/2],
               colors='k', linewidths=0.4, linestyles=['-', '-', '-'])
    cb = fig.colorbar(cf, ax=ax, shrink=0.85)
    cb.set_label('Vorticity  ω', fontsize=10)

    ax.streamplot(xc, yc, U2, V2, color='white', linewidth=1.0,
                  density=1.4, arrowsize=0.8, arrowstyle='->')

    ax.set_xlabel('x', fontsize=12);  ax.set_ylabel('y', fontsize=12)
    ax.set_title(f'Vorticity contours + streamlines  —  Re={Re}  (N=128)', fontsize=12)
    ax.set_xlim(0, 1);  ax.set_ylim(0, 1)
    ax.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(f'plots/contour_Re{Re}.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved plots/contour_Re{Re}.png')

# ==============================================================================
#  Plot 9 — 2x2 field panel: vorticity / speed / pressure / u-velocity
# ==============================================================================
for Re in RE_LIST:
    vort_data = load_vort(Re)
    field     = load_field(Re)

    if vort_data is None or field is None:
        continue

    X_c, Y_c, W = vort_data
    N_cc  = int(round(np.sqrt(len(field['x']))))
    h_cc  = 1.0 / N_cc
    xc    = np.linspace(h_cc/2, 1.0 - h_cc/2, N_cc)
    yc    = np.linspace(h_cc/2, 1.0 - h_cc/2, N_cc)
    XX, YY = np.meshgrid(xc, yc)
    U2    = field['u'].reshape(N_cc, N_cc)
    V2    = field['v'].reshape(N_cc, N_cc)
    P2    = field['p'].reshape(N_cc, N_cc)
    Spd2  = field['speed'].reshape(N_cc, N_cc)

    fig, axes = plt.subplots(2, 2, figsize=(12, 11))

    def add_panel(ax, data, X, Y, title, cmap, symm=False):
        if symm:
            vmax = np.percentile(np.abs(data), 98)
            vmax = max(vmax, 1e-10)
            norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
            cf = ax.contourf(X, Y, data, levels=51, cmap=cmap, norm=norm)
        else:
            cf = ax.contourf(X, Y, data, levels=51, cmap=cmap)
        fig.colorbar(cf, ax=ax, shrink=0.85)
        ax.set_title(title, fontsize=11, fontweight='bold')
        ax.set_xlabel('x', fontsize=10);  ax.set_ylabel('y', fontsize=10)
        ax.set_aspect('equal')
        return cf

    add_panel(axes[0,0], W,    X_c,  Y_c,  'Vorticity  ω',  'RdBu_r',  symm=True)
    add_panel(axes[0,1], Spd2, XX,   YY,   'Speed  |U|',    'viridis')
    add_panel(axes[1,0], P2,   XX,   YY,   'Pressure  p',   'RdYlBu_r',symm=True)
    add_panel(axes[1,1], U2,   XX,   YY,   'u-velocity',    'coolwarm', symm=True)

    fig.suptitle(f'2D flow fields  —  Re={Re}  (N=128, converged)', fontsize=14)
    fig.tight_layout()
    fig.savefig(f'plots/fields_Re{Re}.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved plots/fields_Re{Re}.png')

# ==============================================================================
#  Plot 10 — Animation GIF from vorticity snapshots (Re=1000)
# ==============================================================================
snap_dir  = 'results/vorticity_snapshots'
snap_files = sorted(glob.glob(os.path.join(snap_dir, 'snap_*.dat')))

if snap_files:
    try:
        from matplotlib.animation import FuncAnimation, PillowWriter

        frames = []
        steps  = []
        for sf in snap_files:
            try:
                data = np.loadtxt(sf)
                frames.append(data)
                step_n = int(os.path.basename(sf).replace('snap_','').replace('.dat',''))
                steps.append(step_n)
            except Exception:
                pass

        if frames:
            x_anim = np.linspace(0, 1, frames[0].shape[1])
            y_anim = np.linspace(0, 1, frames[0].shape[0])
            X_anim, Y_anim = np.meshgrid(x_anim, y_anim)

            all_vals = np.concatenate([f.ravel() for f in frames])
            wmax_anim = np.percentile(np.abs(all_vals), 98)
            wmax_anim = max(wmax_anim, 1.0)

            fig_a, ax_a = plt.subplots(figsize=(6, 5.5))
            cf_a = ax_a.contourf(X_anim, Y_anim, frames[0],
                                  levels=np.linspace(-wmax_anim, wmax_anim, 40),
                                  cmap='RdBu_r', extend='both')
            fig_a.colorbar(cf_a, ax=ax_a, shrink=0.85).set_label('Vorticity ω', fontsize=9)
            ax_a.set_xlabel('x', fontsize=11);  ax_a.set_ylabel('y', fontsize=11)
            ax_a.set_aspect('equal')
            ax_a.set_title(f'Re=1000  step={steps[0]:6d}', fontsize=12)

            def update(frame_idx):
                ax_a.clear()
                ax_a.contourf(X_anim, Y_anim, frames[frame_idx],
                               levels=np.linspace(-wmax_anim, wmax_anim, 40),
                               cmap='RdBu_r', extend='both')
                ax_a.set_xlabel('x', fontsize=11);  ax_a.set_ylabel('y', fontsize=11)
                ax_a.set_aspect('equal')
                ax_a.set_title(f'Vorticity — Re=1000   step={steps[frame_idx]:6d}', fontsize=12)
                ax_a.set_xlim(0, 1);  ax_a.set_ylim(0, 1)

            ani = FuncAnimation(fig_a, update, frames=len(frames), interval=120)
            writer = PillowWriter(fps=8)
            ani.save('plots/animation_Re1000.gif', writer=writer, dpi=100)
            plt.close(fig_a)
            print(f'Saved plots/animation_Re1000.gif  ({len(frames)} frames)')
        else:
            print('  [skip] No valid snapshot frames loaded')

    except ImportError as e:
        print(f'  [skip] Animation requires Pillow: {e}')
    except Exception as e:
        print(f'  [warn] Animation failed: {e}')
else:
    print('  [skip] No vorticity snapshots found (run solver_omp first)')

print('\nAll plots saved to plots/')

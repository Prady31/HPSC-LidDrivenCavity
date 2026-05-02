"""
plot_paper.py
=============
Paper-ready figure generator for the Lid-Driven Cavity HPSC project.

All figures are sized for an IEEE 2-column paper using latexify:
  - Single-column : latexify(columns=2)  →  3.39 in wide
  - Double-column : latexify(columns=1)  →  6.90 in wide

Colours : Colorbrewer Set1 (4 hues for Re values), print-safe & colour-blind aware
Output  : figs/fig_*.pdf   (PDF vector figures for Overleaf)

Run from the HPSC Project root:
    python3 plot_paper.py
"""

import os, csv, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.optimize import curve_fit
from math import sqrt

# ── output folder ──────────────────────────────────────────────────────────────
OUT = 'figs'
os.makedirs(OUT, exist_ok=True)

# ── Colorbrewer Set1 (4 colours for Re=100,400,1000,3200) ─────────────────────
CB_SET1  = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3']   # red,blue,green,purple
MARKERS  = ['o', 's', '^', 'D']
RE_LIST  = [100, 400, 1000, 3200]
RE_LABEL = [r'$Re=100$', r'$Re=400$', r'$Re=1000$', r'$Re=3200$']

# ── latexify ──────────────────────────────────────────────────────────────────
def latexify(fig_width=None, fig_height=None, columns=2):
    """Configure matplotlib RC for publication-quality figures.
    columns=2 → single-column (3.39 in)
    columns=1 → full-width   (6.90 in)
    """
    assert columns in [1, 2]
    # Use wider figures so labels have more breathing room
    if fig_width is None:
        fig_width = 4.2 if columns == 2 else 7.5
    if fig_height is None:
        golden = (sqrt(5) - 1.0) / 2.0
        fig_height = fig_width * golden
    fig_height = min(fig_height, 9.0)

    params = {
        'text.usetex'        : False,
        'mathtext.fontset'   : 'cm',          # Computer Modern math
        'font.family'        : 'serif',
        'font.serif'         : ['DejaVu Serif', 'Times New Roman', 'Palatino'],
        # ── font sizes (readable at column width) ─────────────────────
        'font.size'          : 11,
        'axes.labelsize'     : 11,
        'axes.titlesize'     : 11,
        'legend.fontsize'    : 10,
        'legend.handlelength': 1.8,
        'legend.borderpad'   : 0.5,
        'xtick.labelsize'    : 10,
        'ytick.labelsize'    : 10,
        # ── figure / lines ────────────────────────────────────────────
        'figure.figsize'     : [fig_width, fig_height],
        'lines.linewidth'    : 1.6,
        'lines.markersize'   : 6,
        'axes.linewidth'     : 0.8,
        'xtick.major.width'  : 0.8,
        'ytick.major.width'  : 0.8,
        'xtick.major.size'   : 4,
        'ytick.major.size'   : 4,
        'xtick.minor.size'   : 2,
        'ytick.minor.size'   : 2,
        'grid.linewidth'     : 0.5,
        'grid.alpha'         : 0.4,
        # ── padding so axis labels never overlap ticks ─────────────────
        'axes.labelpad'      : 6,
        'xtick.major.pad'    : 4,
        'ytick.major.pad'    : 4,
    }
    matplotlib.rcParams.update(params)

def save(fig, name, tight=True):
    path = os.path.join(OUT, name)
    if tight:
        fig.savefig(path, dpi=300, bbox_inches='tight', format='pdf')
    else:
        fig.savefig(path, dpi=300, format='pdf')
    plt.close(fig)
    print(f'  saved {path}')

# ── data helpers ──────────────────────────────────────────────────────────────
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

def re_fname(base, Re, ext='.csv'):
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
    return d['x'].reshape(N, N), d['y'].reshape(N, N), d['vorticity'].reshape(N, N)

def read_timing_files():
    """Return dict: threads -> {Re: wall_time_s}"""
    thread_map = {}
    for t in [1, 2, 4, 8, 16, 24]:
        p = f'results/timing_omp_t{t}.csv'
        d = read_csv(p)
        if d is not None:
            thread_map[t] = {int(round(re)): wt
                             for re, wt in zip(d['Re'], d['wall_time_s'])}
    return thread_map

def read_serial():
    d = read_csv('results/timing.csv')
    if d is None:
        return {}
    return {int(round(re)): wt for re, wt in zip(d['Re'], d['wall_time_s'])}


# ==============================================================================
# Figure 1 — u-velocity profiles (Ghia validation), 2×2 full-width
# ==============================================================================
print('Figure 1: u-velocity profiles...')
latexify(columns=1, fig_height=6.0)
fig, axes = plt.subplots(2, 2, figsize=matplotlib.rcParams['figure.figsize'])

for idx, (Re, lbl) in enumerate(zip(RE_LIST, RE_LABEL)):
    ax = axes[idx // 2][idx % 2]
    col = CB_SET1[idx]

    ghia = read_csv(re_fname('ghia_compare_Re', Re))
    if ghia is None:
        ghia = read_csv(re_fname('ghia_compare_omp_Re', Re))
    sim  = read_csv(re_fname('u_profile_Re', Re))
    if sim is None:
        sim  = read_csv(re_fname('u_profile_omp_Re', Re))

    if ghia is not None:
        ax.scatter(ghia['u_ghia'], ghia['y'], s=28, facecolors='none',
                   edgecolors='k', linewidths=0.7, zorder=5,
                   label=r'Ghia \textit{et al.} (1982)')
    if sim is not None:
        ax.plot(sim['u_sim'], sim['y'], color=col, lw=1.3,
                label=r'Present ($N{=}128$)')

    ax.set_title(lbl)
    ax.set_xlabel(r'$u$')
    ax.set_ylabel(r'$y$')
    ax.axvline(0, color='grey', lw=0.5, ls='--')
    ax.set_xlim(-0.55, 1.05)
    ax.set_ylim(0, 1)
    ax.grid(True)
    ax.legend(loc='upper left')

fig.tight_layout(pad=0.5)
save(fig, 'fig_u_profiles.pdf')

# ==============================================================================
# Figure 2 — v-velocity profiles, 2×2 full-width
# ==============================================================================
print('Figure 2: v-velocity profiles...')
latexify(columns=1, fig_height=6.0)
fig, axes = plt.subplots(2, 2, figsize=matplotlib.rcParams['figure.figsize'])

for idx, (Re, lbl) in enumerate(zip(RE_LIST, RE_LABEL)):
    ax = axes[idx // 2][idx % 2]
    col = CB_SET1[idx]

    sim = read_csv(re_fname('v_profile_Re', Re))
    if sim is None:
        sim = read_csv(re_fname('v_profile_omp_Re', Re))

    if sim is not None:
        ax.plot(sim['x'], sim['v_sim'], color=col, lw=1.3)

    ax.set_title(lbl)
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$v$')
    ax.axhline(0, color='grey', lw=0.5, ls='--')
    ax.set_xlim(0, 1)
    ax.grid(True)

fig.tight_layout(pad=0.5)
save(fig, 'fig_v_profiles.pdf')

# ==============================================================================
# Figure 3 — Ghia L∞ validation error bar chart, single-column
# ==============================================================================
print('Figure 3: Ghia error bar chart...')
latexify(columns=2)
fig, ax = plt.subplots()

errors = []
for Re in RE_LIST:
    g = read_csv(re_fname('ghia_compare_Re', Re))
    if g is None:
        g = read_csv(re_fname('ghia_compare_omp_Re', Re))
    errors.append(float(np.nanmax(g['error'])) if g is not None else 0.0)

x_pos = np.arange(len(RE_LIST))
bars = ax.bar(x_pos, errors, color=CB_SET1, edgecolor='k', linewidth=0.6,
              width=0.55)
for bar, err in zip(bars, errors):
    ax.text(bar.get_x() + bar.get_width()/2,
            bar.get_height() + max(errors)*0.015,
            f'{err:.4f}', ha='center', va='bottom', fontsize=10)

ax.set_xticks(x_pos)
ax.set_xticklabels([r'$Re{=}100$', r'$Re{=}400$', r'$Re{=}1000$', r'$Re{=}3200$'])
ax.set_ylabel(r'$L_\infty$ error $\,|u_\mathrm{sim} - u_\mathrm{Ghia}|$')
ax.set_xlabel('Reynolds number')
ax.grid(axis='y')
fig.tight_layout(pad=0.4)
save(fig, 'fig_ghia_errors.pdf')

# ==============================================================================
# Figure 4 — Convergence history (all Re), single-column
# ==============================================================================
print('Figure 4: Convergence history...')
latexify(columns=2)
fig, ax = plt.subplots()

for idx, (Re, lbl) in enumerate(zip(RE_LIST, RE_LABEL)):
    d = read_csv(re_fname('convergence_Re', Re))
    if d is None:
        d = read_csv(re_fname('convergence_omp_Re', Re))
    if d is not None:
        mask = d['residual'] > 0
        ax.semilogy(d['step'][mask], d['residual'][mask],
                    color=CB_SET1[idx], lw=1.2, label=lbl)

ax.axhline(1e-6, color='dimgrey', lw=0.8, ls=':', label=r'Tol.\ $10^{-6}$')
ax.set_xlabel(r'Time step')
ax.set_ylabel(r'Max.\ residual')
ax.legend(ncol=2)
ax.grid(True, which='both')
fig.tight_layout(pad=0.4)
save(fig, 'fig_convergence.pdf')

# ==============================================================================
# Figure 5 — Strong-scaling: wall-time + speedup, full-width
# ==============================================================================
print('Figure 5: Strong-scaling...')
latexify(columns=1, fig_height=3.2)
fig, axes = plt.subplots(1, 2, figsize=matplotlib.rcParams['figure.figsize'])

serial  = read_serial()
threads = read_timing_files()
tc = sorted(threads.keys())

ax = axes[0]
for idx, (Re, lbl) in enumerate(zip(RE_LIST, RE_LABEL)):
    t_ser = serial.get(Re)
    if t_ser is None:
        continue
    times = [threads[t].get(Re, np.nan) for t in tc]
    ax.plot(tc, times, marker=MARKERS[idx], color=CB_SET1[idx],
            lw=1.2, markersize=4, label=lbl)
ax.set_xlabel('OpenMP threads $p$')
ax.set_ylabel('Wall time (s)')
ax.set_xticks(tc)
ax.legend(ncol=2, fontsize=10)
ax.grid(True)

ax = axes[1]
for idx, (Re, lbl) in enumerate(zip(RE_LIST, RE_LABEL)):
    t_ser = serial.get(Re)
    if t_ser is None:
        continue
    speedups = [t_ser / threads[t].get(Re, np.nan) for t in tc]
    ax.plot(tc, speedups, marker=MARKERS[idx], color=CB_SET1[idx],
            lw=1.2, markersize=4, label=lbl)

p_ideal = np.array(tc, dtype=float)
ax.plot(tc, p_ideal, 'k--', lw=0.9, label='Ideal')
ax.set_xlabel('OpenMP threads $p$')
ax.set_ylabel(r'Speedup $S(p) = T_1/T_p$')
ax.set_xticks(tc)
ax.set_ylim(bottom=0)
ax.legend(ncol=2, fontsize=10)
ax.grid(True)

fig.tight_layout(pad=0.4)
save(fig, 'fig_speedup.pdf')

# ==============================================================================
# Figure 6 — Amdahl's law fit, single-column
# ==============================================================================
print("Figure 6: Amdahl's law fit...")
latexify(columns=2)
fig, ax = plt.subplots()

def amdahl(p, f):
    return 1.0 / (f + (1.0 - f) / p)

p_fine = np.linspace(1, max(tc) + 2, 300)

for idx, (Re, lbl) in enumerate(zip(RE_LIST, RE_LABEL)):
    t_ser = serial.get(Re)
    if t_ser is None:
        continue
    pts = [(t, t_ser / threads[t].get(Re, np.nan)) for t in tc
           if not np.isnan(threads[t].get(Re, np.nan))]
    if len(pts) < 2:
        continue
    p_pts = np.array([x[0] for x in pts])
    s_pts = np.array([x[1] for x in pts])
    try:
        popt, _ = curve_fit(amdahl, p_pts, s_pts, p0=[0.05], bounds=(0, 1))
        f_s = popt[0]
        ax.scatter(p_pts, s_pts, marker=MARKERS[idx], color=CB_SET1[idx],
                   s=28, zorder=5,
                   label=rf'{lbl} ($f={f_s*100:.1f}\%$)')
        ax.plot(p_fine, amdahl(p_fine, f_s), color=CB_SET1[idx], lw=1.0, ls='--')
    except Exception:
        ax.scatter(p_pts, s_pts, marker=MARKERS[idx], color=CB_SET1[idx],
                   s=28, zorder=5, label=lbl)

ax.plot(p_fine, p_fine, 'k:', lw=0.8, label='Ideal')
ax.set_xlabel('Thread count $p$')
ax.set_ylabel(r'Speedup $S(p)$')
ax.legend(fontsize=10)
ax.set_xticks(tc)
ax.grid(True)
fig.tight_layout(pad=0.4)
save(fig, 'fig_amdahl.pdf')

# ==============================================================================
# Figure 7 — Weak-scaling efficiency, single-column
# ==============================================================================
print('Figure 7: Weak-scaling...')
latexify(columns=2)
fig, ax = plt.subplots()

weak_n64 = read_csv('results/timing_weak_n64.csv')
if weak_n64 is not None and threads:
    t_base = float(weak_n64['wall_time_s'][0])
    t_4t   = threads.get(4, {}).get(1000, None)
    if t_4t is not None:
        eff = (t_base / t_4t) * 100.0
        configs = [r'$N{=}64$, 1 thread', r'$N{=}128$, 4 threads']
        times   = [t_base, t_4t]
        bar_cols = ['#377eb8', '#e41a1c']
        bars = ax.bar(configs, times, color=bar_cols, edgecolor='k',
                      linewidth=0.6, width=0.45)
        for bar, t in zip(bars, times):
            ax.text(bar.get_x() + bar.get_width()/2,
                    bar.get_height() + max(times)*0.015,
                    f'{t:.1f}s', ha='center', va='bottom', fontsize=10,
                    fontweight='bold')
        ax.set_ylabel('Wall time (s)')
        ax.set_ylim(0, max(times) * 1.25)
        ax.set_title(f'Weak-scaling efficiency: {eff:.1f}\\%')
        ax.grid(axis='y')
        ax.text(0.98, 0.96,
                f'$\\eta_w = {eff:.1f}\\%$',
                transform=ax.transAxes,
                ha='right', va='top', fontsize=10,
                bbox=dict(boxstyle='round,pad=0.3', fc='wheat', alpha=0.7))

fig.tight_layout(pad=0.4)
save(fig, 'fig_weak_scaling.pdf')

# ==============================================================================
# Figures 8–11 — Vorticity contours + streamlines, one per Re (single-col)
# ==============================================================================
for Re, lbl in zip(RE_LIST, RE_LABEL):
    print(f'Figure: vorticity contour Re={Re}...')
    vd = load_vort(Re)
    fd = load_field(Re)
    if vd is None or fd is None:
        print(f'  [skip] no data for Re={Re}')
        continue

    X_c, Y_c, W = vd
    N_cc = int(round(np.sqrt(len(fd['x']))))
    h    = 1.0 / N_cc
    xc   = np.linspace(h/2, 1.0 - h/2, N_cc)
    yc   = np.linspace(h/2, 1.0 - h/2, N_cc)
    U2   = fd['u'].reshape(N_cc, N_cc)
    V2   = fd['v'].reshape(N_cc, N_cc)

    latexify(columns=2)
    fig, ax = plt.subplots()

    wmax = max(np.percentile(np.abs(W), 98), 1.0)
    cf = ax.contourf(X_c, Y_c, W,
                     levels=np.linspace(-wmax, wmax, 48),
                     cmap='RdBu_r', extend='both')
    ax.contour(X_c, Y_c, W, levels=[-wmax*0.5, 0, wmax*0.5],
               colors='k', linewidths=0.3)
    cb = fig.colorbar(cf, ax=ax, shrink=0.9, pad=0.02)
    cb.set_label(r'Vorticity $\omega$', labelpad=5)
    cb.ax.tick_params(labelsize=9)

    ax.streamplot(xc, yc, U2, V2, color='white', linewidth=0.6,
                  density=1.2, arrowsize=0.6, arrowstyle='->')
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')
    ax.set_title(lbl)
    ax.set_aspect('equal')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    fig.tight_layout(pad=0.3)
    save(fig, f'fig_contour_Re{Re}.pdf')

# ==============================================================================
# Figure 12 — 2×2 flow fields panel for Re=1000, full-width
# ==============================================================================
print('Figure 12: 2×2 field panel Re=1000...')
Re = 1000
vd = load_vort(Re)
fd = load_field(Re)

if vd is not None and fd is not None:
    X_c, Y_c, W = vd
    N_cc = int(round(np.sqrt(len(fd['x']))))
    h    = 1.0 / N_cc
    xc   = np.linspace(h/2, 1.0 - h/2, N_cc)
    yc   = np.linspace(h/2, 1.0 - h/2, N_cc)
    XX, YY = np.meshgrid(xc, yc)
    U2   = fd['u'].reshape(N_cc, N_cc)
    V2   = fd['v'].reshape(N_cc, N_cc)
    P2   = fd['p'].reshape(N_cc, N_cc)
    Spd  = fd['speed'].reshape(N_cc, N_cc)

    latexify(columns=1, fig_height=6.5)
    fig, axes = plt.subplots(2, 2, figsize=matplotlib.rcParams['figure.figsize'])

    panels = [
        (axes[0,0], W,   X_c, Y_c, r'Vorticity $\omega$',  'RdBu_r',   True),
        (axes[0,1], Spd, XX,  YY,  r'Speed $|\mathbf{U}|$','viridis',   False),
        (axes[1,0], P2,  XX,  YY,  'Pressure p',           'RdYlBu_r',  True),
        (axes[1,1], U2,  XX,  YY,  'u-velocity',          'coolwarm', True),    ]
    for ax, data, X, Y, title, cmap, symm in panels:
        if symm:
            vmax = max(np.percentile(np.abs(data), 98), 1e-10)
            norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
            cf   = ax.contourf(X, Y, data, levels=48, cmap=cmap, norm=norm)
        else:
            cf = ax.contourf(X, Y, data, levels=48, cmap=cmap)
        cb = fig.colorbar(cf, ax=ax, shrink=0.88, pad=0.02)
        cb.ax.tick_params(labelsize=9)
        ax.set_title(title)
        ax.set_xlabel('$x$')
        ax.set_ylabel('$y$')
        ax.set_aspect('equal')

    fig.suptitle('Flow fields -- Re=1000, N=128', y=1.005)
    fig.tight_layout(pad=0.4)
    fig.savefig('figs/fig_fields_Re1000.pdf', dpi=300, bbox_inches='tight', format='pdf')
    plt.close(fig)
    print('  saved figs/fig_fields_Re1000.pdf')

print('\nAll figures saved to figs/')

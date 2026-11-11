"""Experiment 2: fitness trajectories of Mendel's Accountant with 70 new mutations per individual, one assumption varied at a time.
Usage: python3 -I scripts/exp2_analysis.py results/exp2_human figures/"""
import sys, os, json, glob, re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
run_dir, fig_dir = sys.argv[1:3]
SURF, INK, INK2, MUTED, GRID, AXIS = '#fcfcfb', '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c3c2b7'
C_OBS, C_MIT, C_NEU = '#2a78d6', '#eb6834', '#1baf7a'
plt.rcParams.update({'font.family':'sans-serif','font.size':9,'axes.edgecolor':AXIS,'axes.linewidth':0.8,'axes.labelcolor':INK2,'xtick.color':MUTED,'ytick.color':MUTED,
    'xtick.labelcolor':INK2,'ytick.labelcolor':INK2,'axes.facecolor':SURF,'figure.facecolor':SURF,'savefig.facecolor':SURF,'axes.grid':True,'grid.color':GRID,'grid.linewidth':0.8,
    'grid.linestyle':'-','axes.axisbelow':True,'legend.frameon':False,'legend.fontsize':8,'axes.titlecolor':INK,'axes.titlesize':10,'axes.titleweight':'bold','axes.spines.top':False,'axes.spines.right':False})
labels = {'base_Udel5': 'base: U_del = 5 of 70, Sanford DFE, probability selection, h² = 0.2, N = 1,000', 'dfe_all_deleterious': 'all 70 mutations deleterious (Sanford DFE)',
          'dfe_empirical_Udel2': 'U_del = 2, high-impact DFE (Gibson\'s Lynch setting)', 'sel_fulltrunc': 'full truncation selection', 'sel_partialtrunc': 'partial truncation selection (0.5)',
          'heritability_1': 'heritability 1.0 (no environmental noise)', 'fav_1e-3': 'beneficial fraction 1e-3 (0.07 per individual per generation)', 'fav_1e-2': 'beneficial fraction 1e-2 (0.7 per individual per generation)',
          'base_N10000': 'base with N = 10,000'}
styles = {'base_Udel5': (C_MIT, '-', 2.5), 'dfe_all_deleterious': (C_MIT, ':', 2), 'dfe_empirical_Udel2': (C_MIT, '--', 2), 'sel_fulltrunc': ('#8c3d14', '-', 1.5), 'sel_partialtrunc': ('#8c3d14', '--', 1.5),
          'heritability_1': ('#c99a2e', '-', 1.5), 'fav_1e-3': (C_OBS, '-', 1.5), 'fav_1e-2': (C_OBS, '--', 1.5), 'base_N10000': ('#6a4c93', '-', 1.5)}
summary = {}
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.4))
for name in labels:
    f = os.path.join(run_dir, name, 'mendel.fit')
    if not os.path.exists(f): continue
    fit = np.loadtxt(f, comments='#', ndmin=2); hst = np.loadtxt(os.path.join(run_dir, name, 'mendel.hst'), comments='#', ndmin=2)
    g, w, dele, fav, pop = fit[:, 0], fit[:, 3], hst[:, 1], hst[:, 3], fit[:, 1]
    z = np.where((w <= 0.01) | (pop < 0.1*pop[0]))[0]; gext = int(g[z[0]]) if z.size else None   # extinction = fitness <= 0.01 or population collapsed below 10% of its size
    def at(x): return float(np.interp(x, g, w)) if g[-1] >= x else None
    slope = (w[-1] - np.interp(max(g[-1]-5000, g[0]), g, w)) / min(5000, g[-1]-g[0]) if g[-1] > g[0] else None
    summary[name] = dict(last_gen=int(g[-1]), fitness_end=float(w[-1]), fitness_1000=at(1000), fitness_10000=at(10000), fitness_25000=at(25000), extinction_gen=gext,
                         del_per_indiv_end=float(dele[-1]), fav_per_indiv_end=float(fav[-1]), del_accum_per_gen=float(dele[-1]/g[-1]), fitness_slope_per_gen_last5000=slope)
    c, ls, lw = styles[name]; ax1.plot(g, w, color=c, ls=ls, lw=lw); ax2.plot(g, dele, color=c, ls=ls, lw=lw)
    print(f"{name:22s} gen {g[-1]:6.0f}: fitness {w[-1]:.3f} (1,000: {at(1000)}, 10,000: {at(10000)}), extinct at {gext}, deleterious/indiv {dele[-1]:,.0f} ({dele[-1]/g[-1]:.2f}/gen), favorable/indiv {fav[-1]:.1f}")
json.dump(summary, open(os.path.join(run_dir, 'exp2_summary.json'), 'w'), indent=1)
for d, lab, y in [(24074, 'Neanderthal split (24,000)', 0.30), (26316, 'gorilla split (26,000)', 0.12)]:
    ax1.axvline(d, color=C_OBS, lw=1, ls=':'); ax1.text(d*1.06, y, lab, fontsize=7.5, color=C_OBS, va='center')
ax1.set_xscale('log'); ax1.set_ylim(0, 1.05); ax1.set_xlim(1, 1.2e5); ax1.set_xlabel('Generation (log)'); ax1.set_ylabel('Mean fitness (initial = 1)')
ax2.set_xscale('log'); ax2.set_yscale('log'); ax2.set_xlabel('Generation (log)'); ax2.set_ylabel('Deleterious mutations per individual (log)')
ax1.set_title("Mendel's Accountant with 70 new mutations per individual per generation: mean fitness", loc='left', fontsize=9); ax2.set_title('Accumulated deleterious mutations per individual', loc='left', fontsize=9)
handles = [Line2D([],[],color=styles[n][0],ls=styles[n][1],lw=styles[n][2],label=labels[n]) for n in labels if n in summary]
fig.legend(handles=handles, loc='lower left', bbox_to_anchor=(0.04, 0.0), ncol=2, fontsize=7.6); fig.subplots_adjust(left=0.06, right=0.98, top=0.9, bottom=0.3, wspace=0.22)
fig.savefig(os.path.join(fig_dir, 'fig_F4_exp2_human_inputs.png'), dpi=170); plt.close(fig); print('figure written')

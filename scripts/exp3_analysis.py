"""Experiment 3 (reduced): LTEE-like inputs in mendel-go (diploid, no crossover; see results/exp3_ltee/DESIGN.md) vs the observed LTEE
fitness trajectory (Wiser, Ribeck & Lenski 2013: mean fitness relative to the ancestor rose ~70% by 50,000 generations, power-law form
w = (1 + b t)^a; the parameters used here are given in README section 8 with their source status).
Usage: python3 -I scripts/exp3_analysis.py results/exp3_ltee figures/"""
import sys, os, json, glob
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
labels = {'mutator_fav1e-4': 'mutator rate 0.041 per genome, beneficial fraction 1e-4 (model default)', 'mutator_fav1e-2': 'mutator rate 0.041, beneficial fraction 1e-2 (generous)',
          'ancestral_fav1e-2': 'ancestral rate 4.1e-4, beneficial fraction 1e-2'}
styles = {'mutator_fav1e-4': (C_MIT, '-'), 'mutator_fav1e-2': (C_MIT, '--'), 'ancestral_fav1e-2': ('#8c3d14', '-')}
summary = {}
fig, ax = plt.subplots(figsize=(10, 6))
for name in labels:
    f = os.path.join(run_dir, name, 'mendel.fit')
    if not os.path.exists(f): continue
    fit = np.loadtxt(f, comments='#', ndmin=2); hst = np.loadtxt(os.path.join(run_dir, name, 'mendel.hst'), comments='#', ndmin=2)
    g, w = fit[:, 0], fit[:, 3]
    fixed = {}
    for b in sorted(glob.glob(os.path.join(run_dir, name, 'allele-bins', '*.json'))):
        d = json.load(open(b)); fixed[d['generation']] = dict(dele=d['deleterious'][-1], fav=d['favorable'][-1])
    summary[name] = dict(last_gen=int(g[-1]), fitness_end=float(w[-1]), fitness_50000=float(np.interp(50000, g, w)) if g[-1] >= 50000 else None,
                         del_per_indiv_end=float(hst[-1, 1]), fav_per_indiv_end=float(hst[-1, 3]), fixed=fixed)
    c, ls = styles[name]; ax.plot(g, w, color=c, ls=ls, lw=2)
    print(f"{name}: gen {g[-1]:.0f}, fitness {w[-1]:.3f}, deleterious/indiv {hst[-1,1]:.1f}, favorable/indiv {hst[-1,3]:.2f}, fixed at end {fixed.get(int(g[-1]))}")
json.dump(summary, open(os.path.join(run_dir, 'exp3_summary.json'), 'w'), indent=1)
# observed: relative fitness ~1.7 at 50,000 generations (Lenski 2017 ISME J, quoting Wiser, Ribeck & Lenski 2013; downloads/sources/ltee_fitness_notes.md)
ax.plot(50000, 1.70, 's', ms=10, color=C_OBS, mec=SURF, mew=1.5, zorder=5); ax.plot([0, 50000], [1.0, 1.70], color=C_OBS, lw=1, ls=':', alpha=0.6)
ax.text(50000*0.97, 1.75, 'observed LTEE: relative fitness ~1.7 after 50,000 generations\n(Wiser, Ribeck & Lenski 2013, as stated in Lenski 2017); rise is concave (power law)', fontsize=7.8, color=INK2, ha='right')
ax.axhline(1.0, color=GRID, lw=1)
ax.set_xlabel('Generation'); ax.set_ylabel('Mean fitness relative to the ancestor'); ax.set_ylim(0, 1.9); ax.set_xlim(0, 60000)
ax.set_title("LTEE-like inputs in Mendel's Accountant (diploid, no recombination: the closest the program allows) vs the observed E. coli fitness", loc='left', fontsize=9.5)
handles = [Line2D([],[],color=styles[n][0],ls=styles[n][1],lw=2,label="Mendel's Accountant: " + labels[n]) for n in labels if n in summary] + \
          [Line2D([],[],marker='s',color=C_OBS,lw=1,ls=':',ms=8,mec=SURF,label='Observed LTEE grand-mean fitness: ~1.7 at 50,000 generations (dotted line only joins start and end)')]
ax.legend(handles=handles, loc='lower left'); fig.tight_layout(); fig.savefig(os.path.join(fig_dir, 'fig_F6_exp3_ltee.png'), dpi=170); plt.close(fig); print('figure written')

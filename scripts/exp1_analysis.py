"""Experiment 1: replication of Sanford et al. 2008 (Proc. 6th ICC) default human run. Published values (downloads/sources/sanford_2008_icc6.txt):
Fig. 1a: 8,730 deleterious and 1,011 favorable mutations per individual at generation 1,000 with selection (9,000 / 1,000 without);
Fig. 1b: without selection fitness 0 at generation 882; Fig. 5: after 10,000 generations 16,116 deleterious and 3,193 beneficial alleles fixed.
Usage: python3 -I scripts/exp1_analysis.py results/exp1_sanford2008 figures/"""
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
def load(name):
    d = os.path.join(run_dir, name)
    if not os.path.exists(os.path.join(d, 'mendel.fit')): return None
    fit = np.loadtxt(os.path.join(d, 'mendel.fit'), comments='#', ndmin=2); hst = np.loadtxt(os.path.join(d, 'mendel.hst'), comments='#', ndmin=2)
    bins = {}
    for f in sorted(glob.glob(os.path.join(d, 'allele-bins', '*.json'))):
        b = json.load(open(f)); bins[b['generation']] = dict(del_fixed=b['deleterious'][-1], fav_fixed=b['favorable'][-1], neu_fixed=b['neutral'][-1])
    return dict(gen=fit[:, 0], fitness=fit[:, 3], mean_mutns=fit[:, 7], popsize=fit[:, 1], dele=hst[:, 1], fav=hst[:, 3], bins=bins)
out = {}
sel, nosel, long = load('fig1_selection'), load('fig1_noselection_h0'), load('fig5_10k')   # rr1 = reproductive_rate 1.0, i.e. truly no selection (a rate of 2 still culls half the offspring)
if sel is not None:
    i = np.searchsorted(sel['gen'], 1000) if sel['gen'][-1] >= 1000 else len(sel['gen'])-1
    out['fig1_selection'] = dict(gen=int(sel['gen'][i]), deleterious_per_indiv=float(sel['dele'][i]), favorable_per_indiv=float(sel['fav'][i]), fitness=float(sel['fitness'][i]),
                                 published=dict(deleterious=8730, favorable=1011))
    print(f"with selection at gen {sel['gen'][i]:.0f}: deleterious/indiv {sel['dele'][i]:.0f} (published 8,730), favorable/indiv {sel['fav'][i]:.0f} (published 1,011), mean fitness {sel['fitness'][i]:.3f}")
if nosel is not None:
    z = np.where(nosel['fitness'] <= 0)[0]; gext = int(nosel['gen'][z[0]]) if z.size else None
    out['fig1_noselection'] = dict(last_gen=int(nosel['gen'][-1]), extinction_gen=gext, fitness_at_end=float(nosel['fitness'][-1]), deleterious_per_indiv_end=float(nosel['dele'][-1]), published_extinction_gen=882,
                                   fitness_at_500=float(np.interp(500, nosel['gen'], nosel['fitness'])))
    print(f"no selection: last generation {nosel['gen'][-1]:.0f}, fitness at end {nosel['fitness'][-1]:.3f}, fitness reached 0 at {gext} (published 882); deleterious/indiv at end {nosel['dele'][-1]:.0f}")
if long is not None:
    g = max(long['bins']) if long['bins'] else None
    if g: out['fig5_10k'] = dict(gen=g, **long['bins'][g], fitness=float(np.interp(g, long['gen'], long['fitness'])), published=dict(del_fixed=16116, fav_fixed=3193))
    print(f"10,000-generation run: last gen {long['gen'][-1]:.0f}, fitness {long['fitness'][-1]:.3f}, fixed at gen {g}: {long['bins'].get(g)} (published 16,116 deleterious, 3,193 favorable)")
    # neutral expectation for the total number of fixations in this run: U_gametic = 5, Ne = 1.33 N from the mating scheme (Experiment 0)
    Qd = json.load(open(os.path.join(os.path.dirname(run_dir.rstrip('/')), 'neutral_fixation_Q.json'))); tauN = np.array(Qd['tau_over_N']); Qv = np.array(Qd['Q'])
    cum = np.concatenate([[0.], np.cumsum(0.5*(Qv[1:]+Qv[:-1])*np.diff(tauN))]); x = g/(1.33*1000); qbar = float(np.interp(x, tauN, cum)/x)
    tot = long['bins'][g]['del_fixed'] + long['bins'][g]['fav_fixed']; exp_neu = 5.0*g*qbar
    out['fig5_10k']['total_fixed'] = tot; out['fig5_10k']['neutral_expectation_Ne1.33N'] = exp_neu; out['fig5_10k']['published_total'] = 16116+3193
    print(f"   total fixed {tot:,} vs neutral expectation U_gametic x G x Qbar = 5 x {g:,} x {qbar:.3f} = {exp_neu:,.0f} (ratio {tot/exp_neu:.2f}); published total 19,309 (ratio {19309/exp_neu:.2f})")
json.dump(out, open(os.path.join(run_dir, 'exp1_summary.json'), 'w'), indent=1)
# figure: fitness and mutation counts vs generation, with the published points
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
for r, lab, ls in [(sel, 'with selection (6 offspring per female)', '-'), (nosel, 'no selection (random survival: heritability 1e-4, 4 offspring per female)', ':'), (long, '10,000-generation run (selection)', '--')]:
    if r is None: continue
    ax1.plot(r['gen'], r['fitness'], color=C_MIT, lw=2, ls=ls); ax2.plot(r['gen'], r['dele'], color=C_MIT, lw=2, ls=ls); ax2.plot(r['gen'], r['fav'], color=C_OBS, lw=1.5, ls=ls)
ax1.axvline(882, color=INK2, lw=1, ls='--'); ax1.text(900, 0.5, 'published: fitness 0 at\ngeneration 882 without selection', fontsize=7.5, color=INK2)
ax2.plot(1000, 8730, 's', ms=8, color=INK2, mfc='none', mew=1.5); ax2.plot(1000, 1011, 's', ms=8, color=INK2, mfc='none', mew=1.5)
ax2.text(1030, 8730, 'published Fig. 1a: 8,730 deleterious', fontsize=7.5, color=INK2, va='center'); ax2.text(1030, 1011, '1,011 favorable', fontsize=7.5, color=INK2, va='center')
ax1.set_xscale('log'); ax2.set_xscale('log'); ax2.set_yscale('log'); ax1.set_ylim(0, 1.05)
ax1.set_xlabel('Generation (log)'); ax1.set_ylabel('Mean fitness'); ax2.set_xlabel('Generation (log)'); ax2.set_ylabel('Mutations per individual (log)')
ax1.set_title("Sanford et al. 2008 default human run, re-run in mendel-go: mean fitness", loc='left', fontsize=9); ax2.set_title('Mutations per individual (orange deleterious, blue favorable)', loc='left', fontsize=9)
handles = [Line2D([],[],color=C_MIT,lw=2,label='with selection, 1,000 generations'), Line2D([],[],color=C_MIT,lw=2,ls=':',label='no selection'), Line2D([],[],color=C_MIT,lw=2,ls='--',label='with selection, 10,000 generations'),
           Line2D([],[],marker='s',color=INK2,lw=0,ms=7,mfc='none',label='published values (Sanford et al. 2008)')]
fig.legend(handles=handles, loc='lower left', bbox_to_anchor=(0.05, 0.0), ncol=4, fontsize=8); fig.subplots_adjust(left=0.07, right=0.98, top=0.9, bottom=0.2, wspace=0.25)
fig.savefig(os.path.join(fig_dir, 'fig_F3_exp1_sanford2008_replication.png'), dpi=170); plt.close(fig); print('figure written')

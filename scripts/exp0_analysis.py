"""Experiment 0: neutral fixations in Mendel's Accountant vs Wright-Fisher theory.
Reads results/exp0_neutral/<run>/allele-bins/*.json (neutral alleles in the 99-100% bin = cumulative fixations), compares with
U_gametic x G x Qbar(G/Ne) where Qbar comes from results/neutral_fixation_Q.json, for Ne = N and for a fitted Ne.
Usage: python3 -I scripts/exp0_analysis.py results/exp0_neutral results/neutral_fixation_Q.json figures/"""
import sys, os, json, glob, re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
run_dir, q_json, fig_dir = sys.argv[1:4]
SURF, INK, INK2, MUTED, GRID, AXIS = '#fcfcfb', '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c3c2b7'
C_OBS, C_MIT, C_NEU = '#2a78d6', '#eb6834', '#1baf7a'
plt.rcParams.update({'font.family':'sans-serif','font.size':9,'axes.edgecolor':AXIS,'axes.linewidth':0.8,'axes.labelcolor':INK2,'xtick.color':MUTED,'ytick.color':MUTED,
    'xtick.labelcolor':INK2,'ytick.labelcolor':INK2,'axes.facecolor':SURF,'figure.facecolor':SURF,'savefig.facecolor':SURF,'axes.grid':True,'grid.color':GRID,'grid.linewidth':0.8,
    'grid.linestyle':'-','axes.axisbelow':True,'legend.frameon':False,'legend.fontsize':8,'axes.titlecolor':INK,'axes.titlesize':10,'axes.titleweight':'bold','axes.spines.top':False,'axes.spines.right':False})
Qd = json.load(open(q_json)); tauN = np.array(Qd['tau_over_N']); Q = np.array(Qd['Q']); Tfix = Qd['mean_fix_time_over_N']
cum = np.concatenate([[0.], np.cumsum(0.5*(Q[1:]+Q[:-1])*np.diff(tauN))])
def qbar(x):
    x = np.asarray(x, dtype=float); xs = np.maximum(x, 1e-9)
    return np.where(x <= tauN[-1], np.interp(x, tauN, cum)/xs, 1.0 - Tfix/xs)
def expected(G, Ne, Ug): return Ug*np.asarray(G)*qbar(np.asarray(G)/Ne)

runs = []
for ini in sorted(glob.glob(os.path.join(run_dir, '*.ini'))):
    txt = open(ini).read()
    get = lambda k, cast=float: cast(re.search(r'^\s*'+k+r'\s*=\s*([^\s#]+)', txt, re.M).group(1).strip('"'))
    N, G_tot, U = get('pop_size', int), get('num_generations', int), get('mutn_rate')
    name = get('case_id', str).replace('exp0_', '')
    gens, fixed, seg = [], [], []
    for f in sorted(glob.glob(os.path.join(run_dir, name, 'allele-bins', '*.json'))):
        d = json.load(open(f)); gens.append(d['generation']); fixed.append(d['neutral'][-1]); seg.append(sum(d['neutral'][:-1]))
    gens, fixed, seg = np.array(gens), np.array(fixed), np.array(seg)
    Ug = U/2
    # fit Ne (ratio Ne/N) by least squares on the trajectory excluding the first 2N generations (where counts are tiny)
    m = gens >= 2*N
    grid = np.linspace(0.8, 2.5, 171); sse = [np.sum((fixed[m] - expected(gens[m], r*N, Ug))**2) for r in grid]
    r_fit = grid[int(np.argmin(sse))]
    hst = np.loadtxt(os.path.join(run_dir, name, 'mendel.hst'), comments='#')
    runs.append(dict(name=name, N=N, G=G_tot, U=U, gens=gens.tolist(), fixed=fixed.tolist(), segregating=seg.tolist(), Ne_over_N_fit=float(r_fit),
                     final_fixed=int(fixed[-1]), expected_Ne_eq_N=float(expected(G_tot, N, Ug)), expected_Ne_fit=float(expected(G_tot, r_fit*N, Ug)),
                     mean_alleles_per_indiv_final=float(hst[-1, 2])))
    print(f"{name}: N={N}, {G_tot} gen ({G_tot/N:.0f} N), U_gametic={Ug}: fixed {fixed[-1]} | WF expectation Ne=N: {expected(G_tot, N, Ug):.0f} (ratio {fixed[-1]/expected(G_tot, N, Ug):.2f}) | "
          f"fitted Ne/N = {r_fit:.2f}: {expected(G_tot, r_fit*N, Ug):.0f} (ratio {fixed[-1]/expected(G_tot, r_fit*N, Ug):.2f}) | segregating alleles at end {seg[-1]} (4Ne U_dip at Ne=N: {4*N*U:.0f}) | alleles per individual {hst[-1,2]:.0f}")
runs.sort(key=lambda r: r['N'])
json.dump(runs, open(os.path.join(run_dir, 'exp0_summary.json'), 'w'), indent=1)

fig, axes = plt.subplots(1, 3, figsize=(13, 5.2), sharey=False)
for ax, r in zip(axes, runs):
    N, Ug = r['N'], r['U']/2; g = np.linspace(1, r['G'], 400)
    ax.plot(g, expected(g, N, Ug), color=C_NEU, lw=2, ls='--'); ax.plot(g, expected(g, r['Ne_over_N_fit']*N, Ug), color=C_NEU, lw=2)
    ax.plot(g, Ug*g, color=INK2, lw=1, ls=':')
    ax.plot(r['gens'], r['fixed'], 'o', ms=6, color=C_MIT, mec=SURF, mew=1.2, zorder=5)
    ax.set_title(f"N = {N:,}, {r['G']:,} generations ({r['G']/N:.0f} N); fitted Ne = {r['Ne_over_N_fit']:.2f} N", loc='left', fontsize=8.5)
    ax.set_xlabel('Generation'); ax.set_ylabel('Cumulative neutral fixations' if ax is axes[0] else '')
    for k in [2, 4, 8]:
        if k*N <= r['G']: ax.axvline(k*N, color=GRID, lw=1); ax.text(k*N, ax.get_ylim()[1]*0.02 if False else 0, f'{k}N', fontsize=7, color=MUTED, ha='center', va='bottom')
handles = [Line2D([],[],marker='o',color=C_MIT,lw=0,ms=7,mec=SURF,label="Mendel's Accountant (mendel-go 1.2.5), neutral alleles at 99–100% frequency, this analysis"),
           Line2D([],[],color=C_NEU,lw=2,ls='--',label='Wright-Fisher: U × G × Q̄(G/N), Ne = N'),
           Line2D([],[],color=C_NEU,lw=2,label='Wright-Fisher with fitted Ne (fixed-offspring-number mating gives Ne > N)'),
           Line2D([],[],color=INK2,lw=1,ls=':',label='Asymptote k = U (Kimura)')]
fig.legend(handles=handles, loc='lower left', bbox_to_anchor=(0.04, 0.0), ncol=2, fontsize=8)
fig.suptitle("Experiment 0: Mendel's Accountant fixes neutral mutations at Kimura's rate after the 4Ne transient (1 neutral mutation per individual, no selection)", x=0.02, ha='left', fontsize=10, fontweight='bold', color=INK)
fig.subplots_adjust(left=0.06, right=0.99, top=0.88, bottom=0.2, wspace=0.22)
fig.savefig(os.path.join(fig_dir, 'fig_F1b_exp0_own_runs.png'), dpi=170); plt.close(fig); print('figure written')

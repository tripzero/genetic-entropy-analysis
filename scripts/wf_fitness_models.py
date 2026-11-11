"""Standard-model comparison for Experiment 2b: a free-recombination (infinite-sites, independent loci) Wright-Fisher population
with the SAME inputs as the Mendel's Accountant runs, under additive vs multiplicative combination of fitness effects.
Selection is the deterministic diploid change dp = -p q s_eff [h q + (1-h) p] followed by binomial sampling (equivalent to
probability selection with no environmental noise); under additive fitness the marginal selection coefficient is s / W-bar
(relative fitness differences grow as mean fitness falls), under multiplicative fitness it is s.
DFEs: 'sanford' = Weibull of Sanford et al. 2007 (genome 3e9, 0.1% of mutations with effect > 0.1); 'gamma' = gamma(shape 0.19,
mean 0.01) [approximation of the human nonsynonymous DFE of Kim, Huber & Lohmueller 2017 / Eyre-Walker & Keightley 2007;
parameters not verified against those texts in this session].
Usage: python3 -I scripts/wf_fitness_models.py results/wf_models.json figures/"""
import sys, os, json, time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
out_json, fig_dir = sys.argv[1:3]
SURF, INK, INK2, MUTED, GRID, AXIS = '#fcfcfb', '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c3c2b7'
C_OBS, C_MIT, C_NEU = '#2a78d6', '#eb6834', '#1baf7a'
plt.rcParams.update({'font.family':'sans-serif','font.size':9,'axes.edgecolor':AXIS,'axes.linewidth':0.8,'axes.labelcolor':INK2,'xtick.color':MUTED,'ytick.color':MUTED,
    'xtick.labelcolor':INK2,'ytick.labelcolor':INK2,'axes.facecolor':SURF,'figure.facecolor':SURF,'savefig.facecolor':SURF,'axes.grid':True,'grid.color':GRID,'grid.linewidth':0.8,
    'grid.linestyle':'-','axes.axisbelow':True,'legend.frameon':False,'legend.fontsize':8,'axes.titlecolor':INK,'axes.titlesize':10,'axes.titleweight':'bold','axes.spines.top':False,'axes.spines.right':False})

def sanford_dfe(rng, n, G=3e9, q=0.001, theta=0.1):
    a = np.log(G); g = np.log(-np.log(theta)/a)/np.log(q)
    return np.exp(-a*rng.random(n)**g)
def gamma_dfe(rng, n, shape=0.19, mean=0.01):
    return np.minimum(rng.gamma(shape, mean/shape, n), 1.0)

def run(N=1000, gens=10000, U_dip=5.0, dfe='sanford', fitness='additive', h=0.5, seed=1, record_every=10):
    rng = np.random.default_rng(seed); twoN = 2*N
    p = np.zeros(0); s = np.zeros(0)
    fixed_add, fixed_logmult, n_fixed = 0.0, 0.0, 0
    hist = []
    draw = sanford_dfe if dfe == 'sanford' else gamma_dfe
    W = 1.0
    for t in range(1, gens+1):
        n_new = rng.poisson(N*U_dip)            # new mutation copies: 2N gene copies x U_dip/2 per gamete
        if n_new:
            p = np.concatenate([p, np.full(n_new, 1.0/twoN)]); s = np.concatenate([s, draw(rng, n_new)])
        q = 1 - p
        seg_add = np.sum(s*(2*h*p*q + p*p)); seg_logmult = np.sum(np.log1p(-np.minimum(s*(2*h*p*q + p*p), 0.999999)))
        W = (1 - fixed_add - seg_add) if fitness == 'additive' else float(np.exp(fixed_logmult + seg_logmult))
        if t % record_every == 0 or t == 1:
            hist.append((t, W, n_fixed, int(p.size)))
        if fitness == 'additive' and W <= 0:
            hist.append((t, 0.0, n_fixed, int(p.size))); break
        s_eff = s/W if fitness == 'additive' else s
        dp = -p*q*s_eff*(h*q + (1-h)*p)
        p2 = np.clip(p + dp, 0, 1)
        x = rng.binomial(twoN, p2); p = x/twoN
        fx = p >= 1.0
        if fx.any():
            n_fixed += int(fx.sum()); fixed_add += float(s[fx].sum()); fixed_logmult += float(np.log1p(-s[fx]).sum())
        keep = (p > 0) & (p < 1); p = p[keep]; s = s[keep]
    return dict(N=N, gens=gens, U_dip=U_dip, dfe=dfe, fitness=fitness, h=h, hist=hist, final_W=W, n_fixed=n_fixed, extinct_gen=(hist[-1][0] if W <= 0 else None))

cases = [dict(U_dip=5.0, dfe='sanford', fitness='additive'), dict(U_dip=5.0, dfe='sanford', fitness='multiplicative'),
         dict(U_dip=2.0, dfe='gamma', fitness='additive'), dict(U_dip=2.0, dfe='gamma', fitness='multiplicative'),
         dict(U_dip=70.0, dfe='sanford', fitness='additive'), dict(U_dip=70.0, dfe='sanford', fitness='multiplicative', gens=3000)]
results = []
if len(sys.argv) > 3 and sys.argv[3] == '--plot-only':
    results = json.load(open(out_json)); cases = []
for i, c in enumerate(cases):
    t0 = time.time(); r = run(seed=100+i, **c); results.append(r)
    print(f"U_del={c['U_dip']:>4} {c['dfe']:8s} {c['fitness']:14s}: final W {r['final_W']:.3f} at gen {r['hist'][-1][0]}, fixed {r['n_fixed']}, segregating {r['hist'][-1][3]}, extinct {r['extinct_gen']}  ({time.time()-t0:.0f}s)", flush=True)
    json.dump(results, open(out_json, 'w'))
fig, ax = plt.subplots(figsize=(10, 6.2))
style = {('sanford', 5.0): (C_MIT, 2.2), ('sanford', 70.0): ('#8c3d14', 1.8), ('gamma', 2.0): (C_NEU, 2.2)}
for r in results:
    h = np.array(r['hist']); c, lw = style[(r['dfe'], r['U_dip'])]
    ax.plot(h[:, 0], h[:, 1], color=c, lw=lw, ls='-' if r['fitness'] == 'additive' else '--')
ax.set_xscale('log'); ax.set_ylim(0, 1.05); ax.set_xlim(1, 10000)
ax.text(3000*1.05, 0.5, 'U_del = 70 multiplicative\nrun to 3,000 generations only', fontsize=7.5, color=INK2)
ax.set_xlabel('Generation (log)'); ax.set_ylabel('Mean fitness relative to a mutation-free genome')
ax.set_title('Same inputs, additive (solid, as in Mendel\'s Accountant) vs multiplicative (dashed) fitness: the DFE sets the decline', loc='left')
handles = [Line2D([],[],color=C_MIT,lw=2.2,label='Sanford DFE, U_del = 5 per individual (Gibson et al. 2013 setting)'), Line2D([],[],color='#8c3d14',lw=1.8,label='Sanford DFE, U_del = 70 (all new mutations deleterious)'),
           Line2D([],[],color=C_NEU,lw=2.2,label='gamma DFE (shape 0.19, mean 0.01), U_del = 2 (empirical-like; parameters approximate)'),
           Line2D([],[],color=INK2,lw=1.5,label='additive combination (Mendel\'s Accountant, mendel-go: the only mode implemented)'), Line2D([],[],color=INK2,lw=1.5,ls='--',label='multiplicative combination (standard population genetics)')]
ax.legend(handles=handles, loc='lower left'); fig.tight_layout(); fig.savefig(os.path.join(fig_dir, 'fig_F5_additive_vs_multiplicative.png'), dpi=170); plt.close(fig); print('figure written')

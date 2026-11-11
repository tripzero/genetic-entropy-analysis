"""Figures: observed fixed SNV differences per lineage for (1) human vs chimpanzee, (2) eastern vs western gorilla,
(3) modern humans vs Neanderthal, against generations since the population split, with the neutral-theory expectation and
the genetic-entropy model (Mendel's Accountant) as published by its authors (Rupe & Sanford 2013; Sanford et al. 2008;
Gibson et al. 2013). Usage: python3 -I scripts/divergence_figures.py figures/ results/neutral_fixation_Q.json"""
import sys, os, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import params as P

out_dir, q_json = sys.argv[1], sys.argv[2]
SURF, INK, INK2, MUTED, GRID, AXIS = '#fcfcfb', '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c3c2b7'
C_OBS, C_MIT, C_NEU = '#2a78d6', '#eb6834', '#1baf7a'
plt.rcParams.update({'font.family':'sans-serif','font.size':9,'axes.edgecolor':AXIS,'axes.linewidth':0.8,
    'axes.labelcolor':INK2,'xtick.color':MUTED,'ytick.color':MUTED,'xtick.labelcolor':INK2,'ytick.labelcolor':INK2,
    'axes.facecolor':SURF,'figure.facecolor':SURF,'savefig.facecolor':SURF,'axes.grid':True,'grid.color':GRID,
    'grid.linewidth':0.8,'grid.linestyle':'-','axes.axisbelow':True,'legend.frameon':False,'legend.fontsize':8,
    'axes.titlecolor':INK,'axes.titlesize':10,'axes.titleweight':'bold','axes.spines.top':False,'axes.spines.right':False})

# ---------------------------------------------------------------- neutral fixation transient (results/neutral_fixation_Q.json)
Qd = json.load(open(q_json)); tauN = np.array(Qd['tau_over_N']); Q = np.array(Qd['Q']); Tfix = Qd['mean_fix_time_over_N']
cum = np.concatenate([[0.], np.cumsum(0.5*(Q[1:]+Q[:-1])*np.diff(tauN))])
def factor(x):
    """(1/x) integral_0^x Q(u) du: fraction of the asymptotic count U*G that has fixed within G = x*N generations."""
    x = np.asarray(x, dtype=float); xs = np.maximum(x, 1e-9)
    return np.where(x <= tauN[-1], np.interp(x, tauN, cum)/xs, 1.0 - Tfix/xs)
def fixed_by(G, N, U, f):  return f*U*np.asarray(G)*factor(np.asarray(G)/N)   # fixed derived alleles in one lineage
def pairwise(G, Ne, U, f): return f*U*(np.asarray(G) + 2*Ne)                   # pairwise divergence per lineage

MA = P.MA; GB = P.GIBSON2013
ma_total = [r[4] for r in MA['rs2013_table1']]; ma_sel = [r[5] for r in MA['rs2013_table1'] if r[5] > 0]
U_rs = MA['rs2013_U_indiv']/2; Gp = MA['rs2013_gens']
ratios = {}
for rid, N, desc, mpi, tot, sel in MA['rs2013_table1']:
    if desc == 'all neutral':
        exp = U_rs*Gp*factor(Gp/N); ratios[N] = tot/exp
        print(f"Rupe & Sanford {rid}: N={N:,}, neutral fixations {tot}, Kimura expectation with the 4N transient {exp:.0f}, ratio {tot/exp:.3f}")
k_ma = float(np.mean(list(ratios.values())))
G_EXT0 = MA['s2008_no_selection_extinction_gen']                        # 882: U=10, no selection (Sanford et al. 2008)
G_EXT1, G_EXT2 = GB['default_U5']['extinction_gen'], GB['lynch_dfe_U3']['extinction_gen_probability_selection']   # 40,831 and 83,647 with selection
SUMMARY = []

def G_range(d):
    Tc, Tl, Th = d['T_years']; gc, gl, gh = d['gen_years']; return Tc/gc, Tl/gh, Th/gl

def make_figure(d, mode, fname, title, obs_text):
    Gc, Gl, Gh = G_range(d); Fc, Fl, Fh = d['F_obs']; U, Ul, Uh = d['U_hap']; f, fl, fh = d['f']; Ne, Nl, Nh = d['Ne_anc']
    G = np.logspace(3, 6, 400)
    if mode == 'pairwise':
        neu, neu_lo, neu_hi = pairwise(G, Ne, U, f), pairwise(G, Nl, Ul, fl), pairwise(G, Nh, Uh, fh); neu_c = float(pairwise(Gc, Ne, U, f))
    else:
        neu, neu_lo, neu_hi = fixed_by(G, Ne, U, f), fixed_by(G, Nh, Ul, fl), fixed_by(G, Nl, Uh, fh); neu_c = float(fixed_by(Gc, Ne, U, f))
    fig, ax = plt.subplots(figsize=(11, 8.6))
    # --- the entropy model's extinction horizon (deleterious mutations included)
    ax.axvline(G_EXT0, color=C_MIT, lw=1.2, ls=':')
    ax.axvspan(G_EXT1, G_EXT2, facecolor=C_MIT, alpha=0.12, lw=0, hatch='///', edgecolor=C_MIT)
    ax.axvspan(G_EXT2, 1e6, color=C_MIT, alpha=0.07, lw=0)
    ax.text(G_EXT0*1.06, 1.6, f"no selection, 10 mutations per individual:\nfitness 0 at generation {G_EXT0}\n(Sanford et al. 2008)", fontsize=7.2, color=INK2, va='bottom')
    ax.text(G_EXT1*0.96, 1.6, f"with selection (published human-\nparameter runs): fitness 0 at\ngeneration {G_EXT1:,} and {G_EXT2:,}\n(Gibson et al. 2013)", fontsize=7.2, color=INK2, va='bottom', ha='right')
    # --- neutral theory
    ax.fill_between(G, neu_lo, neu_hi, color=C_NEU, alpha=0.15, lw=0); ax.plot(G, neu, color=C_NEU, lw=2)
    # --- Mendel's Accountant neutral-class fixation rate x measured mutation rate
    ax.plot(G, k_ma*neu, color=C_MIT, lw=2, ls='--')
    # --- Mendel's Accountant as published (Rupe & Sanford 2013 Table 1)
    ax.plot([Gp, Gp], [min(ma_total), max(ma_total)], color=C_MIT, lw=2.4, solid_capstyle='round')
    ax.plot(Gp, np.exp(np.mean(np.log(ma_total))), 'o', ms=9, color=C_MIT, mec=SURF, mew=1.5, zorder=5)
    ax.plot([Gp*1.12, Gp*1.12], [min(ma_sel), max(ma_sel)], color=C_MIT, lw=2.4, solid_capstyle='round')
    ax.plot(Gp*1.12, np.exp(np.mean(np.log(ma_sel))), '^', ms=9, color=C_MIT, mec=SURF, mew=1.5, zorder=5)
    ax.text(Gp*0.9, 2.2e3, f"Mendel's Accountant as published\n(Rupe & Sanford 2013, Table 1): 300,000 generations\nat 0.01 mutations per individual per generation;\ntotal fixations {min(ma_total):,}–{max(ma_total):,} (circle),\n'selective' fixations {min(ma_sel):,}–{max(ma_sel):,} (triangle)",
            fontsize=7.2, color=INK2, ha='right', va='center')
    # arrow: the only difference between the published count and the model's own prediction at the measured rate is the input mutation rate
    y_top = k_ma*(pairwise(Gp, Ne, U, f) if mode == 'pairwise' else fixed_by(Gp, Ne, U, f)); y_bot = max(ma_total)*1.3
    ax.annotate('', xy=(Gp, y_top*0.75), xytext=(Gp, y_bot), arrowprops=dict(arrowstyle='->', color=C_MIT, lw=1.2))
    ax.text(Gp*0.93, 3e5, f"× {U/U_rs:,.0f}: the measured gametic mutation rate\n({U:.1f} per genome per generation) over the\n0.005 used in the published runs", fontsize=7.2, color=INK2, va='center', ha='right')
    # --- observed
    ax.add_patch(Rectangle((Gl, Fl), Gh-Gl, Fh-Fl, facecolor=C_OBS, alpha=0.15, lw=0))
    ax.plot([Gl, Gh], [Fc, Fc], color=C_OBS, lw=2); ax.plot(Gc, Fc, 's', ms=9, color=C_OBS, mec=SURF, mew=1.5, zorder=6)
    ax.text(1.08e3, 1.3e8, obs_text, fontsize=7.8, color=INK2, va='top', ha='left')
    ax.annotate('', xy=(Gc, Fc), xytext=(Gc, neu_c), arrowprops=dict(arrowstyle='-', color=INK2, lw=0.9))
    right = Gc > 1e5
    ax.text(Gc*(0.93 if right else 1.07), Fc*0.45 if right else (np.sqrt(Fc*neu_c) if abs(np.log10(Fc/neu_c)) > 0.25 else Fc*0.45),
            f"observed / neutral theory = {Fc/neu_c:.1f}\nobserved / model's own neutral rate at the measured mutation rate = {Fc/(k_ma*neu_c):.1f}", fontsize=7.6, color=INK2, va='center', ha='right' if right else 'left')
    if Gl <= Gp <= Gh or mode == 'pairwise' and d is P.HUMAN_CHIMP:
        ax.text(Gp*0.9, 7e4, f"observed / model as published = {Fc/max(ma_total):,.0f}–{Fc/min(ma_total):,.0f}", fontsize=7.6, color=INK2, ha='right', va='center')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(1e3, 1e6); ax.set_ylim(1, 3e8)
    ax.set_xlabel('Generations since the population split (log)'); ax.set_ylabel('Fixed SNV differences per lineage (log)')
    ax.set_title(title, loc='left')
    handles = [Line2D([],[],marker='s',color=C_OBS,lw=2,ms=8,mec=SURF,label='Observed; box = range of split time, generation time and count'),
               Line2D([],[],color=C_NEU,lw=2,label=('Neutral theory: f × U × (G + 2 Ne) with the measured mutation rate; band = parameter ranges' if mode=='pairwise'
                                                   else 'Neutral theory: f × U × G × (fraction of fixations completed within G generations, Ne 7–15k); band = parameter ranges')),
               Line2D([],[],color=C_MIT,lw=2,ls='--',label=f"Mendel's Accountant neutral-class fixation rate ({k_ma:.2f} × Kimura in its own 300,000-generation runs) × measured mutation rate"),
               Line2D([],[],marker='o',color=C_MIT,lw=0,ms=8,mec=SURF,label="Mendel's Accountant as published: total fixations at 0.01 mutations per individual per generation"),
               Line2D([],[],marker='^',color=C_MIT,lw=0,ms=8,mec=SURF,label="Mendel's Accountant as published: 'selective' fixations (total minus the all-neutral run)"),
               Line2D([],[],color=C_MIT,lw=1.2,ls=':',label="Mendel's Accountant with deleterious mutations, no selection: extinction at generation 882"),
               Patch(facecolor=C_MIT, alpha=0.25, hatch='///', edgecolor=C_MIT, label=f"Mendel's Accountant with deleterious mutations and selection: extinction at generations {G_EXT1:,}–{G_EXT2:,}; lighter shade = beyond")]
    fig.legend(handles=handles, loc='lower left', bbox_to_anchor=(0.06, 0.0), ncol=1, fontsize=8)
    fig.subplots_adjust(left=0.08, right=0.98, top=0.95, bottom=0.26)
    fig.savefig(os.path.join(out_dir, fname), dpi=170); plt.close(fig)
    if mode == 'pairwise': n_lo, n_hi = float(pairwise(Gl, Nl, Ul, fl)), float(pairwise(Gh, Nh, Uh, fh))
    else: n_lo, n_hi = float(fixed_by(Gl, Nh, Ul, fl)), float(fixed_by(Gh, Nl, Uh, fh))
    SUMMARY.append(dict(dataset=d['label'], mode=mode, G=[Gc, Gl, Gh], F_obs=[Fc, Fl, Fh], neutral=[neu_c, n_lo, n_hi], ma_neutral_scaled=[k_ma*neu_c, k_ma*n_lo, k_ma*n_hi],
                        obs_over_neutral=Fc/neu_c, obs_over_ma_neutral=Fc/(k_ma*neu_c), ma_published_total=[min(ma_total), max(ma_total)], ma_published_selective=[min(ma_sel), max(ma_sel)],
                        obs_over_published=[Fc/max(ma_total), Fc/min(ma_total)], U_over_published_U=U/U_rs))
    print(f"{d['label']}: G {Gc:,.0f} ({Gl:,.0f}–{Gh:,.0f}); observed {Fc:,.0f}; neutral {neu_c:,.0f} ({n_lo:,.0f}–{n_hi:,.0f}); obs/neutral {Fc/neu_c:.2f}; "
          f"obs/MA-neutral {Fc/(k_ma*neu_c):.2f}; obs/published {Fc/max(ma_total):,.0f}–{Fc/min(ma_total):,.0f}")

make_figure(P.HUMAN_CHIMP, 'pairwise', 'fig_human_chimp.png',
            'Human vs chimpanzee: 17.5 million fixed SNVs per lineage against the genetic-entropy model and neutral theory',
            'Observed: ~35 M SNV differences = 1.23% of aligned sites (Chimpanzee Sequencing Consortium 2005),\n15–20 M per lineage; split 5.5–7.5 Myr, 20–29-year generations (190,000–375,000 generations)')
make_figure(P.GORILLA, 'pairwise', 'fig_gorilla_gorilla.png',
            'Eastern vs western gorilla: 4.4 million pairwise SNV differences per lineage against the model and neutral theory',
            'Observed: 0.30% (0.23–0.35%) pairwise divergence between eastern and western gorilla genomes (Scally et al. 2012: mean\nsequence divergence 1.75 Myr at 1e-9 per bp per year) = 3.3–5.1 M per lineage; split 150 kyr–1.75 Myr, 19-year generations')
make_figure(P.HUMAN_NEA, 'fixed', 'fig_human_neanderthal.png',
            'Modern humans vs Neanderthal: at least 31,389 fixed derived SNVs against the genetic-entropy model and neutral theory',
            'Observed: 31,389 SNVs fixed derived in 1,094 present-day humans and ancestral in the Altai Neanderthal and Denisovan genomes\n(Prüfer et al. 2014; a lower bound); split 550–765 kyr, 25–29-year generations (19,000–31,000 generations)')

# ---------------------------------------------------------------- F1: neutral fixation in the model vs Kimura (published numbers)
fig, ax = plt.subplots(figsize=(10, 6.8))
G = np.logspace(2, 6, 500)
for N, ls in [(100, '-'), (1000, '-'), (10000, '-')]:
    y = G*factor(G/N); ax.plot(G, y, color=C_NEU, lw=2, ls=ls, alpha=0.5 + 0.25*np.log10(N/100)/2)
    ax.text(1.05e6, y[-1]*(0.55 if N == 100 else 0.7 if N == 1000 else 1.0), f'N = {N:,}', fontsize=8, color=INK2, va='center')
    ax.axvline(4*N, color=GRID, lw=1); ax.text(4*N*1.04, 1.2e2, f'4N = {4*N:,}', fontsize=7, color=MUTED, rotation=90, va='bottom')
ax.plot(G, G, color=INK2, lw=1, ls='--')
ax.text(1.5e2, 2.3e2, 'k = U (Kimura 1968): one neutral fixation per gametic mutation,\nreached after the ~4N-generation fixation transient', fontsize=8, color=INK2, rotation=40, rotation_mode='anchor')
for rid, N, desc, mpi, tot, sel in MA['rs2013_table1']:
    if desc == 'all neutral':
        ax.plot(Gp, tot/U_rs, 'o', ms=9, color=C_MIT, mec=SURF, mew=1.5, zorder=5)
        ax.annotate(f"Rupe & Sanford 2013 run {rid}: N = {N:,}, {tot:,} neutral fixations\nin 300,000 generations at a gametic rate of 0.005 = {tot/U_rs:,.0f} per unit rate\n(Kimura with transient: {U_rs*Gp*factor(Gp/N)/U_rs:,.0f}; ratio {ratios[N]:.2f})",
                    (Gp, tot/U_rs), xytext=(-12, -42 if N == 1000 else 18), textcoords='offset points', fontsize=7.4, color=INK2, ha='right', va='center')
ax.axvline(13500, color=C_MIT, lw=1.2, ls=':'); ax.text(13500*1.05, 1.1e2, "13,500 generations: length of Rupe & Sanford's Figure 1 runs\n(N = 10,000 'not even approaching fixation rate equilibrium')", fontsize=7.4, color=INK2, va='bottom')
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(1e2, 1.3e6); ax.set_ylim(1e2, 2e6)
ax.set_xlabel('Generations since the start of the run (log)'); ax.set_ylabel('Cumulative neutral fixations per unit gametic mutation rate (log)')
ax.set_title("Neutral fixation in Mendel's Accountant: published deep-time counts sit on Kimura's k = U after the 4N transient", loc='left')
handles = [Line2D([],[],color=C_NEU,lw=2,label='Wright-Fisher expectation U × G × (fraction of fixations completed within G generations), from the fixation-time simulation (results/neutral_fixation_Q.json)'),
           Line2D([],[],color=INK2,lw=1,ls='--',label='Asymptote k = U'),
           Line2D([],[],marker='o',color=C_MIT,lw=0,ms=8,mec=SURF,label="Mendel's Accountant all-neutral runs (Rupe & Sanford 2013, Table 1)")]
ax.legend(handles=handles, loc='upper left'); fig.tight_layout(); fig.savefig(os.path.join(out_dir, 'fig_F1_neutral_fixation_transient.png'), dpi=170); plt.close(fig)

# ---------------------------------------------------------------- F2: published fitness trajectories of the model vs the lineages that exist
fig, ax = plt.subplots(figsize=(11, 8.6))
g = np.logspace(1, 6, 600)
# Sanford et al. 2008: no selection, U = 10, fitness 0 at 882 (linear in generations)
ax.plot(g, np.clip(1 - g/G_EXT0, 0, 1), color=C_MIT, lw=2, ls=':')
# Gibson 2013 default DFE, U_del = 5, partial truncation, h2 = 0.4, N = 1,000: 0.79 at 10,000, linear to 0 at 40,831
x1 = np.array([1, 10000, G_EXT1]); y1 = np.array([1, 0.79, 0]); ax.plot(x1, y1, color=C_MIT, lw=2)
ax.plot(10000, 0.84, 'o', ms=8, color=C_MIT, mec=SURF, mew=1.5)
# Gibson 2013 Lynch DFE, U_del = 3, N = 10,000, h2 = 0.4: only the published endpoints are drawn; the path between them is not published
ax.plot(100000, 0.31, 'D', ms=8, color=C_MIT, mec=SURF, mew=1.5, zorder=5); ax.plot([1, 100000], [1, 0.31], color=C_MIT, lw=1, ls=':', alpha=0.6)
ax.plot(G_EXT2, 0.0, 'D', ms=8, color=C_MIT, mec=SURF, mew=1.5, zorder=5, clip_on=False); ax.plot([1, G_EXT2], [1, 0.0], color=C_MIT, lw=1, ls=':', alpha=0.6)
ax.text(100000*0.94, 0.31, "Lynch's DFE, 3 deleterious mutations per individual,\npartial truncation: 0.31 at 100,000", fontsize=7.4, color=INK2, va='center', ha='right')
ax.text(G_EXT2*1.06, 0.02, "Lynch's DFE, probability selection:\nfitness 0 at 83,647", fontsize=7.4, color=INK2, va='bottom')
ax.plot([10000, 10000], [0.22, 0.57], color=C_MIT, lw=2.4, solid_capstyle='round'); ax.plot(10000, 0.22, '^', ms=8, color=C_MIT, mec=SURF, mew=1.5, zorder=5); ax.plot(10000, 0.57, '^', ms=8, color=C_MIT, mec=SURF, mew=1.5, zorder=5)
ax.text(10000*0.94, 0.40, "Lynch's DFE, 12 offspring per female (10 removed):\n0.22 (probability) and 0.57 (partial truncation)\nat 10,000 generations", fontsize=7.4, color=INK2, va='center', ha='right')
ax.text(10000*1.06, 0.80, 'default DFE, N = 10,000: 0.84 at 10,000;\nN = 1,000: 0.79 at 10,000, then linear to 0 at 40,831', fontsize=7.4, color=INK2, va='top')
# standard model: mutation-selection balance keeps mean fitness constant; Lynch's relaxed-selection extrapolation
ax.axhline(0.86, color=C_NEU, lw=2); ax.fill_between(g, 0.82, 0.90, color=C_NEU, alpha=0.15, lw=0)
ax.text(12, 0.93, 'Standard model: mutation–selection balance; mean fitness constant at e^(−U_del) relative to a mutation-free genome\n(0.86 drawn for U_del = 0.15 strongly selected mutations per gamete; the level is a convention, the slope is zero)', fontsize=7.4, color=INK2)
for r in P.LYNCH2016['relaxed_decline_per_gen']:
    ax.plot(g, np.clip(1 - r*g, 0, 1), color=INK2, lw=1, ls='--')
ax.text(12, 0.13, 'Lynch 2016: 0.1–0.5% fitness decline per generation\nIF selection were completely relaxed\n(extrapolated from mouse mutation-accumulation lines)', fontsize=7.4, color=INK2, va='center')
# the lineages that exist
for d, lab, y, ha in [(P.HUMAN_NEA, 'Neanderthal–modern human split', 1.085, 'right'), (P.GORILLA, 'eastern–western gorilla split', 1.03, 'left'), (P.HUMAN_CHIMP, 'human–chimpanzee split', 1.06, 'center')]:
    Gc, Gl, Gh = G_range(d); ax.axvspan(Gl, Gh, color=C_OBS, alpha=0.10, lw=0); ax.axvline(Gc, color=C_OBS, lw=1.5)
    ax.text(Gc*(0.97 if ha == 'right' else 1.03 if ha == 'left' else 1.0), y, f'{lab}: {Gc:,.0f} generations', fontsize=7.4, color=C_OBS, ha=ha, va='bottom')
ax.set_xscale('log'); ax.set_xlim(10, 1e6); ax.set_ylim(0, 1.14)
ax.set_xlabel('Generations (log)'); ax.set_ylabel('Mean population fitness (initial = 1)')
ax.set_title("Published Mendel's Accountant fitness trajectories (human-like parameters) against the generations the real lineages persisted", loc='left')
handles = [Line2D([],[],color=C_MIT,lw=2,ls=':',label='Sanford et al. 2008: 10 mutations per individual, no selection; fitness 0 at generation 882'),
           Line2D([],[],color=C_MIT,lw=2,label='Gibson et al. 2013: default DFE, 5 deleterious mutations per individual, partial truncation, h² = 0.4, N = 1,000; 0.79 at 10,000, fitness 0 at 40,831'),
           Line2D([],[],marker='D',color=C_MIT,lw=1,ls=':',ms=7,mec=SURF,label="Gibson et al. 2013: Lynch's DFE, 3 deleterious mutations per individual, N = 10,000; published endpoints (0.31 at 100,000; fitness 0 at 83,647), path interpolated"),
           Line2D([],[],marker='^',color=C_MIT,lw=0,ms=7,mec=SURF,label="Gibson et al. 2013: Lynch's DFE with 12 offspring per female, fitness at 10,000 generations"),
           Line2D([],[],color=C_NEU,lw=2,label='Standard model: constant mean fitness at mutation–selection balance'),
           Line2D([],[],color=INK2,lw=1,ls='--',label='Lynch 2016: 0.1–0.5% per generation if selection were completely relaxed'),
           Patch(facecolor=C_OBS, alpha=0.3, label='Generations since the three population splits (the lineages persisted)')]
fig.legend(handles=handles, loc='lower left', bbox_to_anchor=(0.06, 0.0), fontsize=7.6)
fig.subplots_adjust(left=0.07, right=0.98, top=0.95, bottom=0.25)
fig.savefig(os.path.join(out_dir, 'fig_F2_published_fitness_trajectories.png'), dpi=170); plt.close(fig)
json.dump(dict(k_ma_neutral_ratio=k_ma, rs2013_neutral_ratios=ratios, extinction_generations=dict(no_selection=G_EXT0, with_selection=[G_EXT1, G_EXT2]), comparisons=SUMMARY), open(os.path.join(os.path.dirname(out_dir.rstrip('/')) or '.', 'results', 'divergence_summary.json'), 'w'), indent=1)
print('done')

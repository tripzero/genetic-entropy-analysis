"""Conditional fixation-time distribution of a new neutral mutation in a Wright-Fisher population.
Q(tau) = fraction of mutations destined to fix that have fixed within tau generations. In the diffusion limit Q depends
on tau only through tau/N (N = diploid population size), so one simulation at N_sim gives Q for any N by rescaling.
Expected neutral fixations completed within G generations of a split = (U/2 per gamete... i.e. U_haploid) x sum_{t<=G} Q(t).
Output: results/neutral_fixation_Q.json with the grid tau/N and Q, plus the mean conditional fixation time."""
import sys, json
import numpy as np
out = sys.argv[1]
N = 1000; twoN = 2*N; n_traj = 4_000_000; max_t = 40*N
rng = np.random.default_rng(12345)
x = np.ones(n_traj, dtype=np.int64)
fix_t = []
for t in range(1, max_t+1):
    if x.size == 0: break
    x = rng.binomial(twoN, x/twoN)
    nf = int((x >= twoN).sum())
    if nf: fix_t.extend([t]*nf)
    x = x[(x > 0) & (x < twoN)]
fix_t = np.array(fix_t)
print(f"N={N}: {n_traj} new mutations, expected eventual fixations {n_traj/twoN:.0f}, observed {fix_t.size}, "
      f"mean fixation time {fix_t.mean()/N:.3f} N generations (theory ~4N for large N), median {np.median(fix_t)/N:.3f} N, "
      f"unfinished at {max_t/N:.0f}N: {x.size}")
grid = np.arange(0, max_t+1)
Q = np.searchsorted(np.sort(fix_t), grid, side='right') / (n_traj/twoN)
# thin the grid for storage
keep = np.unique(np.concatenate([np.arange(0, 2000), np.arange(2000, max_t+1, 10)]))
json.dump({'N_sim': N, 'n_traj': n_traj, 'expected_fixations': n_traj/twoN, 'observed_fixations': int(fix_t.size),
           'mean_fix_time_over_N': float(fix_t.mean()/N), 'tau_over_N': (keep/N).tolist(), 'Q': Q[keep].tolist()}, open(out, 'w'))
for k in [0.5, 1, 1.35, 2, 3, 4, 6, 8, 12, 20]:
    i = int(k*N); print(f"  Q(tau={k:>4}N) = {Q[i]:.3f}   cumulative-sum factor sum_t<=tau Q(t) / tau = {Q[:i+1].sum()/i:.3f}")

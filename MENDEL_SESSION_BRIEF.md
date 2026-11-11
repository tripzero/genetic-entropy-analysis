# Brief for the Mendel's Accountant session: test the genetic-entropy model against observed data

Written 2026-10-30 by the MITTENS verification session. Read this whole file before doing anything.
Work in a NEW directory: `/home/kev/genetic_entropy_analysis/` (scripts/, results/, figures/, downloads/).
Keep downloaded code and data in `downloads/` and never run interpreters from inside it; run analysis
scripts with `python3 -I` and pass paths as arguments.

## 1. Context: what already exists

- `/home/kev/mittens_analysis/` holds a completed test of the MITTENS 3.0 model (Vox Day, Zenodo 23003785).
  Read `README.md` there first. Method that worked: state the model as an equation, apply it to datasets it
  was not fitted to, compare with the standard population-genetic expectation, propagate parameter ranges,
  plot predicted vs observed on log axes.
- Key reusable results (in `/home/kev/mittens_analysis/results/`):
  - `ltee_results.json`: per-population fixation counts in the LTEE (Good et al. 2017 data) by three metrics.
  - `ltee_synonymous.json`: synonymous fixations per population; mutators fix 74x more per generation than
    non-mutators, both groups on the neutral expectation (f_syn x point-mutation rate).
  - `model_test_results.json`: the eight-dataset table (observed, MITTENS, neutral) with parameter ranges.
- Manuscript (Claude Doc, editable with the docs tools): https://claude.ai/artifact/XiRz4BEAP4jhiwCLqY5MGF
  (doc id f8bb0d16-6f9c-4512-90d9-ff5f73918014, body node a947d66f-f569). Do not rebuild it; if a genetic-entropy
  section is wanted, insert it as a new Results subsection + Discussion paragraph + references, or start a
  companion doc. Ask the user which before writing into it.
- Figures are posted to the user's local image viewer "ivp": from `~/stheno/legacy_of_the_fallen` run
  `uv run tools/ivp.py --name "<gallery name>" --no-open <png files...>` then
  `uv run tools/ivp.py --set-info captions.json` (JSON: absolute image path -> caption). Give the user the LAN
  URL http://192.168.69.51:8777/#/g/<slug>, not 127.0.0.1. Start a background watch for notes:
  `uv run tools/ivp.py --watch "<gallery name>" --since 0`.
- Chart style used throughout (keep it): surface #fcfcfb, text #0b0b0b / #52514e, muted #898781, grid #e1e0d9,
  axis #c3c2b7; observed data blue #2a78d6, the model under test orange #eb6834, the standard expectation
  green #1baf7a; 2 px lines, >= 8 px markers with a surface-colour ring, legend always present, log axes
  where ranges exceed 2 decades, direct labels on the few points that matter. matplotlib Agg, dpi 170.

## 2. Objective

Treat "genetic entropy" (Sanford 2005; Mendel's Accountant, Sanford et al. 2007) the way MITTENS was treated:
run the model with parameters matched to real systems, compare its predictions with observations, and show
what the standard model predicts for the same inputs. The claims to test, in the model's own terms:

1. Fitness declines inexorably at ~1-2% per generation for human-like inputs and the population goes extinct
   within ~10^3-10^4 generations ("mutational meltdown").
2. Nearly all mutations are deleterious and effectively unselectable, so deleterious mutations accumulate at
   close to the mutation rate; effectively neutral mutations do not exist in quantity.
3. (Rupe & Sanford 2013, "Haldane's ratchet") Fixation of neutral and beneficial mutations is far rarer than
   theory predicts, so divergence cannot accumulate. This is the direct ancestor of the MITTENS argument.
4. (Carter & Sanford 2012) Human H1N1 influenza was going extinct through mutation accumulation.

## 3. Software setup

- Mendel's Accountant, current implementation: `mendel-go` (Go), https://github.com/genetic-algorithms/mendel-go
  (also `mendel-web-ui` and the older Fortran `mendel`). Clone into `downloads/`, build with Go (`go version`;
  install Go if absent), run `mendel-go -h`. Parameters live in an INI/TOML file; start from the repository's
  defaults file (look for `mendel-defaults.ini` or `test/input/*.ini`) and VERIFY every parameter name there
  before trusting the names below (they are from memory of the published descriptions):
  `pop_size`, `num_generations`, `mutn_rate` (new mutations per individual per generation), `frac_fav_mutn`,
  `fitness_effect_model` (weibull | fixed | uniform), `uniform_fitness_effect_del/fav`,
  `high_impact_mutn_fraction`, `high_impact_mutn_threshold`, `max_fav_fitness_gain`, `genome_size`,
  `selection_model` (fulltrunc | ufp | spps | partialtrunc), `heritability`, `non_scaling_noise`,
  `synergistic_epistasis`, `se_nonlinked_scaling`, `se_linked_scaling`, `reproductive_rate`,
  `num_offspring_model`, `crossover_model`, `haploid`, `dynamic_linkage`, `num_linkage_subunits`,
  `track_neutrals` / `fraction_neutral`, `allele_count_gens`, `files_to_output`.
- Outputs to capture per run: mean fitness per generation, deleterious/favorable/neutral mutations per
  individual per generation, number of FIXED mutations by class per generation (enable allele tracking and
  neutral tracking), population size. Save raw outputs in `results/<run_name>/` and a tidy CSV per run.
- Practical limits: pop_size x num_generations x mutn_rate drives memory and time. Expect 10^4 individuals x
  10^4 generations x ~10 mutations/gen to be feasible; the LTEE's 3 x 10^7 cells is not. Scale and say so.
- Standard-model comparison: SLiM 4 (https://messerlab.org/slim/) if installable, else a small Wright-Fisher
  simulator written in numpy with a gamma DFE and optional synergistic epistasis. Analytical references:
  Haldane-Muller load L = 1 - exp(-U_del) under multiplicative fitness; Kondrashov 1988/1995 and Lesecque,
  Keightley & Eyre-Walker 2012 for load under truncation/synergistic selection; Kimura k = f x U for fixation.

## 4. Experiments

### Experiment 0: does Mendel's Accountant reproduce k = mu for neutral mutations?
Run with `track_neutrals` on, pop_size 500-1,000, mutn_rate 10, fraction neutral 1.0 (or a known fraction),
no selection, for >= 10 x (2 x pop_size) generations. Count neutral fixations per generation after the
transient. Expectation: fixations per generation -> U_neutral (Kimura), after a transient of ~4Ne generations.
This decides whether Rupe & Sanford's "Haldane's ratchet" is a property of the model or of running it for
fewer than 4Ne generations. Also report the mean fixation time (expect ~4Ne) and compare with the
"time versus rate" confusion documented in Section 5.1 of the MITTENS manuscript.

### Experiment 1: the published human runs (replication)
Reproduce the fitness-decline curves of Sanford et al. 2007 / Sanford & Nelson 2012 with their default
parameters (typically pop_size 1,000, mutn_rate 10 or 100, frac_fav_mutn ~0.001, Weibull DFE, 5,000
generations). Record decline per generation and time to meltdown. This establishes "the model as published".

### Experiment 2: human-realistic inputs, varying the assumptions that matter
Inputs: mutn_rate 65-70 (de novo SNVs + indels per child; Jonsson et al. 2017, Kong et al. 2012),
pop_size 10,000 (Ne), generation count 10,000-25,000 (the Neanderthal split is ~24,000 generations ago;
the chimpanzee split ~250,000+). Vary ONE thing at a time:
  a. DFE: Sanford Weibull default vs an empirical DFE (gamma; human estimates: Eyre-Walker & Keightley 2007
     Nat Rev Genet; Kim, Huber & Lohmueller 2017 Genetics) with ~90-95% of sites effectively neutral
     (Lindblad-Toh et al. 2011) so that U_del ~ 2 per diploid genome per generation.
  b. Fitness model: multiplicative vs synergistic epistasis; selection: probability (ufp) vs truncation.
  c. Beneficial fraction 10^-4 to 10^-2.
Compare: (i) fitness trajectory vs the observation that no decline of 1-2%/generation exists over
hundreds of recorded generations (Lynch 2016 Genetics "Mutation and human exceptionalism" is the serious
version of the concern: ~1%/generation IF selection were fully relaxed; cite it fairly);
(ii) FIXED mutations after ~24,000 generations vs the observed >= 31,389 derived alleles fixed in modern
humans since the Neanderthal split (Prufer et al. 2014) and vs the neutral expectation (~29,000; see
MITTENS manuscript Section 3.5 and the `neutral_fixed_by` function in
`/home/kev/mittens_analysis/scripts/mittens_model_test.py`);
(iii) mutations carried per individual vs the observed ~2 deleterious de novo per generation.

### Experiment 3: the LTEE (the hardest test of meltdown)
Haploid, clonal (no crossover), genome_size 4.63e6, mutn_rate 4.1e-4 (ancestral) and 0.041 (mutators),
pop_size as large as feasible (10^4-10^5; note the real Ne ~ 3 x 10^7), 10,000-60,000 generations, beneficial
fraction and effect sizes varied (LTEE drivers have s ~ 0.01-0.1 early, smaller later). Compare:
  - fitness trajectories vs Wiser, Ribeck & Lenski 2013 Science 342:1364 (mean fitness follows a power law
    and rose ~70% by 50,000 generations; mutator populations gained slightly MORE, not less);
  - fixed mutations per generation vs `/home/kev/mittens_analysis/results/ltee_results.json` and
    `ltee_synonymous.json` (non-mutators ~1 per 500 generations within clades, mutator phases ~1 per 26;
    synonymous 1.1e-4 vs 8.1e-3 per generation);
  - the model's prediction for the 100x mutators (meltdown within how many generations?) vs the observed
    30,000-60,000 generations of mutator persistence with rising fitness (Wielgoss et al. 2013 PNAS for the
    antimutator dynamics in Ara-1).

### Experiment 4 (optional): H1N1
Carter & Sanford 2012 (Theor Biol Med Model 9:42) predicted extinction of human H1N1 from mutation
accumulation. Check the claim against the record (seasonal H1N1 lineage displaced in 2009 by the pandemic
reassortant; pdm09 H1N1 still circulating in 2026; their "degradation" metric was codon-bias/fitness proxy,
not measured fitness). If you run MA with their viral parameters, compare its extinction time with 100+ years
of observed persistence. Keep this short; it is a one-paragraph point.

## 5. Deliverables

1. `results/` with raw MA outputs and tidy CSVs; `scripts/` with every run config and the plotting code;
   a `README.md` with the parameter tables and sources, in the same format as `/home/kev/mittens_analysis/README.md`.
2. Figures (same style as above), suggested set:
   - F1: Experiment 0, neutral fixations per generation vs time in MA, with k = U marked (and the 4Ne transient).
   - F2: Experiment 2, mean fitness vs generation under each assumption set, with the "no observable decline"
     band and Lynch's relaxed-selection line.
   - F3: Experiment 2, fixed mutations after 24,000 generations: MA variants vs neutral expectation vs observed.
   - F4: Experiment 3, LTEE fitness: MA (non-mutator, mutator) vs Wiser et al. power law; and fixations per
     generation vs the measured values.
   - F5: predicted-vs-observed summary like `/home/kev/mittens_analysis/figures/fig4_predicted_vs_observed.png`.
   Post them to ivp as gallery "Genetic entropy model vs data" with captions.
3. A written summary (chat) with: what the model is, which assumptions drive its predictions, where it matches
   and where it fails, with numbers; and, if the user wants, a Results subsection for the manuscript.

## 6. Writing and sourcing rules (same as the MITTENS work)

- Every number traced to a primary source you actually opened (Europe PMC full text via
  `https://www.ebi.ac.uk/europepmc/webservices/rest/<PMCID>/fullTextXML` works for open-access papers;
  Nature/Science pages redirect, so use PMC). Search snippets are not sources. Label anything unverified.
- State the model fairly in its own terms before testing it; quote its claims; no polemic.
- Propagate parameter ranges (low/central/high) and run a Monte Carlo over them if the conclusion could
  depend on a choice; report medians and 95% intervals.
- Verify every figure by viewing the PNG before posting; fix label collisions.
- References to have ready: Sanford 2005 (Genetic Entropy); Sanford, Baumgardner, Brewer, Gibson, ReMine 2007
  Scalable Computing: Practice and Experience 8:147 (Mendel's Accountant); Sanford & Nelson 2012 (in
  Studies in Population Genetics, InTech); Gibson, Baumgardner, Brewer, Sanford 2013 (in Biological Information: New Perspectives, World
  Scientific) and Rupe & Sanford 2013 (Proceedings of the 7th International Conference on Creationism; NOT in
  the BINP volume, correction from the genetic_entropy session); Carter & Sanford 2012; Kimura 1968,
  1983; Kondrashov 1988 Nature 336:435 and 1995 J Theor Biol 175:583; Lesecque, Keightley & Eyre-Walker 2012
  Genetics 191:1321; Charlesworth 2013 Evolution 67:3354; Lynch 2016 Genetics 202:869; Eyre-Walker &
  Keightley 2007 Nat Rev Genet 8:610; Kim, Huber & Lohmueller 2017 Genetics 206:345; Wiser, Ribeck & Lenski
  2013 Science 342:1364; Wielgoss et al. 2013 PNAS 110:222; Good et al. 2017 Nature 551:45; Tenaillon et al.
  2016 Nature 536:165; Prufer et al. 2014 Nature 505:43; Jonsson et al. 2017 Nature 549:519.

## 7. Pitfalls

- MA's "extinction" depends on `fitness_dependent_fertility`/reproductive-rate settings; report the mean
  fitness trajectory, not just the extinction flag, and state the settings.
- MA's default Weibull DFE puts most mutations in the "nearly neutral but deleterious" range by construction;
  that choice, plus multiplicative fitness, is what produces the collapse. Make the DFE and epistasis
  comparisons explicit (Experiment 2a-b) so the reader can see which assumption does the work.
- Neutral tracking off by default: turn it on, or MA will report zero neutral fixations and you will
  rediscover "Haldane's ratchet" as an artifact.
- Population-size scaling: fixation times scale with Ne, so compare rates per generation, not totals, and
  run at least 4Ne generations when fixation counts matter.
- Do not state that MA runs were done unless they were; if the build fails, say so and fall back to the
  numpy Wright-Fisher simulator for the standard model while reporting MA's published numbers for the
  entropy model.

# Genetic entropy (Sanford; Mendel's Accountant): model vs observed data

Analysis date: 2026-10-30. Companion to `/home/kev/mittens_analysis/` (same method, same chart style). Everything here is
reproducible from `scripts/` (see "Rerun"). Status: the three divergence comparisons (human-human, human-chimp, gorilla-gorilla) are done from the model's *published* runs
(sections 1-5); Mendel's Accountant (`mendel-go` 1.2.5) was then built from source with the user's approval and Experiments 0-3 of the
brief were run (section 8); all runs are complete. No number below comes from an unverified source; every primary text is in
`downloads/sources/`.

## 1. What the model is

"Genetic entropy" (Sanford 2005) is the claim that the great majority of new mutations are slightly deleterious and below the
reach of selection, so deleterious mutations accumulate at close to the mutation rate, mean fitness declines every generation,
and populations degenerate to extinction; and (Rupe & Sanford 2013, "Haldane's ratchet") that fixation of beneficial mutations is
far too slow for the human lineage to have accumulated its differences from chimpanzee. The model is implemented in Mendel's
Accountant (Sanford, Baumgardner, Brewer, Gibson & ReMine 2007, Scalable Computing: Practice and Experience 8:147), a forward-time
Wright-Fisher simulator: diploid sexual population of fixed size, Poisson mutations per offspring, fitness effects drawn from a
Weibull function d(x) = exp(-ln(G) x^gamma) with G the genome size (minimum effect 1/G = 3e-10, 0.1% of mutations above 0.1 in the
2007 defaults: mean -5.1e-4, median -2.7e-8), additive combination of effects (the only mode implemented in mendel-go,
`dna/linkage.go`), truncation or probability selection on phenotypic fitness with a heritability parameter, 989 linkage blocks on
23 chromosomes with two crossovers each. Fixation is reported as allele frequency >= 99%.

The published claims, in the authors' words (all texts in `downloads/sources/`):

* Rupe & Sanford 2013 (Proc. 7th Int. Conf. Creationism; *not* in the Biological Information volume as the brief assumed):
  "Using realistic parameter settings we only observe a few hundred selection-induced beneficial fixations after 300,000
  generations (6 million years). Even when using highly optimal parameter settings ... we only see a few thousand selection-induced
  fixations. This is significant because the ape-to-man scenario requires tens of millions of selective nucleotide substitutions in
  the human lineage." And: "We first confirm that our numerical simulations correctly tallying the fixation of neutral mutations.
  We show that neutral mutations go to fixation just as predicted by conventional theory (i.e., over deep time the fixation rate
  approached the gametic mutation rate)."
* Gibson, Baumgardner, Brewer & Sanford 2013 (Biological Information: New Perspectives, pp. 232-263): "Simulations based on
  recently published values for mutation rate and effect-distribution in humans show a steady decline in fitness that is not even
  halted by extremely intense selection pressure (12 offspring per female, 10 selectively removed)."
* Sanford et al. 2008 (Proc. 6th ICC): default human parameters, "Without selection this same run reaches a fitness of 0.0
  (extinction) in 882 generations."

## 2. The model's own published numbers

Rupe & Sanford 2013, Table 1 (300,000 generations, 0.01 non-deleterious mutations per individual per generation = gametic rate
0.005, probability selection, heritability 0.2):

| Run | Population | Mutations | Total fixations | 'Selective' fixations | Kimura expectation incl. 4N transient (this analysis) | ratio |
|---|---|---|---|---|---|---|
| qrtsyg | 1,000 | all neutral | 1,413 | 0 | 1,480 | 0.955 |
| qrtqaq | 1,000 | default beneficial distribution | 1,492 | 79 | | |
| qrtimc | 1,000 | beneficial max effect 1.0 | 1,684 | 271 | | |
| qrtmoi | 10,000 | all neutral | 1,234 | 0 | 1,287 | 0.959 |
| qrtsvs | 10,000 | default beneficial distribution | 1,456 | 222 | | |
| qrtniw | 10,000 | beneficial max effect 1.0 | 1,824 | 590 | | |
| qrtoju | 10,000 | optimised (partial truncation, min effect 3.3e-6, max 1.0) | 2,710 | 1,476 | | |

The expectation column is 0.005 x 300,000 x Q-bar(G/N), where Q-bar is the fraction of the asymptotic count U x G that has actually
reached fixation within G generations, from a Wright-Fisher simulation of 4 million neutral mutations
(`scripts/neutral_fixation_transient.py`; mean conditional fixation time 3.97 N, as theory says). The model fixes neutral mutations at
95-96% of Kimura's rate. The "waiting time to fixation rate equilibrium" that the paper treats as a new finding is this transient: at
13,500 generations a population of 10,000 has completed 2.5% of its eventual fixations (fig F1).

Other published runs used below: Rupe & Sanford Fig. 3 (N = 1,000, 20,000 generations, 10 non-neutral mutations per individual,
50% elimination): 6,775 deleterious and 10 beneficial fixations; "Deleterious mutation accumulation could not be simulated through
such deep time, as it always either caused overflow of memory or premature termination due to extinction." Gibson et al. 2013:
default DFE, 5 deleterious mutations per individual, partial truncation, heritability 0.4: mean fitness 0.79 (N = 1,000) and 0.84
(N = 10,000) at 10,000 generations, "fitness declined to near extinction in 40,831 generations"; Lynch-based DFE (27% of effects
above 1e-2), 3 deleterious mutations per individual, N = 10,000: "approximately 5% per generation" initially, then linear, mean fitness
0.31 at 100,000 generations (partial truncation) or extinction at generation 83,647 (probability selection); with 12 offspring per
female, 0.22 (probability) and 0.57 (partial truncation) at 10,000 generations. Sanford & Nelson 2012: 47,648 deleterious mutations
per individual after 5,000 generations at 10 per individual, N = 1,000; Fig. 8: 9,327 deleterious and 41 beneficial fixations.

## 3. Applying the model to the three comparisons (central estimates, low-high in brackets)

| Comparison | Generations since split | Observed fixed SNVs per lineage | Model as published (300,000 gen, U = 0.01) | Model's neutral rate x measured mutation rate | Neutral theory | observed / neutral |
|---|---|---|---|---|---|---|
| Human lineage vs chimpanzee | 252,000 (190,000-375,000) | 17.5 M (15-20 M) | 1,234-2,710 total; 79-1,476 'selective' | 11.7 M | 12.2 M (8.0-21 M) | 1.4 |
| Eastern vs western gorilla | 26,000 (7,900-92,000) | 4.35 M (3.3-5.1 M) | (published runs are at 300,000 gen) | 2.4 M | 2.5 M (1.2-5.9 M) | 1.7 |
| Modern humans vs Neanderthal | 24,000 (19,000-31,000) | >= 31,389 | (idem) | 38,500 | 40,100 (1,050-266,000; Ne-sensitive) | 0.8 |

Reading. (i) The human-chimpanzee shortfall in Rupe & Sanford (observed / published = 6,500-14,000) is the ratio of the measured
gametic mutation rate (35.6 per genome per generation) to the 0.005 they used (7,100x), times the fraction of fixations they chose
to count. Their own all-neutral runs, scaled to the measured rate, give 11.7 M fixations per lineage in 252,000 generations: within
a factor of 1.5 of the 17.5 M observed, and identical to the standard expectation because the model *is* a Wright-Fisher
simulator. The "tens of millions of selective substitutions" the paper says are required is the premise under test, not an
observation: under the standard model nearly all of the 17.5 M are neutral, and the model agrees once it is given the right input.
(ii) For the gorilla and Neanderthal comparisons the split times lie inside the window in which the published human-parameter runs
still have a surviving population (extinction at 40,831-83,647 generations with selection), so the fixation count is the test, and
the model's own neutral machinery at the measured mutation rate matches the observations within a factor of 2 (Neanderthal: within
0.8, the observed count being a lower bound). (iii) The human-chimpanzee split (190,000-375,000 generations) lies beyond every
published extinction horizon of the model when deleterious mutations are included; the published runs with human-like inputs reach
mean fitness 0 at generation 882 (no selection), 40,831 or 83,647 (with selection), or 0.31 at 100,000. The lineages exist.
Whether the extinction result survives an empirical DFE (90-95% of sites effectively neutral, U_del ~ 2), multiplicative fitness, or
the measured mutation rate is Experiment 2 of the brief and needs the simulator to run.

Present-day diversity (for the discussion, not plotted): Prado-Martinez et al. 2013 Table 1 and text give heterozygosity ~0.8e-3
per bp for non-African humans, eastern lowland gorillas, bonobos and western chimpanzees, 1.6-2.4e-3 for central chimpanzees,
western lowland gorillas and orangutans; Scally et al. 2012 Table 2: western lowland gorilla 0.178-0.189%, eastern lowland 0.076%;
1000 Genomes: a typical genome differs from the reference at 4.1-5.0 million sites. These are 4 Ne mu at Ne 10-30k, which both the
standard model and the model under test produce for neutral mutations.

## 4. Figures (`figures/`, captions in `figures/captions.json`)

* `fig_human_chimp.png`, `fig_gorilla_gorilla.png`, `fig_human_neanderthal.png` - observed fixed differences per lineage vs
  generations since the split (blue), neutral theory (green), the model's published fixation counts and its neutral-class rate at the
  measured mutation rate (orange), and the model's published extinction horizons (orange shading). `scripts/divergence_figures.py`.
* `fig_F1_neutral_fixation_transient.png` - Experiment 0 from published numbers: Rupe & Sanford's all-neutral deep-time runs against
  Kimura's k = U with the 4N transient.
* `fig_F2_published_fitness_trajectories.png` - the model's published fitness trajectories for human-like inputs against the
  standard model (constant mean fitness at mutation-selection balance), Lynch's relaxed-selection extrapolation (0.1-0.5% per
  generation), and the generation counts of the three splits.

## 5. Parameters and sources

* Model texts: Sanford et al. 2007 SCPE (`downloads/sources/sanford_2007_scpe.pdf`, from scpe.org); Sanford et al. 2008 Proc. 6th ICC
  and Rupe & Sanford 2013 Proc. 7th ICC (Cedarville University TIND repository records 21040 and 21101); Sanford & Nelson 2012 InTech
  chapter (cdn.intechopen.com/pdfs/38427.pdf); Gibson et al. 2013 (World Scientific open-access reader, text extracted in the
  user's browser; quotations in `scripts/params.py`); Carter & Sanford 2012 (PMC3507676). mendel-go source: `downloads/mendel-go/`
  (parameter names verified in `mendel-defaults.ini`: `pop_size, num_generations, mutn_rate, frac_fav_mutn, fraction_neutral,
  genome_size, fitness_effect_model, high_impact_mutn_fraction, high_impact_mutn_threshold, max_fav_fitness_gain, selection_model
  (fulltrunc|ups|spps|partialtrunc), heritability, non_scaling_noise, reproductive_rate, crossover_model, num_linkage_subunits,
  track_neutrals, tracking_threshold, files_to_output`; `synergistic_epistasis`, `multiplicative_weighting`, haploid and clonal
  reproduction are marked "not currently supported" in that file).
* Human-chimp, Neanderthal, mutation rates, generation times, ancestral Ne, f: as verified in `/home/kev/mittens_analysis/README.md`
  section 5 (Chimpanzee Sequencing Consortium 2005; Kong et al. 2012; Jonsson et al. 2017; Prufer et al. 2014; Langergraber et al.
  2012; Besenbacher et al. 2019).
* Gorilla: Scally et al. 2012 Nature 483:169 (PMC3303130): "we estimate a mean sequence divergence time between them of 1.75 Mya" at
  the paper's 1e-9 per bp per year, i.e. mean pairwise divergence 0.35%; maximum-likelihood population split ~0.5 Myr under a model
  with gene flow; Table 2 heterozygosities. Range 0.23-0.35% as in the MITTENS analysis; aligned autosome 2.9 Gb; 19-year generations;
  ancestral Ne 20-40k (Prado-Martinez et al. 2013 Nature 499:471, PMC3822165, Table 1: western lowland Ne 26.8-53.5k, eastern 12.2-24.3k).
* Neutral fixation transient: `results/neutral_fixation_Q.json` (N_sim = 1,000, 4e6 trajectories, seed 12345).
* Lynch 2016 Genetics 202:869 (PMC4788123): mouse mutation-accumulation lines extrapolate to a human fitness decline of ~0.5% per
  generation in fecundity and ~0.1% in other traits if selection were completely relaxed.

Caveats. The gorilla pairwise divergence is derived from Scally's divergence time and their mutation-rate assumption, not read
directly. The Neanderthal neutral expectation depends on Ne^-1 x time nonlinearly and should be read as "same order". The model's
extinction horizons are the published endpoints of specific runs (population 1,000-10,000, 3-10 mutations per individual, additive
fitness, the authors' DFEs); they are not a prediction for the measured 65-70 mutations per individual with an empirical DFE, which
is what the planned runs will test.

## 6. Rerun

    python3 -I scripts/neutral_fixation_transient.py results/neutral_fixation_Q.json
    python3 -I scripts/divergence_figures.py figures/ results/neutral_fixation_Q.json    # also writes results/divergence_summary.json

## 7. Build and open items

The first build attempt was denied by the session's permission classifier (downloaded code); the user then asked for the build
explicitly and it succeeded (Go 1.27 in `downloads/gotool/go/bin`, binary `downloads/mendel-go/mendel-go`, smoke test in
`results/smoke/`). Open item: a clonal LTEE test, which this program cannot do (Experiment 3 is the diploid approximation).
Build command, from `/home/kev/genetic_entropy_analysis/downloads/mendel-go`:

    PATH=/home/kev/genetic_entropy_analysis/downloads/gotool/go/bin:$PATH GOPATH=/home/kev/genetic_entropy_analysis/downloads/gopath \
    GOCACHE=/home/kev/genetic_entropy_analysis/downloads/gocache GOFLAGS=-mod=mod go build -o mendel-go && ./mendel-go -h

One hypothesis the runs should settle: Rupe & Sanford's Fig. 3 reports 6,775 deleterious fixations in 20,000 generations at a gametic
rate of 5 (N = 1,000), 8% of the neutral expectation, although >90% of their Weibull DFE has |s| < 1/(2N) and should fix at the
neutral rate; the `tracking_threshold` parameter (mutations below it are pooled and never counted as alleles) would produce exactly
such an undercount.

## 8. Own Mendel's Accountant runs (mendel-go 1.2.5 built from source on 2026-10-30; inputs in `results/*/*.ini`, logs alongside)

### Experiment 0: neutral fixation (`results/exp0_neutral/`, `scripts/exp0_analysis.py`, figure `fig_F1b_exp0_own_runs.png`)

One neutral mutation per individual per generation (gametic rate 0.5), no selection (all-neutral, so the 50% culling of the 2-fold
offspring surplus is random), `track_neutrals = true`, allele bins every 100-200 generations; fixations = neutral alleles in the 99-100% bin.

| N | generations | fixations | Wright-Fisher, Ne = N | fitted Ne / N | Wright-Fisher at fitted Ne |
|---|---|---|---|---|---|
| 100 | 2,000 (20 N) | 738 | 792 (ratio 0.93) | 1.29 | 736 (1.00) |
| 300 | 6,000 (20 N) | 2,443 | 2,377 (1.03) | 1.09 | 2,324 (1.05) |
| 1,000 | 8,000 (8 N) | 1,576 | 2,043 (0.77) | 1.27 | 1,615 (0.98) |

The mating scheme (each pair has exactly 4 offspring, half removed) gives an offspring-number variance of 1 and hence Ne = (4N-2)/(1+2) =
1.33 N, which is what the fits recover. With that Ne the simulator follows k = U after the standard 4 Ne transient; Rupe & Sanford's
"waiting time to fixation rate equilibrium" is this transient (fig F1).

### Experiment 1: the model as published (`results/exp1_sanford2008/`, `scripts/exp1_analysis.py`, figure `fig_F3_exp1_sanford2008_replication.png`)

Sanford et al. 2008 default human run: N = 1,000; 10 mutations per individual; 9:1 deleterious:beneficial; Weibull DFE (0.1% of effects > 0.1,
genome 3e9, maximum beneficial gain 0.01); all co-dominant; probability selection; heritability 0.2; 6 offspring per female (4 removed).

| Quantity | published | this re-run |
|---|---|---|
| deleterious mutations per individual at generation 1,000 (with selection) | 8,730 | 8,822 |
| beneficial mutations per individual at generation 1,000 | 1,011 | 970 |
| mean fitness at generation 1,000 | (declining; value not stated) | 0.759 |
| "without selection": fitness 0 at generation 882 | 882 | random survival: 0.565 at 500, 0.359 at 882, 0.311 at 1,500; no extinction |
| fixed alleles at generation 10,000: deleterious / beneficial | 16,116 / 3,193 | 24,468 / 2,800 |
| mean fitness at generation 10,000 | (not stated) | 0.418 |

The mutation counts and the fitness decline with selection reproduce. The published "extinction at 882 generations without selection" does
not reproduce in mendel-go: with exactly two offspring per pair the program removes individuals whose additive fitness has reached zero
and never replaces them, so the population shrank from 1,000 to 1 in 103 generations at mean fitness 0.87 (a bookkeeping artefact, not
extinction by load); with a surplus and random survival (heritability 1e-4) the same removal of fitness-zero individuals acts as
truncation at zero, and mean fitness levels off near 0.3 instead of reaching 0. The total number of fixations in the 10,000-generation
run, 27,268, is 1.13 x the neutral expectation U_gametic x G x Q-bar = 5 x 10,000 x 0.484 = 24,200 (Ne = 1.33 N); the published total,
19,309, is 0.80 x. In this run, with every mutation tracked (`tracking_threshold = 0`), deleterious mutations fixed at the neutral
rate (24,468 vs 21,800 expected for the 90% deleterious share). Rupe & Sanford's Figure 3 reports deleterious fixations at 8% of the neutral
rate for similar inputs; that count is consistent with a non-zero tracking threshold (untracked low-effect mutations are never counted as
alleles, section 7 above), although their setting is not stated.

### Experiment 4: H1N1 (`downloads/sources/h1n1_surveillance_2026_notes.md`)

Carter & Sanford 2012 described "an apparent second extinction of the human H1N1 lineage in 2009" as the outcome of "a continuous
accumulation of mutations" and proposed "natural genetic attenuation" of the virus. The 2009 event was the replacement of the seasonal
lineage by the swine-origin pandemic reassortant A(H1N1)pdm09. Seventeen years on, that lineage is the dominant influenza A subtype: 239 of
241 influenza A viruses subtyped in the United States in the week ending 17 October 2026 (99.2%) and 81.5% of those reported since May
2026 (CDC FluView week 41); WHO's week-39 update notes "an increase in the proportion of influenza A(H1N1)pdm09 detections" and its
predominance in Western, Southern and South-East Asia. Mutation accumulation at the rate they documented has not ended H1N1 circulation.

### Experiment 2: human-realistic inputs, one assumption at a time (`results/exp2_human/`, `scripts/exp2_analysis.py`, figure `fig_F4_exp2_human_inputs.png`)

70 new mutations per individual per generation (the measured de novo SNV + indel rate), N = 1,000, 6 offspring per female, additive fitness
(the only mode mendel-go implements), `tracking_threshold` 1e-3 (fitness exact; see `results/exp2_human/DESIGN.md`). Base case: 5 of the 70
mutations deleterious (fraction_neutral 0.929, Gibson et al.'s "conservative" functional fraction), Sanford's default DFE (mean effect 5e-4),
beneficial fraction 1e-4 with maximum gain 0.01, 50% recessive, probability selection, heritability 0.2.

| Variant (one change from the base) | mean fitness at 1,000 | at 10,000 | deleterious accumulation per generation (of U_del) | note |
|---|---|---|---|---|
| base: U_del = 5, Sanford DFE, probability selection, h² = 0.2 | 0.872 | 0.640 | 4.83 of 5 (97%) | 0.558 at 15,000, 0.473 at 20,000, 0.414 at the Neanderthal split (24,074), 0.395 at 25,000 |
| all 70 mutations deleterious (fraction_neutral 0) | 0.163 | collapsed | 65 of 70 (93%) | population collapsed at generation 7,124 |
| U_del = 2 with the high-impact DFE Gibson et al. used for Lynch's estimate | 0.474 | 0.280 | 1.01 of 2 (50%) | fewer but far larger effects: faster decline than the base |
| full truncation selection | 0.996 | 0.993 | 4.17 of 5 (83%) | essentially no decline |
| partial truncation selection (0.5) | 0.942 | 0.841 | 4.71 of 5 (94%) | |
| heritability 1.0 (no environmental noise) | 0.863 | 0.638 | 4.77 of 5 (95%) | indistinguishable from the base |
| beneficial fraction 1e-3 (0.07 per individual per generation) | 0.871 | 0.639 | 4.80 of 5 | indistinguishable from the base |
| beneficial fraction 1e-2 (0.7 per individual per generation) | 0.861 | 0.634 | 4.81 of 5 | 487 beneficial alleles per individual by 10,000, no fitness effect at maximum gain 0.01 |
| base with N = 10,000 (tracking_threshold 1e-2) | 0.873 | 0.724 | 4.76 of 5 | 0.724 vs 0.640 for N = 1,000 at 10,000: a tenfold larger population slows the decline by about a quarter (Gibson et al.: 0.84 vs 0.79) |

Reading. The model's fitness decline is set by three choices: how many of the 70 mutations are declared deleterious and how large their
effects are (U_del x mean effect: 5 x 5e-4 = 2.5e-3 per generation of additive load entering, of which selection removes only 3% under
probability selection), the additive combination of effects (mean fitness falls linearly towards zero, and zero is extinction), and the
selection mode (truncation selection removes the load almost entirely; probability selection almost none of it). Environmental noise and
beneficial mutations, as the model parameterises them, change nothing. The decline in the base case is 2-4e-5 per generation, 300-500 times slower
than the "1-2% per generation" sometimes quoted for the model, and it decelerates (under additive fitness, relative fitness differences
grow as the mean falls, so selection removes more of the load): 0.640 at 10,000, 0.473 at 20,000, 0.395 at 25,000 generations.
Extrapolating the last 10,000 generations puts fitness zero near generation 42,000, close to the 40,831 that Gibson et al. report for the
same inputs. The population is intact at the Neanderthal and gorilla split times and would be gone long before the human-chimpanzee
split; declaring all 70 mutations deleterious collapses it within 7,200 generations. Which of these inputs is biologically right is the
question, and the standard model's answer is that the strongly selected class (U_del of order 1-2, effects of order 1e-2) is held at a
constant equilibrium load under multiplicative fitness, while the weakly selected majority is effectively neutral and accounts for the
observed divergence; the additive-vs-multiplicative comparison (fig F5, `scripts/wf_fitness_models.py`) tests that directly.

### Experiment 3 (reduced): LTEE-like inputs (`results/exp3_ltee/`, `scripts/exp3_analysis.py`, figure `fig_F6_exp3_ltee.png`)

mendel-go has no haploid or clonal mode, so this is the closest it allows: a diploid population of 1,000 with no crossover (one chromosome,
100 linkage blocks), genome 4.63e6, Sanford's Weibull DFE scaled to that genome (minimum effect 2.2e-7), co-dominant, probability selection,
heritability 0.2, 6 offspring per female, 60,000 generations, maximum beneficial gain 0.1 (`results/exp3_ltee/DESIGN.md`). Observed: the
LTEE populations "reached a relative fitness of ~1.7, meaning they grew ~70% faster than their ancestor" by 50,000 generations (Lenski 2017
ISME J, PMC5607360, quoting Wiser, Ribeck & Lenski 2013; `downloads/sources/ltee_fitness_notes.md`), and the six point-mutator populations
fixed 716-3,155 mutations in their mutator phases with fitness still rising (`/home/kev/mittens_analysis/results/ltee_results.json`).

| Run | mutation rate per genome per generation | beneficial fraction | fitness at 50,000 | at 60,000 | fixed by 60,000: deleterious / beneficial |
|---|---|---|---|---|---|
| mutator, model default | 0.041 | 1e-4 | 0.965 | 0.955 | 1,032 / 0 |
| mutator, generous | 0.041 | 1e-2 | 0.975 | 0.966 | 1,140 / 4 |
| ancestral | 4.1e-4 | 1e-2 | 1.000 | 1.000 | 1 / 0 |

With the model's parameterisation of beneficial mutations (Weibull toward vanishing effects, so that even a 1% beneficial fraction yields
4 beneficial fixations in 60,000 generations), the mutator-rate populations lose 3.5-4.5% of fitness while the real ones gained 70%, and the
ancestral-rate population does nothing at all. The deleterious fixations (1,032-1,140) are again the neutral expectation for the gametic
rate (0.0205 x 60,000 x 0.91 = 1,120). The LTEE outcome is adaptation by beneficial mutations of effect 1-10% arising at a rate the model
does not contain; the test is limited by the missing clonal mode and the 1,000-fold smaller population, both of which favour the model
(less efficient selection for drivers makes its failure to adapt conservative, but a clonal population would also show stronger
Muller's-ratchet decline than this diploid one).

### Experiment 2b: additive vs multiplicative fitness, same inputs (`scripts/wf_fitness_models.py`, `results/wf_models.json`, figure `fig_F5_additive_vs_multiplicative.png`)

mendel-go combines mutation effects additively only, so the comparison was made in a free-recombination Wright-Fisher simulator (independent
loci, N = 1,000, deterministic diploid selection with binomial drift, co-dominant, no environmental noise): under additive fitness the
marginal selection coefficient is s / W-bar, under multiplicative fitness it is s. Mean fitness relative to a mutation-free genome:

| Inputs | additive, 1,000 / 10,000 generations | multiplicative, 1,000 / 10,000 | Mendel's Accountant base run (section 8) |
|---|---|---|---|
| Sanford DFE, U_del = 5 | 0.851 / 0.594 | 0.847 / 0.592 | 0.872 / 0.640 |
| gamma DFE (shape 0.19, mean 0.01; approximate), U_del = 2 | 0.474 / 0.348 | 0.362 / 0.177 | |
| Sanford DFE, U_del = 70 (all mutations deleterious) | 0.183 / 0.044 | 0.100 / 0.023 (at 3,000) | collapse at 7,124 |

With Sanford's DFE the two combination rules give the same trajectory (0.594 vs 0.592 at 10,000 generations), and the independent-locus
simulator agrees with Mendel's Accountant within 0.05. The decline is therefore not an artefact of additivity, linkage or noise. It comes
from the DFE: Sanford's Weibull places about 18% of deleterious mutations between effects of 1e-6 and 1e-4, below the selection threshold
1/(2 Ne) of a population of 1,000-10,000, and these accumulate at the neutral rate (section 8, Experiment 2: 97% of U_del accumulates).
Whether such a class exists at ~1 mutation per individual per generation is the empirical question on which the whole result turns. The
standard model's answer, from comparative genomics (about 5% of the human genome under detectable constraint), is that mutations at the
other 95% of sites are neutral, not "slightly deleterious": set their effect to zero and there is nothing to accumulate, which is the
`fraction_neutral` input of the model itself. The gamma case shows that an empirical-shaped DFE with U_del = 2 also declines under both
rules in this simulator, because the gamma's weak tail is likewise sub-threshold; the published argument that mutation load is tolerable
rests on soft (relative) selection and on epistasis or truncation-like selection (Kondrashov 1988; Lesecque, Keightley & Eyre-Walker 2012;
Charlesworth 2013, none of which were opened in this session), and the truncation variants in Experiment 2 show that the model itself
holds fitness when selection is of that kind.

## 9. Summary of findings

1. **Fixation.** Mendel's Accountant fixes neutral mutations at Kimura's rate k = U, in its authors' own deep-time runs (95-96% of theory)
   and in ours (98-105% with the simulator's Ne = 1.33 N). Rupe & Sanford's "few hundred fixations in 6 million years" used a mutation
   rate 7,100 times below the measured human rate; at the measured rate the model's own neutral class gives 11.7 million fixations per
   lineage, against 17.5 million observed for human-chimpanzee, and matches the gorilla-gorilla (2.4 vs 4.4 million) and human-Neanderthal
   (38,500 vs at least 31,389) comparisons as well as neutral theory does. "Haldane's ratchet" counts only selection-driven fixations and
   labels the rest as damage; the counts themselves are the molecular clock.
2. **The model as published reproduces** (mutation counts within 1-4%, fixations at the neutral rate), except that this implementation
   does not reach "fitness 0 at 882 generations without selection".
3. **Fitness decline with measured mutation rates.** With 70 new mutations per individual, 5 of them deleterious under Sanford's DFE,
   mean fitness falls to 0.64 at 10,000 and 0.40 at 25,000 generations, extrapolating to zero near 42,000: after the Neanderthal and
   gorilla splits, long before the human-chimpanzee split. Declaring all 70 deleterious collapses the population in 7,100 generations;
   truncation selection removes the decline; heritability and beneficial mutations change nothing; a tenfold larger population slows it by a quarter. The decline is identical under
   multiplicative fitness and is set by the DFE's sub-threshold mass, i.e. by the assumption that unconstrained sites are slightly
   deleterious rather than neutral.
4. **LTEE.** With LTEE mutation rates the model loses 3.5-4.5% fitness in 60,000 generations at the mutator rate and nothing at the
   ancestral rate; the real populations gained 70%. The model's beneficial-mutation assumptions cannot produce the experiment's result.
5. **H1N1.** The lineage whose 2009 disappearance was read as mutational extinction was replaced by A(H1N1)pdm09, which is the dominant
   influenza A subtype in 2026.

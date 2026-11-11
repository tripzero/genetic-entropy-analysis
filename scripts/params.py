"""Parameters for the genetic-entropy (Mendel's Accountant) test. Every number carries its source; 'MITTENS' means the value
was verified against the primary text in the companion analysis /home/kev/mittens_analysis (README section 5) and reused.
Ranges are (central, low, high)."""

# ---------------------------------------------------------------- the model as published (primary texts in downloads/sources/)
MA = dict(
    # Rupe & Sanford 2013, Proc. 7th Int. Conf. Creationism (downloads/sources/rupe_sanford_2013_icc.txt), Table 1:
    # 300,000 generations, non-deleterious mutation rate 0.01 per individual per generation (gametic rate 0.005), probability selection.
    rs2013_U_indiv=0.01, rs2013_gens=300_000,
    rs2013_table1=[  # run id, pop size, description, mutations/individual, total fixations, selective fixations
        ('qrtsyg', 1000, 'all neutral', 2870, 1413, 0),
        ('qrtqaq', 1000, 'default beneficial distribution', 3048, 1492, 79),
        ('qrtimc', 1000, 'extended beneficial distribution (max effect 1.0)', 3413, 1684, 271),
        ('qrtmoi', 10000, 'all neutral', 2997, 1234, 0),
        ('qrtsvs', 10000, 'default beneficial distribution', 3469, 1456, 222),
        ('qrtniw', 10000, 'extended beneficial distribution (max effect 1.0)', 4190, 1824, 590),
        ('qrtoju', 10000, 'optimised selection settings (partial truncation, max effect 1.0, min 3.3e-6)', 5977, 2710, 1476)],
    # Rupe & Sanford 2013 Figure 3: pop 1,000, 20,000 generations, 10 non-neutral mutations per individual (9.99 deleterious : 0.01 beneficial),
    # 50% selective elimination: 6,775 deleterious and 10 beneficial fixations.
    rs2013_fig3=dict(N=1000, gens=20000, U_indiv=10.0, del_fix=6775, ben_fix=10),
    # Rupe & Sanford 2013 Figure 1: 1 neutral mutation per individual, 13,500 generations; pop 10,000 'not even approaching fixation rate equilibrium'.
    # Rupe & Sanford 2013, Methods: "Deleterious mutation accumulation could not be simulated through such deep time, as it always either caused
    # overflow of memory or premature termination due to extinction."
    # Sanford et al. 2008, Proc. 6th ICC (downloads/sources/sanford_2008_icc6.txt): default human parameters, pop 1,000, 10 mutations per individual,
    # 6 offspring per female (4 removed), heritability 0.2, 9:1 deleterious:beneficial. Fig. 1b caption: "Without selection this same run reaches
    # a fitness of 0.0 (extinction) in 882 generations." Fig. 5: after 10,000 generations 16,116 deleterious and 3,193 beneficial alleles fixed.
    s2008_no_selection_extinction_gen=882,
    s2008_fig5=dict(N=1000, gens=10000, U_indiv=10.0, del_fix=16116, ben_fix=3193, frac_fav=0.1),
    # Sanford & Nelson 2012 (InTech chapter, downloads/sources/sanford_nelson_2012_intech.txt): Fig. 4, default settings, 5,000 generations:
    # mean deleterious mutation count per individual 47,648, beneficial 9.8. Fig. 9: clonal run extinct 'in just over 5000 generations'.
    sn2012_fig4=dict(N=1000, gens=5000, U_indiv=10.0, del_per_indiv=47648),
    # Weibull DFE (Sanford et al. 2007 SCPE eq. 3.1; mendel-go dna/mutation.go): |s| = exp(-ln(G) x^gamma), x~U(0,1).
    # 2007 defaults theta=0.1, q=0.001 -> mean 5.06e-4, median 2.7e-8 (Sanford et al. 2007 section 4.2.2). mendel-go defaults theta=0.01, q=0.01.
)

# ---------------------------------------------------------------- observed data and standard-model inputs
HUMAN_CHIMP = dict(
    label='Human lineage vs chimpanzee',
    F_obs=(17.5e6, 15e6, 20e6),              # fixed SNV differences per lineage: ~35 M SNVs, 1.23% (Chimpanzee Consortium 2005) [MITTENS]
    T_years=(6.3e6, 5.5e6, 7.5e6),           # split time [MITTENS; Yoo et al. 2025 5.5-6.3]
    gen_years=(25., 20., 29.),               # [MITTENS]
    Ne_anc=(65000, 50000, 100000),           # ancestral effective size, adds 2Ne generations of coalescence [MITTENS]
    U_hap=(1.25e-8*2.85e9, 1.2e-8*2.7e9, 1.29e-8*3.0e9),  # new SNVs per haploid genome per generation (Kong 2012; Jonsson 2017; aligned 2.85 Gb) [MITTENS]
    f=(0.9, 0.85, 0.95),                     # fraction of new mutations not removed by purifying selection [MITTENS]
    N_modern=10000,                          # Ne used for the fixation transient (irrelevant at 250k+ generations)
)
GORILLA = dict(
    label='Eastern vs western gorilla',
    # Scally et al. 2012 Nature (downloads/sources/scally_2012_gorilla.xml): "mean sequence divergence time between them of 1.75 Mya" at the
    # paper's 1e-9 /bp/yr, i.e. mean pairwise divergence 2 x 1.75e6 x 1e-9 = 0.35% per site; ML population split ~0.5 Myr with gene flow.
    # Pairwise divergence per lineage = d/2 x 2.9e9 aligned autosomal sites. Low end 0.23% from the MITTENS analysis range.
    F_obs=(0.0030/2*2.9e9, 0.0023/2*2.9e9, 0.0035/2*2.9e9),
    T_years=(5e5, 1.5e5, 1.75e6),            # population split: ~0.5 Myr (Scally 2012 ML), 150-180 kyr end of gene flow (Pawar 2023), 1.75 Myr mean sequence divergence
    gen_years=(19., 19., 19.),               # Langergraber et al. 2012 [MITTENS]
    Ne_anc=(30000, 20000, 40000),            # Prado-Martinez 2013 Table 1: western lowland Ne 26.8-53.5k, eastern 12.2-24.3k (theta_w based)
    U_hap=(1.13e-8*2.9e9, 1.0e-8*2.9e9, 1.235e-8*2.9e9),   # gorilla pedigree rate (Besenbacher 2019) [MITTENS]
    f=(0.9, 0.85, 0.95),
    N_modern=30000,
)
HUMAN_NEA = dict(
    label='Modern humans vs Neanderthal',
    # Prufer et al. 2014 Nature 505:43: 31,389 single-nucleotide changes fixed derived in 1,094 present-day humans and ancestral in the
    # Altai Neanderthal and Denisovan genomes; a lower bound (callable fraction, filtering) [MITTENS, verified there]
    F_obs=(31389, 31389, 47083),
    T_years=(650e3, 550e3, 765e3),           # population split [MITTENS]
    gen_years=(27., 25., 29.),
    Ne_anc=(10000, 7000, 15000),             # modern-human Ne governs the fixation transient here (fixed derived alleles, not pairwise divergence)
    U_hap=(1.25e-8*2.85e9, 1.2e-8*2.7e9, 1.29e-8*3.0e9),
    f=(0.9, 0.85, 0.95),
    N_modern=10000,
)
# Present-day within-species diversity (for the discussion): Prado-Martinez et al. 2013 Nature 499:471 (PMC3822165, read 2026-10-30):
# non-African humans, eastern lowland gorillas, bonobos, western chimpanzees ~0.8e-3 heterozygotes/bp; central chimpanzees, western lowland
# gorillas, orangutans 1.6-2.4e-3. Scally 2012 Table 2: western lowland gorilla heterozygosity 0.178-0.189%, eastern lowland 0.076%.
# 1000 Genomes (Auton et al. 2015, PMC4750478): a typical genome differs from the reference at 4.1-5.0 million sites.

# Gibson, Baumgardner, Brewer & Sanford 2013, Biological Information: New Perspectives, pp. 232-263 (read 2026-10-30 in the publisher's
# open-access reader). Fixed settings: 6 offspring per female (4 removed), Weibull DFE with 0.1% of mutations > 0.1 (mean -5.4e-4, median
# -1.4e-7), no beneficials, co-dominant, additive, population 10,000 unless stated.
GIBSON2013 = dict(
    default_U5=dict(U_del=5.0, heritability=0.4, selection='partial truncation', fitness_10k={1000: 0.79, 10000: 0.84},
                    extinction_gen=40831, mutations_per_indiv_at_extinction=174890,  # "fitness declined to near extinction in 40,831 generations"
                    accumulation_per_gen=4.5, note='16% fitness decline in the first 10,000 generations; linear afterwards'),
    lynch_dfe_U3=dict(U_del=3.0, N=10000, heritability=0.4, dfe='27% of effects > 1e-2, minimum 1e-6 (genome size set to 1e6)',
                      initial_decline_per_gen=0.05, fitness_100k_partial_trunc=0.31, extinction_gen_probability_selection=83647,
                      offspring12_fitness_10k={'probability': 0.22, 'partial truncation': 0.57}),
)
# Lynch 2016 Genetics 202:869 (PMC4788123): extrapolating mouse mutation-accumulation lines (mutation rate 6.7x human) gives a human fitness
# decline of ~0.5% per generation in fecundity and ~0.1% per generation in other traits IF selection were completely relaxed.
LYNCH2016 = dict(relaxed_decline_per_gen=(0.001, 0.005))

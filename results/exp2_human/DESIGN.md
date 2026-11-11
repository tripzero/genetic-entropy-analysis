# Experiment 2: human-realistic inputs, one assumption varied at a time
Base: N = 1,000; 70 new mutations per individual per generation (de novo SNVs + indels, Jonsson 2017 / Kong 2012); fraction_neutral 0.92857 so that
U_del = 5 per individual (Gibson et al. 2013's "conservative" functional fraction); Sanford default Weibull DFE (0.1% of mutations with effect > 0.1,
mean 5e-4, min 3e-10); beneficial fraction 1e-4 of all mutations (0.007 per individual) with max gain 0.01; 50% recessive (hetero expression 0.1/0.9);
6 offspring per female (reproductive_rate 3), probability selection (spps), heritability 0.2; additive fitness (the only mode in mendel-go);
tracking_threshold 1e-3 (fitness exact; allele bins not used). 25,000 generations for the base (the Neanderthal split), 10,000 for variants.
Variants: dfe_all_deleterious (fraction_neutral 0: all 70 mutations deleterious, the literal genetic-entropy premise); dfe_empirical_Udel2
(U_del = 2, Weibull shifted to high impact as Gibson et al. did for Lynch's DFE: genome_size 1e6 -> min effect 1e-6, 27% of effects > 0.01);
sel_fulltrunc; sel_partialtrunc (0.5); heritability_1; fav_1e-3 and fav_1e-2 (0.07 and 0.7 beneficial per individual per generation);
base_N10000 (N = 10,000, tracking_threshold 0.01).

# Experiment 3 (reduced): LTEE-like inputs in mendel-go
mendel-go has no haploid or clonal reproduction (recombination_model 3 only), so this is a DIPLOID population of 1,000 with no crossover
(crossover_model none, one chromosome, 100 linkage blocks: whole-haplotype linkage, Mendelian segregation still occurs), genome 4.63e6,
all mutations non-neutral with Sanford's Weibull DFE scaled to this genome (minimum effect 2.2e-7), co-dominant, probability selection,
heritability 0.2, 6 offspring per female, 60,000 generations. Mutation rates: mutator 0.041 per genome per generation (Wielgoss 2011 x 100)
and ancestral 4.1e-4; beneficial fraction 1e-4 (model default) or 1e-2 (generous), maximum beneficial gain 0.1 (early LTEE drivers).
Observed (Wiser, Ribeck & Lenski 2013): mean fitness up ~70% by 50,000 generations, mutators slightly higher; no meltdown.

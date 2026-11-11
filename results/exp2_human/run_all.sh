#!/bin/bash
cd /home/kev/genetic_entropy_analysis/results/exp2_human
for n in base_Udel5 dfe_all_deleterious dfe_empirical_Udel2 sel_fulltrunc sel_partialtrunc heritability_1 fav_1e-3 fav_1e-2 base_N10000; do
  ( /home/kev/genetic_entropy_analysis/downloads/mendel-go/mendel-go -f $n.ini -D /home/kev/genetic_entropy_analysis/downloads/mendel-go/mendel-defaults.ini > $n.log 2>&1; echo "$n exit $? $(date +%H:%M:%S)" >> done.txt ) &
done
wait
echo ALL_DONE >> done.txt

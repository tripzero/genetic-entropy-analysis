#!/bin/bash
cd /home/kev/genetic_entropy_analysis/results/exp3_ltee
for n in mutator_fav1e-4 mutator_fav1e-2 ancestral_fav1e-2; do
  ( /home/kev/genetic_entropy_analysis/downloads/mendel-go/mendel-go -f $n.ini -D /home/kev/genetic_entropy_analysis/downloads/mendel-go/mendel-defaults.ini > $n.log 2>&1; echo "$n exit $? $(date +%H:%M:%S)" >> done.txt ) &
done
wait; echo ALL_DONE >> done.txt

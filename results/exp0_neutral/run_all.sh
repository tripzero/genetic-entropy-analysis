#!/bin/bash
cd /home/kev/genetic_entropy_analysis/results/exp0_neutral
for n in N100 N300 N1000; do
  ( /home/kev/genetic_entropy_analysis/downloads/mendel-go/mendel-go -f $n.ini -D /home/kev/genetic_entropy_analysis/downloads/mendel-go/mendel-defaults.ini > $n.log 2>&1; echo "$n exit $?" >> done.txt ) &
done
wait
echo ALL_DONE >> done.txt

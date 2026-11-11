#!/bin/bash
cd /home/kev/genetic_entropy_analysis/results/exp1_sanford2008
for n in fig1_selection fig1_noselection fig5_10k; do
  ( /home/kev/genetic_entropy_analysis/downloads/mendel-go/mendel-go -f $n.ini -D /home/kev/genetic_entropy_analysis/downloads/mendel-go/mendel-defaults.ini > $n.log 2>&1; echo "$n exit $? $(date +%H:%M:%S)" >> done.txt ) &
done
wait
echo ALL_DONE >> done.txt

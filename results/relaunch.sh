#!/bin/bash
cd /home/kev/genetic_entropy_analysis/results/exp1_sanford2008
( /home/kev/genetic_entropy_analysis/downloads/mendel-go/mendel-go -f fig1_noselection_rr1.ini -D /home/kev/genetic_entropy_analysis/downloads/mendel-go/mendel-defaults.ini > fig1_noselection_rr1.log 2>&1; echo "fig1_noselection_rr1 exit $? $(date +%H:%M:%S)" >> done.txt ) &
/home/kev/genetic_entropy_analysis/results/exp2_human/run_all.sh
wait

#!/bin/bash
filepath=$1
sim=$2
scoef=$3
Nw=$4
tick=$5
outpath=$6

for cutoff in 0.01 0.05 0.1 0.2 0.3 0.4 0.5; do
    echo "$filepath" "$sim" "$scoef" "$Nw" "$tick" "$cutoff" $( Rscript power_analysis.R --infile $filepath --cutoff $cutoff) >> $outpath
done
#mv "$filepath" "{$filepath%.txt}_processed.txt"

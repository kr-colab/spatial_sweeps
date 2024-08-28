#!/bin/bash

step=2000000


for chr in 2R 3R; do
    header=`tabix -H ag1000g.agam.n1470.merged_variants.$chr\.ann.vcf.gz | grep "##contig=<ID=$chr,length="`
    length=`echo $header | awk '{sub(/.*=/,"");sub(/>/,"");print}'` 

    # subset by window and run scan
    endwindow=$step
    for startwindow in `seq 1 $step $length`; do
        echo "processing $startwindow to $endwindow"
        sbatch scripts/windowed_genome_scan.batch $chr $startwindow $endwindow
        endwindow=$((endwindow+step))
    done
done

#!/bin/bash
chr=$1
startwindow=$2
stopwindow=$3

# subset chromosome
tabix -h ag1000g.agam.n1470.merged_variants.$chr\.ann.vcf.gz \
    $chr\:$startwindow\-$stopwindow > data/vcf/$chr\.$startwindow\-$stopwindow.ann.vcf

python scripts/spatial_genome_scan.py --vcf data/vcf/$chr\.$startwindow\-$stopwindow.ann.vcf --sample_data data/ag1000g_v3_gambiae.txt --out out/$chr\.$startwindow\-$stopwindow



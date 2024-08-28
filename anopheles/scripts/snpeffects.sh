#!/bin/bash
# run SNPEff and process vcfs
chrom=$1
 
#if [ ! -f data/vcf/ag1000g.agam.n1470.filtered_variants.$chrom\.ann.vcf ]; then
#   # run SNPEff
#   java -Xmx8g -jar ~/snpEff/snpEff.jar Anopheles_gambiae \
#   data/vcf/ag1000g.agam.n1470.filtered_variants.$chrom\.vcf.gz > data/vcf/ag1000g.agam.n1470.filtered_variants.$chrom\.ann.vcf
#fi
# bgzip

#if [ ! -f data/vcf/ag1000g.agam.n1470.filtered_variants.$chrom\.ann.vcf.gz ]; then
bgzip -c -i data/vcf/ag1000g.agam.n1470.filtered_variants.$chrom\.ann.vcf > data/vcf/ag1000g.agam.n1470.filtered_variants.$chrom\.ann.vcf.gz
#fi

# tabix

#if [ ! -f data/vcf/ag1000g.agam.n1470.filtered_variants.$chrom\.ann.vcf.gz.tbi ]; then
tabix -f -p vcf data/vcf/ag1000g.agam.n1470.filtered_variants.$chrom\.ann.vcf.gz
#fi

#!/bin/bash

# filter anopheles VCF

basepath=/sietch_colab/data_share/agp3_trees/vcf/unphased_vcf/gamb

chromosome=$1

zgrep -e PASS -e '#' $basepath\/$chromosome\_sitefilters.vcf.gz > $chromosome\_pass_positions.vcf

bgzip $chromosome\_pass_positions.vcf

tabix -p vcf $chromosome\_pass_positions.vcf.gz

bcftools isec -c none -p $chromosome\_filtered -n=2 -w1 $basepath\/ag1000g.agam.n1470.merged_variants.$chromosome\.vcf.gz $chromosome\_pass_positions.vcf.gz

mv $chromosome\_filtered/0000.vcf ag1000g.agam.n1470.filtered_variants.$chromosome\.vcf

bgzip ag1000g.agam.n1470.filtered_variants.$chromosome\.vcf

tabix -p vcf ag1000g.agam.n1470.filtered_variants.$chromosome\.vcf.gz

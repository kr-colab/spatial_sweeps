#!/bin/bash

vcfpath="/sietch_colab/crehmann/filtered_ag3/"

for chr in 2L 2R 3L 3R; do 
    Rscript scripts/popPCA.R $vcfpath$chr\.bed
done
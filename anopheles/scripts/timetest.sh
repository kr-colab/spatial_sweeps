#!/bin/bash

echo $( date )
python scripts/spatial_genome_scan.py --vcf data/X.0-10000.vcf --sample_data data/ag1000g_v3_gambiae.txt --out test
echo $( date )

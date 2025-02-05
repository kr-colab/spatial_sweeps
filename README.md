# Sweeps in space!

This repo contains code for "Sweeps in space: leveraging geographic data to identify
beneficial alleles in *Anopheles gambiae*".

The main simulation used can be found at `simulation/scripts/spatial_sweep.slim`; functions for processing tree sequences and calculating per-variant frequency and area can be found in `simulation/scripts/frequency_area.py`.

Spatial genome scans on *Anopheles gambiae* data were carried out using `anopheles/scripts/anopheles_SNP_genome_scan.py` and analyzed using `anopheles/anopheles-analysis.Rmd`, which can be visualized and browsed at https://kr-colab.github.io/spatial_sweeps/anopheles/anopheles-analysis.html.

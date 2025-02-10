# Sweeps in space!

This repo contains code for ["Sweeps in space: leveraging geographic data to identify
beneficial alleles in *Anopheles gambiae*"](https://www.biorxiv.org/content/10.1101/2025.02.07.637123v1).

------------------

## Simulation code

The main simulation used can be found at `simulation/scripts/spatial_sweep.slim`; functions for processing tree sequences and calculating per-variant frequency and area can be found in `simulation/scripts/frequency_area.py`.

## Analysis

Spatial genome scans on *Anopheles gambiae* data were carried out using `anopheles/scripts/anopheles_SNP_genome_scan.py` and analyzed using `anopheles/anopheles-analysis.Rmd`.

The entire analysis can be visualized and browsed [here](https://kr-colab.github.io/spatial_sweeps/anopheles/anopheles-analysis.html), including an interactive frequency-area plot of genome-wide WSF outliers and **searchable table** of all SF outliers identified in our analysis.

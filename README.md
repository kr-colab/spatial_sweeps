# Sweeps in space!

This repo contains code for ["Sweeps in space: leveraging geographic data to identify
beneficial alleles in *Anopheles gambiae*"](https://doi.org/10.1093/molbev/msaf141).

## Spatial-frequency genome scan

To run your own spatial-frequency genome scan, use the `spacefreq` library found in this repository.

### Installation

Requires Python ≥ 3.12 and [uv](https://docs.astral.sh/uv/).

```bash
cd spacefreq
uv sync
```

### Required data

- **Genotypes**: VCF or VCF.gz format (bgzipped + tabix-indexed for region queries)
- **Metadata**: tab-separated file with columns `sampleID`, `x` (longitude), `y` (latitude)

```
sampleID    x        y
AG0058-C    -14.917  13.567
AV0043-C    -9.53    8.48
```

### Command-line interface

```bash
uv run spacefreq scan \
    --genotypes path/to/genotypes.vcf.gz \
    --metadata  path/to/metadata.txt \
    --out       path/to/output.txt
```

Optional arguments:

| Flag | Default | Description |
|------|---------|-------------|
| `--start INT` | — | Genomic start position (inclusive) |
| `--stop INT` | — | Genomic stop position (inclusive) |
| `--chromosome STR` | — | Chromosome to scan |
| `--min-locs INT` | 3 | Minimum unique sampling locations to record area |
| `--transect FLOAT` | 1.0 | Transect width (km) for two-location alleles |
| `--sample-area FLOAT` | 1.0 | Default area (km²) for single-location alleles |
| `--no-count-missing-as-ancestral` | — | Exclude missing genotypes from frequency denominator |

The output is a tab-separated file with columns `position`, `alternate`, `frequency`, and `area` (km²), one row per observed alternate allele.

### Python API

```python
import pandas as pd
from spacefreq import scan
from spacefreq.genotypes import read_vcf

genotypes, positions, sample_ids = read_vcf("genotypes.vcf.gz", region="2L:1-500000")
metadata = pd.read_csv("metadata.txt", sep="\t")

result = scan(genotypes, metadata, positions)
print(result.head())
#    position  alternate  frequency        area
# 0      1249          1   0.000374   45.21
# 1      1829          1   0.002381  312.84
```

When the VCF contains more samples than the metadata (e.g. a population subset), `scan()` automatically subsets genotypes to the metadata samples for carrier/area detection while computing frequency globally across all samples.

------------------

## Simulation code

The main simulation used can be found at `simulation/scripts/spatial_sweep.slim`; functions for processing tree sequences and calculating per-variant frequency and area can be found in `simulation/scripts/frequency_area.py`.

## Analysis

Spatial genome scans on *Anopheles gambiae* data were carried out using `anopheles/scripts/anopheles_SNP_genome_scan.py` and analyzed using `anopheles/anopheles-analysis.Rmd`.

The entire analysis can be visualized and browsed [here](https://kr-colab.github.io/spatial_sweeps/anopheles/anopheles-analysis.html), including an interactive frequency-area plot of genome-wide WSF outliers and **searchable table** of all SF outliers identified in our analysis.

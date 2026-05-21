"""Top-level scan() orchestrator."""
from __future__ import annotations

import math
import warnings

import numpy as np
import pandas as pd
import allel

from spacefreq.freq import calculate_frequency
from spacefreq.space import get_carrier_locations, calculate_area


def scan(
    genotypes: allel.GenotypeArray,
    metadata: pd.DataFrame,
    positions: np.ndarray,
    filter: str | None = None,
    start: int | None = None,
    stop: int | None = None,
    chromosome: str | None = None,
    min_locs: int = 3,
    transect: float = 1.0,
    sample_area: float = 1.0,
    count_missing_as_ancestral: bool = True,
    sample_col: str = "sampleID",
    lon_col: str = "x",
    lat_col: str = "y",
) -> pd.DataFrame:
    """Perform a spatial genome scan.

    For every polymorphic variant, records the global frequency and the
    landscape area (convex hull in km²) occupied by carriers of each
    alternate allele.

    Parameters
    ----------
    genotypes:
        GenotypeArray of shape (n_variants, n_samples, ploidy).
    metadata:
        DataFrame containing sample IDs and coordinates.  Column names are
        mapped to the canonical names (``sampleID``, ``x``, ``y``) using the
        *sample_col*, *lon_col*, and *lat_col* parameters.
    positions:
        1-D integer array of genomic positions, one per variant.
    filter:
        Path to a genotype filter file in VCF format (optional, reserved).
    start:
        Inclusive start position for the genomic region to scan.
    stop:
        Inclusive stop position for the genomic region to scan.
    chromosome:
        Chromosome label (informational; not used for filtering here since
        *positions* is already chromosome-specific).
    min_locs:
        Minimum number of unique sampling locations required to record an
        area for a given allele.
    transect:
        Transect width (km) used when exactly two unique locations are found.
    sample_area:
        Default area (km²) for a single sampling location.
    count_missing_as_ancestral:
        Passed through to :func:`spacefreq.freq.calculate_frequency`.
    sample_col:
        Name of the metadata column that holds sample IDs (default ``"sampleID"``).
        Use ``"sample_id"`` for Ag1000G-style metadata.
    lon_col:
        Name of the metadata column that holds longitude (default ``"x"``).
        Use ``"longitude"`` for Ag1000G-style metadata.
    lat_col:
        Name of the metadata column that holds latitude (default ``"y"``).
        Use ``"latitude"`` for Ag1000G-style metadata.

    Returns
    -------
    pd.DataFrame with columns: ``position``, ``alternate``, ``frequency``, ``area``.
    """
    # ------------------------------------------------------------------
    # 0. Normalise metadata column names to canonical form.
    # ------------------------------------------------------------------
    rename_map = {}
    if sample_col != "sampleID" and sample_col in metadata.columns:
        rename_map[sample_col] = "sampleID"
    if lon_col != "x" and lon_col in metadata.columns:
        rename_map[lon_col] = "x"
    if lat_col != "y" and lat_col in metadata.columns:
        rename_map[lat_col] = "y"
    if rename_map:
        metadata = metadata.rename(columns=rename_map)

    # Capture attributes BEFORE any slicing; numpy slices do not preserve
    # custom attributes on ndarray subclasses.
    embedded_ids: list[str] | None = getattr(genotypes, "sample_ids", None)
    chromosomes: np.ndarray | None = getattr(genotypes, "chromosomes", None)

    # ------------------------------------------------------------------
    # 1. Positional filtering
    # ------------------------------------------------------------------
    pos_mask = np.ones(len(positions), dtype=bool)
    if start is not None:
        pos_mask &= positions >= start
    if stop is not None:
        pos_mask &= positions <= stop

    genotypes = genotypes[pos_mask]
    positions = positions[pos_mask]
    if chromosomes is not None:
        chromosomes = chromosomes[pos_mask]

    # ------------------------------------------------------------------
    # 2. Keep only polymorphic variants (at least one alt allele observed)
    #    using ALL samples (for global frequency).
    # ------------------------------------------------------------------
    ac_all = genotypes.count_alleles()
    poly_mask = np.sum(ac_all[:, 1:], axis=1) > 1
    genotypes = genotypes[poly_mask]
    positions = positions[poly_mask]
    ac_all = ac_all[poly_mask]
    if chromosomes is not None:
        chromosomes = chromosomes[poly_mask]

    # ------------------------------------------------------------------
    # 3. Align genotypes to metadata samples.
    #
    #    Frequency is always calculated from ALL samples in the genotype
    #    array (global frequency).  Carrier detection for area calculation
    #    uses only the samples present in metadata.
    #
    #    Alignment strategy (in priority order):
    #      a) If genotypes carries sample_ids (AnnotatedGenotypeArray from
    #         read_vcf), intersect with metadata["sampleID"] and take the
    #         matching columns.  Emit a warning when samples are dropped.
    #      b) If len(metadata) == n_samples, assume rows are already aligned.
    #      c) Otherwise raise ValueError.
    # ------------------------------------------------------------------
    metadata = metadata.reset_index(drop=True)
    n_gt_samples = genotypes.shape[1]

    if embedded_ids is not None and "sampleID" in metadata.columns:
        sid_to_col = {sid: i for i, sid in enumerate(embedded_ids)}
        meta_in_vcf = metadata["sampleID"].isin(sid_to_col)
        if not meta_in_vcf.all():
            missing = list(metadata.loc[~meta_in_vcf, "sampleID"])
            warnings.warn(
                f"{len(missing)} metadata sample(s) not found in the genotype "
                f"file and will be excluded from area calculations: {missing[:5]}"
                f"{'...' if len(missing) > 5 else ''}",
                UserWarning,
                stacklevel=2,
            )
            metadata = metadata[meta_in_vcf].reset_index(drop=True)

        if n_gt_samples != len(metadata):
            warnings.warn(
                f"Genotype file contains {n_gt_samples} samples but metadata "
                f"has {len(metadata)} entries. Subsetting genotypes to the "
                f"{len(metadata)} samples present in metadata for carrier "
                "detection; frequency is still computed globally.",
                UserWarning,
                stacklevel=2,
            )

        col_indices = np.array(
            [sid_to_col[sid] for sid in metadata["sampleID"]], dtype=int
        )
        gt_for_area = genotypes.take(col_indices, axis=1)

    elif n_gt_samples == len(metadata):
        gt_for_area = genotypes

    else:
        raise ValueError(
            f"Genotype array has {n_gt_samples} samples but metadata has "
            f"{len(metadata)} rows and no sample_ids are embedded in the "
            "genotype array.  Pass genotypes from read_vcf() (which embeds "
            "sample IDs) so scan() can align automatically, or pre-align "
            "genotypes and metadata before calling scan()."
        )

    # ------------------------------------------------------------------
    # 4. Calculate per-allele frequency array using ALL samples (global).
    # ------------------------------------------------------------------
    freq_array = calculate_frequency(genotypes, count_missing_as_ancestral)

    # ------------------------------------------------------------------
    # 5. Main loop: for each variant × alt allele, compute freq + area.
    # ------------------------------------------------------------------
    rows = []
    n_variants = len(positions)
    max_alt = ac_all.shape[1] - 1

    for vi in range(n_variants):
        calls_for_area = np.asarray(gt_for_area[vi])
        allele_counts = ac_all[vi]

        for alt in range(1, min(4, max_alt + 1)):  # alleles 1, 2, 3
            if alt >= len(allele_counts) or allele_counts[alt] == 0:
                continue

            freq = float(freq_array[vi, alt])
            locs = get_carrier_locations(calls_for_area, metadata, allele=alt)
            area = calculate_area(locs, min_locs=min_locs, transect=transect,
                                  sample_area=sample_area)

            row: dict = {
                "position": int(positions[vi]),
                "alternate": alt,
                "frequency": freq,
                "area": area if not math.isnan(area) else float("nan"),
            }
            if chromosomes is not None:
                row["chromosome"] = chromosomes[vi]
            rows.append(row)

    if not rows:
        cols = ["position", "alternate", "frequency", "area"]
        if chromosomes is not None:
            cols = ["chromosome"] + cols
        return pd.DataFrame(columns=cols)

    df = pd.DataFrame(rows)
    # Drop rows where area could not be calculated (below min_locs).
    df = df.dropna(subset=["area"]).reset_index(drop=True)
    # Put chromosome first when present.
    if "chromosome" in df.columns:
        cols = ["chromosome", "position", "alternate", "frequency", "area"]
        df = df[cols]
    return df

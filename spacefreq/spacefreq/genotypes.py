"""Reading and parsing genotype files (VCF / ZARR)."""
from __future__ import annotations

import numpy as np
import allel


def read_vcf(
    path: str,
    region: str | None = None,
) -> tuple[allel.GenotypeArray, np.ndarray, list[str]]:
    """Read a VCF (or bgzipped VCF) and return genotypes, positions, and sample IDs.

    The returned :class:`allel.GenotypeArray` has a ``sample_ids`` attribute
    attached so that :func:`spacefreq.scan.scan` can automatically align
    samples to a metadata table without the caller needing to pass sample IDs
    separately.

    Parameters
    ----------
    path:
        Path to the VCF or VCF.gz file.
    region:
        Optional genomic region string in the form ``"chrom:start-stop"``
        (1-based, inclusive) passed directly to ``allel.read_vcf``.

    Returns
    -------
    tuple of:
        genotypes : allel.GenotypeArray, shape (n_variants, n_samples, ploidy)
                    with a ``sample_ids`` attribute carrying the VCF sample names.
        positions : np.ndarray of int, shape (n_variants,)
        sample_ids : list[str]
    """
    kwargs: dict = dict(
        fields=["samples", "calldata/GT", "variants/POS", "variants/CHROM"],
    )
    if region is not None:
        kwargs["region"] = region

    callset = allel.read_vcf(path, **kwargs)

    sample_ids = list(callset["samples"])
    genotypes = allel.GenotypeArray(callset["calldata/GT"])
    # Attach sample IDs and chromosomes directly; GenotypeArray (an ndarray
    # subclass) allows arbitrary attribute assignment.
    genotypes.sample_ids = sample_ids
    genotypes.chromosomes = callset["variants/CHROM"]
    positions = callset["variants/POS"]

    return genotypes, positions, sample_ids

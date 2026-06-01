"""Operations to calculate variant allele frequency."""
from __future__ import annotations

import numpy as np
import allel


def calculate_frequency(
    genotypes: allel.GenotypeArray,
    count_missing_as_ancestral: bool = True,
) -> np.ndarray:
    """Calculate alternate allele frequencies for every variant.

    Parameters
    ----------
    genotypes:
        GenotypeArray of shape (n_variants, n_samples, ploidy).
    count_missing_as_ancestral:
        If True (default), missing calls (-1) are treated as homozygous
        reference and the denominator is ``ploidy × n_samples`` for every
        site.  If False, missing samples are removed from the calculation
        and the denominator is ``ploidy × n_non_missing`` per site.

    Returns
    -------
    np.ndarray of shape (n_variants, n_alleles) containing frequencies.
    Allele index 0 is the reference; index 1, 2, … are alternates.
    """
    n_variants, n_samples, ploidy = genotypes.shape

    # Always report at least alleles 0 (ref) and 1 (first alt) so callers
    # can safely index freq[:, 1] even for monomorphic sites.
    max_allele = max(1, int(genotypes.count_alleles().shape[1]) - 1)
    # Re-compute with explicit max_allele to get a consistently-sized array.
    ac = genotypes.count_alleles(max_allele=max_allele)

    if count_missing_as_ancestral:
        denominator = n_samples * ploidy
        freq = ac / denominator
    else:
        n_missing = genotypes.count_missing(axis=1)  # shape (n_variants,)
        denominator = (n_samples - n_missing) * ploidy  # shape (n_variants,)
        # Avoid division by zero for fully-missing sites
        denominator = np.where(denominator == 0, 1, denominator)
        freq = ac / denominator[:, np.newaxis]

    return np.asarray(freq, dtype=float)

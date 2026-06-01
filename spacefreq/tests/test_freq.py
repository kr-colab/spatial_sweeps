"""Tests for spacefreq.freq: allele frequency calculations."""
import numpy as np
import pytest
import allel

from spacefreq.freq import calculate_frequency


def _make_gt(data):
    """Build a GenotypeArray from a list-of-lists (variants × samples × ploidy)."""
    return allel.GenotypeArray(np.array(data, dtype="i1"))


# ---------------------------------------------------------------------------
# Basic correctness
# ---------------------------------------------------------------------------

class TestCalculateFrequencyBasic:
    def test_all_ref_homozygous(self):
        # 1 variant, 4 samples, all [0, 0] → alt-1 freq = 0.0
        gt = _make_gt([[[0, 0], [0, 0], [0, 0], [0, 0]]])
        freq = calculate_frequency(gt)
        assert freq.shape[0] == 1
        assert freq[0, 1] == pytest.approx(0.0)

    def test_all_alt_homozygous(self):
        # 1 variant, 4 samples, all [1, 1] → alt-1 freq = 1.0
        gt = _make_gt([[[1, 1], [1, 1], [1, 1], [1, 1]]])
        freq = calculate_frequency(gt)
        assert freq[0, 1] == pytest.approx(1.0)

    def test_half_heterozygous(self):
        # 1 variant, 4 samples: 2 het [0,1], 2 hom-ref [0,0] → freq = 0.25
        gt = _make_gt([[[0, 1], [0, 1], [0, 0], [0, 0]]])
        freq = calculate_frequency(gt)
        assert freq[0, 1] == pytest.approx(0.25)

    def test_all_heterozygous(self):
        # 1 variant, 4 samples all [0, 1] → freq = 0.5
        gt = _make_gt([[[0, 1], [0, 1], [0, 1], [0, 1]]])
        freq = calculate_frequency(gt)
        assert freq[0, 1] == pytest.approx(0.5)

    def test_multiple_variants(self):
        gt = _make_gt([
            [[1, 1], [1, 1]],  # variant 0: all hom-alt → 1.0
            [[0, 0], [0, 0]],  # variant 1: all hom-ref → 0.0
            [[0, 1], [0, 1]],  # variant 2: all het → 0.5
        ])
        freq = calculate_frequency(gt)
        assert freq[0, 1] == pytest.approx(1.0)
        assert freq[1, 1] == pytest.approx(0.0)
        assert freq[2, 1] == pytest.approx(0.5)


# ---------------------------------------------------------------------------
# Missing data
# ---------------------------------------------------------------------------

class TestCalculateFrequencyMissingData:
    def test_missing_as_ancestral_default(self):
        # 1 variant, 4 samples: 2 het, 2 missing (-1)
        # count_missing_as_ancestral=True (default): denom = 2 × 4 = 8
        # alt count = 2 → freq = 2/8 = 0.25
        gt = _make_gt([[[0, 1], [0, 1], [-1, -1], [-1, -1]]])
        freq = calculate_frequency(gt, count_missing_as_ancestral=True)
        assert freq[0, 1] == pytest.approx(0.25)

    def test_missing_excluded(self):
        # Same genotypes but count_missing_as_ancestral=False: denom = 2 × 2 = 4
        # alt count = 2 → freq = 2/4 = 0.5
        gt = _make_gt([[[0, 1], [0, 1], [-1, -1], [-1, -1]]])
        freq = calculate_frequency(gt, count_missing_as_ancestral=False)
        assert freq[0, 1] == pytest.approx(0.5)

    def test_all_missing_as_ancestral_returns_zero(self):
        gt = _make_gt([[[- 1, -1], [-1, -1]]])
        freq = calculate_frequency(gt, count_missing_as_ancestral=True)
        assert freq[0, 1] == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# Range constraint: all frequencies must be in [0, 1]
# ---------------------------------------------------------------------------

class TestCalculateFrequencyRange:
    def test_frequencies_in_range(self):
        rng = np.random.default_rng(42)
        raw = rng.integers(0, 2, size=(10, 20, 2))
        gt = allel.GenotypeArray(raw.astype("i1"))
        freq = calculate_frequency(gt)
        valid = freq[:, 1:]
        assert np.all(valid >= 0.0)
        assert np.all(valid <= 1.0)

"""Tests for spacefreq.space: carrier location and area calculations."""
import math
import numpy as np
import pandas as pd
import pytest

from spacefreq.space import get_carrier_locations, calculate_area


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_metadata(sample_ids, xs, ys):
    return pd.DataFrame({"sampleID": sample_ids, "x": xs, "y": ys})


# ---------------------------------------------------------------------------
# get_carrier_locations
# ---------------------------------------------------------------------------

class TestGetCarrierLocations:
    def setup_method(self):
        # 5 samples at 3 distinct locations
        self.metadata = _make_metadata(
            ["s0", "s1", "s2", "s3", "s4"],
            [-10.0, -10.0, -9.0,  -9.0, -8.0],
            [  9.0,   9.0,  8.5,   8.5,  8.0],
        )

    def test_hom_ref_not_carrier(self):
        # All hom-ref → no carriers
        calls = np.array([[0, 0]] * 5, dtype="i1")
        locs = get_carrier_locations(calls, self.metadata, allele=1)
        assert len(locs) == 0

    def test_hom_alt_is_carrier(self):
        # samples 0 and 1 are hom-alt
        calls = np.array([[1, 1], [1, 1], [0, 0], [0, 0], [0, 0]], dtype="i1")
        locs = get_carrier_locations(calls, self.metadata, allele=1)
        assert len(locs) == 1  # s0 and s1 share the same location → 1 unique

    def test_het_is_carrier(self):
        # sample 4 is het [0, 1] → should appear as carrier
        calls = np.array([[0, 0], [0, 0], [0, 0], [0, 0], [0, 1]], dtype="i1")
        locs = get_carrier_locations(calls, self.metadata, allele=1)
        assert len(locs) == 1
        assert locs[0, 0] == pytest.approx(-8.0)
        assert locs[0, 1] == pytest.approx(8.0)

    def test_unique_locations_returned(self):
        # s0 and s1 share same location; both hom-alt → only 1 unique row
        calls = np.array([[1, 1], [1, 1], [0, 0], [0, 0], [0, 0]], dtype="i1")
        locs = get_carrier_locations(calls, self.metadata, allele=1)
        assert locs.shape == (1, 2)

    def test_multiple_carriers_multiple_locations(self):
        # s0 (loc A), s2 (loc B), s4 (loc C) all het
        calls = np.array([[0, 1], [0, 0], [0, 1], [0, 0], [0, 1]], dtype="i1")
        locs = get_carrier_locations(calls, self.metadata, allele=1)
        assert locs.shape == (3, 2)

    def test_multiallelic_other_alt_not_carrier(self):
        # At a tri-allelic site: sample 0 carries allele 2, not allele 1
        calls = np.array([[2, 0], [0, 0], [0, 0], [0, 0], [0, 0]], dtype="i1")
        locs = get_carrier_locations(calls, self.metadata, allele=1)
        assert len(locs) == 0

    def test_multiallelic_het_two_alts_counted_for_each(self):
        # sample 0 has [1, 2]: carrier of allele 1 AND allele 2
        calls = np.array([[1, 2], [0, 0], [0, 0], [0, 0], [0, 0]], dtype="i1")
        locs_1 = get_carrier_locations(calls, self.metadata, allele=1)
        locs_2 = get_carrier_locations(calls, self.metadata, allele=2)
        assert len(locs_1) == 1
        assert len(locs_2) == 1

    def test_multiallelic_non_carrier_excluded(self):
        # sample 0 has [1, 3]: carrier of 1 and 3, NOT 2
        calls = np.array([[1, 3], [0, 0], [0, 0], [0, 0], [0, 0]], dtype="i1")
        locs_2 = get_carrier_locations(calls, self.metadata, allele=2)
        assert len(locs_2) == 0


# ---------------------------------------------------------------------------
# calculate_area
# ---------------------------------------------------------------------------

class TestCalculateArea:
    def test_no_locations_returns_nan(self):
        coords = np.empty((0, 2))
        assert math.isnan(calculate_area(coords))

    def test_below_min_locs_returns_nan(self):
        # 2 locations but min_locs=3 → NaN
        coords = np.array([[-10.0, 9.0], [-9.0, 8.5]])
        assert math.isnan(calculate_area(coords, min_locs=3))

    def test_one_location_returns_sample_area(self):
        coords = np.array([[-10.0, 9.0]])
        area = calculate_area(coords, min_locs=1, sample_area=2.5)
        assert area == pytest.approx(2.5)

    def test_two_locations_returns_nonneg(self):
        coords = np.array([[-10.0, 9.0], [-9.0, 8.5]])
        area = calculate_area(coords, min_locs=2, transect=1.0)
        assert area >= 0.0

    def test_two_locations_transect_scales_area(self):
        coords = np.array([[-10.0, 9.0], [-9.0, 8.5]])
        area1 = calculate_area(coords, min_locs=2, transect=1.0)
        area2 = calculate_area(coords, min_locs=2, transect=2.0)
        assert area2 == pytest.approx(area1 * 2.0, rel=1e-6)

    def test_three_locations_returns_positive_area(self):
        # A triangle in West Africa
        coords = np.array([
            [-10.0,  9.0],
            [ -9.0,  9.0],
            [ -9.5,  8.0],
        ])
        area = calculate_area(coords, min_locs=3)
        assert area > 0.0

    def test_area_is_nonnegative(self):
        coords = np.array([
            [-14.9, 13.6],
            [ -9.5,  8.5],
            [-10.0,  9.3],
            [ -9.0,  8.0],
        ])
        area = calculate_area(coords, min_locs=3)
        assert area >= 0.0

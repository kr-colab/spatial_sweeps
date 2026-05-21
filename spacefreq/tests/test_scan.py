"""Integration tests for spacefreq.scan using real VCF + metadata."""
import numpy as np
import pandas as pd
import pytest
import allel

from spacefreq.genotypes import read_vcf
from spacefreq.scan import scan


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_gt(data, sample_ids=None):
    gt = allel.GenotypeArray(np.array(data, dtype="i1"))
    if sample_ids is not None:
        gt.sample_ids = sample_ids
    return gt


def _synthetic_inputs():
    """6 samples at 3 non-collinear locations; 2 variants that survive the polymorphism filter."""
    # samples: s0/s1 @ (-10, 9), s2/s3 @ (-9, 9), s4/s5 @ (-9.5, 8)
    # variant 0: s0, s2, s4 heterozygous → alt_count=3, carriers at 3 locations
    # variant 1: s4, s5 heterozygous → alt_count=2, carriers at 1 location (→ NaN area)
    raw = np.array([
        [[1, 0], [0, 0], [1, 0], [0, 0], [1, 0], [0, 0]],
        [[0, 0], [0, 0], [0, 0], [0, 0], [0, 1], [0, 1]],
    ], dtype="i1")
    gt = _make_gt(raw, sample_ids=["s0", "s1", "s2", "s3", "s4", "s5"])
    positions = np.array([100, 200])
    meta = pd.DataFrame({
        "sampleID": ["s0", "s1", "s2", "s3", "s4", "s5"],
        "x": [-10.0, -10.0,  -9.0,  -9.0, -9.5, -9.5],
        "y": [  9.0,   9.0,   9.0,   9.0,  8.0,  8.0],
    })
    return gt, positions, meta

VCF = "/sietch_colab/data_share/Ag1000G/Ag3.0/vcf/phased_vcf/gamb/gamb.2L.phased.n1470.derived.vcf.gz"
METADATA = "/home/crehmann/spatial_sweeps/anopheles/data/admixture_k1_metadata.txt"
REGION = "2L:1-500000"


@pytest.fixture(scope="module")
def scan_inputs():
    genotypes, positions, sample_ids = read_vcf(VCF, region=REGION)
    metadata = pd.read_csv(METADATA, sep="\t")
    return genotypes, positions, sample_ids, metadata


@pytest.fixture(scope="module")
def scan_result(scan_inputs):
    genotypes, positions, sample_ids, metadata = scan_inputs
    return scan(genotypes, metadata, positions)


class TestScanOutputShape:
    def test_returns_dataframe(self, scan_result):
        assert isinstance(scan_result, pd.DataFrame)

    def test_has_required_columns(self, scan_result):
        required = {"position", "alternate", "frequency", "area"}
        assert required.issubset(scan_result.columns)

    def test_nonempty(self, scan_result):
        assert len(scan_result) > 0


class TestScanFrequencies:
    def test_frequencies_above_zero(self, scan_result):
        assert (scan_result["frequency"] > 0.0).all()

    def test_frequencies_below_one(self, scan_result):
        assert (scan_result["frequency"] < 1.0).all()


class TestScanAreas:
    def test_areas_nonneg_or_nan(self, scan_result):
        areas = scan_result["area"]
        valid = areas.dropna()
        assert (valid >= 0.0).all()


class TestScanColumnMapping:
    """scan() must accept non-canonical metadata column names."""

    def test_alternate_col_names_give_same_result(self):
        gt, positions, meta_canonical = _synthetic_inputs()
        result_canonical = scan(gt, meta_canonical, positions, min_locs=3)

        meta_alt = meta_canonical.rename(columns={
            "sampleID": "sample_id",
            "x": "longitude",
            "y": "latitude",
        })
        result_alt = scan(
            gt, meta_alt, positions, min_locs=3,
            sample_col="sample_id", lon_col="longitude", lat_col="latitude",
        )

        pd.testing.assert_frame_equal(
            result_canonical.reset_index(drop=True),
            result_alt.reset_index(drop=True),
        )

    def test_canonical_columns_still_work_by_default(self):
        gt, positions, meta = _synthetic_inputs()
        result = scan(gt, meta, positions, min_locs=3)
        assert isinstance(result, pd.DataFrame)
        assert len(result) > 0

    def test_result_contains_area_for_three_location_variant(self):
        gt, positions, meta = _synthetic_inputs()
        result = scan(gt, meta, positions, min_locs=3)
        # variant 0 at position 100 has carriers at 3 locations → area should be present
        assert 100 in result["position"].values
        row = result[result["position"] == 100].iloc[0]
        assert row["area"] > 0.0


class TestScanPositionalFilter:
    def test_start_stop_restricts_positions(self, scan_inputs):
        genotypes, positions, sample_ids, metadata = scan_inputs
        start, stop = 100_000, 300_000
        result = scan(genotypes, metadata, positions, start=start, stop=stop)
        assert (result["position"] >= start).all()
        assert (result["position"] <= stop).all()

    def test_start_stop_fewer_rows_than_full(self, scan_inputs, scan_result):
        genotypes, positions, sample_ids, metadata = scan_inputs
        result_sub = scan(genotypes, metadata, positions, start=100_000, stop=300_000)
        assert len(result_sub) <= len(scan_result)

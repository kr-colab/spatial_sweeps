"""Integration tests for spacefreq.scan using real VCF + metadata."""
import numpy as np
import pandas as pd
import pytest
import allel

from spacefreq.genotypes import read_vcf
from spacefreq.scan import scan

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

"""CLI smoke tests for spacefreq via typer.testing.CliRunner."""
import pandas as pd
import pytest
from typer.testing import CliRunner

from spacefreq.cli import app
from spacefreq.genotypes import read_vcf
from spacefreq.scan import scan

VCF = "/sietch_colab/data_share/Ag1000G/Ag3.0/vcf/phased_vcf/gamb/gamb.2L.phased.n1470.derived.vcf.gz"
METADATA = "/home/crehmann/spatial_sweeps/anopheles/data/admixture_k1_metadata.txt"
REGION = "2L:1-500000"

runner = CliRunner()


@pytest.fixture(scope="module")
def cli_output(tmp_path_factory):
    out = tmp_path_factory.mktemp("cli_out") / "result.txt"
    result = runner.invoke(
        app,
        [
            "scan",
            "--genotypes", VCF,
            "--metadata", METADATA,
            "--out", str(out),
            "--start", "1",
            "--stop", "500000",
        ],
    )
    return result, out


class TestCliExitCode:
    def test_exit_zero(self, cli_output):
        result, _ = cli_output
        assert result.exit_code == 0, result.output


class TestCliOutputFile:
    def test_output_file_created(self, cli_output):
        _, out = cli_output
        assert out.exists()

    def test_output_parseable_as_tsv(self, cli_output):
        _, out = cli_output
        df = pd.read_csv(out, sep="\t", index_col=0)
        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0

    def test_output_has_required_columns(self, cli_output):
        _, out = cli_output
        df = pd.read_csv(out, sep="\t", index_col=0)
        required = {"position", "alternate", "frequency", "area"}
        assert required.issubset(df.columns)


class TestCliColumnMapping:
    """--sample-col / --lon-col / --lat-col must produce output identical to
    the equivalent API call with the same column name arguments."""

    def test_alt_col_names_produce_valid_output(self, tmp_path_factory):
        # Load canonical metadata, rename columns, write to a temp TSV
        canonical_meta = pd.read_csv(METADATA, sep="\t")
        alt_meta = canonical_meta.rename(columns={
            "sampleID": "sample_id",
            "x": "longitude",
            "y": "latitude",
        })
        meta_path = tmp_path_factory.mktemp("meta") / "alt_meta.tsv"
        alt_meta.to_csv(meta_path, sep="\t", index=False)

        out = tmp_path_factory.mktemp("out") / "result_alt.txt"
        result = runner.invoke(
            app,
            [
                "scan",
                "--genotypes", VCF,
                "--metadata", str(meta_path),
                "--out", str(out),
                "--start", "1",
                "--stop", "500000",
                "--sample-col", "sample_id",
                "--lon-col", "longitude",
                "--lat-col", "latitude",
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()
        df = pd.read_csv(out, sep="\t", index_col=0)
        assert len(df) > 0
        assert {"position", "alternate", "frequency", "area"}.issubset(df.columns)


class TestCliMatchesApi:
    def test_cli_matches_api_output(self, cli_output):
        """CLI and Python API must produce identical results for the same inputs."""
        _, out = cli_output
        cli_df = pd.read_csv(out, sep="\t", index_col=0).reset_index(drop=True)

        genotypes, positions, _ = read_vcf(VCF, region=REGION)
        metadata = pd.read_csv(METADATA, sep="\t")
        api_df = scan(genotypes, metadata, positions, start=1, stop=500_000).reset_index(drop=True)

        pd.testing.assert_frame_equal(
            cli_df[["position", "alternate", "frequency", "area"]].reset_index(drop=True),
            api_df[["position", "alternate", "frequency", "area"]].reset_index(drop=True),
            check_exact=False,
            rtol=1e-5,
        )

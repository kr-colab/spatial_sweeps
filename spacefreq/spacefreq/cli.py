"""Typer-based CLI for spacefreq."""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd
import typer
from rich.console import Console

from spacefreq.genotypes import read_vcf
from spacefreq.scan import scan as _scan

app = typer.Typer(help="Spatial genome scan: allele frequency and landscape area.")
console = Console()


@app.callback()
def _main() -> None:
    """spacefreq — spatial genome scanning tools."""


@app.command()
def scan(
    genotypes: Path = typer.Option(..., help="Path to VCF or VCF.gz genotype file."),
    metadata: Path = typer.Option(..., help="Path to sample metadata TSV."),
    out: Path = typer.Option(..., help="Output path for the results TSV."),
    filter: Optional[Path] = typer.Option(None, help="Path to genotype filter VCF (optional)."),
    start: Optional[int] = typer.Option(None, help="Genomic start position (inclusive)."),
    stop: Optional[int] = typer.Option(None, help="Genomic stop position (inclusive)."),
    chromosome: Optional[str] = typer.Option(None, help="Chromosome label."),
    min_locs: int = typer.Option(3, help="Minimum unique locations to record area."),
    transect: float = typer.Option(1.0, help="Transect width (km) for two-location alleles."),
    sample_area: float = typer.Option(1.0, help="Default area (km²) for single-location alleles."),
    count_missing_as_ancestral: bool = typer.Option(
        True,
        "--count-missing-as-ancestral/--no-count-missing-as-ancestral",
        help="Treat missing genotypes as homozygous reference.",
    ),
    sample_col: str = typer.Option("sampleID", help="Metadata column containing sample IDs."),
    lon_col: str = typer.Option("x", help="Metadata column containing longitude."),
    lat_col: str = typer.Option("y", help="Metadata column containing latitude."),
) -> None:
    """Run a spatial genome scan on a VCF genotype file."""
    # Build region string for read_vcf if chromosome + range are given
    region: str | None = None
    if chromosome is not None:
        if start is not None and stop is not None:
            region = f"{chromosome}:{start}-{stop}"
        else:
            region = chromosome

    console.print(f"Reading genotypes from [bold]{genotypes}[/bold] …")
    gt, positions, sample_ids = read_vcf(str(genotypes), region=region)

    console.print(f"Reading metadata from [bold]{metadata}[/bold] …")
    meta_df = pd.read_csv(str(metadata), sep="\t")

    console.print(f"Running scan over {len(positions):,} variants …")
    result = _scan(
        gt,
        meta_df,
        positions,
        filter=str(filter) if filter is not None else None,
        start=start,
        stop=stop,
        chromosome=chromosome,
        min_locs=min_locs,
        transect=transect,
        sample_area=sample_area,
        count_missing_as_ancestral=count_missing_as_ancestral,
        sample_col=sample_col,
        lon_col=lon_col,
        lat_col=lat_col,
    )

    result.to_csv(str(out), sep="\t")
    console.print(f"[green]Done.[/green] {len(result):,} records written to [bold]{out}[/bold].")

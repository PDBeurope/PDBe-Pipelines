"""
Command-line interface for the mmCIF standardisation pipeline.

Usage examples:

    # Standardise a single file
    standard-cif --input model.cif --output standard_model.cif

    # Supply a local CCD file to enrich non-standard ligand _chem_comp data
    standard-cif --input model.cif --output standard_model.cif \\
                 --ccd-path /data/components.cif

    # Enable verbose logging
    standard-cif --input model.cif --output standard_model.cif --verbose
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

import click

from .pipeline import standardise


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--input", "-i",
    "input_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Input mmCIF coordinate file.",
)
@click.option(
    "--output", "-o",
    "output_path",
    required=True,
    type=click.Path(dir_okay=False, path_type=Path),
    help="Output path for the standardised mmCIF file.",
)
@click.option(
    "--ccd-path",
    "ccd_path",
    default=None,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help=(
        "Optional path to a Chemical Component Dictionary (CCD) mmCIF file. "
        "When supplied, _chem_comp is enriched with authoritative CCD data "
        "for any component present in the file."
    ),
)
@click.option(
    "--verbose", "-v",
    is_flag=True,
    default=False,
    help="Enable DEBUG-level logging.",
)
def main(
    input_path: Path,
    output_path: Path,
    ccd_path: Path | None,
    verbose: bool,
) -> None:
    """Standardise an mmCIF coordinate file to OneDep upload format."""
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(levelname)s  %(message)s",
    )
    standardise(input_path, output_path, ccd_path=ccd_path)

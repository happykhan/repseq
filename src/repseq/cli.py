"""CLI entry point for repseq."""

from __future__ import annotations

import click

from repseq.evaluate import run_evaluate
from repseq.select import run_select
from repseq.sweep import run_sweep


@click.group()
@click.version_option()
def cli():
    """repseq -- select representative bacterial isolates for long-read sequencing."""


@cli.command()
@click.option(
    "--assemblies",
    required=True,
    type=click.Path(exists=True, file_okay=False),
    help="Folder of .fasta/.fa/.fna assembly files.",
)
@click.option(
    "--tree",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="Pre-built Newick tree. If absent, mashtree is run.",
)
@click.option(
    "--kleborate",
    "kleborate_tsv",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="Pre-run Kleborate TSV. If absent, kleborate is run.",
)
@click.option(
    "--plasmid-finder",
    "plasmidfinder_tsv",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="Pre-run PlasmidFinder merged TSV. If absent, ABRicate (plasmidfinder db) is run.",
)
@click.option(
    "--hamronization",
    "hamronization_tsv",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="hAMRonization TSV (from AMRFinder+, RGI, ResFinder etc.). Takes priority over --kleborate.",
)
@click.option("--n", "n_select", default=10, type=int, help="Number of isolates to select.")
@click.option(
    "--alpha",
    default=0.5,
    type=click.FloatRange(0.0, 1.0),
    help="Budget split: alpha=1.0 -> pure phylo, alpha=0.0 -> pure AMR/replicon.",
)
@click.option(
    "--cooccurrence",
    is_flag=True,
    default=False,
    help="Add REP+AMR co-occurrence features to capture plasmid-gene linkage.",
)
@click.option(
    "--output-dir",
    default=".",
    type=click.Path(file_okay=False),
    help="Output directory (created if needed).",
)
def select(assemblies, tree, kleborate_tsv, plasmidfinder_tsv, hamronization_tsv, n_select, alpha, cooccurrence, output_dir):
    """Select N representative isolates from an assembly collection."""
    run_select(
        assemblies_dir=assemblies,
        tree_path=tree,
        kleborate_path=kleborate_tsv,
        plasmidfinder_path=plasmidfinder_tsv,
        hamronization_path=hamronization_tsv,
        n=n_select,
        alpha=alpha,
        cooccurrence=cooccurrence,
        output_dir=output_dir,
    )


@cli.command()
@click.option(
    "--selected",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="selected.txt file with one sample ID per line.",
)
@click.option(
    "--ground-truth",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Complete Kleborate TSV for all samples.",
)
@click.option(
    "--tree",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Complete Newick tree for Faith PD calculation.",
)
@click.option(
    "--output-dir",
    default=".",
    type=click.Path(file_okay=False),
    help="Output directory.",
)
def evaluate(selected, ground_truth, tree, output_dir):
    """Evaluate a selection against ground truth."""
    run_evaluate(
        selected_path=selected,
        ground_truth_path=ground_truth,
        tree_path=tree,
        output_dir=output_dir,
    )


@cli.command()
@click.option(
    "--assemblies",
    required=True,
    type=click.Path(exists=True, file_okay=False),
    help="Folder of .fasta/.fa/.fna assembly files.",
)
@click.option(
    "--tree",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="Pre-built Newick tree.",
)
@click.option(
    "--kleborate",
    "kleborate_tsv",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="Pre-run Kleborate TSV.",
)
@click.option(
    "--hamronization",
    "hamronization_tsv",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="hAMRonization TSV. Takes priority over --kleborate.",
)
@click.option("--n", "n_select", default=10, type=int, help="Number of isolates to select.")
@click.option(
    "--ground-truth",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Complete Kleborate TSV for evaluation.",
)
@click.option(
    "--output-dir",
    default=".",
    type=click.Path(file_okay=False),
    help="Output directory.",
)
def sweep(assemblies, tree, kleborate_tsv, hamronization_tsv, n_select, ground_truth, output_dir):
    """Sweep alpha from 0 to 1, generating a Pareto curve."""
    run_sweep(
        assemblies_dir=assemblies,
        tree_path=tree,
        kleborate_path=kleborate_tsv,
        hamronization_path=hamronization_tsv,
        n=n_select,
        ground_truth_path=ground_truth,
        output_dir=output_dir,
    )

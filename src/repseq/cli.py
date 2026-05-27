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
    "--method",
    type=click.Choice(["split", "joint"]),
    default="split",
    help="Selection method: 'split' (default alpha budget) or 'joint' (k-medoids on blended tree+AMR distance).",
)
@click.option(
    "--joint-weight",
    "joint_weight",
    default=0.5,
    type=click.FloatRange(0.0, 1.0),
    help="Weight for AMR distance in joint method (0=pure phylo, 1=pure AMR). Only used with --method joint.",
)
@click.option(
    "--output-dir",
    default=".",
    type=click.Path(file_okay=False),
    help="Output directory (created if needed).",
)
def select(assemblies, tree, kleborate_tsv, plasmidfinder_tsv, hamronization_tsv, n_select, alpha, method, joint_weight, output_dir):
    """Select N representative isolates from an assembly collection."""
    run_select(
        assemblies_dir=assemblies,
        tree_path=tree,
        kleborate_path=kleborate_tsv,
        plasmidfinder_path=plasmidfinder_tsv,
        hamronization_path=hamronization_tsv,
        n=n_select,
        alpha=alpha,
        output_dir=output_dir,
        method=method,
        joint_weight=joint_weight,
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


@cli.command()
@click.option(
    "--assemblies",
    required=True,
    type=click.Path(exists=True, file_okay=False),
    help="Folder of .fasta/.fa/.fna assembly files (used for auto-running tools if needed).",
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
    help="Pre-run Kleborate TSV.",
)
@click.option(
    "--plasmid-finder",
    "plasmidfinder_tsv",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="Pre-run PlasmidFinder merged TSV.",
)
@click.option(
    "--hamronization",
    "hamronization_tsv",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="hAMRonization TSV. Takes priority over --kleborate.",
)
@click.option(
    "--abricate-replicons",
    "abricate_replicons_tsv",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="Pre-run ABRicate (plasmidfinder db) TSV for replicon features.",
)
@click.option("--n", "n_select", default=10, type=int, help="Number of isolates to select.")
@click.option("--pop-size", default=100, type=int, help="NSGA-III population size.")
@click.option("--generations", default=300, type=int, help="Number of NSGA-III generations.")
@click.option("--seed", default=42, type=int, help="Random seed for reproducibility.")
@click.option(
    "--output-dir",
    default=".",
    type=click.Path(file_okay=False),
    help="Output directory (created if needed).",
)
def nsga3(assemblies, tree, kleborate_tsv, plasmidfinder_tsv, hamronization_tsv, abricate_replicons_tsv, n_select, pop_size, generations, seed, output_dir):
    """Multi-objective selection using NSGA-III (phylo distance + AMR + replicon coverage)."""
    from repseq.nsga3 import run_nsga3

    run_nsga3(
        assemblies_dir=assemblies,
        tree_path=tree,
        kleborate_path=kleborate_tsv,
        plasmidfinder_path=plasmidfinder_tsv,
        hamronization_path=hamronization_tsv,
        abricate_replicons_path=abricate_replicons_tsv,
        n=n_select,
        output_dir=output_dir,
        pop_size=pop_size,
        n_gen=generations,
        seed=seed,
    )


@cli.command()
@click.option(
    "--csv",
    "csv_path",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="Flat CSV/TSV input file (columns: isolate_id, site, tier, replicons). "
         "Primary entry point — see README for format details.",
)
@click.option(
    "--metadata",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="Metadata TSV or Excel file with isolate data (legacy format with column mapping).",
)
@click.option(
    "--guarantee-sites/--no-guarantee-sites",
    "guarantee_sites",
    default=True,
    help="Guarantee at least 1 isolate per sentinel site (default: on).",
)
@click.option(
    "--id-col",
    "id_col",
    default="isolate_id",
    help="Column name for isolate ID (--metadata mode only).",
)
@click.option(
    "--lab-col",
    "lab_col",
    default="laboratory",
    help="Column name for laboratory/site (--metadata mode only).",
)
@click.option(
    "--date-col",
    "date_col",
    default="spec_date",
    help="Column name for specimen date (--metadata mode only).",
)
@click.option(
    "--tier-col",
    "tier_col",
    default=None,
    help="Column with pre-assigned tier (Tier1_CR, Tier2_MDR, etc.).",
)
@click.option(
    "--resist-pattern-col",
    "resist_pattern_col",
    default=None,
    help="Column with semicolon-separated resistant drug abbreviations.",
)
@click.option(
    "--r-drugs-list-col",
    "r_drugs_list_col",
    default=None,
    help="Column with comma-separated resistant drugs (raw).",
)
@click.option(
    "--plasmids-col",
    "plasmids_col",
    default=None,
    help="Column with comma-separated plasmid replicons.",
)
@click.option(
    "--carb-nonsus-col",
    "carb_nonsus_col",
    default=None,
    help="Column with boolean for carbapenem non-susceptibility.",
)
@click.option(
    "--colistin-r-col",
    "colistin_r_col",
    default=None,
    help="Column with boolean for colistin resistance.",
)
@click.option(
    "--exclude-drugs",
    "exclude_drugs_str",
    default=None,
    help="Comma-separated drugs to exclude from resist_pattern (e.g. SAM,CTT,TGC).",
)
@click.option(
    "--output-dir",
    default=".",
    type=click.Path(file_okay=False),
    help="Output directory (created if needed).",
)
def june(csv_path, metadata, guarantee_sites, id_col, lab_col, date_col,
         tier_col, resist_pattern_col, r_drugs_list_col, plasmids_col,
         carb_nonsus_col, colistin_r_col, exclude_drugs_str, output_dir):
    """Prioritise isolates for long-read sequencing using June Gayeta's method.

    Use --csv for the simple flat CSV format (recommended), or --metadata for
    the legacy column-mapping format. One of --csv or --metadata is required.
    """
    if csv_path and metadata:
        raise click.UsageError("Use --csv or --metadata, not both.")
    if not csv_path and not metadata:
        raise click.UsageError("One of --csv or --metadata is required.")

    if csv_path:
        from repseq.june import run_june_from_csv

        run_june_from_csv(
            csv_path=csv_path,
            output_dir=output_dir,
            guarantee_sites=guarantee_sites,
        )
    else:
        from repseq.june import run_june_prioritisation

        exclude_drugs = None
        if exclude_drugs_str:
            exclude_drugs = {d.strip() for d in exclude_drugs_str.split(",") if d.strip()}

        run_june_prioritisation(
            metadata_path=metadata,
            output_dir=output_dir,
            id_col=id_col,
            lab_col=lab_col,
            date_col=date_col,
            tier_col=tier_col,
            resist_pattern_col=resist_pattern_col,
            r_drugs_list_col=r_drugs_list_col,
            plasmids_col=plasmids_col,
            carb_nonsus_col=carb_nonsus_col,
            colistin_r_col=colistin_r_col,
            exclude_drugs=exclude_drugs,
            guarantee_sites=guarantee_sites,
        )


@cli.command(name="diversity-curve")
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
    help="Pre-run Kleborate TSV.",
)
@click.option(
    "--plasmid-finder",
    "plasmidfinder_tsv",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="Pre-run PlasmidFinder merged TSV.",
)
@click.option(
    "--hamronization",
    "hamronization_tsv",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="hAMRonization TSV. Takes priority over --kleborate.",
)
@click.option(
    "--abricate-replicons",
    "abricate_replicons_path",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="Pre-run ABRicate (plasmidfinder db) TSV for replicon features.",
)
@click.option("--max-k", default=30, type=int, help="Maximum number of representatives to evaluate.")
@click.option(
    "--output-dir",
    default=".",
    type=click.Path(file_okay=False),
    help="Output directory (created if needed).",
)
def diversity_curve(assemblies, tree, kleborate_tsv, plasmidfinder_tsv, hamronization_tsv, abricate_replicons_path, max_k, output_dir):
    """Plot diversity saturation curves to help choose the number of representatives."""
    from repseq.diversity import run_diversity_curve

    run_diversity_curve(
        assemblies_dir=assemblies,
        tree_path=tree,
        kleborate_path=kleborate_tsv,
        plasmidfinder_path=plasmidfinder_tsv,
        hamronization_path=hamronization_tsv,
        abricate_replicons_path=abricate_replicons_path,
        max_k=max_k,
        output_dir=output_dir,
    )

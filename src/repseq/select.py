"""Main selection logic: budget-split between PARNAS phylo and greedy AMR set cover."""

import os

import click
import pandas as pd

from repseq.phylo import run_mashtree, run_parnas, find_assemblies
from repseq.amr_cover import run_kleborate, parse_kleborate, greedy_set_cover, run_plasmidfinder, parse_plasmidfinder
from repseq.plots import plot_elbow, plot_scatter, plot_tree_heatmap


def run_select(
    assemblies_dir: str,
    tree_path: str | None,
    kleborate_path: str | None,
    plasmidfinder_path: str | None,
    n: int,
    alpha: float,
    output_dir: str,
) -> list[str]:
    """Run the full selection pipeline.

    Returns list of selected sample IDs.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Step 1: build tree if needed
    if tree_path is None:
        tree_path = run_mashtree(assemblies_dir, output_dir)

    # Step 2: run Kleborate if needed
    if kleborate_path is None:
        kleborate_path = run_kleborate(assemblies_dir, output_dir)

    # Step 3: budget split
    n_phylo = round(alpha * n)
    n_amr = n - n_phylo
    click.echo(f"\nBudget split: {n_phylo} phylogenetic + {n_amr} AMR/replicon (alpha={alpha})")

    # Step 4: PARNAS phylogenetic selection
    phylo_selected = run_parnas(tree_path, n_phylo, output_dir)

    # Step 5: parse Kleborate
    binary_matrix, features = parse_kleborate(kleborate_path)
    click.echo(f"Kleborate matrix: {binary_matrix.shape[0]} samples x {binary_matrix.shape[1]} features")

    # Step 5b: run PlasmidFinder and add replicon features
    if plasmidfinder_path is None:
        plasmidfinder_path = run_plasmidfinder(assemblies_dir, output_dir)
    binary_matrix = parse_plasmidfinder(plasmidfinder_path, binary_matrix)
    features = list(binary_matrix.columns)

    # Step 6: greedy set cover for AMR/replicon diversity
    amr_selected = greedy_set_cover(binary_matrix, phylo_selected, n_amr)

    # Step 7: combine
    all_selected = phylo_selected + amr_selected
    click.echo(f"\nFinal selection: {len(all_selected)} samples")
    click.echo(f"  Phylogenetic ({n_phylo}): {phylo_selected}")
    click.echo(f"  AMR/replicon ({n_amr}): {amr_selected}")

    # Write selected.txt
    selected_path = os.path.join(output_dir, "selected.txt")
    with open(selected_path, "w") as fh:
        for sid in all_selected:
            fh.write(sid + "\n")
    click.echo(f"\nSelected samples written to {selected_path}")

    # Write report.tsv
    _write_report(binary_matrix, phylo_selected, amr_selected, output_dir)

    # Write coverage_summary.txt
    _write_coverage_summary(binary_matrix, all_selected, features, output_dir)

    # Generate plots
    diversity_csv = os.path.join(output_dir, "diversity_scores.csv")
    if os.path.exists(diversity_csv):
        plot_elbow(diversity_csv, n, output_dir)

    plot_scatter(tree_path, binary_matrix, all_selected, output_dir)
    plot_tree_heatmap(tree_path, binary_matrix, all_selected, output_dir)

    return all_selected


def _write_report(
    binary_matrix: pd.DataFrame,
    phylo_selected: list[str],
    amr_selected: list[str],
    output_dir: str,
):
    """Write report.tsv with per-sample details."""
    all_selected = set(phylo_selected + amr_selected)
    rows = []
    for sid in binary_matrix.index:
        amr_genes = []
        replicons = []
        for col in binary_matrix.columns:
            if binary_matrix.loc[sid, col] == 1:
                if col.startswith("AMR:"):
                    amr_genes.append(col[4:])
                elif col.startswith("REP:"):
                    replicons.append(col[4:])

        if sid in phylo_selected:
            slot_type = "phylo"
        elif sid in amr_selected:
            slot_type = "amr"
        else:
            slot_type = "not_selected"

        rows.append({
            "sample_id": sid,
            "selected": "yes" if sid in all_selected else "no",
            "slot_type": slot_type,
            "amr_genes": ";".join(amr_genes) if amr_genes else "-",
            "replicons": ";".join(replicons) if replicons else "-",
            "n_amr_genes": len(amr_genes),
            "n_replicons": len(replicons),
        })

    report_df = pd.DataFrame(rows)
    report_path = os.path.join(output_dir, "report.tsv")
    report_df.to_csv(report_path, sep="\t", index=False)
    click.echo(f"Report written to {report_path}")


def _write_coverage_summary(
    binary_matrix: pd.DataFrame,
    selected: list[str],
    features: list[str],
    output_dir: str,
):
    """Write human-readable coverage summary."""
    amr_features = [f for f in features if f.startswith("AMR:")]
    rep_features = [f for f in features if f.startswith("REP:")]

    # All features in collection
    all_amr = set()
    all_rep = set()
    for sid in binary_matrix.index:
        row = binary_matrix.loc[sid]
        all_amr.update(f for f in amr_features if row[f] == 1)
        all_rep.update(f for f in rep_features if row[f] == 1)

    # Features covered by selection
    sel_amr = set()
    sel_rep = set()
    for sid in selected:
        if sid in binary_matrix.index:
            row = binary_matrix.loc[sid]
            sel_amr.update(f for f in amr_features if row[f] == 1)
            sel_rep.update(f for f in rep_features if row[f] == 1)

    pct_amr = (len(sel_amr) / len(all_amr) * 100) if all_amr else 100.0
    pct_rep = (len(sel_rep) / len(all_rep) * 100) if all_rep else 100.0

    summary_path = os.path.join(output_dir, "coverage_summary.txt")
    with open(summary_path, "w") as fh:
        fh.write("repseq Coverage Summary\n")
        fh.write("=" * 40 + "\n\n")
        fh.write(f"Total samples in collection: {len(binary_matrix)}\n")
        fh.write(f"Samples selected: {len(selected)}\n\n")
        fh.write(f"AMR gene features:\n")
        fh.write(f"  Total unique: {len(all_amr)}\n")
        fh.write(f"  Covered by selection: {len(sel_amr)}\n")
        fh.write(f"  Coverage: {pct_amr:.1f}%\n\n")
        fh.write(f"Replicon types:\n")
        fh.write(f"  Total unique: {len(all_rep)}\n")
        fh.write(f"  Covered by selection: {len(sel_rep)}\n")
        fh.write(f"  Coverage: {pct_rep:.1f}%\n\n")
        fh.write(f"Selected samples:\n")
        for sid in selected:
            fh.write(f"  {sid}\n")

    click.echo(f"Coverage summary written to {summary_path}")

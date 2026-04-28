"""Sweep alpha from 0 to 1 and generate Pareto curve."""

from __future__ import annotations

import os

import pandas as pd

from repseq.amr_cover import run_abricate, run_kleborate
from repseq.evaluate import run_evaluate
from repseq.log import print_header, print_message
from repseq.phylo import run_mashtree
from repseq.plots import plot_pareto
from repseq.select import run_select


def run_sweep(
    assemblies_dir: str,
    tree_path: str | None,
    kleborate_path: str | None,
    hamronization_path: str | None,
    n: int,
    ground_truth_path: str,
    output_dir: str,
) -> None:
    """Run selection + evaluation at alpha = 0.0, 0.1, ..., 1.0."""
    os.makedirs(output_dir, exist_ok=True)

    # Pre-build expensive inputs once so they are reused across all alpha sweeps
    if tree_path is None:
        tree_path = run_mashtree(assemblies_dir, output_dir)

    if hamronization_path is None and kleborate_path is None:
        kleborate_path = run_kleborate(assemblies_dir, output_dir)

    # Pre-run replicon typing once (ABRicate plasmidfinder db) — reused across all alphas
    abricate_replicons_path = run_abricate(assemblies_dir, output_dir, db="plasmidfinder")

    alphas = [round(a * 0.1, 1) for a in range(11)]
    results = []

    for alpha in alphas:
        print_header(f"Sweep  alpha = {alpha}")

        sweep_dir = os.path.join(output_dir, f"alpha_{alpha:.1f}")
        os.makedirs(sweep_dir, exist_ok=True)

        run_select(
            assemblies_dir=assemblies_dir,
            tree_path=tree_path,
            kleborate_path=kleborate_path,
            plasmidfinder_path=None,
            hamronization_path=hamronization_path,
            abricate_replicons_path=abricate_replicons_path,
            n=n,
            alpha=alpha,
            output_dir=sweep_dir,
        )

        selected_path = os.path.join(sweep_dir, "selected.txt")
        metrics = run_evaluate(
            selected_path=selected_path,
            ground_truth_path=ground_truth_path,
            tree_path=tree_path,
            output_dir=sweep_dir,
        )
        metrics["alpha"] = alpha
        results.append(metrics)

    pareto_path = os.path.join(output_dir, "pareto.tsv")
    pareto_df = pd.DataFrame(results)
    cols = ["alpha"] + [c for c in pareto_df.columns if c != "alpha"]
    pareto_df = pareto_df[cols]
    pareto_df.to_csv(pareto_path, sep="\t", index=False)
    print_message(f"Pareto table written to {pareto_path}", "success")

    # Write guidance note alongside the Pareto table
    note_path = os.path.join(output_dir, "pareto_note.txt")
    with open(note_path, "w") as fh:
        fh.write(
            "Pareto curve: alpha trade-off guide\n"
            "=====================================\n\n"
            "The Pareto sweep does NOT tell you which alpha to use — that is a\n"
            "scientific decision based on your study goals.\n\n"
            "How to choose alpha:\n"
            "  - alpha = 1.0  → pure phylogenetic diversity (ignore AMR/plasmid profile)\n"
            "  - alpha = 0.0  → pure AMR/replicon diversity (ignore phylogeny)\n"
            "  - alpha = 0.5  → equal weighting (default; good starting point)\n\n"
            "Look for the elbow in pareto_plot.png: the point where adding more\n"
            "phylo budget stops improving Faith PD meaningfully, or vice versa.\n"
            "Use the random_mean_faith_pd_pct and random_mean_amr_pct columns in\n"
            "pareto.tsv to check that your chosen alpha beats random selection on\n"
            "both axes.\n"
        )
    print_message(f"Pareto guidance note written to {note_path}", "success")

    plot_pareto(pareto_path, output_dir)
    print_message("Sweep complete.", "success")

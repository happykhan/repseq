"""Sweep alpha from 0 to 1 and generate Pareto curve."""

import os

import click
import pandas as pd

from repseq.select import run_select
from repseq.evaluate import run_evaluate
from repseq.plots import plot_pareto


def run_sweep(
    assemblies_dir: str,
    tree_path: str | None,
    kleborate_path: str | None,
    n: int,
    ground_truth_path: str,
    output_dir: str,
):
    """Run selection + evaluation at alpha = 0.0, 0.1, ..., 1.0."""
    os.makedirs(output_dir, exist_ok=True)

    # Pre-build tree and Kleborate if not provided, so they are reused across sweeps
    if tree_path is None:
        from repseq.phylo import run_mashtree
        tree_path = run_mashtree(assemblies_dir, output_dir)

    if kleborate_path is None:
        from repseq.amr_cover import run_kleborate
        kleborate_path = run_kleborate(assemblies_dir, output_dir)

    alphas = [round(a * 0.1, 1) for a in range(11)]
    results = []

    for alpha in alphas:
        click.echo(f"\n{'='*50}")
        click.echo(f"Sweep: alpha = {alpha}")
        click.echo(f"{'='*50}")

        sweep_dir = os.path.join(output_dir, f"alpha_{alpha:.1f}")
        os.makedirs(sweep_dir, exist_ok=True)

        # Run selection
        run_select(
            assemblies_dir=assemblies_dir,
            tree_path=tree_path,
            kleborate_path=kleborate_path,
            n=n,
            alpha=alpha,
            output_dir=sweep_dir,
        )

        # Run evaluation
        selected_path = os.path.join(sweep_dir, "selected.txt")
        metrics = run_evaluate(
            selected_path=selected_path,
            ground_truth_path=ground_truth_path,
            tree_path=tree_path,
            output_dir=sweep_dir,
        )
        metrics["alpha"] = alpha
        results.append(metrics)

    # Write Pareto table
    pareto_path = os.path.join(output_dir, "pareto.tsv")
    pareto_df = pd.DataFrame(results)
    # Reorder columns so alpha is first
    cols = ["alpha"] + [c for c in pareto_df.columns if c != "alpha"]
    pareto_df = pareto_df[cols]
    pareto_df.to_csv(pareto_path, sep="\t", index=False)
    click.echo(f"\nPareto table written to {pareto_path}")

    # Plot Pareto curve
    plot_pareto(pareto_path, output_dir)

    click.echo("\nSweep complete.")

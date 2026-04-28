"""Evaluate a selection against ground truth."""

import os

import click
import pandas as pd
from Bio import Phylo

from repseq.amr_cover import parse_kleborate


def faith_pd(tree, selected_names: set[str]) -> float:
    """Calculate Faith's phylogenetic diversity for selected taxa.

    Faith PD = sum of branch lengths connecting the selected taxa in the tree.
    """
    if not selected_names:
        return 0.0

    # Get all terminal names, stripping extensions
    from pathlib import Path
    terminal_map = {}
    for t in tree.get_terminals():
        name = Path(t.name).stem if t.name else ""
        terminal_map[name] = t

    # Find selected terminals
    selected_terminals = set()
    for name in selected_names:
        if name in terminal_map:
            selected_terminals.add(terminal_map[name])

    if len(selected_terminals) < 2:
        return 0.0

    # Calculate total branch length of minimum spanning subtree
    # For each internal node, include its branch length if it is an ancestor
    # of at least one selected terminal
    total_pd = 0.0

    def _has_selected_descendant(clade):
        """Check if clade has any selected descendant and sum PD."""
        nonlocal total_pd
        if clade.is_terminal():
            return clade in selected_terminals
        has_sel = False
        for child in clade.clades:
            child_has = _has_selected_descendant(child)
            if child_has:
                has_sel = True
                bl = child.branch_length if child.branch_length else 0.0
                total_pd += bl
        return has_sel

    _has_selected_descendant(tree.root)
    return total_pd


def run_evaluate(
    selected_path: str,
    ground_truth_path: str,
    tree_path: str,
    output_dir: str,
) -> dict:
    """Evaluate selection coverage metrics.

    Returns dict with pct_amr_covered, pct_replicons_covered, pct_faith_pd.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Read selected sample IDs
    with open(selected_path) as fh:
        selected = {line.strip() for line in fh if line.strip()}
    click.echo(f"Evaluating {len(selected)} selected samples...")

    # Parse ground truth Kleborate
    binary_matrix, features = parse_kleborate(ground_truth_path)

    amr_features = [f for f in features if f.startswith("AMR:")]
    rep_features = [f for f in features if f.startswith("REP:")]

    # Features present in the full collection
    all_amr = set()
    all_rep = set()
    for sid in binary_matrix.index:
        row = binary_matrix.loc[sid]
        all_amr.update(f for f in amr_features if row.get(f, 0) == 1)
        all_rep.update(f for f in rep_features if row.get(f, 0) == 1)

    # Features covered by selected samples
    sel_amr = set()
    sel_rep = set()
    for sid in selected:
        if sid in binary_matrix.index:
            row = binary_matrix.loc[sid]
            sel_amr.update(f for f in amr_features if row.get(f, 0) == 1)
            sel_rep.update(f for f in rep_features if row.get(f, 0) == 1)

    pct_amr = (len(sel_amr) / len(all_amr) * 100) if all_amr else 100.0
    pct_rep = (len(sel_rep) / len(all_rep) * 100) if all_rep else 100.0

    # Faith PD
    tree = Phylo.read(tree_path, "newick")
    all_names = set()
    for t in tree.get_terminals():
        from pathlib import Path
        name = Path(t.name).stem if t.name else ""
        all_names.add(name)

    total_faith = faith_pd(tree, all_names)
    selected_faith = faith_pd(tree, selected)
    pct_faith = (selected_faith / total_faith * 100) if total_faith > 0 else 100.0

    # Write coverage metrics
    metrics = {
        "pct_amr_covered": round(pct_amr, 2),
        "pct_replicons_covered": round(pct_rep, 2),
        "pct_faith_pd": round(pct_faith, 2),
        "n_selected": len(selected),
        "total_amr_features": len(all_amr),
        "covered_amr_features": len(sel_amr),
        "total_replicon_types": len(all_rep),
        "covered_replicon_types": len(sel_rep),
        "total_faith_pd": round(total_faith, 6),
        "selected_faith_pd": round(selected_faith, 6),
    }

    metrics_path = os.path.join(output_dir, "coverage_metrics.tsv")
    metrics_df = pd.DataFrame([metrics])
    metrics_df.to_csv(metrics_path, sep="\t", index=False)
    click.echo(f"Coverage metrics written to {metrics_path}")

    # Human-readable summary
    click.echo(f"\n--- Coverage Summary ---")
    click.echo(f"AMR gene coverage:    {pct_amr:.1f}% ({len(sel_amr)}/{len(all_amr)} features)")
    click.echo(f"Replicon coverage:    {pct_rep:.1f}% ({len(sel_rep)}/{len(all_rep)} types)")
    click.echo(f"Faith PD coverage:    {pct_faith:.1f}%")

    return metrics

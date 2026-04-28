"""Evaluate a selection against ground truth."""

from __future__ import annotations

import os
import random
from pathlib import Path

import pandas as pd
from Bio import Phylo

from repseq.amr_cover import parse_kleborate
from repseq.log import print_header, print_message

_RANDOM_REPS = 100


def faith_pd(tree, selected_names: set[str]) -> float:
    """Sum of branch lengths in the minimum spanning subtree of *selected_names*."""
    if not selected_names:
        return 0.0

    terminal_map = {}
    for t in tree.get_terminals():
        name = Path(t.name).stem if t.name else ""
        terminal_map[name] = t

    selected_terminals = {terminal_map[n] for n in selected_names if n in terminal_map}
    if len(selected_terminals) < 2:
        return 0.0

    total_pd = 0.0

    def _has_selected_descendant(clade) -> bool:
        nonlocal total_pd
        if clade.is_terminal():
            return clade in selected_terminals
        has_sel = False
        for child in clade.clades:
            if _has_selected_descendant(child):
                has_sel = True
                total_pd += child.branch_length or 0.0
        return has_sel

    _has_selected_descendant(tree.root)
    return total_pd


def _minimax_metrics(tree, selected_names: set[str], all_names: set[str]) -> dict[str, float]:
    """For each unselected sample compute distance to its nearest selected neighbour.

    Returns max and mean of those distances (the minimax criterion PARNAS optimises).
    """
    terminal_map = {Path(t.name).stem: t for t in tree.get_terminals() if t.name}
    selected_terminals = [terminal_map[n] for n in selected_names if n in terminal_map]
    unselected = [n for n in all_names - selected_names if n in terminal_map]

    if not selected_terminals or not unselected:
        return {"max_minimax_dist": 0.0, "mean_minimax_dist": 0.0}

    min_dists = []
    for name in unselected:
        t = terminal_map[name]
        min_dists.append(min(tree.distance(t, s) for s in selected_terminals))

    return {
        "max_minimax_dist":  round(max(min_dists), 6),
        "mean_minimax_dist": round(sum(min_dists) / len(min_dists), 6),
    }


def _random_baseline(
    tree,
    binary_matrix: pd.DataFrame,
    all_names: set[str],
    n_selected: int,
    total_faith: float,
    n_reps: int = _RANDOM_REPS,
) -> dict[str, float]:
    """Compare repseq selection against random draws of the same size.

    Returns mean Faith PD % and mean AMR coverage % over *n_reps* random selections.
    """
    amr_cols = [c for c in binary_matrix.columns if c.startswith("AMR:")]
    all_amr = {c for c in amr_cols if binary_matrix[c].sum() > 0}
    pool = [n for n in all_names if n in binary_matrix.index]
    k = min(n_selected, len(pool))

    faith_scores: list[float] = []
    amr_scores: list[float] = []

    for _ in range(n_reps):
        rand_sel = set(random.sample(pool, k))
        fd = faith_pd(tree, rand_sel)
        faith_scores.append(fd / total_faith * 100 if total_faith > 0 else 100.0)

        rand_amr: set[str] = set()
        for sid in rand_sel:
            row = binary_matrix.loc[sid]
            rand_amr.update(c for c in amr_cols if row.get(c, 0) == 1)
        amr_scores.append(len(rand_amr) / len(all_amr) * 100 if all_amr else 100.0)

    return {
        "random_mean_faith_pd_pct": round(sum(faith_scores) / n_reps, 2),
        "random_mean_amr_pct":      round(sum(amr_scores) / n_reps, 2),
    }


def run_evaluate(
    selected_path: str,
    ground_truth_path: str,
    tree_path: str,
    output_dir: str,
) -> dict:
    """Evaluate a selection: Faith PD, minimax distance, AMR coverage, and random baseline.

    Returns dict suitable for use in sweep Pareto table.
    """
    os.makedirs(output_dir, exist_ok=True)

    with open(selected_path) as fh:
        selected = {line.strip() for line in fh if line.strip()}
    print_message(f"Evaluating {len(selected)} selected samples...", "info")

    binary_matrix, features = parse_kleborate(ground_truth_path)
    amr_features = [f for f in features if f.startswith("AMR:")]
    rep_features = [f for f in features if f.startswith("REP:")]

    # Vectorised feature-coverage computation
    all_amr = set(binary_matrix[amr_features].columns[binary_matrix[amr_features].any()]) if amr_features else set()
    all_rep = set(binary_matrix[rep_features].columns[binary_matrix[rep_features].any()]) if rep_features else set()

    sel_in = [s for s in selected if s in binary_matrix.index]
    sel_amr = set(binary_matrix.loc[sel_in, amr_features].columns[binary_matrix.loc[sel_in, amr_features].any()]) if sel_in and amr_features else set()
    sel_rep = set(binary_matrix.loc[sel_in, rep_features].columns[binary_matrix.loc[sel_in, rep_features].any()]) if sel_in and rep_features else set()

    pct_amr = len(sel_amr) / len(all_amr) * 100 if all_amr else 100.0
    pct_rep = len(sel_rep) / len(all_rep) * 100 if all_rep else 100.0

    # Faith PD
    tree = Phylo.read(tree_path, "newick")
    all_names = {Path(t.name).stem for t in tree.get_terminals() if t.name}
    total_faith = faith_pd(tree, all_names)
    selected_faith = faith_pd(tree, selected)
    pct_faith = selected_faith / total_faith * 100 if total_faith > 0 else 100.0

    # Minimax distance (what PARNAS actually optimises)
    minimax = _minimax_metrics(tree, selected, all_names)

    # Random baseline (100 draws of the same N)
    print_message(f"Computing random baseline ({_RANDOM_REPS} draws)...", "info")
    baseline = _random_baseline(tree, binary_matrix, all_names, len(selected), total_faith)

    metrics = {
        "n_selected":               len(selected),
        "pct_faith_pd":             round(pct_faith, 2),
        "total_faith_pd":           round(total_faith, 6),
        "selected_faith_pd":        round(selected_faith, 6),
        "max_minimax_dist":         minimax["max_minimax_dist"],
        "mean_minimax_dist":        minimax["mean_minimax_dist"],
        "pct_amr_covered":          round(pct_amr, 2),
        "pct_replicons_covered":    round(pct_rep, 2),
        "total_amr_features":       len(all_amr),
        "covered_amr_features":     len(sel_amr),
        "total_replicon_types":     len(all_rep),
        "covered_replicon_types":   len(sel_rep),
        **baseline,
    }

    metrics_path = os.path.join(output_dir, "coverage_metrics.tsv")
    pd.DataFrame([metrics]).to_csv(metrics_path, sep="\t", index=False)
    print_message(f"Coverage metrics written to {metrics_path}", "success")

    print_header("Coverage Summary")
    print_message(f"Faith PD:             {pct_faith:.1f}%  (random baseline: {baseline['random_mean_faith_pd_pct']:.1f}%)", "info")
    print_message(f"AMR gene coverage:    {pct_amr:.1f}%  (random baseline: {baseline['random_mean_amr_pct']:.1f}%)", "info")
    print_message(f"Replicon coverage:    {pct_rep:.1f}%", "info")
    print_message(f"Minimax dist (max):   {minimax['max_minimax_dist']:.4f}", "info")
    print_message(f"Minimax dist (mean):  {minimax['mean_minimax_dist']:.4f}", "info")

    return metrics

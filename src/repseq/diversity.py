"""Diversity-vs-k curves for choosing the number of representatives."""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo
from scipy.spatial.distance import pdist, squareform

from repseq.amr_cover import build_feature_matrix
from repseq.log import print_header, print_message
from repseq.phylo import run_mashtree, run_parnas
from repseq.plots import plot_diversity_curve


def greedy_jaccard_curve(jac_dm: np.ndarray, max_k: int) -> dict[int, float]:
    """Return {k: mean_pairwise_jaccard} for k=2..max_k using a greedy algorithm.

    At each step, the candidate that maximises the mean pairwise Jaccard
    distance among the selected set is added.

    Parameters
    ----------
    jac_dm : np.ndarray
        Symmetric Jaccard distance matrix (n x n), with 0 on the diagonal.
    max_k : int
        Maximum subset size to evaluate.

    Returns
    -------
    dict[int, float]
        Mapping from k to mean pairwise Jaccard distance.
    """
    n = jac_dm.shape[0]
    max_k = min(max_k, n)

    # Edge case: no diversity at all
    if jac_dm.max() == 0.0:
        return {k: 0.0 for k in range(2, max_k + 1)}

    # Initialise with the pair having the highest Jaccard distance
    best_i, best_j = 0, 1
    best_val = jac_dm[0, 1]
    for i in range(n):
        for j in range(i + 1, n):
            if jac_dm[i, j] > best_val:
                best_val = jac_dm[i, j]
                best_i, best_j = i, j

    selected = [best_i, best_j]
    # current_sum tracks sum of all ordered pairs (i.e. jac_dm[sel, sel].sum())
    # For a symmetric matrix with 0 diagonal, sum = 2 * jac_dm[i, j] for a pair
    current_sum = 2.0 * jac_dm[best_i, best_j]

    result: dict[int, float] = {}
    # k=2: mean = sum / (2*1) = jac_dm[i,j]
    result[2] = current_sum / (2 * 1)

    all_indices = set(range(n))

    for k in range(3, max_k + 1):
        candidates = all_indices - set(selected)
        if not candidates:
            break

        best_candidate = -1
        best_new_mean = -1.0
        n_ordered_pairs = k * (k - 1)

        for c in candidates:
            # Adding c: new_sum = current_sum + 2 * sum(jac_dm[c, s] for s in selected)
            contrib = 2.0 * jac_dm[c, selected].sum()
            new_sum = current_sum + contrib
            new_mean = new_sum / n_ordered_pairs
            if new_mean > best_new_mean:
                best_new_mean = new_mean
                best_candidate = c

        selected.append(best_candidate)
        current_sum += 2.0 * jac_dm[best_candidate, selected[:-1]].sum()
        result[k] = best_new_mean

    return result


def run_diversity_curve(
    assemblies_dir: str,
    tree_path: str | None,
    kleborate_path: str | None,
    plasmidfinder_path: str | None,
    hamronization_path: str | None,
    abricate_replicons_path: str | None,
    max_k: int,
    output_dir: str,
) -> str:
    """Compute and plot diversity saturation curves for phylo, AMR, and replicon diversity.

    For each k from 2 to max_k, measures how much diversity is captured by k
    representatives under three different criteria. Outputs a TSV and a plot.

    Returns path to diversity_curve.png.
    """
    os.makedirs(output_dir, exist_ok=True)
    print_header("Diversity Saturation Curve")

    # --- Build tree if needed ---
    if tree_path is None:
        tree_path = run_mashtree(assemblies_dir, output_dir)

    # --- Build feature matrix ---
    has_amr_source = hamronization_path is not None or kleborate_path is not None
    matrix_assemblies_dir = None if has_amr_source else assemblies_dir
    binary_matrix, kleborate_path = build_feature_matrix(
        assemblies_dir=matrix_assemblies_dir,
        output_dir=output_dir,
        hamronization_path=hamronization_path,
        kleborate_path=kleborate_path,
        plasmidfinder_path=plasmidfinder_path,
        abricate_replicons_path=abricate_replicons_path,
    )

    # --- Intersect tree terminals with matrix ---
    tree = Phylo.read(tree_path, "newick")
    tree_terminals = {Path(t.name).stem for t in tree.get_terminals() if t.name}
    names = sorted(set(binary_matrix.index) & tree_terminals)

    if len(names) < 3:
        print_message(
            f"Only {len(names)} samples in both tree and matrix -- need at least 3.",
            "error",
        )
        return os.path.join(output_dir, "diversity_curve.png")

    # Cap max_k to number of samples minus 1
    max_k = min(max_k, len(names) - 1)
    print_message(
        f"Computing diversity curves for k=2..{max_k} across {len(names)} samples",
        "info",
    )

    # --- Phylo curve from PARNAS ---
    diversity_csv = os.path.join(output_dir, "diversity_scores.csv")
    if not os.path.exists(diversity_csv):
        print_message(
            "diversity_scores.csv not found -- running PARNAS to generate it",
            "info",
        )
        run_parnas(tree_path, max_k, output_dir, n_total=max_k)

    phylo_df = pd.read_csv(diversity_csv)
    phylo_df.columns = phylo_df.columns.str.strip()
    phylo_curve: dict[int, float] = dict(
        zip(
            phylo_df["Representatives"].astype(int),
            phylo_df["Diversity_covered"].astype(float),
            strict=False,
        )
    )

    # --- AMR Jaccard curve ---
    amr_cols = [c for c in binary_matrix.columns if c.startswith("AMR:")]
    if amr_cols:
        amr_arr = binary_matrix.loc[names, amr_cols].values.astype(float)
        amr_jac = squareform(pdist(amr_arr, metric="jaccard"))
        amr_jac = np.nan_to_num(amr_jac, nan=0.0)
        if amr_jac.max() == 0.0:
            print_message(
                "AMR Jaccard matrix is all zeros -- all samples have identical AMR profiles. "
                "AMR curve will be flat.",
                "warning",
            )
        amr_curve = greedy_jaccard_curve(amr_jac, max_k)
    else:
        print_message("No AMR features found -- skipping AMR curve", "warning")
        amr_curve = {k: np.nan for k in range(2, max_k + 1)}

    # --- Replicon Jaccard curve ---
    rep_cols = [c for c in binary_matrix.columns if c.startswith("REP:")]
    if rep_cols:
        rep_arr = binary_matrix.loc[names, rep_cols].values.astype(float)
        rep_jac = squareform(pdist(rep_arr, metric="jaccard"))
        rep_jac = np.nan_to_num(rep_jac, nan=0.0)
        if rep_jac.max() == 0.0:
            print_message(
                "Replicon Jaccard matrix is all zeros -- all samples have identical replicon profiles. "
                "Replicon curve will be flat.",
                "warning",
            )
        rep_curve = greedy_jaccard_curve(rep_jac, max_k)
    else:
        print_message("No replicon features found -- skipping replicon curve", "warning")
        rep_curve = {k: np.nan for k in range(2, max_k + 1)}

    # --- Build output DataFrame ---
    rows = []
    for k in range(2, max_k + 1):
        rows.append({
            "k": k,
            "phylo_pct": phylo_curve.get(k, np.nan),
            "amr_jaccard_div": amr_curve.get(k, np.nan),
            "rep_jaccard_div": rep_curve.get(k, np.nan),
        })

    curve_df = pd.DataFrame(rows)

    tsv_path = os.path.join(output_dir, "diversity_curve.tsv")
    curve_df.to_csv(tsv_path, sep="\t", index=False)
    print_message(f"Diversity curve data written to {tsv_path}", "success")

    # --- Plot ---
    png_path = plot_diversity_curve(curve_df, output_dir)

    # --- Summary: k at which each curve reaches 80%, 90%, 95% of its max ---
    print_header("Diversity Curve Summary")
    for col_label, col_name in [
        ("Phylogenetic", "phylo_pct"),
        ("AMR (Jaccard)", "amr_jaccard_div"),
        ("Replicon (Jaccard)", "rep_jaccard_div"),
    ]:
        vals = curve_df[col_name].dropna()
        if vals.empty or vals.max() == 0:
            print_message(f"{col_label}: no data or all zeros", "warning")
            continue
        col_max = vals.max()
        for threshold_pct in [80, 90, 95]:
            threshold = col_max * threshold_pct / 100.0
            hits = curve_df.loc[curve_df[col_name] >= threshold, "k"]
            if not hits.empty:
                k_at = int(hits.iloc[0])
                print_message(
                    f"{col_label}: {threshold_pct}% of max reached at k={k_at}",
                    "info",
                )
            else:
                print_message(
                    f"{col_label}: {threshold_pct}% of max not reached in k=2..{max_k}",
                    "warning",
                )

    return png_path

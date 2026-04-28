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
    """Return {k: pct_of_total_pairwise_diversity} for k=2..max_k.

    Greedy: at each step add the sample with the highest sum of distances to
    the current selection (maximises cumulative pairwise diversity).

    Value is normalised to % of the total pairwise diversity in the full
    collection, so the curve starts low at k=2 and grows toward 100% at k=N —
    directly comparable to PARNAS's Diversity_covered column.
    """
    n = jac_dm.shape[0]
    max_k = min(max_k, n)

    total_sum = jac_dm.sum() / 2.0  # sum of all unordered pairs
    if total_sum == 0.0:
        return {k: 0.0 for k in range(2, max_k + 1)}

    # Start with the pair with highest Jaccard distance
    upper = np.triu(jac_dm, k=1)
    best_i, best_j = (int(x) for x in np.unravel_index(np.argmax(upper), upper.shape))
    selected = [best_i, best_j]
    pair_sum = float(jac_dm[best_i, best_j])

    result: dict[int, float] = {2: pair_sum / total_sum * 100}

    # dist_to_sel[i] = sum of distances from sample i to all currently selected samples
    dist_to_sel = jac_dm[:, selected].sum(axis=1).astype(float)

    all_indices = set(range(n))
    for _ in range(3, max_k + 1):
        candidates = list(all_indices - set(selected))
        if not candidates:
            break
        gains = dist_to_sel[candidates]
        best_c = candidates[int(np.argmax(gains))]
        pair_sum += dist_to_sel[best_c]
        selected.append(best_c)
        dist_to_sel += jac_dm[:, best_c]
        result[len(selected)] = pair_sum / total_sum * 100

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
        normed = vals / col_max * 100.0
        for threshold_pct in [80, 90, 95]:
            hits = curve_df.loc[normed >= threshold_pct, "k"]
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

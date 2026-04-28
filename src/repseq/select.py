"""Main selection logic: budget-split between PARNAS phylo and greedy AMR set cover."""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo
from scipy.spatial.distance import pdist, squareform

from repseq.amr_cover import (
    add_cooccurrence_features,
    greedy_set_cover,
    parse_abricate_replicons,
    parse_hamronization,
    parse_kleborate,
    parse_plasmidfinder,
    read_st_map,
    run_abricate,
    run_kleborate,
)
from repseq.joint import compute_joint_dist, run_kmedoids, tree_to_dist_matrix
from repseq.log import print_message
from repseq.phylo import find_assemblies, run_mashtree, run_parnas
from repseq.plots import plot_elbow, plot_scatter, plot_tree_heatmap


def run_select(
    assemblies_dir: str,
    tree_path: str | None,
    kleborate_path: str | None,
    plasmidfinder_path: str | None,
    hamronization_path: str | None,
    n: int,
    alpha: float,
    output_dir: str,
    cooccurrence: bool = False,
    abricate_replicons_path: str | None = None,
    method: str = "split",
    joint_weight: float = 0.5,
) -> list[str]:
    """Run the full selection pipeline.

    AMR feature priority: hAMRonization > Kleborate > ABRicate (auto-run).
    Replicon priority: --plasmid-finder > pre-run ABRicate > auto-run ABRicate.

    Returns list of selected sample IDs.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Step 1: build tree if needed
    if tree_path is None:
        tree_path = run_mashtree(assemblies_dir, output_dir)

    # Step 2: build AMR binary matrix from best available source
    #         (needed for both split and joint methods)
    if hamronization_path:
        print_message("Using hAMRonization input for AMR features", "info")
        binary_matrix, _ = parse_hamronization(hamronization_path)
    elif kleborate_path:
        print_message("Using Kleborate input for AMR features", "info")
        binary_matrix, _ = parse_kleborate(kleborate_path)
    else:
        print_message("No AMR input provided — running ABRicate (ncbi db) automatically", "info")
        kleborate_path = run_kleborate(assemblies_dir, output_dir)
        binary_matrix, _ = parse_kleborate(kleborate_path)

    print_message(
        f"AMR matrix: {binary_matrix.shape[0]} samples x {binary_matrix.shape[1]} features", "info"
    )

    # Step 4a: pad matrix so every assembly is represented (even if it had no hits)
    all_assembly_stems = {Path(p).stem for p in find_assemblies(assemblies_dir)}
    missing_stems = all_assembly_stems - set(binary_matrix.index)
    if missing_stems:
        print_message(f"Padding {len(missing_stems)} assemblies with no AMR hits into matrix", "info")
        zero_rows = pd.DataFrame(
            0,
            index=sorted(missing_stems),
            columns=binary_matrix.columns if len(binary_matrix.columns) else pd.Index([]),
        )
        binary_matrix = pd.concat([binary_matrix, zero_rows])

    # Step 5: add replicon features from best available source
    if plasmidfinder_path:
        print_message("Using pre-run PlasmidFinder input for replicon features", "info")
        binary_matrix = parse_plasmidfinder(plasmidfinder_path, binary_matrix)
    elif abricate_replicons_path:
        print_message("Using pre-run ABRicate replicon output", "info")
        binary_matrix = parse_abricate_replicons(abricate_replicons_path, binary_matrix)
    else:
        print_message("Running ABRicate (plasmidfinder db) for replicon features", "info")
        abricate_rep_path = run_abricate(assemblies_dir, output_dir, db="plasmidfinder")
        binary_matrix = parse_abricate_replicons(abricate_rep_path, binary_matrix)

    # Step 5b: optionally enrich with co-occurrence (REP+AMR) features
    if cooccurrence:
        binary_matrix = add_cooccurrence_features(binary_matrix)

    features = list(binary_matrix.columns)

    # Step 6: select samples (method-dependent)
    if method == "joint":
        # Joint k-medoids on blended tree + AMR distance
        tree = Phylo.read(tree_path, "newick")
        tree_terminals = {Path(t.name).stem for t in tree.get_terminals() if t.name}
        names = [n_id for n_id in binary_matrix.index if n_id in tree_terminals]

        if len(names) < n:
            print_message(
                f"Only {len(names)} samples present in both tree and AMR matrix "
                f"(requested {n}); selecting all of them.",
                "warning",
            )

        # Tree distance matrix
        tree_dm = tree_to_dist_matrix(tree, names)

        # AMR Jaccard distance matrix (handle all-zero rows gracefully)
        amr_arr = binary_matrix.loc[names].values.astype(float)
        amr_dm = squareform(pdist(amr_arr, metric="jaccard"))
        amr_dm = np.nan_to_num(amr_dm, nan=0.0)

        # Blend and run k-medoids
        combined = compute_joint_dist(tree_dm, amr_dm, joint_weight)
        indices = run_kmedoids(combined, n)
        all_selected = [names[i] for i in indices]
        phylo_selected: list[str] = []
        amr_selected = all_selected

        print_message(
            f"Joint k-medoids selected {len(all_selected)} samples "
            f"(joint_weight={joint_weight})",
            "success",
        )
        slot_label = "joint"
    else:
        # Default: budget-split between PARNAS phylo and greedy AMR set cover
        n_phylo = round(alpha * n)
        n_amr = n - n_phylo
        print_message(
            f"Budget split: {n_phylo} phylogenetic + {n_amr} AMR/replicon (alpha={alpha})",
            "info",
        )
        phylo_selected = run_parnas(tree_path, n_phylo, output_dir, n_total=n)
        amr_selected = greedy_set_cover(binary_matrix, phylo_selected, n_amr)
        all_selected = phylo_selected + amr_selected
        print_message(f"Final selection: {len(all_selected)} samples", "success")
        print_message(f"  Phylogenetic ({n_phylo}): {phylo_selected}", "info")
        print_message(f"  AMR/replicon ({n_amr}): {amr_selected}", "info")
        slot_label = "amr"

    # Write selected.txt
    selected_path = os.path.join(output_dir, "selected.txt")
    with open(selected_path, "w") as fh:
        for sid in all_selected:
            fh.write(sid + "\n")
    print_message(f"Selected samples written to {selected_path}", "success")

    # Write report.tsv
    _write_report(binary_matrix, phylo_selected, amr_selected, output_dir, slot_label=slot_label)

    # Write coverage_summary.txt
    st_map = read_st_map(kleborate_path) if kleborate_path else {}
    _write_coverage_summary(binary_matrix, all_selected, features, output_dir, st_map=st_map)

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
    slot_label: str = "amr",
) -> None:
    """Write report.tsv with per-sample details.

    Parameters
    ----------
    slot_label : str
        Label used for AMR-selected samples in the slot_type column.
        Pass "joint" when using the joint k-medoids method.
    """
    all_selected = set(phylo_selected + amr_selected)
    phylo_set = set(phylo_selected)
    amr_set = set(amr_selected)

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

        if sid in phylo_set:
            slot_type = "phylo"
        elif sid in amr_set:
            slot_type = slot_label
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
    print_message(f"Report written to {report_path}", "success")


def _write_coverage_summary(
    binary_matrix: pd.DataFrame,
    selected: list[str],
    features: list[str],
    output_dir: str,
    st_map: dict[str, str] | None = None,
) -> None:
    """Write human-readable coverage summary."""
    amr_features = [f for f in features if f.startswith("AMR:")]
    rep_features = [f for f in features if f.startswith("REP:")]

    # Vectorised: find features present anywhere in the collection
    if amr_features:
        all_amr = set(binary_matrix[amr_features].columns[binary_matrix[amr_features].any(axis=0)])
    else:
        all_amr: set[str] = set()

    if rep_features:
        all_rep = set(binary_matrix[rep_features].columns[binary_matrix[rep_features].any(axis=0)])
    else:
        all_rep: set[str] = set()

    # Features covered by selection (vectorised)
    sel_in_matrix = [s for s in selected if s in binary_matrix.index]
    if sel_in_matrix and amr_features:
        sel_amr = set(
            binary_matrix.loc[sel_in_matrix, amr_features]
            .columns[binary_matrix.loc[sel_in_matrix, amr_features].any(axis=0)]
        )
    else:
        sel_amr: set[str] = set()

    if sel_in_matrix and rep_features:
        sel_rep = set(
            binary_matrix.loc[sel_in_matrix, rep_features]
            .columns[binary_matrix.loc[sel_in_matrix, rep_features].any(axis=0)]
        )
    else:
        sel_rep: set[str] = set()

    pct_amr = (len(sel_amr) / len(all_amr) * 100) if all_amr else 100.0
    pct_rep = (len(sel_rep) / len(all_rep) * 100) if all_rep else 100.0

    # ST coverage
    if st_map:
        all_sts = set(st_map.values())
        sel_sts = {st_map[s] for s in selected if s in st_map}
        pct_st = len(sel_sts) / len(all_sts) * 100 if all_sts else 100.0
    else:
        all_sts: set[str] = set()
        sel_sts: set[str] = set()
        pct_st = 100.0

    summary_path = os.path.join(output_dir, "coverage_summary.txt")
    with open(summary_path, "w") as fh:
        fh.write("repseq Coverage Summary\n")
        fh.write("=" * 40 + "\n\n")
        fh.write(f"Total samples in collection: {len(binary_matrix)}\n")
        fh.write(f"Samples selected: {len(selected)}\n\n")
        fh.write("AMR gene features:\n")
        fh.write(f"  Total unique: {len(all_amr)}\n")
        fh.write(f"  Covered by selection: {len(sel_amr)}\n")
        fh.write(f"  Coverage: {pct_amr:.1f}%\n\n")
        fh.write("Replicon types:\n")
        fh.write(f"  Total unique: {len(all_rep)}\n")
        fh.write(f"  Covered by selection: {len(sel_rep)}\n")
        fh.write(f"  Coverage: {pct_rep:.1f}%\n\n")
        if all_sts:
            fh.write("Sequence types (MLST):\n")
            fh.write(f"  Total unique: {len(all_sts)}\n")
            fh.write(f"  Covered by selection: {len(sel_sts)}\n")
            fh.write(f"  Coverage: {pct_st:.1f}%\n\n")
        fh.write("Selected samples:\n")
        for sid in selected:
            fh.write(f"  {sid}\n")

    print_message(f"Coverage summary written to {summary_path}", "success")

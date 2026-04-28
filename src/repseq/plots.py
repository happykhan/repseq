"""All plotting functions for repseq."""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib
import matplotlib.colors
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from Bio import Phylo
from sklearn.manifold import MDS

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from repseq.log import print_message

# Drug class -> colour mapping (matches reference figure palette)
DRUG_CLASS_COLORS = {
    "AGly":          "#4575b4",
    "Bla":           "#fee090",
    "Bla_Carb":      "#fdae61",
    "Bla_ESBL":      "#ffffbf",
    "Bla_ESBL_inhR": "#e0f3f8",
    "Bla_inhR":      "#abd9e9",
    "Col":           "#74add1",
    "Fcyn":          "#9e0142",
    "Flq":           "#d73027",
    "Gly":           "#3288bd",
    "MLS":           "#762a83",
    "Phe":           "#e66101",
    "Rif":           "#f46d43",
    "Sul":           "#d9ef8b",
    "Tet":           "#1a9850",
    "Tgc":           "#66bd63",
    "Tmt":           "#313695",
    "REP":           "#b2abd2",
    "Other":         "#cccccc",
}

DRUG_CLASS_NAMES = {
    "AGly":          "Aminoglycoside",
    "Bla":           "Beta-lactam",
    "Bla_Carb":      "Carbapenem",
    "Bla_ESBL":      "ESBL",
    "Bla_ESBL_inhR": "ESBL + inhR",
    "Bla_inhR":      "Beta-lactam inhR",
    "Col":           "Colistin",
    "Fcyn":          "Fosfomycin",
    "Flq":           "Fluoroquinolone",
    "Gly":           "Glycopeptide",
    "MLS":           "Macrolide / MLS",
    "Phe":           "Phenicol",
    "Rif":           "Rifampicin",
    "Sul":           "Sulfonamide",
    "Tet":           "Tetracycline",
    "Tgc":           "Tigecycline",
    "Tmt":           "Trimethoprim",
    "REP":           "Plasmid replicon",
    "Other":         "Other",
}


def _col_drug_class(col_name: str) -> str:
    """Extract drug class key from an AMR: or REP: feature column name."""
    if col_name.startswith("REP:"):
        return "REP"
    if col_name.startswith("AMR:"):
        parts = col_name.split(":")
        if len(parts) >= 2:
            return parts[1].replace("_acquired", "")
    return "Other"


def _get_leaf_order(tree) -> list[str]:
    """Get leaf names in tree traversal order (topology order)."""
    names = []
    for clade in tree.get_terminals():
        name = clade.name if clade.name else ""
        name = Path(name).stem
        names.append(name)
    return names


def _draw_tree_axes(
    tree,
    ax,
    leaf_order: list[str],
    selected_ids: list[str],
    n_leaves: int,
) -> None:
    """Draw a dendrogram of *tree* on *ax*, rows aligned with *leaf_order*."""
    leaf_y = {name: i for i, name in enumerate(leaf_order)}
    positions: dict[int, tuple[float, float]] = {}

    def _assign(clade, cumx: float = 0.0):
        x = cumx + (clade.branch_length or 0.0)
        if clade.is_terminal():
            name = Path(clade.name or "").stem
            positions[id(clade)] = (x, leaf_y.get(name, 0))
        else:
            for child in clade.clades:
                _assign(child, x)
            ys = [positions[id(c)][1] for c in clade.clades]
            positions[id(clade)] = (x, (min(ys) + max(ys)) / 2.0)

    _assign(tree.root)

    max_x = max(v[0] for v in positions.values()) if positions else 1.0
    if max_x == 0:
        max_x = 1.0

    def _draw(clade):
        px, py = positions[id(clade)]
        pxn = px / max_x
        for child in clade.clades:
            cx, cy = positions[id(child)]
            cxn = cx / max_x
            ax.plot([pxn, cxn], [cy, cy], color="black", linewidth=0.7, solid_capstyle="butt")
            ax.plot([pxn, pxn], [py, cy], color="black", linewidth=0.7)
            _draw(child)

    _draw(tree.root)

    for clade in tree.get_terminals():
        name = Path(clade.name or "").stem
        if name in selected_ids and name in leaf_y:
            px, py = positions[id(clade)]
            ax.plot(px / max_x, py, "o", color="red", markersize=4, zorder=10)

    ax.set_xlim(-0.02, 1.05)
    ax.set_ylim(n_leaves - 0.5, -0.5)
    ax.axis("off")


def _draw_drug_class_bar(ax, drug_classes: list[str], n_features: int) -> None:
    """Draw a coloured rectangle per column indicating drug class."""
    for i, dc in enumerate(drug_classes):
        color = DRUG_CLASS_COLORS.get(dc, "#cccccc")
        ax.add_patch(plt.Rectangle(
            (i - 0.5, 0), 1, 1,
            color=color, ec="white", lw=0.4,
            transform=ax.transData,
        ))
    ax.set_xlim(-0.5, n_features - 0.5)
    ax.set_ylim(0, 1)
    ax.axis("off")


def _draw_legend(ax, drug_classes: list[str]) -> None:
    """Draw drug-class legend + present/absent key in *ax*."""
    ax.axis("off")
    y = 0.98

    ax.text(0.05, y, "Drug class", fontsize=8, fontweight="bold",
            transform=ax.transAxes, va="top")
    y -= 0.04

    seen: list[str] = []
    for dc in drug_classes:
        if dc not in seen:
            seen.append(dc)

    for dc in seen:
        color = DRUG_CLASS_COLORS.get(dc, "#cccccc")
        label = DRUG_CLASS_NAMES.get(dc, dc)
        ax.add_patch(plt.Rectangle(
            (0.02, y - 0.015), 0.1, 0.025,
            transform=ax.transAxes, color=color, ec="black", lw=0.5,
        ))
        ax.text(0.16, y - 0.003, label, fontsize=7, transform=ax.transAxes, va="center")
        y -= 0.038

    y -= 0.02
    ax.text(0.05, y, "Key", fontsize=8, fontweight="bold",
            transform=ax.transAxes, va="top")
    y -= 0.04

    for label, color in [
        ("Present",    "#f46d43"),
        ("Absent",     "#d0d0d0"),
        ("Not typed",  "#f5f5f5"),
        ("Selected",   "red"),
    ]:
        if label == "Selected":
            ax.plot(
                [0.07], [y - 0.003], "o",
                color=color, markersize=6, transform=ax.transAxes, zorder=10,
                clip_on=False,
            )
        else:
            ax.add_patch(plt.Rectangle(
                (0.02, y - 0.015), 0.1, 0.025,
                transform=ax.transAxes, color=color, ec="black", lw=0.5,
            ))
        ax.text(0.16, y - 0.003, label, fontsize=7, transform=ax.transAxes, va="center")
        y -= 0.035


def plot_elbow(diversity_csv: str, n_chosen: int, output_dir: str) -> str:
    """Plot diversity covered vs n, with vertical line at chosen n."""
    out_path = os.path.join(output_dir, "elbow_plot.png")
    try:
        df = pd.read_csv(diversity_csv)
    except Exception as e:
        print_message(f"Cannot read diversity CSV for elbow plot: {e}", "warning")
        return out_path

    cols = df.columns.tolist()
    if len(cols) >= 2:
        x_col, y_col = cols[0], cols[1]
    else:
        print_message("Diversity CSV has fewer than 2 columns, skipping elbow plot.", "warning")
        return out_path

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(df[x_col], df[y_col], "o-", color="#2c7bb6", linewidth=2, markersize=6)
    ax.axvline(x=n_chosen, color="red", linestyle="--", linewidth=1.5, label=f"n = {n_chosen}")
    ax.set_xlabel("Number of representatives")
    ax.set_ylabel("Diversity captured")
    ax.set_title("Phylogenetic diversity vs number of representatives")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print_message(f"Elbow plot saved to {out_path}", "success")
    return out_path


def plot_scatter(
    tree_path: str,
    binary_matrix: pd.DataFrame,
    selected_ids: list[str],
    output_dir: str,
) -> str:
    """PCoA of Mash distances, coloured by AMR gene count, selected starred."""
    out_path = os.path.join(output_dir, "scatter_plot.png")

    try:
        tree = Phylo.read(tree_path, "newick")
    except Exception as e:
        print_message(f"Cannot parse tree for scatter plot: {e}", "warning")
        return out_path

    terminals = list(tree.get_terminals())
    names = []
    for t in terminals:
        name = t.name
        if name:
            name = Path(name).stem
        names.append(name)

    n = len(names)
    if n < 3:
        print_message("Fewer than 3 samples, skipping scatter plot.", "warning")
        return out_path

    dist_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            try:
                d = tree.distance(terminals[i], terminals[j])
            except Exception:
                d = 0.0
            dist_matrix[i, j] = d
            dist_matrix[j, i] = d

    mds = MDS(
        n_components=2,
        metric="precomputed",
        random_state=42,
        normalized_stress="auto",
        n_init=4,
        init="random",
    )
    coords = mds.fit_transform(dist_matrix)

    amr_cols = [c for c in binary_matrix.columns if c.startswith("AMR:")]
    counts = []
    for name in names:
        if name in binary_matrix.index and amr_cols:
            counts.append(int(binary_matrix.loc[name, amr_cols].sum()))
        else:
            counts.append(0)

    max_count = max(counts) if counts and max(counts) > 0 else 1

    fig, ax = plt.subplots(figsize=(10, 8))
    cmap = plt.cm.YlOrRd
    norm = matplotlib.colors.Normalize(vmin=0, vmax=max_count)

    for i, name in enumerate(names):
        c = cmap(norm(counts[i]))
        is_sel = name in selected_ids
        ax.scatter(
            coords[i, 0], coords[i, 1],
            c=[c], marker="*" if is_sel else "o",
            s=200 if is_sel else 60,
            zorder=10 if is_sel else 5,
            edgecolors="black" if is_sel else "grey",
            linewidths=1.5 if is_sel else 0.5,
        )
        if is_sel:
            ax.annotate(name, (coords[i, 0], coords[i, 1]),
                        fontsize=7, ha="left", va="bottom",
                        xytext=(5, 5), textcoords="offset points")

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.7)
    cbar.set_label("AMR gene count")

    ax.set_xlabel("MDS dimension 1")
    ax.set_ylabel("MDS dimension 2")
    ax.set_title("PCoA of Mash distances (selected samples starred)")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print_message(f"Scatter plot saved to {out_path}", "success")
    return out_path


def plot_tree_heatmap(
    tree_path: str,
    binary_matrix: pd.DataFrame,
    selected_ids: list[str],
    output_dir: str,
) -> str:
    """Phylogenetic tree (left) + AMR/replicon heatmap (right)."""
    out_path = os.path.join(output_dir, "tree_heatmap.png")

    try:
        tree = Phylo.read(tree_path, "newick")
    except Exception as e:
        print_message(f"Cannot parse tree for tree+heatmap plot: {e}", "warning")
        return out_path

    leaf_order = _get_leaf_order(tree)
    n_leaves = len(leaf_order)

    if n_leaves == 0:
        print_message("No leaves in tree, skipping tree+heatmap plot.", "warning")
        return out_path

    missing = [s for s in leaf_order if s not in binary_matrix.index]
    if missing:
        print_message(f"{len(missing)} samples in tree not found in AMR/replicon matrix.", "warning")

    if binary_matrix.empty or len(binary_matrix.columns) == 0:
        print_message("No features to plot in heatmap.", "warning")
        return out_path

    # Sort columns: AMR classes (alphabetical by class then gene), then REP
    amr_cols = sorted([c for c in binary_matrix.columns if c.startswith("AMR:")])
    rep_cols = sorted([c for c in binary_matrix.columns if c.startswith("REP:")])
    ordered_cols = amr_cols + rep_cols
    n_features = len(ordered_cols)

    # Expand matrix to cover all tree leaves (NaN for leaves not in binary_matrix)
    heatmap_data = pd.DataFrame(np.nan, index=leaf_order, columns=ordered_cols, dtype=float)
    for sid in binary_matrix.index:
        if sid in heatmap_data.index:
            heatmap_data.loc[sid, ordered_cols] = (
                binary_matrix.loc[sid, ordered_cols].values.astype(float)
            )

    # Figure dimensions
    row_h = 0.22
    col_w = 0.22
    tree_w = 2.8
    bar_h = 0.45
    legend_w = 2.2
    heat_w = max(4.0, n_features * col_w)
    heat_h = max(4.0, n_leaves * row_h)

    fig_w = tree_w + heat_w + legend_w + 0.8
    fig_h = bar_h + heat_h + 0.6

    fig = plt.figure(figsize=(fig_w, fig_h))

    gs = gridspec.GridSpec(
        2, 3,
        height_ratios=[bar_h, heat_h],
        width_ratios=[tree_w, heat_w, legend_w],
        hspace=0.005,
        wspace=0.12,
        figure=fig,
    )

    ax_tree = fig.add_subplot(gs[:, 0])
    ax_bar = fig.add_subplot(gs[0, 1])
    ax_heat = fig.add_subplot(gs[1, 1])
    ax_leg = fig.add_subplot(gs[:, 2])

    # Tree
    _draw_tree_axes(tree, ax_tree, leaf_order, selected_ids, n_leaves)

    # Drug class colour bar
    drug_classes = [_col_drug_class(c) for c in ordered_cols]
    _draw_drug_class_bar(ax_bar, drug_classes, n_features)

    # Heatmap
    # Encode: NaN -> -1 (no data), 0 -> 0 (absent), 1 -> 1 (present)
    plot_vals = np.where(np.isnan(heatmap_data.values), -1.0, heatmap_data.values)

    cmap_heat = matplotlib.colors.ListedColormap(["#f5f5f5", "#d0d0d0", "#f46d43"])
    norm_heat = matplotlib.colors.BoundaryNorm([-1.5, -0.5, 0.5, 1.5], cmap_heat.N)

    ax_heat.imshow(plot_vals, cmap=cmap_heat, norm=norm_heat,
                   aspect="auto", interpolation="none")

    # Grid lines
    ax_heat.set_xticks(np.arange(n_features + 1) - 0.5, minor=True)
    ax_heat.set_yticks(np.arange(n_leaves + 1) - 0.5, minor=True)
    ax_heat.grid(which="minor", color="white", linewidth=0.3)
    ax_heat.tick_params(which="minor", bottom=False, left=False)

    # Column labels
    display_cols = [c.split(":")[-1] for c in ordered_cols]
    ax_heat.set_xticks(range(n_features))
    ax_heat.set_xticklabels(display_cols, rotation=90, fontsize=6, ha="center")
    ax_heat.tick_params(axis="x", bottom=False, top=True, labelbottom=False, labeltop=True)
    ax_heat.xaxis.set_label_position("top")

    # Row labels
    ax_heat.set_yticks(range(n_leaves))
    ax_heat.set_yticklabels(leaf_order, fontsize=6)
    ax_heat.yaxis.tick_left()
    for i, tick_label in enumerate(ax_heat.get_yticklabels()):
        name = leaf_order[i]
        if name in selected_ids:
            tick_label.set_color("red")
            tick_label.set_fontweight("bold")
        elif name not in binary_matrix.index:
            tick_label.set_color("#aaaaaa")

    # Legend
    _draw_legend(ax_leg, drug_classes)

    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print_message(f"Tree+heatmap plot saved to {out_path}", "success")
    return out_path


def plot_pareto(pareto_tsv: str, output_dir: str) -> str:
    """Plot Pareto curve from sweep results."""
    out_path = os.path.join(output_dir, "pareto_plot.png")
    df = pd.read_csv(pareto_tsv, sep="\t")

    fig, ax = plt.subplots(figsize=(10, 6))

    if "pct_amr_covered" in df.columns and "pct_replicons_covered" in df.columns:
        ax.plot(df["alpha"], df["pct_amr_covered"], "o-", label="AMR coverage %", color="#d7191c")
        ax.plot(
            df["alpha"], df["pct_replicons_covered"], "s-",
            label="Replicon coverage %", color="#2c7bb6",
        )
    if "pct_faith_pd" in df.columns:
        ax.plot(df["alpha"], df["pct_faith_pd"], "^-", label="Faith PD %", color="#1a9641")

    ax.set_xlabel("Alpha (phylo budget fraction)")
    ax.set_ylabel("Coverage (%)")
    ax.set_title("Pareto curve: coverage vs alpha")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(0, 105)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print_message(f"Pareto plot saved to {out_path}", "success")
    return out_path


def plot_nsga3_front(pareto_df: pd.DataFrame, output_dir: str) -> str:
    """3-panel pairwise scatter of the NSGA-III Pareto front."""
    out_path = os.path.join(output_dir, "pareto_front.png")

    pairs = [
        ("minimax_dist", "amr_profile_div", "Minimax distance", "AMR profile diversity (Jaccard)"),
        ("minimax_dist", "rep_profile_div", "Minimax distance", "Replicon profile diversity (Jaccard)"),
        ("amr_profile_div", "rep_profile_div", "AMR profile diversity (Jaccard)", "Replicon profile diversity (Jaccard)"),
    ]

    rec_mask = pareto_df["is_recommended"].astype(bool)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    for ax, (xcol, ycol, xlabel, ylabel) in zip(axes, pairs, strict=True):
        ax.scatter(
            pareto_df[xcol][~rec_mask],
            pareto_df[ycol][~rec_mask],
            c="grey", s=30, alpha=0.6, edgecolors="none",
        )
        if rec_mask.any():
            ax.scatter(
                pareto_df[xcol][rec_mask],
                pareto_df[ycol][rec_mask],
                c="red", marker="*", s=200, zorder=10,
            )
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)

    fig.suptitle("NSGA-III Pareto Front", fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print_message(f"Pareto front plot saved to {out_path}", "success")
    return out_path


def plot_nsga3_parallel(pareto_df: pd.DataFrame, output_dir: str) -> str:
    """Parallel coordinates plot of the NSGA-III Pareto front."""
    out_path = os.path.join(output_dir, "pareto_parallel.png")

    obj_cols = ["minimax_dist", "amr_profile_div", "rep_profile_div"]
    labels = ["Minimax dist", "AMR diversity", "Replicon diversity"]

    vals = pareto_df[obj_cols].values.copy()
    col_min = vals.min(axis=0)
    col_max = vals.max(axis=0)
    col_range = col_max - col_min
    col_range[col_range == 0] = 1.0
    normed = (vals - col_min) / col_range

    rec_mask = pareto_df["is_recommended"].values.astype(bool)
    x_positions = np.arange(len(obj_cols))

    fig, ax = plt.subplots(figsize=(7, 5))

    # Non-recommended solutions
    for i in range(len(normed)):
        if not rec_mask[i]:
            ax.plot(x_positions, normed[i], color="lightgrey", lw=0.8, alpha=0.5)

    # Recommended solution on top
    for i in range(len(normed)):
        if rec_mask[i]:
            ax.plot(x_positions, normed[i], color="red", lw=2.5, zorder=10)

    ax.set_xticks(x_positions)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Normalised value")
    ax.set_ylim(-0.05, 1.05)
    ax.set_title("NSGA-III Pareto Front -- Parallel Coordinates", fontweight="bold")
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print_message(f"Parallel coordinates plot saved to {out_path}", "success")
    return out_path


def plot_diversity_curve(curve_df: pd.DataFrame, output_dir: str) -> str:
    """Plot a 3-curve diversity saturation (elbow) plot.

    Each curve is normalised to 0--100% of its own maximum so all three
    criteria are visually comparable on the same y-axis.

    Parameters
    ----------
    curve_df : pd.DataFrame
        Must have columns: k, phylo_pct, amr_jaccard_div, rep_jaccard_div.
    output_dir : str
        Directory in which to save diversity_curve.png.

    Returns
    -------
    str
        Path to the saved PNG.
    """
    out_path = os.path.join(output_dir, "diversity_curve.png")

    fig, ax = plt.subplots(figsize=(9, 5.5))

    curves = [
        ("phylo_pct", "Phylogenetic (PARNAS)", "#2c7bb6", "o"),
        ("amr_jaccard_div", "AMR profile (Jaccard)", "#d7191c", "s"),
        ("rep_jaccard_div", "Replicon profile (Jaccard)", "#1a9641", "^"),
    ]

    for col, label, colour, marker in curves:
        vals = curve_df[col].copy()
        if vals.isna().all():
            continue
        col_max = vals.max()
        normed = vals / col_max * 100.0 if col_max > 0 else vals * 0.0
        ax.plot(
            curve_df["k"],
            normed,
            marker=marker,
            color=colour,
            linewidth=2,
            markersize=5,
            label=label,
        )

    ax.axhline(y=90, color="grey", linestyle="--", linewidth=1, alpha=0.7, label="90% threshold")

    ax.set_xlabel("k (number of representatives)")
    ax.set_ylabel("% of diversity at k=max")
    ax.set_ylim(0, 105)
    ax.set_title("Diversity saturation curve — how many representatives are needed?")
    ax.legend(loc="lower right")
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print_message(f"Diversity curve plot saved to {out_path}", "success")
    return out_path

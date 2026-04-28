"""All plotting functions for dustmeselecta."""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import Phylo
from sklearn.manifold import MDS

import click


def plot_elbow(diversity_csv: str, n_chosen: int, output_dir: str) -> str:
    """Plot diversity covered vs n, with vertical line at chosen n."""
    out_path = os.path.join(output_dir, "elbow_plot.png")
    try:
        df = pd.read_csv(diversity_csv)
    except Exception as e:
        click.echo(f"Warning: cannot read diversity CSV for elbow plot: {e}")
        return out_path

    # PARNAS diversity CSV typically has columns: k, diversity (or similar)
    # Try to detect the columns
    cols = df.columns.tolist()
    if len(cols) >= 2:
        x_col, y_col = cols[0], cols[1]
    else:
        click.echo("Warning: diversity CSV has fewer than 2 columns, skipping elbow plot.")
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
    click.echo(f"Elbow plot saved to {out_path}")
    return out_path


def plot_scatter(
    tree_path: str,
    binary_matrix: pd.DataFrame,
    selected_ids: list[str],
    output_dir: str,
) -> str:
    """PCoA of Mash distances, coloured by dominant resistance class, selected starred."""
    out_path = os.path.join(output_dir, "scatter_plot.png")

    # Parse tree to get distance matrix
    try:
        tree = Phylo.read(tree_path, "newick")
    except Exception as e:
        click.echo(f"Warning: cannot parse tree for scatter plot: {e}")
        return out_path

    # Get all terminal names
    terminals = [t for t in tree.get_terminals()]
    names = []
    for t in terminals:
        name = t.name
        if name:
            # Strip file extension if present
            from pathlib import Path
            name = Path(name).stem
        names.append(name)

    n = len(names)
    if n < 3:
        click.echo("Warning: fewer than 3 samples, skipping scatter plot.")
        return out_path

    # Build distance matrix from tree
    dist_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            try:
                d = tree.distance(terminals[i], terminals[j])
            except Exception:
                d = 0.0
            dist_matrix[i, j] = d
            dist_matrix[j, i] = d

    # MDS for 2D embedding (PCoA-like)
    mds = MDS(
        n_components=2,
        metric="precomputed",
        random_state=42,
        normalized_stress="auto",
        n_init=4,
        init="random",
    )
    coords = mds.fit_transform(dist_matrix)

    # Determine dominant resistance class per sample
    amr_cols = [c for c in binary_matrix.columns if c.startswith("AMR:")]
    # Group by resistance class prefix (the gene column name)
    class_counts = {}
    for sid in names:
        if sid in binary_matrix.index:
            row = binary_matrix.loc[sid]
            class_counts[sid] = sum(row[amr_cols]) if amr_cols else 0
        else:
            class_counts[sid] = 0

    # Colour by total AMR gene count (binned)
    counts = [class_counts.get(name, 0) for name in names]
    max_count = max(counts) if counts and max(counts) > 0 else 1

    fig, ax = plt.subplots(figsize=(10, 8))
    cmap = plt.cm.YlOrRd
    norm = matplotlib.colors.Normalize(vmin=0, vmax=max_count)

    for i, name in enumerate(names):
        c = cmap(norm(counts[i]))
        is_selected = name in selected_ids
        marker = "*" if is_selected else "o"
        size = 200 if is_selected else 60
        zorder = 10 if is_selected else 5
        ax.scatter(
            coords[i, 0], coords[i, 1],
            c=[c], marker=marker, s=size, zorder=zorder,
            edgecolors="black" if is_selected else "grey",
            linewidths=1.5 if is_selected else 0.5,
        )
        if is_selected:
            ax.annotate(
                name, (coords[i, 0], coords[i, 1]),
                fontsize=7, ha="left", va="bottom",
                xytext=(5, 5), textcoords="offset points",
            )

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
    click.echo(f"Scatter plot saved to {out_path}")
    return out_path


def _get_leaf_order(tree) -> list[str]:
    """Get leaf names in tree traversal order (topology order)."""
    from pathlib import Path
    names = []
    for clade in tree.get_terminals():
        name = clade.name if clade.name else ""
        name = Path(name).stem
        names.append(name)
    return names


def plot_tree_heatmap(
    tree_path: str,
    binary_matrix: pd.DataFrame,
    selected_ids: list[str],
    output_dir: str,
) -> str:
    """Phylogenetic tree (left) + AMR gene/replicon heatmap (right), selected highlighted."""
    out_path = os.path.join(output_dir, "tree_heatmap.png")

    try:
        tree = Phylo.read(tree_path, "newick")
    except Exception as e:
        click.echo(f"Warning: cannot parse tree for tree+heatmap plot: {e}")
        return out_path

    leaf_order = _get_leaf_order(tree)
    n_leaves = len(leaf_order)

    if n_leaves == 0:
        click.echo("Warning: no leaves in tree, skipping tree+heatmap plot.")
        return out_path

    # Reorder binary matrix to match tree leaf order
    available = [s for s in leaf_order if s in binary_matrix.index]
    missing = [s for s in leaf_order if s not in binary_matrix.index]
    if missing:
        click.echo(f"Warning: {len(missing)} samples in tree not found in Kleborate matrix.")

    # Only plot samples present in both tree and matrix
    plot_order = [s for s in leaf_order if s in binary_matrix.index]
    if not plot_order:
        click.echo("Warning: no overlapping samples between tree and matrix.")
        return out_path

    heatmap_data = binary_matrix.loc[plot_order]

    # Determine figure dimensions
    n_features = len(heatmap_data.columns)
    fig_width = max(12, 6 + n_features * 0.3)
    fig_height = max(6, n_leaves * 0.3 + 2)

    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, max(1, n_features * 0.15)], wspace=0.02)

    # Left panel: tree as dendrogram
    ax_tree = fig.add_subplot(gs[0])
    ax_tree.set_xlim(0, 1)
    ax_tree.set_ylim(0, n_leaves)

    # Draw tree using Bio.Phylo
    Phylo.draw(
        tree, do_show=False, axes=ax_tree,
        label_func=lambda c: "",  # We add labels manually
    )

    # Adjust y-axis to align with heatmap
    ax_tree.set_ylabel("")
    ax_tree.set_xlabel("")
    ax_tree.set_yticks([])

    # Get y positions of terminals from the drawn tree
    # Bio.Phylo.draw assigns y positions 0..n-1 to terminals
    y_positions = {}
    for i, name in enumerate(leaf_order):
        y_positions[name] = i

    # Add labels with highlighting for selected samples
    for name in plot_order:
        y_pos = y_positions.get(name, 0)
        colour = "red" if name in selected_ids else "black"
        weight = "bold" if name in selected_ids else "normal"
        # Add a red marker for selected samples
        if name in selected_ids:
            ax_tree.plot(
                ax_tree.get_xlim()[1] * 0.95, y_pos + 1,
                "o", color="red", markersize=5, zorder=10,
            )

    ax_tree.spines["top"].set_visible(False)
    ax_tree.spines["right"].set_visible(False)

    # Right panel: heatmap
    ax_heat = fig.add_subplot(gs[1])

    if n_features > 0:
        # Clean up feature names for display
        display_cols = []
        for c in heatmap_data.columns:
            if c.startswith("AMR:"):
                display_cols.append(c[4:])
            elif c.startswith("REP:"):
                display_cols.append(c[4:])
            else:
                display_cols.append(c)

        sns.heatmap(
            heatmap_data.values,
            ax=ax_heat,
            cmap="YlOrRd",
            cbar=False,
            xticklabels=display_cols,
            yticklabels=False,
            linewidths=0.5,
            linecolor="white",
        )
        ax_heat.set_xticklabels(ax_heat.get_xticklabels(), rotation=90, fontsize=7)

        # Add sample labels on right side with highlighting
        ax_heat.set_yticks(np.arange(len(plot_order)) + 0.5)
        labels = []
        for name in plot_order:
            labels.append(name)
        ax_heat.set_yticklabels(labels, fontsize=7)
        ax_heat.yaxis.tick_right()

        # Bold + red for selected samples
        for i, label in enumerate(ax_heat.get_yticklabels()):
            if plot_order[i] in selected_ids:
                label.set_color("red")
                label.set_fontweight("bold")
    else:
        ax_heat.text(0.5, 0.5, "No AMR/replicon features detected",
                     ha="center", va="center", transform=ax_heat.transAxes)

    ax_heat.set_title("AMR genes & replicon types")

    fig.suptitle("Phylogenetic tree with AMR/replicon heatmap", fontsize=12, y=1.01)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    click.echo(f"Tree+heatmap plot saved to {out_path}")
    return out_path


def plot_pareto(pareto_tsv: str, output_dir: str) -> str:
    """Plot Pareto curve from sweep results."""
    out_path = os.path.join(output_dir, "pareto_plot.png")
    df = pd.read_csv(pareto_tsv, sep="\t")

    fig, ax = plt.subplots(figsize=(10, 6))

    if "pct_amr_covered" in df.columns and "pct_replicons_covered" in df.columns:
        ax.plot(df["alpha"], df["pct_amr_covered"], "o-", label="AMR coverage %", color="#d7191c")
        ax.plot(df["alpha"], df["pct_replicons_covered"], "s-", label="Replicon coverage %", color="#2c7bb6")
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
    click.echo(f"Pareto plot saved to {out_path}")
    return out_path

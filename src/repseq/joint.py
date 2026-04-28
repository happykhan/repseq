"""Joint k-medoids selection on a blended phylogenetic + AMR distance matrix."""

from __future__ import annotations

from pathlib import Path

import numpy as np


def tree_to_dist_matrix(tree, names: list[str]) -> np.ndarray:  # type: ignore[type-arg]
    """Build an (n x n) pairwise tip-to-tip distance matrix from a BioPython Phylo tree.

    Parameters
    ----------
    tree : Bio.Phylo tree object
        Parsed Newick/Nexus tree.
    names : list[str]
        Sample names (stems) defining row/column order.

    Returns
    -------
    np.ndarray
        Symmetric distance matrix with shape (len(names), len(names)).
    """
    terminals = list(tree.get_terminals())
    stem_to_terminal = {Path(t.name).stem: t for t in terminals if t.name}

    n = len(names)
    dist = np.zeros((n, n))
    for i in range(n):
        ti = stem_to_terminal.get(names[i])
        if ti is None:
            continue
        for j in range(i + 1, n):
            tj = stem_to_terminal.get(names[j])
            if tj is None:
                continue
            try:
                d = tree.distance(ti, tj)
            except Exception:
                d = 0.0
            dist[i, j] = d
            dist[j, i] = d
    return dist


def compute_joint_dist(tree_dm: np.ndarray, amr_dm: np.ndarray, w: float) -> np.ndarray:
    """Blend normalised phylogenetic and AMR distance matrices.

    Parameters
    ----------
    tree_dm : np.ndarray
        Phylogenetic distance matrix.
    amr_dm : np.ndarray
        AMR (Jaccard) distance matrix.
    w : float
        Weight for AMR distance (0 = pure phylo, 1 = pure AMR).

    Returns
    -------
    np.ndarray
        Blended distance matrix with values in [0, 1].
    """
    tree_max = tree_dm.max()
    amr_max = amr_dm.max()

    tree_norm = tree_dm / tree_max if tree_max > 0 else tree_dm.copy()
    amr_norm = amr_dm / amr_max if amr_max > 0 else amr_dm.copy()

    return (1 - w) * tree_norm + w * amr_norm


def run_kmedoids(dist: np.ndarray, k: int, seed: int = 42) -> list[int]:
    """K-medoids clustering with k-means++ initialisation and PAM-style iteration.

    Parameters
    ----------
    dist : np.ndarray
        Symmetric (n x n) distance/dissimilarity matrix.
    k : int
        Number of medoids to select.
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    list[int]
        Indices (row positions) of the k selected medoids.
    """
    n = dist.shape[0]
    if k >= n:
        return list(range(n))

    rng = np.random.default_rng(seed)

    # k-means++ initialisation
    medoids = [int(rng.integers(0, n))]
    for _ in range(k - 1):
        # Squared distance to nearest existing medoid
        min_dists = np.min(dist[:, medoids], axis=1) ** 2
        # Zero out already-chosen medoids
        min_dists[medoids] = 0.0
        total = min_dists.sum()
        if total == 0:
            # All remaining points are equidistant (or identical); pick uniformly
            remaining = [i for i in range(n) if i not in medoids]
            medoids.append(int(rng.choice(remaining)))
        else:
            probs = min_dists / total
            chosen = int(rng.choice(n, p=probs))
            medoids.append(chosen)

    # PAM-style iteration (max 100 rounds)
    for _ in range(100):
        # Assign each point to nearest medoid
        assignments = np.argmin(dist[:, medoids], axis=1)

        # For each cluster, find the member minimising sum of intra-cluster distances
        new_medoids = []
        for ci in range(k):
            cluster_members = np.where(assignments == ci)[0]
            if len(cluster_members) == 0:
                # Empty cluster: keep old medoid
                new_medoids.append(medoids[ci])
                continue
            # Sum of distances from each member to all other members in the cluster
            intra_dists = dist[np.ix_(cluster_members, cluster_members)]
            total_dists = intra_dists.sum(axis=1)
            best_local = int(cluster_members[np.argmin(total_dists)])
            new_medoids.append(best_local)

        if sorted(new_medoids) == sorted(medoids):
            break
        medoids = new_medoids

    return medoids

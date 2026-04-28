"""NSGA-III multi-objective optimisation for representative isolate selection."""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo
from pymoo.algorithms.moo.nsga3 import NSGA3
from pymoo.core.crossover import Crossover
from pymoo.core.mutation import Mutation
from pymoo.core.problem import ElementwiseProblem
from pymoo.core.sampling import Sampling
from pymoo.optimize import minimize
from pymoo.util.ref_dirs import get_reference_directions

from repseq.amr_cover import build_feature_matrix
from repseq.joint import tree_to_dist_matrix
from repseq.log import print_header, print_message
from repseq.phylo import run_mashtree
from repseq.plots import plot_nsga3_front, plot_nsga3_parallel


class _RepseqProblem(ElementwiseProblem):
    def __init__(self, dist_matrix: np.ndarray, amr_matrix: np.ndarray, rep_matrix: np.ndarray, k: int):
        self.D = dist_matrix
        self.A = amr_matrix
        self.R = rep_matrix
        self.k = k
        self.n_total = dist_matrix.shape[0]
        super().__init__(
            n_var=k,
            n_obj=3,
            xl=np.zeros(k, dtype=int),
            xu=np.full(k, self.n_total - 1, dtype=int),
            vtype=int,
        )

    def _evaluate(self, x, out, *args, **kwargs):
        idx = list(x)
        # f1: minimax distance (minimise directly)
        unselected = list(set(range(self.n_total)) - set(idx))
        if unselected:
            min_dists = self.D[np.ix_(unselected, idx)].min(axis=1)
            f1 = float(min_dists.max())
        else:
            f1 = 0.0
        # f2: AMR coverage (negate for minimisation)
        f2 = -float((self.A[idx, :].sum(axis=0) > 0).mean()) if self.A.shape[1] > 0 else -1.0
        # f3: replicon coverage (negate for minimisation)
        f3 = -float((self.R[idx, :].sum(axis=0) > 0).mean()) if self.R.shape[1] > 0 else -1.0
        out["F"] = [f1, f2, f3]


class _SubsetSampling(Sampling):
    def _do(self, problem, n_samples, **kwargs):
        return np.array([
            np.random.choice(problem.n_total, problem.k, replace=False)
            for _ in range(n_samples)
        ])


class _SubsetCrossover(Crossover):
    def __init__(self):
        super().__init__(2, 2)

    def _do(self, problem, X, **kwargs):
        n_matings = X.shape[1]
        k, n_total = problem.k, problem.n_total
        Y = np.zeros_like(X)
        for i in range(n_matings):
            pool = list(set(X[0, i].tolist()) | set(X[1, i].tolist()))
            if len(pool) < k:
                extra = list(set(range(n_total)) - set(pool))
                np.random.shuffle(extra)
                pool += extra[: k - len(pool)]
            np.random.shuffle(pool)
            Y[0, i] = pool[:k]
            np.random.shuffle(pool)
            Y[1, i] = pool[:k]
        return Y


class _SubsetMutation(Mutation):
    def __init__(self, prob: float = 0.3):
        super().__init__()
        self.prob = prob

    def _do(self, problem, X, **kwargs):
        n_total = problem.n_total
        for i in range(X.shape[0]):
            if np.random.random() < self.prob:
                selected = list(X[i])
                unselected = list(set(range(n_total)) - set(selected))
                if unselected:
                    j = np.random.randint(len(selected))
                    selected[j] = int(np.random.choice(unselected))
                    X[i] = selected
        return X


def select_display_solutions(pareto_df: pd.DataFrame, rec_idx: int) -> list[int]:
    """Select 5 anchor solutions from the Pareto front for display.

    Returns solution indices for:
    1. Recommended (closest to ideal)
    2. Best AMR coverage
    3. Best replicon coverage
    4. Best phylo (lowest minimax_dist)
    5. Gap filler (max-min distance to the other 4 in normalised space)
    """
    chosen: list[int] = [rec_idx]

    best_amr = int(pareto_df["pct_amr_covered"].idxmax())
    if best_amr not in chosen:
        chosen.append(best_amr)

    best_rep = int(pareto_df["pct_rep_covered"].idxmax())
    if best_rep not in chosen:
        chosen.append(best_rep)

    best_phylo = int(pareto_df["minimax_dist"].idxmin())
    if best_phylo not in chosen:
        chosen.append(best_phylo)

    if len(chosen) < 5:
        obj = pareto_df[["minimax_dist", "pct_amr_covered", "pct_rep_covered"]].values.copy()
        for col_idx in range(obj.shape[1]):
            col_min, col_max = obj[:, col_idx].min(), obj[:, col_idx].max()
            rng = col_max - col_min
            if rng > 0:
                obj[:, col_idx] = (obj[:, col_idx] - col_min) / rng
            else:
                obj[:, col_idx] = 0.5
        # Invert minimax so higher = better
        obj[:, 0] = 1.0 - obj[:, 0]

        remaining = [i for i in range(len(pareto_df)) if i not in chosen]
        if remaining:
            chosen_obj = obj[chosen]
            best_gap_idx = None
            best_gap_val = -1.0
            for ri in remaining:
                dists_to_chosen = np.linalg.norm(obj[ri] - chosen_obj, axis=1)
                min_dist = float(dists_to_chosen.min())
                if min_dist > best_gap_val:
                    best_gap_val = min_dist
                    best_gap_idx = ri
            if best_gap_idx is not None:
                chosen.append(best_gap_idx)

    # Fill remaining slots if duplicates collapsed the list below 5
    if len(chosen) < 5:
        remaining = [i for i in range(len(pareto_df)) if i not in chosen]
        for ri in remaining:
            chosen.append(ri)
            if len(chosen) >= 5:
                break

    return chosen[:5]


def run_nsga3(
    assemblies_dir: str,
    tree_path: str | None,
    kleborate_path: str | None,
    plasmidfinder_path: str | None,
    hamronization_path: str | None,
    abricate_replicons_path: str | None,
    n: int,
    output_dir: str,
    pop_size: int = 100,
    n_gen: int = 300,
    seed: int = 42,
) -> list[str]:
    """Run NSGA-III multi-objective selection.

    Optimises three objectives simultaneously: minimax phylogenetic
    distance, AMR gene coverage, and replicon type coverage.

    Returns list of sample IDs from the recommended Pareto-optimal solution.
    """
    os.makedirs(output_dir, exist_ok=True)

    print_header("NSGA-III Multi-Objective Selection")

    if tree_path is None:
        tree_path = run_mashtree(assemblies_dir, output_dir)

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

    tree = Phylo.read(tree_path, "newick")
    tree_terminals = {Path(t.name).stem for t in tree.get_terminals() if t.name}
    names = sorted(set(binary_matrix.index) & tree_terminals)

    if len(names) < n:
        print_message(
            f"Only {len(names)} samples in both tree and matrix (requested {n}); "
            f"using all {len(names)}.",
            "warning",
        )
        n = min(n, len(names))

    if len(names) < 2:
        print_message("Fewer than 2 samples -- cannot run NSGA-III.", "error")
        return []

    if len(names) > 200:
        print_message(
            f"Building {len(names)}x{len(names)} tree distance matrix -- this may be slow",
            "warning",
        )
    tree_dm = tree_to_dist_matrix(tree, names)

    amr_cols = [c for c in binary_matrix.columns if c.startswith("AMR:")]
    rep_cols = [c for c in binary_matrix.columns if c.startswith("REP:")]
    amr_matrix = binary_matrix.loc[names, amr_cols].values if amr_cols else np.empty((len(names), 0))
    rep_matrix = binary_matrix.loc[names, rep_cols].values if rep_cols else np.empty((len(names), 0))

    print_message(
        f"Problem: select {n} from {len(names)} samples | "
        f"{len(amr_cols)} AMR features, {len(rep_cols)} replicon features",
        "info",
    )

    problem = _RepseqProblem(tree_dm, amr_matrix, rep_matrix, n)
    ref_dirs = get_reference_directions("das-dennis", 3, n_partitions=12)

    algorithm = NSGA3(
        ref_dirs=ref_dirs,
        pop_size=pop_size,
        sampling=_SubsetSampling(),
        crossover=_SubsetCrossover(),
        mutation=_SubsetMutation(prob=0.3),
    )

    print_message(
        f"Running NSGA-III: pop_size={pop_size}, generations={n_gen}, seed={seed}",
        "info",
    )
    result = minimize(
        problem,
        algorithm,
        ("n_gen", n_gen),
        seed=seed,
        verbose=False,
    )

    opt = result.opt
    if opt is None or len(opt) == 0:
        print_message("NSGA-III returned no optimal solutions -- falling back to final population", "warning")
        opt = result.pop
        if opt is None or len(opt) == 0:
            print_message("No solutions found.", "error")
            return []

    front_X = opt.get("X")
    raw_F = opt.get("F")

    # f1 (minimax) is already a minimisation target (un-negated)
    # f2, f3 (AMR, replicon coverage) were negated -- un-negate them
    front_minimax = raw_F[:, 0]
    front_amr = -raw_F[:, 1]
    front_rep = -raw_F[:, 2]

    print_message(f"Pareto front size: {len(raw_F)} solutions", "success")

    # Normalise objectives for closest-to-ideal recommendation
    front_obj = np.column_stack([front_minimax, front_amr * 100, front_rep * 100])
    obj = front_obj.copy()

    # minimax: lower is better -- invert after normalising
    minimax_vals = obj[:, 0]
    mini_min, mini_max = minimax_vals.min(), minimax_vals.max()
    if mini_max > mini_min:
        obj[:, 0] = 1.0 - (minimax_vals - mini_min) / (mini_max - mini_min)
    else:
        obj[:, 0] = 1.0
    # AMR and replicon are 0-100; normalise to 0-1
    obj[:, 1] = (obj[:, 1] - obj[:, 1].min()) / max(obj[:, 1].max() - obj[:, 1].min(), 1e-9)
    obj[:, 2] = (obj[:, 2] - obj[:, 2].min()) / max(obj[:, 2].max() - obj[:, 2].min(), 1e-9)

    ideal = np.ones(3)
    dists_to_ideal = np.linalg.norm(obj - ideal, axis=1)
    rec_idx = int(np.argmin(dists_to_ideal))

    # Build pareto_front.tsv
    rows = []
    for si in range(len(raw_F)):
        sample_indices = [int(idx) for idx in front_X[si]]
        sample_names = [names[idx] for idx in sample_indices]
        rows.append({
            "solution_id": si,
            "minimax_dist": float(front_minimax[si]),
            "pct_amr_covered": round(float(front_amr[si]) * 100, 2),
            "pct_rep_covered": round(float(front_rep[si]) * 100, 2),
            "is_recommended": si == rec_idx,
            "selected_samples": ";".join(sorted(sample_names)),
        })

    pareto_df = pd.DataFrame(rows)
    pareto_path = os.path.join(output_dir, "pareto_front.tsv")
    pareto_df.to_csv(pareto_path, sep="\t", index=False)
    print_message(f"Pareto front written to {pareto_path}", "success")

    rec_samples = sorted([names[int(idx)] for idx in front_X[rec_idx]])
    rec_path = os.path.join(output_dir, "recommended.txt")
    with open(rec_path, "w") as fh:
        for sid in rec_samples:
            fh.write(sid + "\n")
    print_message(f"Recommended solution written to {rec_path}", "success")

    plot_nsga3_front(pareto_df, output_dir)
    plot_nsga3_parallel(pareto_df, output_dir)

    print_header("NSGA-III Results Summary")
    rec_row = pareto_df.iloc[rec_idx]
    print_message(f"Pareto front: {len(pareto_df)} solutions", "info")
    print_message(
        f"Recommended solution (closest to ideal): "
        f"minimax dist={rec_row['minimax_dist']:.5f}, "
        f"AMR={rec_row['pct_amr_covered']:.1f}%, "
        f"replicons={rec_row['pct_rep_covered']:.1f}%",
        "success",
    )
    print_message(f"Selected {len(rec_samples)} samples: {rec_samples}", "info")

    return rec_samples

# EuSCAPE NSGA-III Example

Multi-objective selection of 20 representative *K. pneumoniae* isolates from the
EuSCAPE collection (58 samples) using NSGA-III.

## Data

- **Tree**: `../euscape/inputs/tree.nwk` (Mash-distance tree, 58 tips)
- **AMR features**: `../euscape/inputs/kleborate.tsv` (Kleborate v3, 87 AMR gene features)
- **Replicon features**: none available (no assemblies in repo; replicon coverage is trivially 100%)

## Command

```bash
repseq nsga2 \
  --assemblies <assemblies_dir> \
  --tree examples/euscape/inputs/tree.nwk \
  --kleborate examples/euscape/inputs/kleborate.tsv \
  --n 20 \
  --pop-size 100 \
  --generations 300 \
  --seed 42 \
  --output-dir examples/euscape_nsga2
```

## Results

| Metric | Value |
|--------|-------|
| Pareto front size | 10 solutions |
| Recommended: mean pairwise distance | 0.0157 |
| Recommended: AMR coverage | 94.3% |
| Recommended: replicon coverage | 100.0% (no REP features) |

## Comparison with `repseq select` (split method, n=20)

The `repseq select` run (alpha=0.5, split method) on the same collection achieved:

| Metric | `select` (split) | `nsga2` (recommended) |
|--------|-------------------|-----------------------|
| AMR coverage | 89.7% | 94.3% |
| Replicon coverage | 91.7% | 100.0%* |

\*No replicon features in NSGA-III run (kleborate-only input); the `select` run
had ABRicate replicon data from assemblies.

The NSGA-III approach finds a Pareto front of 10 non-dominated solutions,
letting the user inspect trade-offs between phylogenetic diversity, AMR gene
coverage, and replicon coverage. The recommended solution (closest to ideal
point) achieves higher AMR coverage than the split method, though a direct
comparison is limited by the absence of replicon data in this example.

## Output files

- `pareto_front.tsv` -- all Pareto-optimal solutions with objectives and sample lists
- `recommended.txt` -- sample IDs from the recommended solution
- `pareto_front.png` -- 3-panel pairwise scatter of Pareto front
- `pareto_parallel.png` -- parallel coordinates plot of Pareto front

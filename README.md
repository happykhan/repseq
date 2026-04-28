# dustmeselecta

Select representative bacterial isolates from a large short-read collection for follow-up long-read sequencing.

Given a folder of assemblies and a budget of N samples, `dms` splits slots between:
- **Phylogenetic diversity** — PARNAS medoid selection on a Mash distance tree
- **AMR/plasmid diversity** — greedy set cover on resistance gene and replicon profiles

The balance between the two is controlled by a single `--alpha` parameter.

## Install

Requires [pixi](https://pixi.sh).

```bash
git clone git@github.com:happykhan/dustmeselecta.git
cd dustmeselecta
pixi install
```

## Quick start

```bash
# Select 96 representatives from a folder of assemblies
pixi run dms select --assemblies assemblies/ --n 96

# With a pre-built tree and Kleborate output
pixi run dms select --assemblies assemblies/ --tree iqtree.nwk --kleborate kleborate.tsv --n 96

# Favour AMR/plasmid diversity (alpha=0.3 means 30% phylo, 70% AMR)
pixi run dms select --assemblies assemblies/ --n 96 --alpha 0.3 --output-dir results/
```

## Subcommands

### `dms select`

Main selection command.

```
Options:
  --assemblies DIR     Folder of .fasta/.fa/.fna assemblies  [required]
  --tree FILE          Pre-built Newick tree (skip mashtree)
  --kleborate FILE     Pre-run Kleborate TSV (skip auto-run)
  --n INT              Number of samples to select  [default: 10]
  --alpha FLOAT        Fraction of slots for phylo diversity, 0–1  [default: 0.5]
  --output-dir DIR     Output directory  [default: current dir]
```

**Algorithm:**
1. Build Mash distance tree with `mashtree` (skipped if `--tree` provided)
2. Run Kleborate on assemblies (skipped if `--kleborate` provided)
3. `n_phylo = round(alpha × n)` slots filled by PARNAS phylogenetic medoids
4. `n_amr = n - n_phylo` slots filled by greedy set cover on AMR gene + replicon profiles
5. Final selection = union of both sets

**Outputs (all in `--output-dir`):**

| File | Description |
|------|-------------|
| `selected.txt` | One sample ID per line — ready for sequencing request |
| `report.tsv` | All samples: selected, cluster, slot type (phylo/amr/not_selected) |
| `diversity_scores.csv` | PARNAS diversity covered at each n from 2 to n |
| `coverage_summary.txt` | Human-readable: % AMR gene combinations + replicon types covered |
| `elbow_plot.png` | Phylo diversity covered vs n — find the elbow |
| `scatter_plot.png` | PCoA of Mash distances, selected samples starred |
| `tree_heatmap.png` | Phylogenetic tree + AMR gene heatmap, selected highlighted in red |

### `dms evaluate`

Score a selection against a complete ground-truth dataset.

```bash
pixi run dms evaluate \
  --selected selected.txt \
  --ground-truth complete_kleborate.tsv \
  --tree complete.treefile \
  --output-dir eval/
```

Outputs `coverage_metrics.tsv` with:
- % AMR gene combinations covered
- % replicon types covered
- % phylogenetic diversity (Faith PD) covered

### `dms sweep`

Run `select` + `evaluate` across α = 0.0, 0.1, … 1.0 and generate a Pareto curve.

```bash
pixi run dms sweep \
  --assemblies assemblies/ \
  --n 96 \
  --ground-truth complete_kleborate.tsv \
  --output-dir sweep/
```

Outputs `pareto.tsv` and `pareto_plot.png` showing the trade-off between phylogenetic and AMR/plasmid coverage at each alpha.

## Test data

```bash
cd test_data && bash download.sh
pixi run dms select --assemblies test_data/ --n 5
```

## Choosing alpha

Run with the default (α = 0.5) first and check the `coverage_summary.txt`. If your research question is primarily about AMR gene context and plasmid architecture, lower alpha (0.2–0.3). If phylogenetic breadth and outbreak reconstruction matter more, raise it (0.7–0.8). Use `dms sweep` with a ground-truth dataset to find the empirically optimal value.

The elbow plot from `dms select` shows where adding more phylogenetic representatives hits diminishing returns — useful for choosing n.

## Roadmap

- [ ] PlasmidFinder integration for replicon typing (`--plasmid-finder`)
- [ ] hAMRonization support for AMRFinderPlus, ResFinder, RGI (`--hamronization`)
- [ ] Ancestral state reconstruction mode (PastML) for state-change-based selection
- [ ] bactscout-style logging

## Background

See: Arredondo-Alonso et al., *GigaScience* 2021 (pangenome-based selection) and PARNAS, *Systematic Biology* 2023 (exact k-medoids on phylogenetic trees).

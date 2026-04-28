# repseq

Select representative bacterial isolates from a large short-read collection for follow-up long-read sequencing.

Given a folder of assemblies and a budget of N samples, `repseq` splits slots between:
- **Phylogenetic diversity** — PARNAS exact k-medoids on a Mash distance tree
- **AMR/plasmid diversity** — greedy set cover on resistance gene and replicon profiles

The balance between the two is controlled by a single `--alpha` parameter.

## Install

Requires [pixi](https://pixi.sh).

```bash
git clone git@github.com:happykhan/repseq.git
cd repseq
pixi install
```

## Quick start

```bash
# Select 20 representatives — fully automatic (runs mashtree + ABRicate)
pixi run repseq select --assemblies assemblies/ --n 20

# With a pre-built tree and Kleborate output (Klebsiella collections)
pixi run repseq select --assemblies assemblies/ --tree tree.nwk --kleborate kleborate.tsv --n 20

# With hAMRonization output (any species — AMRFinder+, RGI, ResFinder etc.)
pixi run repseq select --assemblies assemblies/ --hamronization amrfinder.tsv --n 20

# Favour AMR/plasmid diversity (30% phylo slots, 70% AMR slots)
pixi run repseq select --assemblies assemblies/ --n 20 --alpha 0.3 --output-dir results/

# Add co-occurrence features to capture plasmid-gene linkage
pixi run repseq select --assemblies assemblies/ --n 20 --cooccurrence
```

## Subcommands

### `repseq select`

Main selection command.

```
Options:
  --assemblies DIR       Folder of .fasta/.fa/.fna assemblies              [required]

  AMR gene features (pick one, or omit to auto-run ABRicate):
  --hamronization FILE   hAMRonization TSV — any tool (AMRFinder+, RGI, ResFinder)
  --kleborate FILE       Kleborate TSV — Klebsiella only, also adds ST/virulence metadata

  Replicon/plasmid features (omit to auto-run ABRicate plasmidfinder db):
  --plasmid-finder FILE  Pre-run PlasmidFinder merged TSV

  Selection parameters:
  --tree FILE            Pre-built Newick tree (skips mashtree)
  --n INT                Number of samples to select             [default: 10]
  --alpha FLOAT          Fraction of slots for phylo diversity, 0–1  [default: 0.5]
                           alpha=1.0 → all slots filled by phylogeny only
                           alpha=0.0 → all slots filled by AMR/plasmid diversity only
  --cooccurrence         Add REP+AMR co-occurrence features (off by default)
  --output-dir DIR       Output directory                        [default: .]
```

**If no AMR flags are given**, `repseq` runs `abricate --db ncbi` on all assemblies automatically. This works on any bacterial species.

**If no replicon flag is given**, `repseq` runs `abricate --db plasmidfinder` automatically.

**AMR input priority:** `--hamronization` > `--kleborate` > ABRicate auto-run.

**Algorithm:**
1. Build Mash distance tree with `mashtree` (skipped if `--tree` provided)
2. Build binary AMR presence-absence matrix from best available source (hAMRonization → Kleborate → ABRicate ncbi)
3. Add replicon features from best available source (`--plasmid-finder` → ABRicate plasmidfinder)
4. Pad matrix so every assembly appears as a row, even if it had no AMR or replicon hits
5. `n_phylo = round(alpha × n)` slots filled by PARNAS phylogenetic medoids
6. `n_amr = n − n_phylo` slots filled by greedy set cover on the feature matrix
7. Final selection = union of both sets

**Outputs (all in `--output-dir`):**

| File | Description |
|------|-------------|
| `selected.txt` | One sample ID per line — ready for sequencing request |
| `report.tsv` | All samples: selected flag, slot type (phylo/amr/not_selected), AMR genes, replicons |
| `coverage_summary.txt` | Human-readable % AMR gene and replicon coverage |
| `abricate_ncbi.tsv` | Raw ABRicate AMR output (present only if auto-run) |
| `abricate_plasmidfinder.tsv` | Raw ABRicate replicon output (present only if auto-run) |
| `diversity_scores.csv` | PARNAS diversity covered at each n from 2 to n |
| `elbow_plot.png` | Phylo diversity vs n — use to choose n |
| `scatter_plot.png` | PCoA of Mash distances, selected samples starred |
| `tree_heatmap.png` | Tree + AMR gene/replicon heatmap, selected samples in red |

### `repseq evaluate`

Score a selection against a complete ground-truth dataset.

```bash
pixi run repseq evaluate \
  --selected selected.txt \
  --ground-truth complete_kleborate.tsv \
  --tree complete.nwk \
  --output-dir eval/
```

Outputs `coverage_metrics.tsv` with:

| Metric | Description |
|--------|-------------|
| `pct_faith_pd` | Faith's Phylogenetic Diversity covered (%) |
| `max_minimax_dist` | Max distance from any unselected sample to its nearest selected neighbour — the quantity PARNAS directly optimises |
| `mean_minimax_dist` | Mean of the above across all unselected samples |
| `pct_amr_covered` | % of AMR gene features in the collection covered by selection |
| `pct_replicons_covered` | % of replicon types covered |
| `random_mean_faith_pd_pct` | Mean Faith PD % of 100 random draws of the same N |
| `random_mean_amr_pct` | Mean AMR coverage % of 100 random draws of the same N |

The random baseline columns let you immediately quantify how much better repseq does than chance, which is the first thing a reviewer will ask for.

**Faith's Phylogenetic Diversity (PD):** total branch length of the tree spanned by the selected samples as a fraction of the total branch length of the full collection. A selection covering 90% Faith PD represents 90% of the evolutionary history in the dataset.

**Minimax distance:** the worst-case gap between any unselected sample and its nearest representative. A lower max minimax distance means more uniform coverage of the tree. This is what PARNAS directly minimises, so it is the most appropriate metric for evaluating the phylogenetic component of the selection.

### `repseq sweep`

Run `select` + `evaluate` across α = 0.0, 0.1, … 1.0 and generate a Pareto curve.

```bash
pixi run repseq sweep \
  --assemblies assemblies/ \
  --n 20 \
  --ground-truth complete_kleborate.tsv \
  --output-dir sweep/
```

Outputs `pareto.tsv`, `pareto_plot.png`, and `pareto_note.txt`.

**Important:** the sweep is a visualisation, not an optimisation. It does not tell you which α to use — that is a scientific decision based on your study goals. Read `pareto_note.txt` for guidance. As a rule, look for the elbow in `pareto_plot.png` — the point where adding more phylo budget stops improving Faith PD meaningfully. Use the `random_mean_faith_pd_pct` column to verify the chosen α beats random selection on both axes.

## Choosing alpha

| Goal | Suggested alpha |
|------|----------------|
| Outbreak reconstruction / phylogenetic breadth | 0.7 – 1.0 |
| Balanced (default) | 0.5 |
| AMR gene context and plasmid architecture | 0.2 – 0.4 |
| Pure AMR/replicon diversity | 0.0 |

Run `repseq sweep` with a ground-truth dataset to find the empirically best value for your collection. Check that it beats the random baseline on both axes before committing.

## Co-occurrence features (`--cooccurrence`)

By default, the feature matrix treats AMR gene presence and replicon presence independently. This means a sample with *IncFII* and *CTX-M-15* is indistinguishable from one where those two happen to be on different plasmids.

The `--cooccurrence` flag adds combined features (e.g. `CO:IncFII+CTX-M-15`) that are 1 only when both are present in the same sample. This captures plasmid–gene linkage and rewards selecting samples that carry novel gene-on-plasmid combinations, which is often the point of long-read sequencing. It is off by default because it can expand the feature matrix substantially on large diverse collections.

## Test data

```bash
cd test_data && bash download.sh
pixi run repseq select --assemblies test_data/ --n 5
```

## Scientific notes

**Is the greedy set cover optimal?** No — set cover is NP-hard. The greedy algorithm is guaranteed to find a solution within a factor of (1 − 1/e) ≈ 63% of the optimum (Nemhauser et al. 1978), and in practice performs much better on biological datasets where features cluster by lineage.

**Why report both Faith PD and minimax distance?** PARNAS minimises the maximum distance from any sample to its nearest representative (minimax criterion), but Faith PD — the standard metric in the literature — measures something different (total branch length covered). Both are reported so you can evaluate the selection on the criterion it was optimised for *and* the criterion reviewers expect to see.

**Validating the tree:** `mashtree` produces a neighbour-joining tree from Mash sketches, which is fast but approximate. For publication, verify that the topology agrees with a core-genome SNP tree using your preferred aligner/phylogenetic tool, and use the SNP tree via `--tree` if they disagree substantially.

## Roadmap

- [x] PlasmidFinder integration for replicon typing (`--plasmid-finder`)
- [x] hAMRonization support for AMRFinder+, ResFinder, RGI (`--hamronization`)
- [ ] Ancestral state reconstruction mode (PastML) for state-change-based selection
- [ ] Deduplication of near-identical assemblies before selection
- [ ] ST-aware selection constraint ("must cover all STs with ≥2 samples")

## Background

**PARNAS** — phylogenetic medoid selection:
Trost et al. "PARNAS: Objectively Selecting the Most Representative Taxa on a Phylogeny." *Systematic Biology* 72(5):1052–1063, 2023. https://doi.org/10.1093/sysbio/syad028

**ABRicate** — fast AMR and replicon screening against multiple databases:
https://github.com/tseemann/abricate

**Kleborate** — AMR and virulence typing for Klebsiella (optional, enriches report):
Lam et al. "A genomic surveillance framework and genotyping tool for Klebsiella pneumoniae and its related species complex." *Nature Communications* 12:4188, 2021. https://doi.org/10.1038/s41467-021-24448-3

**PlasmidFinder** — plasmid replicon typing:
Carattoli et al. "In Silico Detection and Typing of Plasmids using PlasmidFinder." *Antimicrobial Agents and Chemotherapy* 58(7):3895–3903, 2014. https://doi.org/10.1128/AAC.02412-14

**hAMRonization** — standardised AMR tool output parsing:
https://github.com/pha4ge/hAMRonization

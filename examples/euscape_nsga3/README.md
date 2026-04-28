# EuSCAPE NSGA-III Worked Example

This walkthrough demonstrates the NSGA-III multi-objective selection workflow on 79 *Klebsiella pneumoniae* isolates from the **EuSCAPE** (European Survey of Carbapenemase-Producing Enterobacteriaceae) collection. Unlike `repseq select`, which splits a fixed budget between phylogenetic and AMR objectives, NSGA-III optimises all three objectives simultaneously and returns a Pareto front of non-dominated solutions.

The goal: select 20 representative isolates for follow-up long-read sequencing, balancing phylogenetic breadth, AMR profile diversity, and plasmid replicon profile diversity without pre-specifying a trade-off parameter (alpha).

## Dataset

| Property | Value |
|----------|-------|
| Species | *Klebsiella pneumoniae* complex |
| Isolates | 79 assemblies |
| Countries | Austria, Belgium, Czech Republic, France, Germany, Greece, Hungary, Ireland, Israel, Italy, Latvia, Luxembourg, Romania, Spain, Turkey, UK |
| Dominant STs | ST512 (n=21), ST258 (n=9), ST11 (n=8), ST147 (n=4) |
| AMR gene features (Kleborate) | 87 unique acquired genes |
| Replicon types (ABRicate plasmidfinder) | 40 unique replicon types (among tree-intersected samples) |

The EuSCAPE collection is a well-characterised surveillance dataset with known epidemiological context. The collection is ST-structured: ST512 and ST258 are closely related NDM/KPC-producing lineages that dominate European carbapenem-resistant *Kpn* surveillance.

Note: 58 of the 79 assemblies have matching entries in the Kleborate output and tree. NSGA-III optimises over this intersection.

## Method

**NSGA-III** (Non-dominated Sorting Genetic Algorithm III) optimises three objectives simultaneously:

1. **Minimax phylogenetic distance** -- minimise the worst-case gap between any unselected sample and its nearest selected sample. This is the same criterion used by PARNAS: the minimax distance measures how well the selected set covers the full phylogeny, with lower values meaning no sample is far from a representative. Unlike mean pairwise distance, minimax directly penalises large gaps in coverage.
2. **AMR profile diversity** -- maximise the mean pairwise Jaccard distance between AMR gene profiles of selected samples. A value of 1.0 means all selected pairs have completely non-overlapping resistance profiles; 0.0 means all selected samples have identical profiles. This rewards selecting samples that are spread across AMR profile space rather than simply ticking every gene box (as set-cover coverage does).
3. **Replicon profile diversity** -- maximise the mean pairwise Jaccard distance between plasmid replicon profiles of selected samples. Interpretation is the same as for AMR: higher values mean more diverse replicon repertoires across the selected set.

Unlike `repseq select` (the split method), NSGA-III does not require an alpha parameter to balance phylogenetic vs. AMR objectives. Instead, it uses Das-Dennis reference directions (12 partitions = 91 reference points for 3 objectives) to maintain diversity across the Pareto front. The algorithm returns all non-dominated solutions, and the user chooses from the front.

**Solution representation:** Each candidate solution is a fixed-size set of k=20 sample indices. Custom crossover (pool-and-resample) and mutation (swap one selected sample for an unselected one) operators ensure valid subsets throughout the search. The population size is 100 and the algorithm runs for 300 generations.

**Recommended solution:** The solution closest to the ideal point in normalised objective space, where minimax distance is inverted (lower = better) and AMR/replicon diversities are normalised to [0, 1]. This is the solution that best balances all three objectives without strongly sacrificing any one.

## Input files

| File | Description |
|------|-------------|
| `../euscape/inputs/tree.nwk` | Mash NJ tree built with `mashtree` (79 tips) |
| `../euscape/inputs/kleborate.tsv` | Kleborate v3 output (AMR gene typing, ST, virulence) |
| `abricate_plasmidfinder.tsv` | ABRicate plasmidfinder output (replicon typing, 464 hits across 79 assemblies) |

The tree and Kleborate results were pre-computed and provided as inputs. The ABRicate plasmidfinder output was generated as the first step of this example (see below).

## Commands

```bash
# Step 1: Run ABRicate plasmidfinder for replicon typing
pixi run repseq abricate --assemblies assemblies/ --db plasmidfinder --output-dir .

# Step 2: Run NSGA-III
pixi run repseq nsga3 \
  --assemblies assemblies/ \
  --tree inputs/tree.nwk \
  --kleborate inputs/kleborate.tsv \
  --abricate-replicons abricate_plasmidfinder.tsv \
  --n 20 \
  --pop-size 100 \
  --generations 300 \
  --seed 42 \
  --output-dir .
```

## Pareto front

NSGA-III returned **31 non-dominated solutions**, each representing a different trade-off between minimax distance, AMR profile diversity, and replicon profile diversity. The recommended solution (closest to the ideal point) achieves **AMR diversity 0.9306**, **replicon diversity 0.8419**, and a minimax distance of **0.01037**.

The table below shows a representative subset of solutions spanning the range of trade-offs:

| Solution | Minimax dist | AMR diversity | Replicon diversity | Note |
|----------|-------------|---------------|-------------------|------|
| 4 | 0.00326 | 0.9285 | 0.8027 | Lowest minimax (best phylo coverage) |
| 16 | 0.01037 | 0.9379 | 0.7968 | Highest AMR diversity |
| 25 | 0.01037 | 0.9164 | 0.8602 | Highest replicon diversity |
| **22** | **0.01037** | **0.9306** | **0.8419** | **Recommended (closest to ideal)** |
| 17 | 0.01037 | 0.9273 | 0.8542 | Strong AMR, good replicons |
| 7 | 0.01030 | 0.9131 | 0.8599 | High replicon, moderate AMR |
| 29 | 0.04902 | 0.9366 | 0.7977 | Gap filler -- high AMR diversity, poor phylo |

The full front ranges from 0.9131--0.9379 AMR diversity, 0.7957--0.8602 replicon diversity, and 0.00326--0.04902 minimax distance. The AMR and replicon diversity metrics are mean pairwise Jaccard distances (0--1 scale), so these values indicate high overall diversity across solutions: even the worst solutions select samples with substantially different AMR profiles from each other.

Note: the AMR diversity metric here is not directly comparable to the AMR coverage percentage reported by `repseq select`. Coverage measures the fraction of unique genes represented at least once; Jaccard diversity measures how different the selected profiles are from each other. A solution with 0.93 Jaccard diversity selects samples whose AMR repertoires overlap on average only 7% of their genes.

![Pareto front scatter](pareto_front.png)

The 3-panel scatter shows pairwise relationships between objectives. The red star marks the recommended solution. The AMR vs. replicon diversity panel shows a clear trade-off: solutions with highest AMR diversity tend to have lower replicon diversity and vice versa.

![Parallel coordinates](pareto_parallel.png)

The parallel coordinates plot shows the same trade-offs. The recommended solution (red line) runs through the upper middle of all three axes, confirming it does not strongly sacrifice any single objective.

## Anchor solutions

Five solutions are selected as anchor points across the Pareto front, each representing an extreme or balanced trade-off. These are chosen deterministically: best on each individual objective, the recommended balanced solution, and a gap-filler that maximises diversity from the other four in normalised objective space.

### Recommended (balanced) -- solution 22

| Metric | Value |
|--------|-------|
| Minimax distance | 0.01037 |
| AMR profile diversity | 0.9306 |
| Replicon profile diversity | 0.8419 |

The recommended solution balances all three objectives. It achieves strong AMR diversity (0.9306) and replicon diversity (0.8419) while maintaining reasonable phylogenetic coverage. This is the best starting point when no single objective dominates the study design.

![Recommended solution heatmap](recommended_heatmap.png)

### Best AMR diversity -- solution 16

| Metric | Value |
|--------|-------|
| Minimax distance | 0.01037 |
| AMR profile diversity | 0.9379 |
| Replicon profile diversity | 0.7968 |

This solution maximises AMR profile diversity, selecting samples with the most dissimilar resistance repertoires. Replicon diversity is lower (0.7968), indicating more overlap in plasmid replicon profiles among the selected samples. Suitable when capturing diverse resistance mechanisms is the primary goal.

![Best AMR heatmap](solution_amr_heatmap.png)

### Best replicon diversity -- solution 25

| Metric | Value |
|--------|-------|
| Minimax distance | 0.01037 |
| AMR profile diversity | 0.9164 |
| Replicon profile diversity | 0.8602 |

This solution maximises replicon profile diversity, selecting samples with the most dissimilar plasmid replicon profiles. AMR diversity is slightly lower (0.9164). This is the best choice when plasmid diversity is the priority, such as for mobilome studies.

![Best replicon heatmap](solution_rep_heatmap.png)

### Best phylo -- solution 4

| Metric | Value |
|--------|-------|
| Minimax distance | 0.00326 |
| AMR profile diversity | 0.9285 |
| Replicon profile diversity | 0.8027 |

The lowest minimax distance means every unselected sample is close to at least one selected representative. AMR diversity is still high (0.9285), making this a strong choice when phylogenetic representativeness is the priority.

![Best phylo heatmap](solution_phylo_heatmap.png)

### Gap filler -- solution 29

| Metric | Value |
|--------|-------|
| Minimax distance | 0.04902 |
| AMR profile diversity | 0.9366 |
| Replicon profile diversity | 0.7977 |

Selected as the solution most distant from the other four anchors in normalised objective space. It has the highest minimax distance (worst phylogenetic coverage) but strong AMR diversity, filling a gap in the displayed trade-off landscape that the extremes and the recommended solution leave uncovered.

![Gap filler heatmap](solution_balanced_heatmap.png)

## Output files

| File | Description |
|------|-------------|
| `pareto_front.tsv` | All 31 Pareto-optimal solutions with objective values and sample lists |
| `recommended.txt` | Sample IDs from the recommended solution |
| `pareto_front.png` | 3-panel pairwise scatter of Pareto front |
| `pareto_parallel.png` | Parallel coordinates plot of Pareto front |
| `abricate_plasmidfinder.tsv` | ABRicate plasmidfinder output (replicon hits) |
| `recommended_heatmap.png` | Tree + AMR/replicon heatmap for recommended solution |
| `solution_amr_heatmap.png` | Tree + heatmap for best AMR diversity solution |
| `solution_rep_heatmap.png` | Tree + heatmap for best replicon diversity solution |
| `solution_phylo_heatmap.png` | Tree + heatmap for best phylo solution |
| `solution_balanced_heatmap.png` | Tree + heatmap for gap-filler solution |

## Notes

- **Replicon objective:** The replicon diversity objective is only meaningful when ABRicate plasmidfinder (or PlasmidFinder) output is provided. With Kleborate-only input and no replicon data, replicon diversity is trivially zero for all solutions and the problem reduces to a 2-objective optimisation. Always provide replicon data for a real 3-objective run.

- **Tree-matrix intersection:** NSGA-III operates on the intersection of tree tips and feature matrix rows. If some assemblies are missing from the Kleborate output (e.g., non-*Kpn* samples filtered out), the effective sample pool is smaller than the number of assemblies.

- **EuSCAPE assembly quality:** Assemblies are short-read Illumina assemblies. The Mash distance tree approximates the true SNP phylogeny; for publication, validate topology against a core-genome SNP tree.

- **Reproducibility:** Results are deterministic for a given seed (default 42). Changing population size or generation count will produce different Pareto fronts.

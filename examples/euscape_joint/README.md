# EuSCAPE joint selection (n=20, joint-weight=0.5)

## Method

The `--method joint` flag replaces the default alpha budget-split with a single
k-medoids optimisation on a blended distance matrix. Two distance matrices are
computed: phylogenetic (tip-to-tip branch lengths from the Newick tree) and AMR
(Jaccard distances on the binary AMR/replicon feature matrix). Both are
normalised to [0, 1] and combined as:

    D_combined = (1 - w) * D_phylo + w * D_amr

where `w` is `--joint-weight` (0.5 by default, equal weighting). K-medoids with
k-means++ initialisation and PAM-style refinement selects k representative
samples that are collectively as spread out as possible in this blended space.

## Command

```bash
repseq select \
  --assemblies /tmp/euscape_fastas/ \
  --tree /tmp/euscape_out/tree.nwk \
  --kleborate /tmp/euscape_out/kleborate.tsv \
  --n 20 \
  --method joint \
  --joint-weight 0.5 \
  --output-dir /tmp/euscape_example/euscape_joint
```

## Selected samples

| sample_id | slot_type | n_amr_genes | n_replicons |
|-----------|-----------|-------------|-------------|
| EuSCAPE_AT023 | joint | 14 | 9 |
| EuSCAPE_DE073 | joint | 5 | 9 |
| EuSCAPE_ES046 | joint | 13 | 4 |
| EuSCAPE_ES094 | joint | 12 | 6 |
| EuSCAPE_GR153 | joint | 1 | 5 |
| EuSCAPE_HU009 | joint | 7 | 6 |
| EuSCAPE_IE008 | joint | 6 | 2 |
| EuSCAPE_IE024 | joint | 8 | 5 |
| EuSCAPE_IL011 | joint | 4 | 3 |
| EuSCAPE_IL046 | joint | 6 | 6 |
| EuSCAPE_IT062 | joint | 9 | 4 |
| EuSCAPE_IT258 | joint | 10 | 6 |
| EuSCAPE_IT278 | joint | 9 | 6 |
| EuSCAPE_RO052 | joint | 0 | 7 |
| EuSCAPE_RS002 | joint | 0 | 7 |
| EuSCAPE_RS105 | joint | 0 | 4 |
| EuSCAPE_TR009 | joint | 0 | 6 |
| EuSCAPE_TR057 | joint | 0 | 6 |
| EuSCAPE_TR203 | joint | 0 | 5 |
| EuSCAPE_UK048 | joint | 0 | 6 |

## Coverage summary

| Metric | Joint (this run) | Split (alpha=0.5) |
|--------|------------------|-------------------|
| AMR genes covered | 45 / 87 (51.7%) | 78 / 87 (89.7%) |
| Replicon types covered | 40 / 48 (83.3%) | 44 / 48 (91.7%) |

## Comparison with split method

Of the 20 samples selected, **10 overlap** with the split (alpha=0.5) selection
and 10 differ.

**Shared (10):** EuSCAPE_DE073, EuSCAPE_ES046, EuSCAPE_ES094, EuSCAPE_HU009,
EuSCAPE_IE008, EuSCAPE_IL046, EuSCAPE_IT062, EuSCAPE_IT278, EuSCAPE_TR009,
EuSCAPE_TR057

**Only in split:** EuSCAPE_BE092, EuSCAPE_CZ007, EuSCAPE_ES085, EuSCAPE_FR056,
EuSCAPE_GR075, EuSCAPE_IT266, EuSCAPE_LU008, EuSCAPE_LV006, EuSCAPE_RO094,
EuSCAPE_TR083

**Only in joint:** EuSCAPE_AT023, EuSCAPE_GR153, EuSCAPE_IE024, EuSCAPE_IL011,
EuSCAPE_IT258, EuSCAPE_RO052, EuSCAPE_RS002, EuSCAPE_RS105, EuSCAPE_TR203,
EuSCAPE_UK048

The split method explicitly maximises each axis separately: PARNAS picks
phylogenetically diverse tips, then greedy set cover fills AMR/replicon gaps.
This produces high coverage on both axes (89.7% AMR, 91.7% replicons) but the
two selection criteria operate independently.

The joint method selects "generalist" samples that score moderately well on both
axes simultaneously. It finds samples that are spread out in the blended
phylo+AMR space, which tends to favour structurally distinct isolates (different
tree positions *and* different resistance profiles). The trade-off is lower raw
AMR feature coverage (51.7%) because the method does not explicitly chase rare
AMR genes the way greedy set cover does. Joint selection is better suited to
scenarios where you want a representative panel that balances both dimensions
rather than maximising coverage of a specific feature type.

# repseq

Select representative bacterial isolates for long-read sequencing from surveillance collections.

`repseq` has two main workflows:

1. **`repseq ghru`** -- Hierarchical priority ranking for sentinel site surveillance (e.g. the ARSRL Kleb survey). This is the primary workflow documented below.
2. **`repseq select`** -- Phylogenetic + AMR diversity selection from assembly collections (research use).

This README covers the `ghru` workflow. It will be reviewed by David and Aruka.

---

## Install

Requires [pixi](https://pixi.sh).

```bash
git clone git@github.com:happykhan/repseq.git
cd repseq
pixi install
```

Verify the install:

```bash
pixi run repseq ghru --help
```

---

## Preparing your input CSV

The simplest way to use `repseq ghru` is with a single flat CSV file. Each row is one isolate.

### Required columns

| Column | Description | Example |
|--------|-------------|---------|
| `isolate_id` | Unique identifier for the isolate (accession number) | `KPN-2024-0042` |
| `site` | Sentinel site name | `Manila` |
| `tier` | Resistance tier (see table below) | `Tier1_CR` |
| `replicons` | Comma-separated PlasmidFinder replicon list. Leave blank if none | `IncFIB(K), IncFII(K)` |

### Tier values

| Tier | Meaning | How to assign |
|------|---------|---------------|
| `Tier1_CR` | Carbapenem-resistant (R or I to any carbapenem) | ETP, IPM, or MEM non-susceptible |
| `Tier2_MDR` | Multi-drug resistant | Resistant to 3+ antimicrobial classes (excluding AMP/PEN), or colistin-resistant |
| `Tier3_nonMDR_resistant` | Resistant but not MDR | Resistant to 1-2 non-PEN classes |
| `Tier4_pansusceptible` | Pan-susceptible | No resistance detected (AMP-only counts as pan-susceptible) |

If you do not have pre-assigned tiers, you can omit the `tier` column and instead provide:
- `resist_pattern` -- semicolon-separated resistant drug abbreviations (e.g. `AMP;CIP;GEN;SXT`)
- `drug_classes` -- number of resistant drug classes
- `carb_nonsus` -- `True` or `False` for carbapenem non-susceptibility

repseq will derive the tier automatically from these columns.

### Optional columns

| Column | Description | Default if absent |
|--------|-------------|-------------------|
| `spec_date` | Specimen collection date (YYYY-MM-DD). Used as a tiebreaker when two isolates have the same profile at the same site | All dates set to 2000-01-01 (no preference) |
| `rp_code` | Pre-computed resistance profile code (e.g. `RP1`). Both `rp_code` and `pp_code` must be present to skip derivation | Derived from `resist_pattern` frequency |
| `pp_code` | Pre-computed plasmid profile code (e.g. `PP1`) | Derived from `replicons` frequency |
| `resist_pattern` | Semicolon-separated resistant drug abbreviations. Used for RP code derivation and tier derivation | `pan-susceptible` |

### Example CSV

```
isolate_id,site,tier,replicons,spec_date,resist_pattern
KPN-001,Manila,Tier1_CR,"IncFIB(K), IncFII(K)",2024-06-15,AMP;ETP;CIP;GEN
KPN-002,Manila,Tier2_MDR,IncX4,2024-05-20,AMP;CIP;GEN;SXT
KPN-003,Cebu,Tier1_CR,"IncFIB(K), IncFII(K)",2024-07-01,AMP;IPM;MEM
KPN-004,Cebu,Tier4_pansusceptible,,2024-04-10,pan-susceptible
KPN-005,Davao,Tier4_pansusceptible,,2024-03-22,pan-susceptible
KPN-006,Manila,Tier1_CR,IncFIB(K),2024-08-01,AMP;ETP;CIP;GEN
KPN-007,Davao,Tier3_nonMDR_resistant,,2024-06-30,AMP;CIP
KPN-008,Cebu,Tier2_MDR,"IncFIB(K), IncX4",2024-05-15,AMP;CIP;GEN;SXT
```

Save this as `input.csv`.

---

## Running the prioritisation

### Step 1: Run the command

```bash
pixi run repseq ghru --csv input.csv --output-dir results/
```

This takes a few seconds. You will see a summary printed to the terminal.

### Step 2: Review the outputs

The `results/` folder will contain:

| File | What it contains |
|------|------------------|
| `priority_full.tsv` | All isolates ranked from 1 to N, with tier, batch, RP/PP codes, and priority rank |
| `primary_batch.tsv` | Only the isolates selected for sequencing (Site_Guarantee + Primary batch) |
| `selected.txt` | One isolate ID per line -- hand this to the sequencing team |
| `profile_summary.tsv` | Summary per resistance profile: how many isolates, how many sites, how many selected |
| `site_tier_summary.tsv` | Cross-tabulation of sites vs tiers |
| `plasmid_profiles.tsv` | Summary of plasmid replicon combinations |

### Step 3: Check the selected list

Open `results/selected.txt`. These are the isolates to sequence. The file contains one accession number per line.

Open `results/primary_batch.tsv` for the full details of each selected isolate, including why it was selected (the `selection_batch` column).

---

## How the selection works

### Site guarantee (hard floor)

Every sentinel site contributes at least one isolate, regardless of resistance tier. This is a hard floor: even if a site has only Tier 4 (pan-susceptible) isolates, one will be selected.

The rationale (from the 22 May TGIBF meeting): absence of resistance at a sentinel site is surveillance data. You cannot distinguish "no resistance detected" from "not sampled" if the site is missing entirely.

For each site, repseq picks the highest-ranked isolate (best tier, then most frequent profile combination, then most recent specimen date). These isolates are tagged `selection_batch = "Site_Guarantee"` in the output.

To disable the site guarantee:

```bash
pixi run repseq ghru --csv input.csv --no-guarantee-sites --output-dir results/
```

### Primary batch (1 per site per profile combination)

After the site guarantees, repseq selects one isolate per site per RP x PP combination. RP is the resistance profile code (isolates with the same resistance drug pattern get the same RP code). PP is the plasmid profile code (isolates with the same replicon combination get the same PP code).

Within each site and RP x PP combination, the most recent isolate (by specimen date) is selected. These are tagged `selection_batch = "Primary"`.

### Secondary batch

All remaining isolates are tagged `selection_batch = "Secondary"`. They appear in `priority_full.tsv` but not in `primary_batch.tsv` or `selected.txt`.

### Priority ranking

All isolates (Site_Guarantee, Primary, and Secondary) are assigned a priority rank from 1 (highest priority) to N (lowest). The ranking follows this order:

1. **Site guarantees first** -- sorted by tier within the guarantee block
2. **Then tier** -- Tier1_CR before Tier2_MDR before Tier3 before Tier4
3. **Then batch** -- Primary before Secondary
4. **Then combo frequency** -- the most widespread RP x PP combination (most isolates across the collection) ranks first
5. **Then site** -- alphabetical within the same combo
6. **Then date** -- most recent specimen date first
7. **Then accession** -- alphabetical for complete determinism

---

## Worked example

Using the 8-isolate CSV from above:

```bash
pixi run repseq ghru --csv input.csv --output-dir example_results/
```

Expected output in `priority_full.tsv` (simplified):

| priority_rank | isolate_id | site | tier | selection_batch | rp_code | pp_code |
|---|---|---|---|---|---|---|
| 1 | KPN-006 | Manila | Tier1_CR | Site_Guarantee | RP1 | PP2 |
| 2 | KPN-003 | Cebu | Tier1_CR | Site_Guarantee | RP2 | PP1 |
| 3 | KPN-007 | Davao | Tier3_nonMDR_resistant | Site_Guarantee | ... | ... |
| 4 | KPN-001 | Manila | Tier1_CR | Primary | RP1 | PP1 |
| ... | ... | ... | ... | ... | ... | ... |

What happened:
- **3 site guarantees** (Manila, Cebu, Davao -- one each). Manila picked KPN-006 (Tier1, most recent). Cebu picked KPN-003 (Tier1). Davao picked KPN-007 (Tier3 -- the best available at that site).
- **Primary batch** fills in the remaining unique site x RP x PP slots not already covered by guarantees.
- **Secondary batch** contains duplicates (same site, same profile combination).

The `selected.txt` file contains the Site_Guarantee and Primary isolate IDs only. The Secondary isolates are backup candidates if a Primary isolate cannot be sequenced.

---

## Legacy metadata format

If your data is already in the column-mapped TSV/Excel format used by RITM, you can use the `--metadata` flag instead:

```bash
pixi run repseq ghru \
  --metadata isolate_metadata.tsv \
  --tier-col tier \
  --resist-pattern-col resist_pattern \
  --plasmids-col Plasmids \
  --output-dir results/
```

The `--metadata` flag supports custom column names via `--id-col`, `--lab-col`, `--date-col`, etc. See `pixi run repseq ghru --help` for all options.

---

## Other subcommands

### `repseq select`

Select N representative isolates from an assembly collection using phylogenetic + AMR diversity. This is a different algorithm from `ghru` -- it uses PARNAS k-medoids on a Mash distance tree combined with greedy set cover on AMR/replicon profiles.

```bash
pixi run repseq select --assemblies assemblies/ --n 20
```

#### Balancing AMR vs plasmid diversity (`--rep-weight`)

By default, AMR gene features and replicon (inc type) features are weighted equally in the greedy set cover. If your collection has many more AMR gene features than inc types — which is typical — the algorithm will naturally fill its budget with AMR-rich isolates and may never pick up rare inc types.

Use `--rep-weight` to shift the balance:

```bash
# Equal weight (default)
pixi run repseq select --assemblies assemblies/ --n 20 --rep-weight 1.0

# Favour plasmid diversity: each new inc type counts as 5 AMR genes
pixi run repseq select --assemblies assemblies/ --n 20 --rep-weight 5.0

# Pure AMR coverage: ignore plasmid features entirely
pixi run repseq select --assemblies assemblies/ --n 20 --rep-weight 0.0
```

After each run, check `coverage_summary.txt` in the output directory to see AMR and replicon coverage percentages. If REP coverage is low, increase `--rep-weight`.

See [`examples/rep_weight/compare_weights.py`](examples/rep_weight/compare_weights.py) for a worked example with simulated data showing the coverage trade-off across five weight values.

### `repseq evaluate`

Score a selection against a complete ground-truth dataset.

```bash
pixi run repseq evaluate \
  --selected selected.txt \
  --ground-truth complete_kleborate.tsv \
  --tree complete.nwk
```

### `repseq sweep`

Run `select` + `evaluate` across alpha values from 0 to 1 and generate a Pareto curve.

```bash
pixi run repseq sweep \
  --assemblies assemblies/ \
  --n 20 \
  --ground-truth complete_kleborate.tsv
```

### `repseq nsga3`

Multi-objective selection using NSGA-III (phylogenetic distance + AMR coverage + replicon coverage).

```bash
pixi run repseq nsga3 --assemblies assemblies/ --n 20
```

### `repseq diversity-curve`

Plot diversity saturation curves to help choose the number of representatives.

```bash
pixi run repseq diversity-curve --assemblies assemblies/
```

---

## Drug abbreviations

The following drug abbreviations are recognised for tier assignment and resistance profile grouping:

| Abbreviation | Drug | Class |
|---|---|---|
| AMP | Ampicillin | PEN (excluded from MDR count) |
| AMC | Amoxicillin-clavulanate | BLI |
| CZO | Cefazolin | 1GC |
| CXA | Cefuroxime | 2GC |
| FOX | Cefoxitin | 2GC |
| CTT | Cefotetan | 2GC |
| CRO | Ceftriaxone | 3GC |
| CAZ | Ceftazidime | 3GC |
| CTX | Cefotaxime | 3GC |
| FEP | Cefepime | 4GC |
| ETP | Ertapenem | CAR |
| IPM | Imipenem | CAR |
| MEM | Meropenem | CAR |
| ATM | Aztreonam | MON |
| GEN | Gentamicin | AMG |
| TOB | Tobramycin | AMG |
| AMK | Amikacin | AMG |
| CIP | Ciprofloxacin | FQN |
| LEV | Levofloxacin | FQN |
| SXT | Trimethoprim-sulfamethoxazole | FOL |
| TCY | Tetracycline | TET |
| TGC | Tigecycline | TET |
| COL | Colistin | COL |
| TZP | Piperacillin-tazobactam | BLI |
| SAM | Sulbactam-ampicillin | BLI |

Carbapenem non-susceptibility (R or I to ETP, IPM, or MEM) triggers Tier 1. Colistin resistance triggers Tier 2. AMP resistance alone does not contribute to MDR counting (intrinsic in KPN).

---

## References

- Gayeta J, RITM Philippines, 2026. Hierarchical sampling prioritisation for long-read sequencing.
- Magiorakos et al. "Multidrug-resistant, extensively drug-resistant and pandrug-resistant bacteria." *Clin Microbiol Infect* 18(3):268-281, 2012.
- Trost et al. "PARNAS: Objectively Selecting the Most Representative Taxa on a Phylogeny." *Systematic Biology* 72(5):1052-1063, 2023.
- Carattoli et al. "In Silico Detection and Typing of Plasmids using PlasmidFinder." *Antimicrobial Agents and Chemotherapy* 58(7):3895-3903, 2014.

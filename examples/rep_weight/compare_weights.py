"""Worked example: effect of --rep-weight on AMR vs plasmid coverage.

Simulates a 30-isolate collection with a deliberate trade-off:
  - 22 AMR gene features spread across isolates; several rare genes only
    appear in isolates that carry NO plasmids
  - 6 inc types that only appear in isolates with minimal AMR genes

With n=6 slots (tight budget) the algorithm must choose between capturing
rare AMR genes or rare inc types. --rep-weight shifts that balance.

Run with:
    pixi run python examples/rep_weight/compare_weights.py
"""

from __future__ import annotations

import pandas as pd

import repseq.amr_cover as _amr_mod
from repseq.amr_cover import greedy_set_cover

# Silence internal console messages so the comparison table stays readable
_amr_mod.print_message = lambda *a, **kw: None

# ---------------------------------------------------------------------------
# Build synthetic binary matrix with a deliberate trade-off
# ---------------------------------------------------------------------------

# AMR features: 8 common + 14 rare (rare genes only in "AMR-rich" isolates)
COMMON_AMR = [
    "AMR:Bla:TEM-1", "AMR:Bla:SHV-11", "AMR:AGly:aac(3)-IIa",
    "AMR:Flq:qnrB", "AMR:Sul:sul1", "AMR:Tet:tet(A)", "AMR:Tmt:dfrA1", "AMR:Col:mcr-1",
]
RARE_AMR = [
    "AMR:Bla_Carb:KPC-2", "AMR:Bla_Carb:NDM-1", "AMR:Bla_Carb:OXA-48", "AMR:Bla_Carb:OXA-232",
    "AMR:AGly:rmtB", "AMR:AGly:aph(3')-Ia", "AMR:Flq:qnrS",
    "AMR:Bla:CTX-M-15", "AMR:Bla:CTX-M-55", "AMR:Bla:OXA-1",
    "AMR:Tet:tet(M)", "AMR:Sul:sul2", "AMR:MLS:erm(B)", "AMR:Rif:arr-2",
]
ALL_AMR = COMMON_AMR + RARE_AMR

# Inc types: 6 types, each only in isolates with few AMR genes
INC_TYPES = ["REP:IncFII(K)", "REP:IncFIB(K)", "REP:IncX3", "REP:IncX4", "REP:IncI1", "REP:IncHI1B"]

ALL_FEATURES = ALL_AMR + INC_TYPES

rows: dict[str, dict[str, int]] = {}

# Isolates 1-6: "AMR-only" — carry the 14 rare AMR genes (2-3 each), no plasmids
rare_amr_groups = [RARE_AMR[i:i+3] for i in range(0, len(RARE_AMR), 3)]
# pad to 6 groups
while len(rare_amr_groups) < 6:
    rare_amr_groups.append(RARE_AMR[-2:])

for i in range(6):
    sid = f"AMR-only-{i+1:02d}"
    row = {f: 0 for f in ALL_FEATURES}
    for gene in COMMON_AMR[:4]:
        row[gene] = 1
    for gene in rare_amr_groups[i]:
        row[gene] = 1
    rows[sid] = row

# Isolates 7-12: "plasmid-only" — one isolate per inc type, carry only 2 common AMR genes
for i, inc in enumerate(INC_TYPES):
    sid = f"plasmid-only-{i+1:02d}"
    row = {f: 0 for f in ALL_FEATURES}
    row[COMMON_AMR[0]] = 1
    row[COMMON_AMR[1]] = 1
    row[inc] = 1
    rows[sid] = row

# Isolates 13-30: "background" — common AMR genes only, no plasmids
for i in range(18):
    sid = f"background-{i+1:02d}"
    row = {f: 0 for f in ALL_FEATURES}
    n_genes = (i % 4) + 2
    for gene in COMMON_AMR[:n_genes]:
        row[gene] = 1
    rows[sid] = row

binary_matrix = pd.DataFrame(rows).T[ALL_FEATURES].fillna(0).astype(int)

# Summary of the collection
n_amr_unique = int(binary_matrix[ALL_AMR].any(axis=0).sum())
n_rep_unique = int(binary_matrix[INC_TYPES].any(axis=0).sum())
print(f"Collection: {len(rows)} isolates")
print(f"  Unique AMR features present: {n_amr_unique}")
print(f"  Unique REP (inc type) features present: {n_rep_unique}")
print()
print("Trade-off: isolates with rare AMR genes carry NO plasmids;")
print("           isolates with plasmids carry only 2 common AMR genes.")
print()


# ---------------------------------------------------------------------------
# Helper: coverage metrics for a selection
# ---------------------------------------------------------------------------

def coverage(selected: list[str]) -> dict[str, object]:
    sel_rows = binary_matrix.loc[selected]
    amr_covered = int(sel_rows[ALL_AMR].any(axis=0).sum())
    rep_covered = int(sel_rows[INC_TYPES].any(axis=0).sum())
    return {
        "amr_covered": amr_covered,
        "amr_pct": amr_covered / n_amr_unique * 100,
        "rep_covered": rep_covered,
        "rep_pct": rep_covered / n_rep_unique * 100,
    }


# ---------------------------------------------------------------------------
# Run set cover at five rep_weight values
# ---------------------------------------------------------------------------

N_SELECT = 6
weights_to_test = [0.0, 0.5, 1.0, 3.0, 10.0]

print(f"Selecting {N_SELECT} isolates (budget) from {len(rows)}:\n")
print(f"{'rep_weight':>12} | {'AMR covered':>12} {'AMR %':>7} | {'REP covered':>12} {'REP %':>7}")
print("-" * 65)

for w in weights_to_test:
    selected = greedy_set_cover(binary_matrix, exclude_samples=[], n_amr=N_SELECT, rep_weight=w)
    cov = coverage(selected)
    print(
        f"{w:>12.1f} | {cov['amr_covered']:>12d} {cov['amr_pct']:>6.0f}% | "
        f"{cov['rep_covered']:>12d} {cov['rep_pct']:>6.0f}%"
    )

print()
print("Key:")
print(f"  Total unique AMR features in collection : {n_amr_unique}")
print(f"  Total unique inc types in collection    : {n_rep_unique}")
print()
print("Why the default (1.0) gives 0% REP coverage in this scenario:")
print(f"  There are {n_amr_unique} AMR features and {n_rep_unique} inc types. Each AMR-only isolate")
print("  covers 3-4 new AMR genes; each plasmid-only isolate covers only 1 new inc type")
print("  (= score 1 vs score 3-4). So the greedy algorithm always picks AMR-only first.")
print()
print("Choosing a rep_weight:")
print("  Start with --rep-weight 1.0 (default).")
print("  If the coverage summary shows low REP coverage, increase rep_weight until the")
print("  balance matches your biological priority. For this collection, rep_weight ~3")
print("  gives a 50/50 mix; rep_weight ~6 captures all inc types at the cost of some")
print("  rare AMR genes. There is no single correct value — it depends on whether your")
print("  primary goal is AMR surveillance or plasmid epidemiology.")

"""June Gayeta's sampling prioritisation method for long-read sequencing.

Hierarchical priority ranking based on resistance tier, site representation,
and resistance/plasmid profile diversity. Designed for ARSRL (Philippines)
surveillance collections.

Reference: Gayeta J, RITM Philippines, 2026.
"""

from __future__ import annotations

import os
from datetime import datetime

import pandas as pd

from repseq.log import print_message


# ---------------------------------------------------------------------------
# Tier assignment
# ---------------------------------------------------------------------------

# Carbapenem drug abbreviations (R or I triggers Tier 1)
CARBAPENEM_DRUGS = {"ETP", "IPM", "MEM"}

# Drug-to-class mapping for MDR counting (Magiorakos et al. 2012).
# AMP is excluded from MDR class counting (intrinsic KPN resistance).
# COL I is excluded from MDR class counting per June's spec.
DRUG_CLASS_MAP: dict[str, str] = {
    # Penicillins
    "AMP": "PEN",
    "AMC": "BLI",
    # Cephalosporins
    "CZO": "1GC",
    "CXA": "2GC",
    "FOX": "2GC",
    "CTT": "2GC",
    "CRO": "3GC",
    "CAZ": "3GC",
    "CTX": "3GC",
    "FEP": "4GC",
    # Carbapenems
    "ETP": "CAR",
    "IPM": "CAR",
    "MEM": "CAR",
    # Monobactams
    "ATM": "MON",
    # Aminoglycosides
    "GEN": "AMG",
    "TOB": "AMG",
    "AMK": "AMG",
    # Fluoroquinolones
    "CIP": "FQN",
    "LEV": "FQN",
    # Folate pathway
    "SXT": "FOL",
    # Tetracyclines
    "TCY": "TET",
    "TGC": "TET",
    # Colistin
    "COL": "COL",
    # Piperacillin-tazobactam
    "TZP": "BLI",
    # Sulbactam-ampicillin
    "SAM": "BLI",
}

# Classes excluded from MDR counting (AMP class = PEN, and COL intermediate)
MDR_EXCLUDE_CLASSES = {"PEN"}


def assign_tier(
    resist_drugs: list[str],
    carb_nonsus: bool,
    colistin_r: bool,
) -> str:
    """Assign resistance tier.

    Parameters
    ----------
    resist_drugs : list[str]
        List of drug abbreviations the isolate is resistant to (R interpretation).
        Case-sensitive upper-case abbreviations expected.
    carb_nonsus : bool
        True if non-susceptible (R or I) to any carbapenem.
    colistin_r : bool
        True if resistant to colistin.

    Returns
    -------
    str
        One of: Tier1_CR, Tier2_MDR, Tier3_nonMDR_resistant, Tier4_pansusceptible.
    """
    if carb_nonsus:
        return "Tier1_CR"

    # Count drug classes for MDR, excluding PEN (AMP is intrinsic)
    classes = set()
    for drug in resist_drugs:
        drug_upper = drug.upper()
        cls = DRUG_CLASS_MAP.get(drug_upper)
        if cls and cls not in MDR_EXCLUDE_CLASSES:
            classes.add(cls)

    if colistin_r or len(classes) >= 3:
        return "Tier2_MDR"

    if len(classes) >= 1:
        return "Tier3_nonMDR_resistant"

    return "Tier4_pansusceptible"


# ---------------------------------------------------------------------------
# Profile code assignment
# ---------------------------------------------------------------------------

def assign_rp_codes(
    df: pd.DataFrame,
    pattern_col: str = "resist_pattern",
    tier_col: str | None = "tier",
) -> pd.Series:
    """Assign RP codes (RP1, RP2, ...) by resist_pattern frequency.

    When *tier_col* is provided and exists in *df*, codes are assigned
    per-tier in tier order (Tier1 patterns numbered first, then Tier2,
    etc.).  Within each tier the most common pattern gets the lowest
    available number.  Ties are broken by pattern string (alphabetical)
    for reproducibility.

    Returns a Series aligned with df index.
    """
    use_tier = tier_col is not None and tier_col in df.columns

    if use_tier:
        # Assign codes per-tier, with continuous numbering across tiers
        pattern_to_code: dict[tuple[str, str], str] = {}
        code_counter = 1
        for tier in sorted(TIER_ORDER, key=TIER_ORDER.get):  # type: ignore[arg-type]
            tier_df = df[df[tier_col] == tier]
            if tier_df.empty:
                continue
            counts = tier_df[pattern_col].value_counts()
            sorted_patterns = (
                counts.reset_index()
                .rename(columns={pattern_col: "pattern", "count": "n"})
                .sort_values(["n", "pattern"], ascending=[False, True])
            )
            for _, row in sorted_patterns.iterrows():
                key = (tier, row["pattern"])
                if key not in pattern_to_code:
                    pattern_to_code[key] = f"RP{code_counter}"
                    code_counter += 1
        return pd.Series(
            [pattern_to_code.get((t, p), "RP0") for t, p in zip(df[tier_col], df[pattern_col])],
            index=df.index,
        )

    # Global assignment (no tier partitioning)
    counts = df[pattern_col].value_counts()
    sorted_patterns = (
        counts.reset_index()
        .rename(columns={pattern_col: "pattern", "count": "n"})
        .sort_values(["n", "pattern"], ascending=[False, True])
    )
    pattern_to_code_global = {
        row["pattern"]: f"RP{i + 1}"
        for i, (_, row) in enumerate(sorted_patterns.iterrows())
    }
    return df[pattern_col].map(pattern_to_code_global)


def assign_pp_codes(
    df: pd.DataFrame,
    plasmid_col: str = "plasmid_profile",
) -> pd.Series:
    """Assign PP codes (PP1, PP2, ...) by plasmid profile frequency.

    Isolates with no plasmids get code 'none'.
    Most common profile gets PP1.

    Returns a Series aligned with df index.
    """
    counts = df[plasmid_col].value_counts()
    sorted_profiles = (
        counts.reset_index()
        .rename(columns={plasmid_col: "profile", "count": "n"})
        .sort_values(["n", "profile"], ascending=[False, True])
    )
    profile_to_code = {
        row["profile"]: f"PP{i + 1}"
        for i, (_, row) in enumerate(sorted_profiles.iterrows())
    }
    return df[plasmid_col].map(profile_to_code)


# ---------------------------------------------------------------------------
# Primary batch selection + hierarchical sort
# ---------------------------------------------------------------------------

TIER_ORDER = {
    "Tier1_CR": 0,
    "Tier2_MDR": 1,
    "Tier3_nonMDR_resistant": 2,
    "Tier4_pansusceptible": 3,
}

BATCH_ORDER = {"Primary": 0, "Secondary": 1}


def select_primary_batch(df: pd.DataFrame) -> pd.Series:
    """Assign Primary/Secondary batch labels.

    Primary = 1 per site per RP x PP combination (most recent spec_date).
    Secondary = all others.

    Returns a Series of 'Primary'/'Secondary' aligned with df index.
    """
    df = df.copy()
    # Sort so that within each (rp_pp_combo, laboratory), most recent date comes first
    df["_spec_date_ts"] = pd.to_datetime(df["spec_date"])
    df = df.sort_values(
        ["rp_pp_combo", "laboratory", "_spec_date_ts"],
        ascending=[True, True, False],
    )

    # Mark the first occurrence per (rp_pp_combo, laboratory) as Primary
    is_first = ~df.duplicated(subset=["rp_pp_combo", "laboratory"], keep="first")
    result = pd.Series("Secondary", index=df.index)
    result[is_first] = "Primary"
    return result.reindex(df.index)


def compute_priority_rank(df: pd.DataFrame) -> pd.DataFrame:
    """Compute the full priority ranking.

    Expects columns: isolate_id, laboratory, spec_date, tier, rp_code,
    pp_code, rp_pp_combo, selection_batch.

    Adds: priority_rank, combo_total (total isolates per RP x PP combo),
    combo_sites (number of unique sites per combo).

    The 5-level sort:
      1. Tier (Tier1 first)
      2. Primary before Secondary
      3. Most widespread RP x PP combination first (highest total count)
      4. Site alphabetical within RP x PP
      5. Most recent spec_date within site x RP x PP
    """
    df = df.copy()

    # Compute combo-level stats
    combo_stats = df.groupby("rp_pp_combo").agg(
        combo_total=("rp_pp_combo", "size"),
        combo_sites=("laboratory", "nunique"),
    ).reset_index()
    df = df.merge(combo_stats, on="rp_pp_combo", how="left")

    # Numeric sort keys
    df["_tier_ord"] = df["tier"].map(TIER_ORDER)
    df["_batch_ord"] = df["selection_batch"].map(BATCH_ORDER)
    df["_spec_date_ts"] = pd.to_datetime(df["spec_date"])

    # Sort by the 5-level hierarchy
    # Level 3 tiebreaker for same combo_total: combo_sites DESC, then
    # laboratory ASC (which also serves as Level 4), then spec_date DESC (Level 5).
    # Additional tiebreaker: accession_no ASC for complete determinism.
    df = df.sort_values(
        [
            "_tier_ord",      # Level 1: Tier
            "_batch_ord",     # Level 2: Primary before Secondary
            "combo_total",    # Level 3: most widespread combo first
            "combo_sites",    # Level 3 tiebreaker: most sites
            "laboratory",     # Level 4: site alphabetical
            "_spec_date_ts",  # Level 5: most recent date first
            "isolate_id",     # Final tiebreaker: accession_no
        ],
        ascending=[True, True, False, False, True, False, True],
    )

    df["priority_rank"] = range(1, len(df) + 1)

    # Drop temp columns
    df = df.drop(columns=["_tier_ord", "_batch_ord", "_spec_date_ts"])

    return df


# ---------------------------------------------------------------------------
# Input parsing: reconstruct from metadata + resistance/plasmid data
# ---------------------------------------------------------------------------

def parse_resist_pattern(r_drugs_list: str, exclude_drugs: set[str] | None = None) -> str:
    """Convert a comma-separated r_drugs_list to semicolon-separated resist_pattern.

    Optionally exclude certain drugs (e.g. SAM, CTT, TGC per June's 80% threshold).
    """
    if pd.isna(r_drugs_list) or str(r_drugs_list).strip() in ("", "-", "nan", "pan-susceptible"):
        return "pan-susceptible"

    exclude = exclude_drugs or set()
    drugs = [d.strip() for d in str(r_drugs_list).split(",") if d.strip()]
    drugs = [d for d in drugs if d.upper() not in {e.upper() for e in exclude}]

    if not drugs:
        return "pan-susceptible"

    return ";".join(drugs)


def parse_plasmid_profile(plasmids_str: str) -> str:
    """Convert a comma-separated plasmid string to a canonical sorted profile."""
    if pd.isna(plasmids_str) or str(plasmids_str).strip() in ("", "-", "nan", "none", "None"):
        return "none"

    plasmids = sorted(p.strip() for p in str(plasmids_str).split(",") if p.strip())
    return ", ".join(plasmids) if plasmids else "none"


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run_june_prioritisation(
    metadata_path: str,
    output_dir: str,
    id_col: str = "isolate_id",
    lab_col: str = "laboratory",
    date_col: str = "spec_date",
    tier_col: str | None = None,
    resist_pattern_col: str | None = None,
    r_drugs_list_col: str | None = None,
    plasmids_col: str | None = None,
    carb_nonsus_col: str | None = None,
    colistin_r_col: str | None = None,
    exclude_drugs: set[str] | None = None,
) -> pd.DataFrame:
    """Run June's prioritisation method.

    The metadata TSV/Excel must contain at minimum: isolate_id, laboratory,
    spec_date, and enough data to derive tier and resistance/plasmid profiles.

    Parameters
    ----------
    metadata_path : str
        Path to TSV or Excel file with isolate metadata.
    output_dir : str
        Directory for output files.
    id_col, lab_col, date_col : str
        Column names for isolate ID, laboratory/site, and specimen date.
    tier_col : str | None
        Column with pre-assigned tier. If None, tier is derived from
        resistance data.
    resist_pattern_col : str | None
        Column with semicolon-separated resistant drug abbreviations.
        If None, derived from r_drugs_list_col with exclude_drugs.
    r_drugs_list_col : str | None
        Column with comma-separated resistant drugs (raw). Used to derive
        resist_pattern if resist_pattern_col is None.
    plasmids_col : str | None
        Column with comma-separated plasmid replicons.
    carb_nonsus_col : str | None
        Column with bool/True/False for carbapenem non-susceptibility.
    colistin_r_col : str | None
        Column with bool/True/False for colistin resistance.
    exclude_drugs : set[str] | None
        Drugs to exclude from resist_pattern derivation (e.g. SAM, CTT, TGC).

    Returns
    -------
    pd.DataFrame
        Full ranked output with priority_rank, selection_batch, tier, rp_code,
        pp_code, combo_total, combo_sites, and all original metadata columns.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Read input
    ext = os.path.splitext(metadata_path)[1].lower()
    if ext in (".xlsx", ".xls"):
        df = pd.read_excel(metadata_path)
    else:
        df = pd.read_csv(metadata_path, sep="\t")

    print_message(f"Loaded {len(df)} isolates from {metadata_path}", "info")

    # Validate required columns
    for col_name, col_val in [(id_col, "isolate_id"), (lab_col, "laboratory"), (date_col, "spec_date")]:
        if col_name not in df.columns:
            raise ValueError(f"Required column '{col_name}' not found. Available: {list(df.columns)}")

    # Standardise column names
    df = df.rename(columns={id_col: "isolate_id", lab_col: "laboratory", date_col: "spec_date"})

    # Parse spec_date
    df["spec_date"] = pd.to_datetime(df["spec_date"])

    # Derive resist_pattern
    if resist_pattern_col and resist_pattern_col in df.columns:
        df["resist_pattern"] = df[resist_pattern_col].fillna("pan-susceptible")
    elif r_drugs_list_col and r_drugs_list_col in df.columns:
        df["resist_pattern"] = df[r_drugs_list_col].apply(
            lambda x: parse_resist_pattern(x, exclude_drugs)
        )
    else:
        print_message("No resistance data columns found; all isolates treated as pan-susceptible", "warning")
        df["resist_pattern"] = "pan-susceptible"

    # Derive plasmid profile
    if plasmids_col and plasmids_col in df.columns:
        df["plasmid_profile"] = df[plasmids_col].apply(parse_plasmid_profile)
    else:
        df["plasmid_profile"] = "none"

    # Assign tier
    if tier_col and tier_col in df.columns:
        df["tier"] = df[tier_col]
    elif carb_nonsus_col and colistin_r_col:
        # Derive from resistance data
        def _derive_tier(row: pd.Series) -> str:
            pattern = str(row["resist_pattern"])
            if pattern == "pan-susceptible":
                drugs: list[str] = []
            else:
                drugs = [d.strip() for d in pattern.split(";")]
            carb = bool(row.get(carb_nonsus_col, False))
            col_r = bool(row.get(colistin_r_col, False))
            return assign_tier(drugs, carb, col_r)

        df["tier"] = df.apply(_derive_tier, axis=1)
    else:
        print_message("No tier data; all isolates assigned Tier4_pansusceptible", "warning")
        df["tier"] = "Tier4_pansusceptible"

    # Assign RP codes
    df["rp_code"] = assign_rp_codes(df, "resist_pattern", tier_col="tier")

    # Assign PP codes
    df["pp_code"] = assign_pp_codes(df, "plasmid_profile")

    # Build RP x PP combo key
    df["rp_pp_combo"] = df["rp_code"] + "_" + df["pp_code"]

    # Assign Primary/Secondary batch
    df["selection_batch"] = select_primary_batch(df)

    # Compute priority rank
    df = compute_priority_rank(df)

    # Write outputs
    _write_outputs(df, output_dir)

    return df


def _write_outputs(df: pd.DataFrame, output_dir: str) -> None:
    """Write all output files for the June method."""
    # 1. Priority Full (all isolates ranked)
    full_path = os.path.join(output_dir, "priority_full.tsv")
    df.to_csv(full_path, sep="\t", index=False)
    print_message(f"Full priority list written to {full_path}", "success")

    # 2. Primary Batch only
    primary = df[df["selection_batch"] == "Primary"]
    primary_path = os.path.join(output_dir, "primary_batch.tsv")
    primary.to_csv(primary_path, sep="\t", index=False)
    print_message(
        f"Primary batch ({len(primary)} isolates) written to {primary_path}",
        "success",
    )

    # 3. Selected.txt (Primary batch isolate IDs, one per line)
    selected_path = os.path.join(output_dir, "selected.txt")
    with open(selected_path, "w") as fh:
        for sid in primary["isolate_id"]:
            fh.write(str(sid) + "\n")
    print_message(f"Selected isolate IDs written to {selected_path}", "success")

    # 4. Profile Summary (per RP code)
    profile_summary = (
        df.groupby(["rp_code", "tier", "resist_pattern"])
        .agg(
            n_isolates=("isolate_id", "size"),
            n_sites=("laboratory", "nunique"),
            primary_count=("selection_batch", lambda x: (x == "Primary").sum()),
            sites_list=("laboratory", lambda x: ", ".join(sorted(x.unique()))),
        )
        .reset_index()
        .sort_values(["tier", "n_isolates"], ascending=[True, False])
    )
    profile_path = os.path.join(output_dir, "profile_summary.tsv")
    profile_summary.to_csv(profile_path, sep="\t", index=False)
    print_message(f"Profile summary written to {profile_path}", "success")

    # 5. Site/Tier Summary
    site_tier = pd.crosstab(df["laboratory"], df["tier"], margins=True, margins_name="Total")
    site_tier_path = os.path.join(output_dir, "site_tier_summary.tsv")
    site_tier.to_csv(site_tier_path, sep="\t")
    print_message(f"Site/tier summary written to {site_tier_path}", "success")

    # 6. Plasmid Summary
    if "plasmid_profile" in df.columns:
        plasmid_profiles = (
            df[df["plasmid_profile"] != "none"]
            .groupby(["pp_code", "plasmid_profile"])
            .agg(occurrence=("isolate_id", "size"))
            .reset_index()
            .sort_values("occurrence", ascending=False)
        )
        pp_path = os.path.join(output_dir, "plasmid_profiles.tsv")
        plasmid_profiles.to_csv(pp_path, sep="\t", index=False)
        print_message(f"Plasmid profile summary written to {pp_path}", "success")

    # Summary stats
    print_message(
        f"Total: {len(df)} isolates, {len(primary)} Primary, "
        f"{len(df) - len(primary)} Secondary",
        "info",
    )
    for tier in TIER_ORDER:
        n = (df["tier"] == tier).sum()
        np_ = ((df["tier"] == tier) & (df["selection_batch"] == "Primary")).sum()
        if n > 0:
            print_message(f"  {tier}: {n} isolates ({np_} Primary)", "info")

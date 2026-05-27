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

BATCH_ORDER = {"Site_Guarantee": 0, "Primary": 1, "Secondary": 2}


def select_site_guarantees(df: pd.DataFrame) -> set:
    """Select the highest-ranked isolate from each unique site.

    For each site, the "best" isolate is the one with the best tier (lowest
    TIER_ORDER value), then the most frequent RP x PP combo (highest
    combo_total), then the most recent spec_date.

    Returns a set of DataFrame index values for the guaranteed isolates.
    """
    df = df.copy()
    df["_tier_ord"] = df["tier"].map(TIER_ORDER)
    df["_spec_date_ts"] = pd.to_datetime(df["spec_date"])

    # Compute combo frequency to use as a tiebreaker
    combo_counts = df["rp_pp_combo"].value_counts().to_dict()
    df["_combo_freq"] = df["rp_pp_combo"].map(combo_counts)

    # Sort: best tier first, then most frequent combo, then most recent date,
    # then isolate_id for determinism
    df = df.sort_values(
        ["_tier_ord", "_combo_freq", "_spec_date_ts", "isolate_id"],
        ascending=[True, False, False, True],
    )

    # Pick the first (best) isolate per site
    guaranteed_idx = set()
    for _site, group in df.groupby("laboratory"):
        guaranteed_idx.add(group.index[0])

    df.drop(columns=["_tier_ord", "_spec_date_ts", "_combo_freq"], inplace=True)
    return guaranteed_idx


def select_primary_batch(
    df: pd.DataFrame,
    guarantee_sites: bool = True,
) -> pd.Series:
    """Assign Site_Guarantee/Primary/Secondary batch labels.

    When *guarantee_sites* is True (default):
      - Site_Guarantee = 1 per sentinel site (highest-ranked isolate from
        each unique site). This ensures every site contributes at least one
        sequence, even if it only has Tier 4 isolates.
      - Primary = 1 per site per RP x PP combination (most recent spec_date),
        excluding any isolates already tagged as Site_Guarantee.
      - Secondary = all others.

    When *guarantee_sites* is False, Site_Guarantee is skipped and the
    original Primary/Secondary logic applies.

    Returns a Series of batch labels aligned with df index.
    """
    df = df.copy()

    # Start with all Secondary
    result = pd.Series("Secondary", index=df.index)

    # Compute combo frequency for site guarantee selection
    if guarantee_sites:
        combo_counts = df["rp_pp_combo"].value_counts().to_dict()
        df["_combo_freq"] = df["rp_pp_combo"].map(combo_counts)
        guaranteed_idx = select_site_guarantees(df)
        result[result.index.isin(guaranteed_idx)] = "Site_Guarantee"
        if "_combo_freq" in df.columns:
            df.drop(columns=["_combo_freq"], inplace=True)

    # Primary: 1 per site per RP x PP combo (most recent date), excluding guarantees
    df["_spec_date_ts"] = pd.to_datetime(df["spec_date"])
    df = df.sort_values(
        ["rp_pp_combo", "laboratory", "_spec_date_ts"],
        ascending=[True, True, False],
    )

    # Only consider non-guarantee isolates for Primary assignment
    non_guarantee_mask = result != "Site_Guarantee"
    df_for_primary = df[non_guarantee_mask.reindex(df.index)]

    is_first = ~df_for_primary.duplicated(
        subset=["rp_pp_combo", "laboratory"], keep="first"
    )
    result[is_first[is_first].index] = "Primary"

    return result.reindex(df.index)


def compute_priority_rank(df: pd.DataFrame) -> pd.DataFrame:
    """Compute the full priority ranking.

    Expects columns: isolate_id, laboratory, spec_date, tier, rp_code,
    pp_code, rp_pp_combo, selection_batch.

    Adds: priority_rank, combo_total (total isolates per RP x PP combo),
    combo_sites (number of unique sites per combo).

    The sort hierarchy:
      - Site_Guarantee isolates rank first (sorted by tier within them)
      - Then the standard 5-level sort for Primary and Secondary:
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

    # _is_guarantee: 0 for Site_Guarantee, 1 for everything else.
    # This ensures guarantees come first, then within the non-guarantee
    # block the standard tier > batch > combo sort applies.
    df["_is_guarantee"] = (df["selection_batch"] != "Site_Guarantee").astype(int)

    df = df.sort_values(
        [
            "_is_guarantee",  # Level 0: Site_Guarantee block first
            "_tier_ord",      # Level 1: Tier
            "_batch_ord",     # Level 2: Primary before Secondary (within non-guarantee)
            "combo_total",    # Level 3: most widespread combo first
            "combo_sites",    # Level 3 tiebreaker: most sites
            "laboratory",     # Level 4: site alphabetical
            "_spec_date_ts",  # Level 5: most recent date first
            "isolate_id",     # Final tiebreaker: accession_no
        ],
        ascending=[True, True, True, False, False, True, False, True],
    )

    df["priority_rank"] = range(1, len(df) + 1)

    # Drop temp columns
    df = df.drop(columns=["_tier_ord", "_batch_ord", "_spec_date_ts", "_is_guarantee"])

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
    guarantee_sites: bool = True,
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
    guarantee_sites : bool
        If True (default), guarantee at least 1 isolate from every sentinel
        site before applying the standard Primary batch selection. Isolates
        selected this way are tagged selection_batch = "Site_Guarantee".

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

    # Assign Site_Guarantee/Primary/Secondary batch
    df["selection_batch"] = select_primary_batch(df, guarantee_sites=guarantee_sites)

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

    # 2. Primary Batch (Site_Guarantee + Primary)
    primary = df[df["selection_batch"].isin(("Site_Guarantee", "Primary"))]
    primary_path = os.path.join(output_dir, "primary_batch.tsv")
    primary.to_csv(primary_path, sep="\t", index=False)
    n_guarantee = (primary["selection_batch"] == "Site_Guarantee").sum()
    n_primary = (primary["selection_batch"] == "Primary").sum()
    print_message(
        f"Primary batch ({len(primary)} isolates: {n_guarantee} site guarantees, "
        f"{n_primary} primary) written to {primary_path}",
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
    n_secondary = (df["selection_batch"] == "Secondary").sum()
    print_message(
        f"Total: {len(df)} isolates, {len(primary)} selected "
        f"({n_guarantee} site guarantees, {n_primary} primary), "
        f"{n_secondary} secondary",
        "info",
    )
    for tier in TIER_ORDER:
        n = (df["tier"] == tier).sum()
        ng = ((df["tier"] == tier) & (df["selection_batch"] == "Site_Guarantee")).sum()
        np_ = ((df["tier"] == tier) & (df["selection_batch"] == "Primary")).sum()
        if n > 0:
            print_message(
                f"  {tier}: {n} isolates ({ng} site guarantees, {np_} primary)",
                "info",
            )


# ---------------------------------------------------------------------------
# Flat CSV input
# ---------------------------------------------------------------------------

_CSV_REQUIRED_COLS = {"isolate_id", "site"}
_CSV_TIER_DERIVATION_COLS = {"resist_pattern", "drug_classes", "carb_nonsus"}


def run_june_from_csv(
    csv_path: str,
    output_dir: str,
    guarantee_sites: bool = True,
) -> pd.DataFrame:
    """Run June's prioritisation from a single flat CSV file.

    The CSV must contain:
      - isolate_id
      - site (sentinel site name)
      - tier (pre-assigned) OR resist_pattern + drug_classes + carb_nonsus
      - replicons (comma-separated PlasmidFinder replicon list)

    Optional columns:
      - spec_date (used as a tiebreaker; defaults to 2000-01-01 if absent)
      - rp_code / pp_code (pre-computed; skip derivation if both present)

    Parameters
    ----------
    csv_path : str
        Path to the flat CSV file.
    output_dir : str
        Directory for output files.
    guarantee_sites : bool
        If True (default), guarantee at least 1 isolate per sentinel site.

    Returns
    -------
    pd.DataFrame
        Full ranked output.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Read input — support both CSV and TSV based on extension
    ext = os.path.splitext(csv_path)[1].lower()
    if ext in (".tsv",):
        df = pd.read_csv(csv_path, sep="\t")
    elif ext in (".xlsx", ".xls"):
        df = pd.read_excel(csv_path)
    else:
        df = pd.read_csv(csv_path)

    print_message(f"Loaded {len(df)} isolates from {csv_path}", "info")

    # Validate required columns
    missing = _CSV_REQUIRED_COLS - set(df.columns)
    if missing:
        raise ValueError(
            f"Required columns missing: {missing}. "
            f"Available: {list(df.columns)}"
        )

    # Standardise: rename site -> laboratory
    df = df.rename(columns={"site": "laboratory"})

    # Handle spec_date (optional)
    if "spec_date" in df.columns:
        df["spec_date"] = pd.to_datetime(df["spec_date"])
    else:
        df["spec_date"] = pd.Timestamp("2000-01-01")

    # Derive or use tier
    if "tier" in df.columns:
        df["tier"] = df["tier"].fillna("Tier4_pansusceptible")
    elif _CSV_TIER_DERIVATION_COLS.issubset(set(df.columns)):
        # Derive tier from resist_pattern, drug_classes, carb_nonsus
        def _derive_tier_csv(row: pd.Series) -> str:
            pattern = str(row["resist_pattern"])
            if pattern in ("pan-susceptible", "", "nan"):
                drugs: list[str] = []
            else:
                drugs = [d.strip() for d in pattern.split(";")]
            carb = _parse_bool(row["carb_nonsus"])
            # Colistin resistance: check if COL is in the drug list
            colistin_r = "COL" in [d.upper() for d in drugs]
            return assign_tier(drugs, carb, colistin_r)

        df["tier"] = df.apply(_derive_tier_csv, axis=1)
    else:
        print_message(
            "No tier column and insufficient columns to derive it; "
            "all isolates assigned Tier4_pansusceptible",
            "warning",
        )
        df["tier"] = "Tier4_pansusceptible"

    # Derive resist_pattern if not present
    if "resist_pattern" not in df.columns:
        df["resist_pattern"] = "pan-susceptible"
    else:
        df["resist_pattern"] = df["resist_pattern"].fillna("pan-susceptible")

    # Derive plasmid_profile from replicons column
    if "replicons" in df.columns:
        df["plasmid_profile"] = df["replicons"].apply(parse_plasmid_profile)
    else:
        df["plasmid_profile"] = "none"

    # Use pre-computed RP/PP codes or derive them
    if "rp_code" in df.columns and "pp_code" in df.columns:
        print_message("Using pre-computed rp_code and pp_code columns", "info")
    else:
        if "rp_code" not in df.columns:
            df["rp_code"] = assign_rp_codes(df, "resist_pattern", tier_col="tier")
        if "pp_code" not in df.columns:
            df["pp_code"] = assign_pp_codes(df, "plasmid_profile")

    # Build RP x PP combo key
    df["rp_pp_combo"] = df["rp_code"] + "_" + df["pp_code"]

    # Assign batch
    df["selection_batch"] = select_primary_batch(df, guarantee_sites=guarantee_sites)

    # Compute priority rank
    df = compute_priority_rank(df)

    # Write outputs
    _write_outputs(df, output_dir)

    return df


def _parse_bool(val: object) -> bool:
    """Parse a value to bool, handling common CSV representations."""
    if isinstance(val, bool):
        return val
    if isinstance(val, (int, float)):
        return bool(val)
    s = str(val).strip().lower()
    return s in ("true", "1", "yes", "y")

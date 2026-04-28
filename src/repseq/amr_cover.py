"""Kleborate runner and greedy set-cover for AMR gene + plasmid replicon diversity."""

from __future__ import annotations

import glob
import json
import os
import re
import shutil
import subprocess
from pathlib import Path

import click
import pandas as pd

from repseq.log import print_message

# Kleborate v3 acquired resistance gene column suffixes (under klebsiella_pneumo_complex__amr__)
AMR_ACQUIRED_SUFFIXES = [
    "AGly_acquired", "Bla_acquired", "Bla_Carb_acquired", "Bla_ESBL_acquired",
    "Bla_ESBL_inhR_acquired", "Bla_inhR_acquired", "Col_acquired",
    "Fcyn_acquired", "Flq_acquired", "Gly_acquired", "MLS_acquired",
    "Phe_acquired", "Rif_acquired", "Sul_acquired", "Tet_acquired",
    "Tgc_acquired", "Tmt_acquired",
]

# Legacy Kleborate v2 column names (fallback)
LEGACY_GENE_COLS = [
    "AGly_genes", "Bla_genes", "Col_genes", "Fcyn_genes", "Flq_genes",
    "Gly_genes", "MLS_genes", "Ntmdz_genes", "Phe_genes", "Rif_genes",
    "Sul_genes", "Tmt_genes", "Tet_genes", "Tgc_genes",
]

# Possible replicon column names
REPLICON_COLS = ["plasmid_replicons"]


def run_kleborate(assemblies_dir: str, output_dir: str) -> str:
    """Run Kleborate on assemblies, return path to the merged output TSV."""
    kleb_outdir = os.path.join(output_dir, "kleborate_output")
    os.makedirs(kleb_outdir, exist_ok=True)

    # Gather assembly files
    patterns = ["*.fasta", "*.fa", "*.fna"]
    assembly_files: list[str] = []
    for pat in patterns:
        assembly_files.extend(glob.glob(os.path.join(assemblies_dir, pat)))
    assembly_files = sorted(set(assembly_files))

    if not assembly_files:
        raise click.ClickException(
            f"No assembly files (.fasta/.fa/.fna) found in {assemblies_dir}"
        )

    cmd = [
        "kleborate",
        "-a",
    ] + assembly_files + [
        "-o", kleb_outdir,
        "-p", "kpsc",
    ]

    print_message(f"Running Kleborate on {len(assembly_files)} assemblies...", "info")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise click.ClickException(f"Kleborate failed:\n{result.stderr}")

    # Find the output file -- Kleborate v3 creates *_output.txt in the outdir
    output_candidates = glob.glob(os.path.join(kleb_outdir, "*_output.txt"))
    if not output_candidates:
        raise click.ClickException(
            f"Kleborate produced no output files in {kleb_outdir}"
        )

    # Use the first (and typically only) output file
    kleb_output = output_candidates[0]

    # Copy to a convenient location
    final_path = os.path.join(output_dir, "kleborate.tsv")
    shutil.copy2(kleb_output, final_path)
    print_message(f"Kleborate output written to {final_path}", "success")
    return final_path


def _find_amr_columns(df: pd.DataFrame) -> list[str]:
    """Find AMR acquired gene columns in the Kleborate output."""
    cols = df.columns.tolist()

    # Try Kleborate v3 prefixed columns first
    v3_cols = [c for c in cols if "__amr__" in c and "_acquired" in c]
    if v3_cols:
        return v3_cols

    # Fallback to legacy v2 columns
    v2_cols = [c for c in LEGACY_GENE_COLS if c in cols]
    if v2_cols:
        return v2_cols

    # Try any column with 'acquired' or 'genes' in the name
    fallback = [c for c in cols if "acquired" in c.lower() or "genes" in c.lower()]
    return fallback


def _find_replicon_columns(df: pd.DataFrame) -> list[str]:
    """Find replicon/plasmid columns in the Kleborate output."""
    cols = df.columns.tolist()
    return [c for c in REPLICON_COLS if c in cols]


def parse_kleborate(kleborate_path: str) -> tuple[pd.DataFrame, list[str]]:
    """Parse Kleborate TSV into a binary presence/absence matrix.

    Returns (binary_matrix, feature_names) where binary_matrix is indexed by sample ID.
    """
    df = pd.read_csv(kleborate_path, sep="\t")

    # Determine sample ID column -- Kleborate v3 uses 'strain'
    id_col = None
    for candidate in ["strain", "sample", "assembly", "name"]:
        if candidate in df.columns:
            id_col = candidate
            break
    if id_col is None:
        raise click.ClickException(
            f"Cannot find sample ID column in Kleborate output. Columns: {list(df.columns)}"
        )

    # Clean sample IDs: strip path prefix and extension
    df[id_col] = df[id_col].apply(lambda x: Path(str(x)).stem)
    df = df.set_index(id_col)

    amr_columns = _find_amr_columns(df)
    replicon_columns = _find_replicon_columns(df)

    print_message(f"Found {len(amr_columns)} AMR columns, {len(replicon_columns)} replicon columns", "info")

    features: dict[str, dict[str, int]] = {}

    # Process AMR gene columns
    for col in amr_columns:
        # Use a cleaner column label for feature names
        short_name = col.split("__")[-1] if "__" in col else col
        for sample_id, val in df[col].items():
            if sample_id not in features:
                features[sample_id] = {}
            if pd.isna(val) or str(val).strip() in ("", "-", "none", "None"):
                continue
            # Genes are semicolon-separated; also handle comma-separated
            raw = str(val)
            genes = [
                g.strip() for g in raw.replace(",", ";").split(";")
                if g.strip() and g.strip() != "-"
            ]
            for gene in genes:
                feat_name = f"AMR:{short_name}:{gene}"
                features[sample_id][feat_name] = 1

    # Process replicon columns
    for col in replicon_columns:
        for sample_id, val in df[col].items():
            if sample_id not in features:
                features[sample_id] = {}
            if pd.isna(val) or str(val).strip() in ("", "-", "none", "None"):
                continue
            replicons = [
                r.strip() for r in str(val).split(";")
                if r.strip() and r.strip() != "-"
            ]
            for rep in replicons:
                feat_name = f"REP:{rep}"
                features[sample_id][feat_name] = 1

    # Include any sample that had no features at all
    for sample_id in df.index:
        if sample_id not in features:
            features[sample_id] = {}

    # Build binary matrix
    all_features = sorted({f for fs in features.values() for f in fs})
    if not all_features:
        print_message("No AMR/replicon features found in Kleborate output.", "warning")
        binary_df = pd.DataFrame(index=sorted(features.keys()))
        return binary_df, []

    rows = []
    sample_ids = sorted(features.keys())
    for sid in sample_ids:
        row = {f: features[sid].get(f, 0) for f in all_features}
        rows.append(row)

    binary_df = pd.DataFrame(rows, index=sample_ids, columns=all_features)
    return binary_df, all_features


def _find_plasmidfinder_db() -> str:
    """Locate PlasmidFinder database directory, searching conda/pixi share paths."""
    pf_bin = shutil.which("plasmidfinder.py")
    if pf_bin:
        # Walk up from bin/ to find share/plasmidfinder*/database
        bin_dir = Path(pf_bin).parent
        prefix = bin_dir.parent
        for candidate in sorted(prefix.glob("share/plasmidfinder*/database"), reverse=True):
            if candidate.is_dir():
                return str(candidate)
    # Fallback: try common system paths
    for fallback in [
        "/usr/share/plasmidfinder/database",
        "/opt/conda/share/plasmidfinder/database",
    ]:
        if Path(fallback).is_dir():
            return fallback
    # Last resort: return directory name only, let PlasmidFinder resolve it
    return "."


def run_plasmidfinder(assemblies_dir: str, output_dir: str) -> str:
    """Run PlasmidFinder on all assemblies, return path to merged results TSV."""
    pf_outdir = os.path.join(output_dir, "plasmidfinder_output")
    os.makedirs(pf_outdir, exist_ok=True)

    patterns = ["*.fasta", "*.fa", "*.fna"]
    assembly_files: list[str] = []
    for pat in patterns:
        assembly_files.extend(glob.glob(os.path.join(assemblies_dir, pat)))
    assembly_files = sorted(set(assembly_files))

    if not assembly_files:
        raise click.ClickException(f"No assembly files found in {assemblies_dir}")

    print_message(f"Running PlasmidFinder on {len(assembly_files)} assemblies...", "info")
    merged_rows = []
    for asm in assembly_files:
        sample_id = Path(asm).stem
        sample_outdir = os.path.join(pf_outdir, sample_id)
        os.makedirs(sample_outdir, exist_ok=True)
        db_path = _find_plasmidfinder_db()
        cmd = [
            "plasmidfinder.py", "-i", asm, "-o", sample_outdir,
            "-p", db_path, "-d", "enterobacteriaceae",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print_message(
                f"PlasmidFinder failed for {sample_id}: {result.stderr.strip()[:200]}",
                "warning",
            )
        plasmids_found: list[str] = []

        json_path = os.path.join(sample_outdir, "data.json")
        tsv_path = os.path.join(sample_outdir, "results_tab.tsv")

        if os.path.exists(json_path):
            try:
                with open(json_path) as jfh:
                    data = json.load(jfh)
                for db_hits in data.get("plasmidfinder", {}).get("results", {}).values():
                    for subdb_hits in db_hits.values():
                        for hit in subdb_hits.values():
                            p = hit.get("plasmid", "")
                            if p and p not in ("-", "nan", ""):
                                plasmids_found.append(p)
            except Exception:
                pass
        elif os.path.exists(tsv_path):
            try:
                pf_df = pd.read_csv(tsv_path, sep="\t")
                plasmid_col = None
                for candidate in ["Plasmid", "PlasmidFinder", "plasmid"]:
                    if candidate in pf_df.columns:
                        plasmid_col = candidate
                        break
                if plasmid_col is None:
                    for col in pf_df.columns:
                        if "plasmid" in col.lower():
                            plasmid_col = col
                            break
                if plasmid_col and not pf_df.empty:
                    plasmids_found = [
                        str(v) for v in pf_df[plasmid_col]
                        if str(v) not in ("-", "nan", "")
                    ]
            except Exception:
                pass

        for p in plasmids_found:
            merged_rows.append({"sample_id": sample_id, "Plasmid": p})

    merged_path = os.path.join(output_dir, "plasmidfinder.tsv")
    if merged_rows:
        pd.DataFrame(merged_rows).to_csv(merged_path, sep="\t", index=False)
    else:
        pd.DataFrame(columns=["sample_id", "Plasmid"]).to_csv(merged_path, sep="\t", index=False)
    n_hits = len(merged_rows)
    print_message(f"PlasmidFinder output written to {merged_path} ({n_hits} replicon hits)", "success")
    return merged_path


def parse_plasmidfinder(pf_path: str, binary_matrix: pd.DataFrame) -> pd.DataFrame:
    """Parse PlasmidFinder TSV and add replicon columns to binary_matrix."""
    try:
        df = pd.read_csv(pf_path, sep="\t")
    except Exception:
        print_message("Could not parse PlasmidFinder output.", "warning")
        return binary_matrix

    if df.empty or "Plasmid" not in df.columns or "sample_id" not in df.columns:
        print_message(
            f"PlasmidFinder output is empty or missing expected columns (found: {list(df.columns)}).",
            "warning",
        )
        return binary_matrix

    for _, row in df.iterrows():
        sample_id = str(row["sample_id"])
        plasmid = str(row["Plasmid"]).strip()
        if not plasmid or plasmid in ("-", "nan"):
            continue
        feat_name = f"REP:{plasmid}"
        if sample_id in binary_matrix.index:
            if feat_name not in binary_matrix.columns:
                binary_matrix[feat_name] = 0
            binary_matrix.loc[sample_id, feat_name] = 1

    rep_cols = [c for c in binary_matrix.columns if c.startswith("REP:")]
    print_message(f"PlasmidFinder added {len(rep_cols)} replicon features", "info")
    return binary_matrix


def greedy_set_cover(
    binary_matrix: pd.DataFrame,
    exclude_samples: list[str],
    n_amr: int,
) -> list[str]:
    """Greedy set-cover selection for AMR + replicon diversity.

    Picks the sample that adds the most uncovered features at each step,
    excluding samples already selected by PARNAS.
    """
    if n_amr <= 0:
        return []

    if binary_matrix.empty or len(binary_matrix.columns) == 0:
        available = [s for s in binary_matrix.index if s not in exclude_samples]
        selected = available[:n_amr]
        print_message(f"No AMR features to cover; selected {len(selected)} samples arbitrarily.", "info")
        return selected

    available = {s for s in binary_matrix.index if s not in exclude_samples}
    if not available:
        print_message("No samples available for AMR set cover.", "warning")
        return []

    # Track which features are already covered by excluded (PARNAS) samples
    covered: set[str] = set()
    for sid in exclude_samples:
        if sid in binary_matrix.index:
            row = binary_matrix.loc[sid]
            covered.update(row[row == 1].index.tolist())

    selected: list[str] = []
    for _ in range(min(n_amr, len(available))):
        best_sample = None
        best_new = -1
        for sid in available:
            row = binary_matrix.loc[sid]
            present = set(row[row == 1].index.tolist())
            new_features = len(present - covered)
            if new_features > best_new:
                best_new = new_features
                best_sample = sid
        if best_sample is None:
            break
        selected.append(best_sample)
        row = binary_matrix.loc[best_sample]
        covered.update(row[row == 1].index.tolist())
        available.discard(best_sample)

    print_message(f"Greedy set cover selected {len(selected)} samples: {selected}", "success")
    return selected


def add_cooccurrence_features(binary_matrix: pd.DataFrame) -> pd.DataFrame:
    """Add REP+AMR co-occurrence columns to binary_matrix.

    A co-occurrence feature CO:IncFII+TEM-1 is 1 only when both IncFII and
    TEM-1 are present in the same sample, capturing plasmid-gene linkage that
    binary presence/absence alone misses.  Only features present in ≥1 sample
    are retained so the matrix does not explode in size.
    """
    amr_cols = [c for c in binary_matrix.columns if c.startswith("AMR:")]
    rep_cols = [c for c in binary_matrix.columns if c.startswith("REP:")]

    if not amr_cols or not rep_cols:
        print_message("No REP or AMR columns — skipping co-occurrence features.", "warning")
        return binary_matrix

    co_data: dict[str, pd.Series] = {}
    for rep in rep_cols:
        for amr in amr_cols:
            col = f"CO:{rep[4:]}+{amr[4:]}"
            vals = (binary_matrix[rep] & binary_matrix[amr]).astype(int)
            if vals.sum() > 0:
                co_data[col] = vals

    if co_data:
        co_df = pd.DataFrame(co_data, index=binary_matrix.index)
        binary_matrix = pd.concat([binary_matrix, co_df], axis=1)
        print_message(f"Added {len(co_df.columns)} co-occurrence (REP+AMR) features", "info")
    else:
        print_message("No co-occurring REP+AMR pairs found.", "info")

    return binary_matrix


# ---------------------------------------------------------------------------
# Drug-class normalisation (ABRicate / hAMRonization → our short codes)
# ---------------------------------------------------------------------------

_CLASS_NORM: dict[str, str] = {
    "AMINOGLYCOSIDE": "AGly",
    "BETA-LACTAM": "Bla",
    "CARBAPENEM": "Bla_Carb",
    "CEPHALOSPORIN": "Bla",
    "COLISTIN": "Col",
    "FLUOROQUINOLONE": "Flq",
    "FOSFOMYCIN": "Fcyn",
    "GLYCOPEPTIDE": "Gly",
    "LINCOSAMIDE": "MLS",
    "MACROLIDE": "MLS",
    "MACROLIDE/LINCOSAMIDE/STREPTOGRAMIN": "MLS",
    "PHENICOL": "Phe",
    "RIFAMYCIN": "Rif",
    "SULFONAMIDE": "Sul",
    "TETRACYCLINE": "Tet",
    "TIGECYCLINE": "Tgc",
    "TRIMETHOPRIM": "Tmt",
}


def _normalise_class(raw: str) -> str:
    """Map a raw drug-class string (any case) to a repseq short code."""
    return _CLASS_NORM.get(raw.upper().strip(), raw.strip())


def _clean_replicon_name(gene: str) -> str:
    """Strip PlasmidFinder database suffixes from a replicon GENE name.

    e.g. 'IncFII(K)_1_Kpn3_JN233704' -> 'IncFII(K)'
    """
    return re.sub(r"_\d+[_.].*$", "", gene)


# ---------------------------------------------------------------------------
# ABRicate runner + parsers
# ---------------------------------------------------------------------------

def run_abricate(assemblies_dir: str, output_dir: str, db: str) -> str:
    """Run ABRicate on all assemblies with *db*, return path to merged TSV."""
    patterns = ["*.fasta", "*.fa", "*.fna"]
    assembly_files: list[str] = []
    for pat in patterns:
        assembly_files.extend(glob.glob(os.path.join(assemblies_dir, pat)))
    assembly_files = sorted(set(assembly_files))

    if not assembly_files:
        raise click.ClickException(f"No assembly files found in {assemblies_dir}")

    out_path = os.path.join(output_dir, f"abricate_{db}.tsv")
    print_message(f"Running ABRicate (db={db}) on {len(assembly_files)} assemblies...", "info")

    cmd = ["abricate", "--db", db, "--quiet"] + assembly_files
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise click.ClickException(f"ABRicate failed:\n{result.stderr}")

    with open(out_path, "w") as fh:
        fh.write(result.stdout)
    print_message(f"ABRicate ({db}) output written to {out_path}", "success")
    return out_path


def _abricate_sample_id(file_value: str) -> str:
    """Extract sample stem from ABRicate FILE column value."""
    return Path(file_value).stem


def parse_abricate_amr(abricate_path: str, binary_matrix: pd.DataFrame) -> pd.DataFrame:
    """Parse ABRicate AMR TSV and add AMR: features to *binary_matrix*.

    Expects columns: FILE, GENE, RESISTANCE, %IDENTITY, %COVERAGE.
    Uses RESISTANCE for drug class and GENE for gene name.
    """
    try:
        df = pd.read_csv(abricate_path, sep="\t")
    except Exception as e:
        print_message(f"Cannot parse ABRicate AMR output: {e}", "warning")
        return binary_matrix

    if df.empty or "GENE" not in df.columns:
        print_message("ABRicate AMR output is empty or missing GENE column.", "warning")
        return binary_matrix

    file_col = "#FILE" if "#FILE" in df.columns else "FILE"
    resistance_col = "RESISTANCE" if "RESISTANCE" in df.columns else None

    new_features: dict[str, dict[str, int]] = {}
    for _, row in df.iterrows():
        sample_id = _abricate_sample_id(str(row.get(file_col, "")))
        if sample_id not in binary_matrix.index:
            continue
        gene = str(row["GENE"]).strip()
        raw_class = str(row[resistance_col]).strip() if resistance_col else "Other"
        drug_class = _normalise_class(raw_class)
        feat = f"AMR:{drug_class}:{gene}"
        new_features.setdefault(sample_id, {})[feat] = 1

    for sample_id, feats in new_features.items():
        for feat, val in feats.items():
            if feat not in binary_matrix.columns:
                binary_matrix[feat] = 0
            binary_matrix.loc[sample_id, feat] = val

    amr_cols = [c for c in binary_matrix.columns if c.startswith("AMR:")]
    print_message(f"ABRicate AMR added features; matrix now has {len(amr_cols)} AMR columns", "info")
    return binary_matrix


def parse_abricate_replicons(abricate_path: str, binary_matrix: pd.DataFrame) -> pd.DataFrame:
    """Parse ABRicate plasmidfinder TSV and add REP: features to *binary_matrix*."""
    try:
        df = pd.read_csv(abricate_path, sep="\t")
    except Exception as e:
        print_message(f"Cannot parse ABRicate replicon output: {e}", "warning")
        return binary_matrix

    if df.empty or "GENE" not in df.columns:
        print_message("ABRicate replicon output is empty.", "warning")
        return binary_matrix

    file_col = "#FILE" if "#FILE" in df.columns else "FILE"

    for _, row in df.iterrows():
        sample_id = _abricate_sample_id(str(row.get(file_col, "")))
        if sample_id not in binary_matrix.index:
            continue
        replicon = _clean_replicon_name(str(row["GENE"]).strip())
        if not replicon or replicon in ("-", "nan"):
            continue
        feat = f"REP:{replicon}"
        if feat not in binary_matrix.columns:
            binary_matrix[feat] = 0
        binary_matrix.loc[sample_id, feat] = 1

    rep_cols = [c for c in binary_matrix.columns if c.startswith("REP:")]
    print_message(f"ABRicate replicons added; matrix now has {len(rep_cols)} REP columns", "info")
    return binary_matrix


# ---------------------------------------------------------------------------
# hAMRonization parser
# ---------------------------------------------------------------------------

def parse_hamronization(hamronization_path: str) -> tuple[pd.DataFrame, list[str]]:
    """Parse a hAMRonization TSV into a binary presence/absence matrix.

    hAMRonization standardises output from AMRFinder+, RGI, ResFinder etc.
    Expected columns: input_file_name, gene_symbol, drug_class.
    Returns (binary_matrix, feature_names).
    """
    try:
        df = pd.read_csv(hamronization_path, sep="\t")
    except Exception as e:
        raise click.ClickException(f"Cannot read hAMRonization file: {e}") from e

    # Locate required columns (names vary slightly across hAMRonization versions)
    file_col = next((c for c in df.columns if "file" in c.lower()), None)
    gene_col = next((c for c in df.columns if "gene_symbol" in c.lower() or c == "gene_name"), None)
    class_col = next((c for c in df.columns if "drug_class" in c.lower()), None)

    if file_col is None or gene_col is None:
        raise click.ClickException(
            f"Cannot find required columns in hAMRonization output. "
            f"Found: {list(df.columns)}"
        )

    features: dict[str, dict[str, int]] = {}
    for _, row in df.iterrows():
        sample_id = Path(str(row[file_col])).stem
        gene = str(row[gene_col]).strip()
        raw_class = str(row[class_col]).strip() if class_col else "Other"
        drug_class = _normalise_class(raw_class)
        if not gene or gene in ("-", "nan"):
            continue
        feat = f"AMR:{drug_class}:{gene}"
        features.setdefault(sample_id, {})[feat] = 1

    all_features = sorted({f for fs in features.values() for f in fs})
    if not all_features:
        print_message("No features found in hAMRonization output.", "warning")
        sample_ids = sorted({Path(str(r[file_col])).stem for _, r in df.iterrows()})
        return pd.DataFrame(index=sample_ids), []

    sample_ids = sorted(features.keys())
    rows = [{f: features[sid].get(f, 0) for f in all_features} for sid in sample_ids]
    binary_df = pd.DataFrame(rows, index=sample_ids, columns=all_features)
    print_message(
        f"hAMRonization: {binary_df.shape[0]} samples x {binary_df.shape[1]} AMR features", "info"
    )
    return binary_df, all_features

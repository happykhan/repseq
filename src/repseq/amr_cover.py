"""Kleborate runner and greedy set-cover for AMR gene + plasmid replicon diversity."""

import os
import subprocess
import glob
from pathlib import Path

import click
import pandas as pd


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
    """Run Kleborate on assemblies, return path to the merged output TSV.

    Kleborate v3 takes individual assembly files via -a and writes output
    to a directory (-o). The main output is
    klebsiella_pneumo_complex_output.txt inside that directory.
    """
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

    click.echo(f"Running Kleborate on {len(assembly_files)} assemblies...")
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
    import shutil
    shutil.copy2(kleb_output, final_path)
    click.echo(f"Kleborate output written to {final_path}")
    return final_path


def _find_amr_columns(df: pd.DataFrame) -> list[str]:
    """Find AMR acquired gene columns in the Kleborate output.

    Handles both Kleborate v3 (prefixed) and v2 (unprefixed) column names.
    """
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

    Returns:
        (binary_matrix, feature_names) where binary_matrix is indexed by sample ID.
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

    click.echo(f"Found {len(amr_columns)} AMR columns, {len(replicon_columns)} replicon columns")

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
            genes = [g.strip() for g in raw.replace(",", ";").split(";")
                     if g.strip() and g.strip() != "-"]
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
            replicons = [r.strip() for r in str(val).split(";")
                         if r.strip() and r.strip() != "-"]
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
        click.echo("Warning: no AMR/replicon features found in Kleborate output.")
        # Create an empty matrix with all samples
        binary_df = pd.DataFrame(index=sorted(features.keys()))
        return binary_df, []

    rows = []
    sample_ids = sorted(features.keys())
    for sid in sample_ids:
        row = {f: features[sid].get(f, 0) for f in all_features}
        rows.append(row)

    binary_df = pd.DataFrame(rows, index=sample_ids, columns=all_features)
    return binary_df, all_features


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
        # No features to cover; just pick samples not yet selected
        available = [s for s in binary_matrix.index if s not in exclude_samples]
        selected = available[:n_amr]
        click.echo(f"No AMR features to cover; selected {len(selected)} samples arbitrarily.")
        return selected

    available = [s for s in binary_matrix.index if s not in exclude_samples]
    if not available:
        click.echo("Warning: no samples available for AMR set cover.")
        return []

    # Track which features are already covered by excluded (PARNAS) samples
    covered = set()
    for sid in exclude_samples:
        if sid in binary_matrix.index:
            row = binary_matrix.loc[sid]
            covered.update(row[row == 1].index.tolist())

    selected: list[str] = []
    for _ in range(min(n_amr, len(available))):
        best_sample = None
        best_new = -1
        for sid in available:
            if sid in selected:
                continue
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
        available = [s for s in available if s != best_sample]

    click.echo(f"Greedy set cover selected {len(selected)} samples: {selected}")
    return selected

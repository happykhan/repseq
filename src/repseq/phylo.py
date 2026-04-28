"""Wrappers for mashtree and PARNAS."""

from __future__ import annotations

import glob
import os
import subprocess
from pathlib import Path

import click

from repseq.log import print_message


def find_assemblies(assemblies_dir: str) -> list[str]:
    """Return sorted list of assembly file paths in the given directory."""
    patterns = ["*.fasta", "*.fa", "*.fna"]
    files: list[str] = []
    for pat in patterns:
        files.extend(glob.glob(os.path.join(assemblies_dir, pat)))
    files = sorted(set(files))
    if not files:
        raise click.ClickException(
            f"No assembly files (.fasta/.fa/.fna) found in {assemblies_dir}"
        )
    return files


def run_mashtree(assemblies_dir: str, output_dir: str) -> str:
    """Run mashtree on assemblies, return path to tree file."""
    tree_path = os.path.join(output_dir, "tree.nwk")
    assembly_files = find_assemblies(assemblies_dir)
    cmd = ["mashtree"] + assembly_files
    print_message(f"Running mashtree on {len(assembly_files)} assemblies...", "info")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise click.ClickException(f"mashtree failed:\n{result.stderr}")
    # mashtree writes the tree to stdout
    with open(tree_path, "w") as fh:
        fh.write(result.stdout)
    print_message(f"Tree written to {tree_path}", "success")
    return tree_path


def run_parnas(tree_path: str, n_phylo: int, output_dir: str) -> list[str]:
    """Run PARNAS to select n_phylo representative samples. Returns list of sample IDs."""
    if n_phylo <= 0:
        return []

    diversity_csv = os.path.join(output_dir, "diversity_scores.csv")
    cmd = [
        "parnas",
        "-t", tree_path,
        "-n", str(n_phylo),
        "--diversity", diversity_csv,
    ]
    print_message(f"Running PARNAS to select {n_phylo} phylogenetically diverse samples...", "info")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise click.ClickException(f"PARNAS failed:\n{result.stderr}")

    # PARNAS outputs selected taxa to stdout, one per line
    selected = []
    for line in result.stdout.strip().splitlines():
        line = line.strip()
        if line:
            sample_id = Path(line).stem
            selected.append(sample_id)
    print_message(f"PARNAS selected {len(selected)} samples: {selected}", "success")
    return selected

#!/usr/bin/env bash
# Download test assemblies for dustmeselecta
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

BASE_URL="https://pub-e661bf7ded744bd79c156d3a4f4323ef.r2.dev"

SAMPLES=(
    Sample_25067ceb
    Sample_213aab83
    Sample_454f0c67
    Sample_eb789daf
    Sample_6d8e0e28
    Sample_58e78605
    Sample_f7e53a6b
    Sample_a30abdfc
    Sample_f937a7ba
    Sample_e236b6bd
    Sample_6d72ba33
    Sample_b54727ae
    Sample_da0a7978
    Sample_9922a0e9
    Sample_4d9ab532
    Sample_256a041e
    Sample_283376d3
    Sample_21984b5f
    Sample_6d9ceb52
    Sample_3cc326fa
)

for sample in "${SAMPLES[@]}"; do
    if [ ! -f "${sample}.fasta" ]; then
        echo "Downloading ${sample}.fasta..."
        wget -q "${BASE_URL}/${sample}.fasta"
    else
        echo "Skipping ${sample}.fasta (already exists)"
    fi
done

# Download sample sheet
if [ ! -f "real_typing_sample_sheet_f4f139b204.csv" ]; then
    echo "Downloading sample sheet..."
    wget -q "${BASE_URL}/real_typing_sample_sheet_f4f139b204.csv"
else
    echo "Skipping sample sheet (already exists)"
fi

echo "Done. $(ls *.fasta 2>/dev/null | wc -l) assemblies downloaded."

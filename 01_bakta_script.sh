#!/usr/bin/env bash
#SBATCH --job-name=bakta_array
#SBATCH --output=bakta_array_%A_%a.out
#SBATCH --error=bakta_array_%A_%a.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

# ============================================================
# Bakta annotation array job
# Annotates multiple genome FASTA files in parallel
# Requires: Bakta >= 1.9, full database downloaded
# Set BAKTA_DB to your local database path before running
# ============================================================

set -euo pipefail

# --- USER SETTINGS: edit before running ---
BAKTA_DB=/path/to/bakta_db          # path to full Bakta database
CONDA_ENV=bio                        # name of your conda environment

# --- Conda activation ---
if [ -f "${HOME}/miniconda3/etc/profile.d/conda.sh" ]; then
    source "${HOME}/miniconda3/etc/profile.d/conda.sh"
elif [ -f "${HOME}/anaconda3/etc/profile.d/conda.sh" ]; then
    source "${HOME}/anaconda3/etc/profile.d/conda.sh"
else
    echo "Could not find conda.sh — adjust path"
    exit 1
fi
conda activate "${CONDA_ENV}"

# --- Output directory ---
OUTDIR=./bakta_outputs
mkdir -p "${OUTDIR}"

# --- Input files ---
shopt -s nullglob
FILES=(*.fa *.fna *.fas)

if [ ${#FILES[@]} -eq 0 ]; then
    echo "No FASTA files found in current directory"
    exit 1
fi

# --- Select file for this array task ---
input="${FILES[$SLURM_ARRAY_TASK_ID]}"
base=$(basename "${input}")
prefix="${base%.*}"
sample_out="${OUTDIR}/${prefix}"
mkdir -p "${sample_out}"

echo "[$(date)] Running Bakta on ${input}"

bakta \
    --db "${BAKTA_DB}" \
    --output "${sample_out}" \
    --prefix "${prefix}" \
    --genus Escherichia \
    --species coli \
    --gram - \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --force \
    --keep-contig-headers \
    --translation-table 11 \
    --verbose \
    "${input}"

echo "[$(date)] Finished ${input}"

#!/usr/bin/env bash
#SBATCH --job-name=panaroo
#SBATCH --output=panaroo_%j.out
#SBATCH --error=panaroo_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

# ============================================================
# Panaroo pangenome analysis
# Constructs pangenome from Bakta-annotated GFF files
# Run from directory containing all .gff files
# Requires: Panaroo >= 1.3
# ============================================================

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"

# --- Collect all GFF files ---
ls *.gff > gff_list.txt
echo "Found $(wc -l < gff_list.txt) GFF files"

# --- Run Panaroo ---
panaroo \
    -i gff_list.txt \
    -o panaroo_output \
    --clean-mode strict \
    --remove-invalid-genes \
    -t "${SLURM_CPUS_PER_TASK}"

echo "[$(date)] Panaroo finished"

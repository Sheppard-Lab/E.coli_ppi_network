#!/usr/bin/env bash

# ============================================================
# 04_string_mapping.sh
# Maps Bakta-annotated K-12 proteome to STRING protein IDs
# using 100% sequence identity, filters interactions by
# confidence score, and integrates pangenome classifications
# to produce node and edge files for network visualisation.
#
# Usage: bash 04_string_mapping.sh
#
# Requirements: MMseqs2, awk, grep, cut, comm
#
# Input files required:
#   - bakta.faa                        : Bakta-annotated K-12 proteome
#   - string.faa                       : STRING K-12 proteome (511145)
#   - 511145.protein.links.v12.0.txt   : STRING interaction file
#   - gene_presence_absence_roary.csv  : Panaroo pangenome output
#
# ============================================================

set -euo pipefail

# --- USER SETTINGS ---
BAKTA_FAA=bakta.faa
STRING_FAA=511145.protein.sequences.v12.0.fa
STRING_LINKS=511145.protein.links.v12.0.txt
PANGENOME=gene_presence_absence_roary.csv
STRING_SCORE=700                        # minimum combined score threshold
K12_ACCESSION="GCA_000005845.2"         # K-12 reference accession in pangenome
CORE_THRESHOLD=95                       # % strains to be classified as core
THREADS=8

OUTDIR=./results
mkdir -p "${OUTDIR}"
TMP="${OUTDIR}/tmp_mmseqs"
mkdir -p "${TMP}"


# ============================================================
# STEP 1: MMseqs2 sequence identity search
# Bakta proteome vs STRING proteome
# ============================================================
echo ""
echo "[STEP 1] MMseqs2 search: Bakta vs STRING"

mmseqs easy-search \
    "${BAKTA_FAA}" \
    "${STRING_FAA}" \
    "${OUTDIR}/results.tsv" \
    "${TMP}" \
    --min-seq-id 1.0 \
    -c 1.0 \
    --cov-mode 0 \
    --format-output "query,target,pident,qcov,tcov,qlen,tlen,evalue" \
    --threads "${THREADS}"

# --- Verification ---
TOTAL_HITS=$(wc -l < "${OUTDIR}/results.tsv")
UNIQUE_BAKTA=$(cut -f1 "${OUTDIR}/results.tsv" | sort -u | wc -l)
UNIQUE_STRING=$(cut -f2 "${OUTDIR}/results.tsv" | sort -u | wc -l)
BAKTA_SIZE=$(grep -c "^>" "${BAKTA_FAA}")
STRING_SIZE=$(grep -c "^>" "${STRING_FAA}")


# ============================================================
# STEP 2: Filter to unambiguous one-to-one mappings
# Remove Bakta tags that hit multiple STRING entries
# ============================================================
echo ""
echo "[STEP 2] Filtering to one-to-one mappings"

cut -f1 "${OUTDIR}/results.tsv" | sort | uniq -c | \
    awk '$1==1 {print $2}' > "${TMP}/one_to_one_bakta.txt"

grep -F -f "${TMP}/one_to_one_bakta.txt" \
    "${OUTDIR}/results.tsv" > "${OUTDIR}/unambiguous_pairs.tsv"

# --- Verification ---
UNAMBIG=$(wc -l < "${OUTDIR}/unambiguous_pairs.tsv")
DROPPED=$((UNIQUE_BAKTA - UNAMBIG))


# ============================================================
# STEP 3: Filter STRING interactions by confidence score
# ============================================================
echo ""
echo "[STEP 3] Filtering STRING interactions (score >= ${STRING_SCORE})"

awk -v score="${STRING_SCORE}" \
    'NR==1 || $3>=score' \
    "${STRING_LINKS}" > "${OUTDIR}/string_${STRING_SCORE}.txt"


# ============================================================
# STEP 4: Replace STRING IDs with Bakta locus tags
# Only retain edges where BOTH partners are in mapped set
# ============================================================
echo ""
echo "[STEP 4] Replacing STRING IDs with Bakta locus tags"

awk 'NR==FNR{map[$2]=$1; next}
     NR==1{next}
     ($1 in map) && ($2 in map){
         printf "%s\t%s\t%s\n", map[$1], map[$2], $3
     }' \
    "${OUTDIR}/unambiguous_pairs.tsv" \
    "${OUTDIR}/string_${STRING_SCORE}.txt" > "${OUTDIR}/bakta_interactions_${STRING_SCORE}.txt"



# ============================================================
# STEP 5: Extract K-12 locus tags and pangenome classification
# Core = >= CORE_THRESHOLD% of strains
# Accessory = < CORE_THRESHOLD% of strains
# ============================================================
echo ""
echo "[STEP 5] Extracting K-12 gene classifications from pangenome"

# Find K-12 column dynamically
K12_COL=$(head -n1 "${PANGENOME}" | tr ',' '\n' | \
    grep -n "${K12_ACCESSION}" | cut -d: -f1)
echo "  K-12 column: ${K12_COL}  (expected: 15)"

awk -F',' -v col="${K12_COL}" -v threshold="${CORE_THRESHOLD}" \
    'NR==1{next}
     {
         locus=$col
         if(locus=="") next
         n=$4+0
         if(n>=threshold) class="core"
         else class="accessory"
         print locus, class, $1, $3
     }' "${PANGENOME}" > "${OUTDIR}/k12_gene_classes.txt"

# --- Verification ---
TOTAL_K12=$(wc -l < "${OUTDIR}/k12_gene_classes.txt")
CORE_COUNT=$(awk '$2=="core"' "${OUTDIR}/k12_gene_classes.txt" | wc -l)
ACC_COUNT=$(awk '$2=="accessory"' "${OUTDIR}/k12_gene_classes.txt" | wc -l)

echo "  Total K-12 genes:    ${TOTAL_K12}  (expected: 4328)"
echo "  Core genes:          ${CORE_COUNT}  (expected: 3322)"
echo "  Accessory genes:     ${ACC_COUNT}  (expected: 1006)"

# ============================================================
# STEP 6: Join pangenome classifications with STRING mappings
# ============================================================
echo ""
echo "[STEP 6] Joining pangenome classes with STRING mappings"

awk 'NR==FNR{class[$1]=$2; gene[$1]=$3; annot[$1]=$4; next}
     {
         locus=$1
         string_id=$2
         if(locus in class)
             print locus, string_id, class[locus], gene[locus], annot[locus]
     }' \
    "${OUTDIR}/k12_gene_classes.txt" \
    "${OUTDIR}/unambiguous_pairs.tsv" > "${OUTDIR}/k12_network_nodes.txt"

# --- Verification ---
TOTAL_NODES=$(wc -l < "${OUTDIR}/k12_network_nodes.txt")
CORE_NODES=$(awk '$3=="core"' "${OUTDIR}/k12_network_nodes.txt" | wc -l)
ACC_NODES=$(awk '$3=="accessory"' "${OUTDIR}/k12_network_nodes.txt" | wc -l)

echo "  Total nodes:         ${TOTAL_NODES}  (expected: 3805)"
echo "  Core nodes:          ${CORE_NODES}  (expected: 3068)"
echo "  Accessory nodes:     ${ACC_NODES}  (expected: 737)"

# ============================================================
# STEP 7: Build final network files for Cytoscape
# nodes.tsv  — node attributes with pangenome classification
# edges_typed.tsv — edges with interaction type annotation
# ============================================================
echo ""
echo "[STEP 7] Building final network files"

# Node file
echo -e "locus_tag\tstring_id\tgene_class\tgene_name\tannotation" \
    > "${OUTDIR}/nodes.tsv"
sed 's/ /\t/g' "${OUTDIR}/k12_network_nodes.txt" >> "${OUTDIR}/nodes.tsv"

# Edge file with interaction type
awk 'NR==FNR{class[$1]=$3; next}
     NR==1{print "protein1\tprotein2\tscore\tinteraction_type"; next}
     {
         c1=class[$1]; c2=class[$2]
         if(c1=="core" && c2=="core") type="core-core"
         else if(c1=="accessory" && c2=="accessory") type="accessory-accessory"
         else type="mixed"
         printf "%s\t%s\t%s\t%s\n", $1, $2, $3, type
     }' \
    "${OUTDIR}/nodes.tsv" \
    "${OUTDIR}/bakta_interactions_${STRING_SCORE}.txt" \
    > "${OUTDIR}/edges_typed.tsv"

# --- Verification ---
FINAL_EDGES=$(awk 'NR>1' "${OUTDIR}/edges_typed.tsv" | wc -l)
CORE_CORE=$(awk '$4=="core-core"' "${OUTDIR}/edges_typed.tsv" | wc -l)
ACC_ACC=$(awk '$4=="accessory-accessory"' "${OUTDIR}/edges_typed.tsv" | wc -l)
MIXED=$(awk '$4=="mixed"' "${OUTDIR}/edges_typed.tsv" | wc -l)

# Connected nodes
cut -f1 "${OUTDIR}/edges_typed.tsv" | grep -v protein1 | sort -u \
    > "${TMP}/edge_sources.txt"
cut -f2 "${OUTDIR}/edges_typed.tsv" | grep -v protein2 | sort -u \
    > "${TMP}/edge_targets.txt"
cat "${TMP}/edge_sources.txt" "${TMP}/edge_targets.txt" | sort -u \
    > "${TMP}/edge_all_proteins.txt"

grep "	core	" "${OUTDIR}/nodes.tsv" | cut -f1 | sort \
    > "${TMP}/core_nodes.txt"
grep "	accessory	" "${OUTDIR}/nodes.tsv" | cut -f1 | sort \
    > "${TMP}/accessory_nodes.txt"

CORE_CONNECTED=$(comm -12 "${TMP}/core_nodes.txt" \
    "${TMP}/edge_all_proteins.txt" | wc -l)
ACC_CONNECTED=$(comm -12 "${TMP}/accessory_nodes.txt" \
    "${TMP}/edge_all_proteins.txt" | wc -l)
CORE_ISOLATED=$(comm -23 "${TMP}/core_nodes.txt" \
    "${TMP}/edge_all_proteins.txt" | wc -l)
ACC_ISOLATED=$(comm -23 "${TMP}/accessory_nodes.txt" \
    "${TMP}/edge_all_proteins.txt" | wc -l)



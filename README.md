# E. coli K-12 Protein-Protein Interaction Network

Scripts and data for mapping *E. coli* K-12 pangenome gene classifications
(core vs accessory) onto STRING protein interaction data, and generating
network figures for visualisation.

---

## Overview

This repository accompanies the analysis of the *E. coli* K-12 protein
interaction network in the context of pangenome gene classification. The
pipeline:

1. Downloads 100 *E. coli* genome assemblies from NCBI
2. Annotates genomes with Bakta
3. Constructs a pangenome with Panaroo
4. Maps Bakta locus tags to STRING protein identifiers using MMseqs2
5. Integrates pangenome classifications (core/accessory) with STRING
   interaction data
6. Produces network files for visualisation and analysis
   with igraph

---

## Repository Structure

```
E.coli_ppi_network/
├── README.md
├── Supplementary Table 1.tsv        # Metadata for 100 E. coli genomes used
├── 01_bakta_script.sh               # Bakta annotation
├── 02_panaroo_script.sh             # Panaroo pangenome construction
├── 03_string_mapping.sh             # STRING mapping + network file generation
├── combined_BCD.py                  # Degree distribution and interaction analysis
├── louvain_hubs.py                  # Louvain community detection and hub annotation
├── gene_presence_absence_roary.csv  # Panaroo output: gene presence/absence
├── unambiguous_pairs.tsv            # Bakta ↔ STRING one-to-one mappings
├── k12_gene_classes.txt             # K-12 locus tags with core/accessory labels
├── nodes.tsv                        # Final node attribute file for Cytoscape
└── edges_typed.tsv                  # Final edge file with interaction type labels
```

## Input Data

### Genome Assemblies
100 *E. coli* genome assemblies were downloaded manually from NCBI.
Accessions and metadata for all 100 genomes are provided in
`Supplementary Table 1.tsv`.


### STRING Database Files
Download STRING v12.0 files for *E. coli* K-12 (taxonomy ID: 511145):
```bash
wget https://stringdb-downloads.org/download/protein.sequences.v12.0/511145.protein.sequences.v12.0.fa.gz
wget https://stringdb-downloads.org/download/protein.links.v12.0/511145.protein.links.v12.0.txt.gz
gunzip *.gz
```

---

## Pipeline

### Step 1 — Genome download
Genomes were downloaded manually from NCBI using accessions in
`Supplementary Table 1.tsv`.

### Step 2 — Bakta annotation
```bash
sbatch --array=0-99 01_bakta_script.sh
```
Annotates all genome assemblies in parallel using a SLURM array job.
Requires the full Bakta database — set `BAKTA_DB` path at the top of
the script before running.

Output: `bakta_outputs/<prefix>/<prefix>.faa` for each genome.
The K-12 proteome (`GCA_000005845.2`) is used as `bakta.faa` in Step 4.

### Step 3 — Panaroo pangenome
```bash
bash 02_panaroo_script.sh
```
Constructs the pangenome from all Bakta GFF outputs.

Output: `gene_presence_absence_roary.csv` (provided in this repository).

### Step 4 — STRING mapping and network preparation
```bash
bash 03_string_mapping.sh
```
This is the core mapping script. It:
- Searches Bakta K-12 proteome against STRING proteome at 100% identity
  and 100% bidirectional coverage (both query and target) using MMseqs2
- Filters to unambiguous one-to-one Bakta ↔ STRING mappings
- Filters STRING interactions to combined score ≥ 700 (high confidence)
- Extracts K-12 gene classifications from the pangenome output
- Produces final node and edge files for network visualisation

Key parameters (adjustable at top of script):
```
STRING_SCORE=700        # interaction confidence threshold
CORE_THRESHOLD=95       # % strains to classify as core
K12_ACCESSION=GCA_000005845.2
```

### Step 5 — Network analysis
```bash
python combined_BCD.py
python louvain_hubs.py
```
Input: `nodes.tsv`, `edges_typed.tsv`

### References 

Schwengers, O., Jelonek, L., Dieckmann, M. A., Beyvers, S., Blom, J., & Goesmann, A. (2021). Bakta: Rapid and standardized annotation of bacterial genomes via alignment-free sequence identification: Find out more about Bakta, the motivation, challenges and applications, here. Microbial Genomics, 7(11). https://doi.org/10.1099/mgen.0.000685

Steinegger, M., & Söding, J. (2017). MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature Biotechnology, 35(11), 1026–1028. https://doi.org/10.1038/nbt.3988

Szklarczyk, D., Kirsch, R., Koutrouli, M., Nastou, K., Mehryary, F., Hachilif, R., Gable, A. L., Fang, T., Doncheva, N. T., Pyysalo, S., Bork, P., Jensen, L. J., & von Mering, C. (2023). The STRING database in 2023: Protein–protein association networks and functional enrichment analyses for any sequenced genome of interest. Nucleic Acids Research, 51(D1), D638–D646. https://doi.org/10.1093/nar/gkac1000

Tonkin-Hill, G., MacAlasdair, N., Ruis, C., Weimann, A., Horesh, G., Lees, J. A., Gladstone, R. A., Lo, S., Beaudoin, C., Floto, R. A., Frost, S. D. W., Corander, J., Bentley, S. D., & Parkhill, J. (2020). Producing polished prokaryotic pangenomes with the Panaroo pipeline. Genome Biology, 21(1), 180. https://doi.org/10.1186/s13059-020-02090-4





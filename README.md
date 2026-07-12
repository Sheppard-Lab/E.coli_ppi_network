# E. coli K-12 PPI network figures

Scripts for generating protein–protein interaction network figures for the
E. coli K-12 pangenome perspective, using STRING v12.0 interaction data with
genes classified as core or accessory.

## Scripts
- `combined_BCD.py` — three stacked panels: interaction degree by gene class
  (core vs accessory), total interactions by type (core–core,
  accessory–accessory, core–accessory), and the number of accessory genes
  interacting with core vs other accessory genes.
- `louvain_hubs.py` — full network layout with Louvain community detection,
  nodes coloured by pangenome class and sized by degree, communities annotated
  with functional labels anchored to marker genes.

## Inputs
- `nodes.tsv` — columns: locus_tag, gene_class, gene_name
- `edges_typed.tsv` — columns: protein1, protein2, score, interaction_type

Update the file paths at the top of each script before running.

## Requirements
python-igraph, pandas, numpy, matplotlib

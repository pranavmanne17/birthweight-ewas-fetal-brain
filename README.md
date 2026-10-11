# Growth-associated EWAS signatures and fetal neurodevelopment

Code and output tables for:

> Manne P. *Maternal Micronutrient Deficiency, DNA Methylation, and Fetal Neurodevelopment in Low-Resource Settings: A Secondary Analysis of Growth-Associated EWAS Signatures.* bioRxiv (2026). doi: [add after posting]

The study re-annotates published DNA methylation signals for two nutrition-sensitive growth phenotypes: 914 birthweight-associated CpGs (729 genes) from the PACE consortium meta-analysis (Küpers et al., 2019) and three SOCS3 CpGs associated with child height in low- and middle-income country cohorts (Issarapu et al., 2023). No new methylation data were generated and no individual-level data were used.

## Contents

| Path | Description |
|---|---|
| `classify_genes.py` | Assigns each gene to one of five functional categories from its GO Biological Process terms (Methods 2.5) |
| `requirements.txt` | Python packages needed to run the script |
| `data/genes_729.txt` | The 729 nearest genes of the 914 birthweight CpGs, one per line |
| `data/cpgs_914.txt` | The 914 CpG probe IDs used as eFORGE input |
| `results/gene_classification.csv` | Keyword category, all matching categories, final category after manual review, and whether review changed it, for all 729 genes |
| `results/supplementary_tables_S1_S2.xlsx` | Supplementary Tables S1 (keyword sets) and S2 (manual reassignments) |
| `results/eforge/` | eFORGE v2.0 output for the full CpG set and the neurodevelopment subset (E081, E082, E091) |
| `results/brainspan/` | BrainSpan RNA-Seq expression values (RPKM) for the neurodevelopment-category genes and SOCS3 |

## Running the classification

```
pip install -r requirements.txt
python classify_genes.py
```

This writes `results/go_gene_classification_keyword.csv`, the automated keyword assignment before manual review. GO annotations change over time; the paper used annotations retrieved from MyGene.info on 21 June 2026, so a rerun may differ slightly. Manual corrections to 30 genes are listed in Supplementary Table S2.

Three gene symbols in the original list had been converted to dates by spreadsheet software. The script maps them back to MARCHF3, MARCHF9 and SEPTIN11.

## Data sources

- Birthweight EWAS summary statistics: Küpers LK et al. *Nat Commun* 2019;10:1893. doi:10.1038/s41467-019-09671-3 (Supplementary Data)
- SOCS3 height CpGs: Issarapu P et al. *Nat Commun* 2023;14:5200. doi:10.1038/s41467-023-40607-0
- Gene Ontology annotations: MyGene.info (Xin J et al. *Genome Biol* 2016;17:91)
- Chromatin states: Roadmap Epigenomics 15-state models via eFORGE v2.0 (Breeze CE et al. *Bioinformatics* 2019)
- Fetal brain expression: BrainSpan Atlas of the Developing Human Brain (www.brainspan.org)

## Licence

Code: MIT (see `LICENSE`). Data derived from the sources above remain subject to their original terms; cite the original studies when reusing them.

## Contact

Pranav Manne · pranavmanne17@gmail.com

# LPL-lactylation-CRLM-analysis
Reproducible analysis code for single-cell RNA-seq, pseudotime, cell-cell communication, spatial transcriptomics, and metabolic scoring used in the CRLM study.
# Analysis Code for the Lactate–LPL CRLM Study

This repository contains the major computational analysis scripts used in this study to improve transparency and reproducibility.

## Analysis Scripts

- `01_scRNA_analysis.R`: scRNA-seq preprocessing, clustering, cell-type annotation, inferCNV analysis, and metabolic scoring
- `02_monocle2_pseudotime.R`: Monocle2 pseudotime analysis of macrophages
- `03_cellchat_analysis.R`: CellChat-based cell–cell communication analysis
- `04_spatial_analysis.R`: Spatial transcriptomic scoring and spatial distribution analysis

## Data Availability

Processed scRNA-seq expression matrices, cell-level metadata, and cell-type annotations are available in the NGDC OMIX database.


## Gene Sets

The complete gene lists used for lactate metabolism, cholesterol metabolism, and macrophage signatures are provided in the Supplementary Materials of the manuscript.

## Metabolic Scoring

Metabolic scores were calculated using AUCell, UCell, singscore, and AddModuleScore. The four algorithm-derived scores were directly summed without additional weighting coefficients.

## Reproducibility

Fixed random seeds were used where applicable, and `sessionInfo()` was recorded to document the software environment.

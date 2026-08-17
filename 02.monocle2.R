#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle)
  library(dplyr)
  library(ggplot2)
  library(ClusterGVis)
  library(org.Hs.eg.db)
  library(clusterProfiler)
})

set.seed(2025)

# ------------------------------------------------------------------------------
# 1. Configuration
# ------------------------------------------------------------------------------

INPUT_RDS <- "data/Macro.rds"
OUT_DIR   <- "results03"
DATA_DIR  <- "data"

SAMPLE_COLUMN <- "orig.ident"
GROUP_COLUMN  <- "ScoringGroup"
SCORE_COLUMN  <- "Scoring"

N_CORES <- as.integer(Sys.getenv("N_CORES", "8"))

N_HVG <- 2000
DDR_NUM_DIM <- 5

PSEUDOTIME_QVAL <- 0.01
N_HEATMAP_GENES <- 100
HEATMAP_CLUSTERS <- 3

# Leave NULL to use Monocle 2 automatic trajectory orientation.
# If a root state is explicitly specified in the final manuscript,
# set the corresponding state number here and report it in Methods.
ROOT_STATE <- NULL

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(INPUT_RDS)) {
  stop("Input file not found: ", INPUT_RDS)
}

# ------------------------------------------------------------------------------
# 2. Load and validate macrophage Seurat object
# ------------------------------------------------------------------------------

scRNA <- readRDS(INPUT_RDS)
DefaultAssay(scRNA) <- "RNA"

required_meta <- c(SAMPLE_COLUMN, GROUP_COLUMN, SCORE_COLUMN)
missing_meta <- setdiff(required_meta, colnames(scRNA@meta.data))

if (length(missing_meta) > 0) {
  stop(
    "Missing required metadata column(s): ",
    paste(missing_meta, collapse = ", ")
  )
}

message("Cells: ", ncol(scRNA))
message("Genes: ", nrow(scRNA))
message("Samples: ", length(unique(scRNA@meta.data[[SAMPLE_COLUMN]])))

# ------------------------------------------------------------------------------
# 3. Unsupervised ordering-gene selection using Seurat HVGs
# ------------------------------------------------------------------------------

# Recompute HVGs within the macrophage subset so trajectory construction
# is independent of ScoringGroup.
scRNA <- FindVariableFeatures(
  scRNA,
  assay = "RNA",
  selection.method = "vst",
  nfeatures = N_HVG,
  verbose = FALSE
)

ordering_genes <- VariableFeatures(scRNA)

# Remove technical/housekeeping gene families commonly excluded from trajectory
# construction in this analysis.
ordering_genes <- ordering_genes[
  !grepl("^(HIST|MT-|HBA|HBB|RPL|RPS)", ordering_genes)
]

if (length(ordering_genes) < 10) {
  stop("Too few HVGs were retained for trajectory inference.")
}

write.table(
  ordering_genes,
  file.path(OUT_DIR, "monocle_ordering_genes_HVG.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# ------------------------------------------------------------------------------
# 4. Convert Seurat object to Monocle 2 CellDataSet
# ------------------------------------------------------------------------------

counts <- GetAssayData(scRNA, assay = "RNA", layer = "counts")

gene_ann <- data.frame(
  gene_short_name = rownames(counts),
  row.names = rownames(counts),
  stringsAsFactors = FALSE
)

cell_ann <- scRNA@meta.data

fd <- new("AnnotatedDataFrame", data = gene_ann)
pd <- new("AnnotatedDataFrame", data = cell_ann)

sc_cds <- newCellDataSet(
  counts,
  phenoData = pd,
  featureData = fd,
  expressionFamily = negbinomial.size(),
  lowerDetectionLimit = 1
)

sc_cds <- estimateSizeFactors(sc_cds)
sc_cds <- estimateDispersions(sc_cds)
sc_cds <- detectGenes(sc_cds, min_expr = 0.1)

ordering_genes <- intersect(ordering_genes, rownames(sc_cds))
sc_cds <- setOrderingFilter(sc_cds, ordering_genes)

pdf(file.path(OUT_DIR, "monocle_ordering_genes_HVG.pdf"), width = 6, height = 5)
print(plot_ordering_genes(sc_cds))
dev.off()

# ------------------------------------------------------------------------------
# 5. DDRTree trajectory inference
# ------------------------------------------------------------------------------

# Sample identity is included as a residual term to reduce sample-level effects.
sc_cds <- reduceDimension(
  sc_cds,
  max_components = 2,
  num_dim = DDR_NUM_DIM,
  reduction_method = "DDRTree",
  cores = N_CORES,
  residualModelFormulaStr = paste0("~", SAMPLE_COLUMN),
  verbose = TRUE
)

if (is.null(ROOT_STATE)) {
  sc_cds <- orderCells(sc_cds)
} else {
  sc_cds <- orderCells(sc_cds, root_state = ROOT_STATE)
}

# ------------------------------------------------------------------------------
# 6. Post hoc visualization of metabolic states along the trajectory
# ------------------------------------------------------------------------------

pdf(file.path(OUT_DIR, "monocle_state.pdf"), width = 7, height = 6)
print(plot_cell_trajectory(sc_cds, color_by = "State"))
dev.off()

pdf(file.path(OUT_DIR, "monocle_pseudotime.pdf"), width = 7, height = 6)
print(plot_cell_trajectory(sc_cds, color_by = "Pseudotime"))
dev.off()

pdf(file.path(OUT_DIR, "monocle_scoringgroup.pdf"), width = 7, height = 6)
print(
  plot_cell_trajectory(sc_cds, color_by = GROUP_COLUMN) +
    scale_color_manual(
      values = c(
        "Low" = "steelblue",
        "Medium" = "grey",
        "High" = "firebrick"
      )
    )
)
dev.off()

pdf(file.path(OUT_DIR, "monocle_scoring.pdf"), width = 7, height = 6)
print(
  plot_cell_trajectory(sc_cds, color_by = SCORE_COLUMN) +
    scale_color_gradientn(colors = c("grey90", "pink", "firebrick"))
)
dev.off()

flow <- pData(sc_cds)

density_plot <- ggplot(
  flow,
  aes(x = Pseudotime, fill = .data[[GROUP_COLUMN]])
) +
  geom_density(bw = 0.5, alpha = 0.5, colour = NA) +
  scale_fill_manual(
    values = c(
      "Low" = "steelblue",
      "Medium" = "grey",
      "High" = "firebrick"
    ),
    name = "Scoring Group"
  ) +
  labs(x = "Pseudotime", y = "Density") +
  theme_classic()

ggsave(
  file.path(OUT_DIR, "monocle_scoringgroup_density.pdf"),
  density_plot,
  width = 7,
  height = 6
)

# ------------------------------------------------------------------------------
# 7. Identify genes associated with pseudotime
# ------------------------------------------------------------------------------

# Pseudotime-associated genes are identified AFTER trajectory construction.
# This avoids using ScoringGroup to define the trajectory.
pseudotime_deg <- differentialGeneTest(
  sc_cds[ordering_genes, ],
  fullModelFormulaStr = "~sm.ns(Pseudotime)",
  cores = N_CORES
)

pseudotime_deg <- pseudotime_deg %>%
  tibble::rownames_to_column("gene") %>%
  arrange(qval)

write.csv(
  pseudotime_deg,
  file.path(OUT_DIR, "monocle_pseudotime_DEGs.csv"),
  row.names = FALSE
)

heatmap_genes <- pseudotime_deg %>%
  filter(!is.na(qval), qval < PSEUDOTIME_QVAL) %>%
  slice_head(n = N_HEATMAP_GENES) %>%
  pull(gene)

if (length(heatmap_genes) < 10) {
  warning(
    "Fewer than 10 genes passed q < ", PSEUDOTIME_QVAL,
    ". Using the top ", N_HEATMAP_GENES, " genes ranked by q-value."
  )
  heatmap_genes <- pseudotime_deg %>%
    filter(!is.na(qval)) %>%
    slice_head(n = N_HEATMAP_GENES) %>%
    pull(gene)
}

heatmap_genes <- intersect(heatmap_genes, rownames(sc_cds))

# ------------------------------------------------------------------------------
# 8. Pseudotime heatmap
# ------------------------------------------------------------------------------

heatmap_data <- plot_pseudotime_heatmap2(
  sc_cds[heatmap_genes, ],
  num_clusters = HEATMAP_CLUSTERS,
  cores = N_CORES
)

mark_genes <- c(
  "TNF", "IL1B", "CD80", "CCL3", "CCL4", "NFKB1",
  "CHI3L1", "EGR2", "CCL18", "APOE", "TREM2", "CYP27A1",
  "LPL", "APOC1", "SPP1", "ACOT2", "FABP5"
)
mark_genes <- intersect(mark_genes, heatmap_genes)

# ------------------------------------------------------------------------------
# 9. Functional enrichment of pseudotime gene clusters
# ------------------------------------------------------------------------------

enrich_types <- c("BP", "CC", "MF", "KEGG")

type_colors <- c(
  BP   = "#E41A1C",
  CC   = "#377EB8",
  MF   = "#4DAF4A",
  KEGG = "#984EA3"
)

enrich_list <- lapply(enrich_types, function(etype) {
  x <- enrichCluster(
    object = heatmap_data,
    OrgDb = org.Hs.eg.db,
    type = etype,
    pvalueCutoff = 0.05,
    topn = 4,
    seed = 5201314
  )

  if (is.null(x) || nrow(x) == 0) {
    return(NULL)
  }

  x$Type <- etype
  x
})

names(enrich_list) <- enrich_types
enrich_list <- Filter(Negate(is.null), enrich_list)

if (length(enrich_list) > 0) {
  enrich <- do.call(rbind, enrich_list)
  anno_term_data <- enrich[, 1:4, drop = FALSE]
  go_col <- unname(type_colors[as.character(enrich$Type)])

  write.csv(
    enrich,
    file.path(OUT_DIR, "monocle_cluster_enrichment.csv"),
    row.names = FALSE
  )
} else {
  enrich <- NULL
  anno_term_data <- NULL
  go_col <- NULL
}

# ------------------------------------------------------------------------------
# 10. Heatmap visualization
# ------------------------------------------------------------------------------

pdf(
  file.path(OUT_DIR, "monocle_heatmap.pdf"),
  height = 9,
  width = 16,
  onefile = FALSE
)

if (!is.null(anno_term_data)) {
  visCluster(
    object = heatmap_data,
    plot.type = "both",
    column_names_rot = 45,
    cluster_rows = FALSE,
    show_row_dend = FALSE,
    markGenes = mark_genes,
    markGenes.side = "left",
    annoTerm.data = anno_term_data,
    go.col = go_col,
    by.go = "anno_link",
    add.bar = TRUE,
    line.side = "left"
  )
} else {
  visCluster(
    object = heatmap_data,
    plot.type = "both",
    column_names_rot = 45,
    cluster_rows = FALSE,
    show_row_dend = FALSE,
    markGenes = mark_genes,
    markGenes.side = "left",
    add.bar = TRUE,
    line.side = "left"
  )
}

dev.off()

# ------------------------------------------------------------------------------
# 11. Save reproducible analysis objects
# ------------------------------------------------------------------------------

saveRDS(sc_cds, file.path(DATA_DIR, "monocle2_cds.rds"))
saveRDS(heatmap_data, file.path(DATA_DIR, "monocle2_heatmap_data.rds"))

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

message("Monocle 2 pseudotime analysis completed.")

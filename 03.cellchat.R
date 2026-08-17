#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(ggplot2)
})

set.seed(2025)

# ------------------------------------------------------------------------------
# 1. Configuration
# ------------------------------------------------------------------------------

INPUT_RDS <- "data/af.rds"
OUT_DIR <- "cellchat"

TUMOR_LABEL <- "Tumor"

SENDER_GROUPS <- c("High", "Medium", "Low")
RECEIVER_GROUPS <- c("C1", "C2", "C3", "C4")
SELECTED_GROUPS <- c(RECEIVER_GROUPS, SENDER_GROUPS)

MIN_CELLS <- 10
NBOOT <- 100
POPULATION_SIZE <- FALSE

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(INPUT_RDS)) {
  stop("Input file not found: ", INPUT_RDS)
}

# ------------------------------------------------------------------------------
# 2. Load and validate the Seurat object
# ------------------------------------------------------------------------------

af <- readRDS(INPUT_RDS)
DefaultAssay(af) <- "RNA"

required_meta <- c("Group", "Malignant_celltype", "orig.ident")
missing_meta <- setdiff(required_meta, colnames(af@meta.data))

if (length(missing_meta) > 0) {
  stop(
    "Missing required metadata column(s): ",
    paste(missing_meta, collapse = ", ")
  )
}

# Select tumor cells first.
tumor_cells <- colnames(af)[af$Group == TUMOR_LABEL]
scRNA <- subset(af, cells = tumor_cells)

# Use a neutral name for CellChat groups because the original metadata column
# contains both malignant epithelial subclusters and macrophage metabolic states.
scRNA$CellChat_group <- as.character(scRNA$Malignant_celltype)

# Retain only the predefined groups used in the analysis.
scRNA <- subset(scRNA, subset = CellChat_group %in% SELECTED_GROUPS)
scRNA$CellChat_group <- factor(scRNA$CellChat_group, levels = SELECTED_GROUPS)
scRNA$CellChat_group <- droplevels(scRNA$CellChat_group)

# ------------------------------------------------------------------------------
# 3. Quality checks before CellChat inference
# ------------------------------------------------------------------------------

group_counts <- as.data.frame(table(scRNA$CellChat_group))
colnames(group_counts) <- c("CellChat_group", "n_cells")
write.csv(
  group_counts,
  file.path(OUT_DIR, "cell_counts_by_group.csv"),
  row.names = FALSE
)

sample_group_counts <- as.data.frame(
  table(scRNA$orig.ident, scRNA$CellChat_group)
)
colnames(sample_group_counts) <- c(
  "sample_id", "CellChat_group", "n_cells"
)
write.csv(
  sample_group_counts,
  file.path(OUT_DIR, "cell_counts_by_sample_and_group.csv"),
  row.names = FALSE
)

missing_groups <- setdiff(SELECTED_GROUPS, unique(as.character(scRNA$CellChat_group)))
if (length(missing_groups) > 0) {
  stop(
    "The following predefined CellChat groups are absent: ",
    paste(missing_groups, collapse = ", ")
  )
}

small_groups <- group_counts$CellChat_group[group_counts$n_cells < MIN_CELLS]
if (length(small_groups) > 0) {
  stop(
    "The following groups contain fewer than ", MIN_CELLS,
    " cells and should not be analyzed: ",
    paste(small_groups, collapse = ", ")
  )
}

# ------------------------------------------------------------------------------
# 4. Create CellChat object
# ------------------------------------------------------------------------------

# CellChat requires normalized expression values for single-cell analysis.
data_input <- GetAssayData(
  scRNA,
  assay = "RNA",
  layer = "data"
)

meta_input <- scRNA@meta.data

cellchat <- createCellChat(
  object = data_input,
  meta = meta_input,
  group.by = "CellChat_group"
)

cellchat@DB <- CellChatDB.human

# ------------------------------------------------------------------------------
# 5. Infer ligand-receptor communication
# ------------------------------------------------------------------------------

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# raw.use = TRUE uses the normalized expression matrix rather than
# PPI-projected/smoothed expression.
cellchat <- computeCommunProb(
  cellchat,
  raw.use = TRUE,
  population.size = POPULATION_SIZE,
  nboot = NBOOT,
  seed.use = 2025
)

cellchat <- filterCommunication(
  cellchat,
  min.cells = MIN_CELLS
)

# Export ALL significant ligand-receptor interactions before any
# pathway-specific visualization.
lr_all <- subsetCommunication(cellchat)
write.csv(
  lr_all,
  file.path(OUT_DIR, "all_significant_LR_interactions.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 6. Pathway-level communication
# ------------------------------------------------------------------------------

cellchat <- computeCommunProbPathway(cellchat)

pathway_all <- subsetCommunication(
  cellchat,
  slot.name = "netP"
)

write.csv(
  pathway_all,
  file.path(OUT_DIR, "all_significant_pathways.csv"),
  row.names = FALSE
)

cellchat <- aggregateNet(cellchat)

# ------------------------------------------------------------------------------
# 7. Global network visualization
# ------------------------------------------------------------------------------

group_size <- as.numeric(table(cellchat@idents))

pdf(
  file.path(OUT_DIR, "global_interaction_count.pdf"),
  width = 8,
  height = 8
)
netVisual_circle(
  cellchat@net$count,
  vertex.weight = group_size,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Number of interactions"
)
dev.off()

pdf(
  file.path(OUT_DIR, "global_interaction_strength.pdf"),
  width = 8,
  height = 8
)
netVisual_circle(
  cellchat@net$weight,
  vertex.weight = group_size,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Interaction strength"
)
dev.off()

pdf(
  file.path(OUT_DIR, "global_heatmap_strength.pdf"),
  width = 8,
  height = 6
)
print(
  netVisual_heatmap(
    cellchat,
    measure = "weight",
    color.heatmap = "Reds"
  )
)
dev.off()

# ------------------------------------------------------------------------------
# 8. Prespecified macrophage -> malignant epithelial communication
# ------------------------------------------------------------------------------

# Export all significant interactions from metabolic-state macrophages
# to malignant epithelial subclusters.
lr_selected <- subsetCommunication(
  cellchat,
  sources.use = SENDER_GROUPS,
  targets.use = RECEIVER_GROUPS
)

write.csv(
  lr_selected,
  file.path(
    OUT_DIR,
    "macrophage_to_malignant_all_LR_interactions.csv"
  ),
  row.names = FALSE
)

pdf(
  file.path(OUT_DIR, "macrophage_to_malignant_strength.pdf"),
  width = 8,
  height = 8
)
netVisual_circle(
  cellchat@net$weight,
  vertex.weight = group_size,
  sources.use = SENDER_GROUPS,
  targets.use = RECEIVER_GROUPS,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Macrophage-to-malignant communication"
)
dev.off()

pdf(
  file.path(OUT_DIR, "macrophage_to_malignant_heatmap.pdf"),
  width = 8,
  height = 6
)
print(
  netVisual_heatmap(
    cellchat,
    measure = "weight",
    sources.use = SENDER_GROUPS,
    targets.use = RECEIVER_GROUPS,
    color.heatmap = "Reds"
  )
)
dev.off()

# Show all significant LR interactions between the predefined sender/receiver
# groups. No pathway is selected at this stage.
bubble_all <- netVisual_bubble(
  cellchat,
  sources.use = SENDER_GROUPS,
  targets.use = RECEIVER_GROUPS,
  remove.isolate = TRUE
)

ggsave(
  file.path(OUT_DIR, "macrophage_to_malignant_all_LR_bubble.pdf"),
  bubble_all,
  width = 8,
  height = 10
)

# ------------------------------------------------------------------------------
# 9. Prespecified SPP1 pathway analysis
# ------------------------------------------------------------------------------

# SPP1 is evaluated because it is a prespecified mechanistic pathway in the
# study, not because it was selected as the strongest pathway after CellChat.
if ("SPP1" %in% cellchat@netP$pathways) {

  spp1_lr <- subsetCommunication(
    cellchat,
    sources.use = SENDER_GROUPS,
    targets.use = RECEIVER_GROUPS,
    signaling = "SPP1"
  )

  write.csv(
    spp1_lr,
    file.path(OUT_DIR, "SPP1_macrophage_to_malignant_LR.csv"),
    row.names = FALSE
  )

  spp1_bubble <- netVisual_bubble(
    cellchat,
    sources.use = SENDER_GROUPS,
    targets.use = RECEIVER_GROUPS,
    signaling = "SPP1",
    remove.isolate = TRUE
  )

  ggsave(
    file.path(OUT_DIR, "SPP1_macrophage_to_malignant_bubble.pdf"),
    spp1_bubble,
    width = 7,
    height = 6
  )

  pdf(
    file.path(OUT_DIR, "SPP1_pathway_network.pdf"),
    width = 8,
    height = 7
  )
  netVisual_aggregate(
    cellchat,
    signaling = "SPP1",
    sources.use = SENDER_GROUPS,
    targets.use = RECEIVER_GROUPS,
    layout = "circle"
  )
  dev.off()

  # Export the contribution of all SPP1-associated LR pairs.
  spp1_pairs <- extractEnrichedLR(
    cellchat,
    signaling = "SPP1",
    geneLR.return = FALSE
  )

  write.csv(
    spp1_pairs,
    file.path(OUT_DIR, "SPP1_enriched_LR_pairs.csv"),
    row.names = FALSE
  )

} else {
  message("SPP1 was not inferred as a significant pathway in this dataset.")
}

# ------------------------------------------------------------------------------
# 10. Save reproducible objects and software information
# ------------------------------------------------------------------------------

saveRDS(
  cellchat,
  file.path(OUT_DIR, "cellchat_object.rds")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(OUT_DIR, "sessionInfo.txt")
)

message("CellChat analysis completed.")

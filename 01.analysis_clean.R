#!/usr/bin/env Rscript

set.seed(2025)
options(stringsAsFactors = FALSE)

# ----------------------------- 1. Configuration -------------------------------
PROJECT_DIR <- normalizePath(Sys.getenv("PROJECT_DIR", "."), mustWork = FALSE)
DATA_DIR    <- file.path(PROJECT_DIR, "data", "10x")
META_FILE   <- file.path(PROJECT_DIR, "metadata", "sample_information.tsv")
CELL_ANNO_FILE <- file.path(PROJECT_DIR, "metadata", "cell_annotation.tsv")
EPI_ANNO_FILE  <- file.path(PROJECT_DIR, "metadata", "epithelial_annotation.tsv")
GENE_ORDER_FILE <- file.path(PROJECT_DIR, "reference", "hg38_gencode_v27.txt")
CHOLESTEROL_GENESET_FILE <- file.path(PROJECT_DIR, "genesets", "cholesterol_genes.txt")
OUT_DIR <- file.path(PROJECT_DIR, "results")

N_CORES <- as.integer(Sys.getenv("N_CORES", "8"))
MAIN_NPCS <- 40
MAIN_RESOLUTION <- 0.5
EPI_NPCS <- 10
EPI_RESOLUTION <- 0.1

QC_MIN_COUNTS <- 1000
QC_MAX_FEATURES <- 6000
QC_MAX_MT <- 10

# inferCNV reference groups used in the study. Update only if the annotation
# table uses different group names.
INFERCNV_REFERENCE_GROUPS <- c(
  "Normal_0", "Normal_1", "Normal_2", "Normal_3", "Normal_4",
  "Normal_5", "Normal_6", "Normal_7", "Normal_8"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "01_global"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "02_epithelial"), recursive = TRUE, showWarnings = FALSE)

# ------------------------------ 2. Packages -----------------------------------
required_packages <- c(
  "Seurat", "dplyr", "tibble", "ggplot2", "patchwork", "harmony",
  "SingleR", "celldex", "infercnv", "readr", "irGSEA"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(harmony)
  library(SingleR)
  library(celldex)
  library(infercnv)
  library(readr)
  library(irGSEA)
})

# -------------------------- 3. Read and merge 10x data ------------------------
if (!dir.exists(DATA_DIR)) stop("10x data directory not found: ", DATA_DIR)

sample_dirs <- list.dirs(DATA_DIR, recursive = FALSE, full.names = TRUE)
if (length(sample_dirs) == 0) stop("No sample directories found in: ", DATA_DIR)

sample_ids <- basename(sample_dirs)
message("Reading ", length(sample_ids), " de-identified samples...")

sc_list <- lapply(seq_along(sample_dirs), function(i) {
  sample_id <- sample_ids[i]
  counts <- Read10X(data.dir = sample_dirs[i])

  obj <- CreateSeuratObject(
    counts = counts,
    project = sample_id,
    min.cells = 3,
    min.features = 200
  )
  obj$orig.ident <- sample_id
  RenameCells(obj, add.cell.id = sample_id)
})
names(sc_list) <- sample_ids

af <- merge(x = sc_list[[1]], y = sc_list[-1])
if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
  af <- JoinLayers(af)
}
rm(sc_list)

# ----------------------------- 4. Quality control -----------------------------
af[["percent.mt"]] <- PercentageFeatureSet(af, pattern = "^MT-")
af[["percent.rb"]] <- PercentageFeatureSet(af, pattern = "^RP[SL]")

af <- subset(
  af,
  subset = nCount_RNA >= QC_MIN_COUNTS &
    nFeature_RNA <= QC_MAX_FEATURES &
    percent.mt <= QC_MAX_MT
)

# Sample-level metadata must contain de-identified sample IDs only.
# Required columns: orig.ident, sample_type, Group, dataset
if (!file.exists(META_FILE)) stop("Metadata file not found: ", META_FILE)
sample_meta <- read.delim(META_FILE, check.names = FALSE)
required_meta_cols <- c("orig.ident", "sample_type", "Group", "dataset")
if (!all(required_meta_cols %in% colnames(sample_meta))) {
  stop("sample_information.tsv must contain: ", paste(required_meta_cols, collapse = ", "))
}

idx <- match(af$orig.ident, sample_meta$orig.ident)
if (anyNA(idx)) stop("Some sample IDs in the Seurat object are absent from sample_information.tsv")
af$sample_type <- sample_meta$sample_type[idx]
af$Group <- sample_meta$Group[idx]
af$dataset <- sample_meta$dataset[idx]

saveRDS(af, file.path(OUT_DIR, "01_global", "01_qc_filtered.rds"))

# Optional QC figure
pdf(file.path(OUT_DIR, "01_global", "01_qc_violin.pdf"), width = 12, height = 8)
print(VlnPlot(
  af,
  group.by = "orig.ident",
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"),
  ncol = 1,
  raster = TRUE
))
dev.off()

# ----------------------- 5. Normalization and clustering ----------------------
af <- NormalizeData(af, normalization.method = "LogNormalize", scale.factor = 10000)
af <- FindVariableFeatures(af, selection.method = "vst", nfeatures = 2000)
af <- ScaleData(af, vars.to.regress = "percent.mt")
af <- RunPCA(af, features = VariableFeatures(af))
af <- RunHarmony(af, group.by.vars = "orig.ident")
af <- FindNeighbors(af, reduction = "harmony", dims = 1:MAIN_NPCS)
af <- FindClusters(af, resolution = MAIN_RESOLUTION)
af <- RunUMAP(af, reduction = "harmony", dims = 1:MAIN_NPCS)

main_markers <- FindAllMarkers(
  af,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
write.csv(
  main_markers,
  file.path(OUT_DIR, "01_global", "02_cluster_markers.csv"),
  row.names = FALSE
)

pdf(file.path(OUT_DIR, "01_global", "03_cluster_umap.pdf"), width = 7, height = 6)
print(DimPlot(af, reduction = "umap", group.by = "seurat_clusters", label = TRUE, raster = TRUE))
dev.off()

# ------------------------- 6. Reference-based annotation ----------------------
hpca_ref <- celldex::HumanPrimaryCellAtlasData()

pred <- SingleR(
  test = GetAssayData(af, assay = "RNA", layer = "data"),
  ref = hpca_ref,
  labels = hpca_ref$label.fine,
  clusters = as.character(af$seurat_clusters)
)

singleR_annotation <- data.frame(
  seurat_clusters = rownames(pred),
  cell_type = pred$labels,
  stringsAsFactors = FALSE
)
write.table(
  singleR_annotation,
  file.path(OUT_DIR, "01_global", "04_SingleR_annotation.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Final manual annotation table used in the study.
# Required columns: seurat_clusters, cell_type
if (!file.exists(CELL_ANNO_FILE)) stop("Cell annotation file not found: ", CELL_ANNO_FILE)
cell_anno <- read.delim(CELL_ANNO_FILE, check.names = FALSE)
if (!all(c("seurat_clusters", "cell_type") %in% colnames(cell_anno))) {
  stop("cell_annotation.tsv must contain: seurat_clusters, cell_type")
}

annotation_map <- setNames(as.character(cell_anno$cell_type), as.character(cell_anno$seurat_clusters))
af$cell_type <- annotation_map[as.character(af$seurat_clusters)]
if (anyNA(af$cell_type)) warning("Some clusters are missing from cell_annotation.tsv")

DefaultAssay(af) <- "RNA"
Idents(af) <- "cell_type"

pdf(file.path(OUT_DIR, "01_global", "05_celltype_umap.pdf"), width = 8, height = 6)
print(DimPlot(af, reduction = "umap", group.by = "cell_type", label = TRUE, raster = TRUE))
dev.off()

saveRDS(af, file.path(OUT_DIR, "01_global", "06_annotated_all_cells.rds"))

# ----------------------- 7. Epithelial-cell subclustering ---------------------
Idents(af) <- "cell_type"
Epi <- subset(af, idents = "Epithelial cell")

Epi <- NormalizeData(Epi, normalization.method = "LogNormalize", scale.factor = 10000)
Epi <- FindVariableFeatures(Epi, selection.method = "vst", nfeatures = 2000)
Epi <- ScaleData(Epi, vars.to.regress = "percent.mt")
Epi <- RunPCA(Epi, features = VariableFeatures(Epi))
Epi <- RunHarmony(Epi, group.by.vars = "orig.ident")
Epi <- FindNeighbors(Epi, reduction = "harmony", dims = 1:EPI_NPCS)
Epi <- FindClusters(Epi, resolution = EPI_RESOLUTION)
Epi <- RunUMAP(Epi, reduction = "harmony", dims = 1:EPI_NPCS)

pdf(file.path(OUT_DIR, "02_epithelial", "01_epithelial_cluster_umap.pdf"), width = 7, height = 6)
print(DimPlot(Epi, reduction = "umap", group.by = "seurat_clusters", label = TRUE, raster = TRUE))
dev.off()

saveRDS(Epi, file.path(OUT_DIR, "02_epithelial", "02_epithelial_clustered.rds"))

# ------------------------------- 8. inferCNV ----------------------------------
# The gene-order file is expected to contain four tab-separated columns:
# gene, chromosome, start, end.
if (!file.exists(GENE_ORDER_FILE)) stop("Gene-order file not found: ", GENE_ORDER_FILE)

gene_order <- readr::read_tsv(
  GENE_ORDER_FILE,
  col_names = c("gene", "chr", "start", "end"),
  show_col_types = FALSE
)
gene_order <- gene_order[!duplicated(gene_order$gene), ]

Epi$Group_celltype <- paste0(Epi$Group, "_", Epi$seurat_clusters)
annotation_file <- file.path(OUT_DIR, "02_epithelial", "infercnv_cell_annotations.tsv")
write.table(
  Epi$Group_celltype,
  annotation_file,
  sep = "\t", quote = FALSE, col.names = FALSE
)

counts_mat <- GetAssayData(Epi, assay = "RNA", layer = "counts")
common_genes <- intersect(gene_order$gene, rownames(counts_mat))

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = counts_mat[common_genes, ],
  annotations_file = annotation_file,
  delim = "\t",
  gene_order_file = GENE_ORDER_FILE,
  ref_group_names = INFERCNV_REFERENCE_GROUPS,
  chr_exclude = c("chrX", "chrY", "chrM")
)

infercnv_out <- file.path(OUT_DIR, "02_epithelial", "infercnv_output")
infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,
  out_dir = infercnv_out,
  no_prelim_plot = TRUE,
  cluster_by_groups = TRUE,
  denoise = TRUE,
  HMM = FALSE,
  min_cells_per_gene = 10,
  num_threads = N_CORES,
  write_expr_matrix = TRUE,
  useRaster = TRUE,
  inspect_subclusters = FALSE
)
saveRDS(infercnv_obj, file.path(OUT_DIR, "02_epithelial", "03_infercnv_object.rds"))

# ----------------------------- 9. CNV score -----------------------------------
obs_file <- file.path(infercnv_out, "infercnv.observations.txt")
ref_file <- file.path(infercnv_out, "infercnv.references.txt")

if (file.exists(obs_file) && file.exists(ref_file)) {
  obs <- read.table(obs_file, header = TRUE, check.names = FALSE)
  ref <- read.table(ref_file, header = TRUE, check.names = FALSE)
  expr <- cbind(obs, ref)

  expr_scale <- scale(t(expr))
  shifted <- sweep(expr_scale, 2, apply(expr_scale, 2, min), "-")
  ranges <- apply(expr_scale, 2, max) - apply(expr_scale, 2, min)
  ranges[ranges == 0] <- 1
  expr_rescaled <- t(2 * sweep(shifted, 2, ranges, "/") - 1)

  cnv_score <- data.frame(
    cell = colnames(expr_rescaled),
    cnv_score = colSums(expr_rescaled^2),
    stringsAsFactors = FALSE
  )
  cnv_score$cell <- gsub("\\.", "-", cnv_score$cell)
  write.csv(
    cnv_score,
    file.path(OUT_DIR, "02_epithelial", "04_cnv_scores.csv"),
    row.names = FALSE
  )
}

# ---------------------- 10. Final epithelial annotation -----------------------
# Required columns: newcelltype, cell_type
if (!file.exists(EPI_ANNO_FILE)) stop("Epithelial annotation file not found: ", EPI_ANNO_FILE)
epi_anno <- read.delim(EPI_ANNO_FILE, check.names = FALSE)
if (!all(c("newcelltype", "cell_type") %in% colnames(epi_anno))) {
  stop("epithelial_annotation.tsv must contain: newcelltype, cell_type")
}

epi_map <- setNames(as.character(epi_anno$cell_type), as.character(epi_anno$newcelltype))
Epi$anno <- epi_map[as.character(Epi$Group_celltype)]
if (anyNA(Epi$anno)) warning("Some epithelial groups are missing from epithelial_annotation.tsv")

# Transfer refined epithelial labels back to the complete object.
af$Epi_celltype <- as.character(af$cell_type)
af@meta.data[colnames(Epi), "Epi_celltype"] <- as.character(Epi$anno)

saveRDS(Epi, file.path(OUT_DIR, "02_epithelial", "05_epithelial_annotated.rds"))
saveRDS(af, file.path(OUT_DIR, "01_global", "07_all_cells_with_epithelial_labels.rds"))

# ---------------------- 11. Multi-algorithm gene-set scoring ------------------
# This reproduces the study's scoring framework:
# AddModuleScore + singscore + AUCell + UCell, followed by direct summation.

read_geneset <- function(file, object) {
  if (!file.exists(file)) stop("Gene-set file not found: ", file)
  genes <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)[, 1]
  genes <- unique(as.character(genes))
  genes <- genes[genes %in% rownames(object)]
  if (length(genes) == 0) stop("No genes from the gene set were found in the Seurat object")
  list(genes)
}

score_gene_set <- function(object, geneset_file, prefix, ncores = 8, seed = 2025) {
  geneset <- read_geneset(geneset_file, object)

  # Remove temporary irGSEA assays from a previous scoring run, if present.
  for (assay_name in c("singscore", "AUCell", "UCell")) {
    if (assay_name %in% Assays(object)) object[[assay_name]] <- NULL
  }

  module_name <- paste0(prefix, "_AddModuleScore")
  object <- AddModuleScore(
    object,
    features = geneset,
    name = module_name,
    seed = seed
  )
  module_col <- paste0(module_name, "1")

  object <- irGSEA.score(
    object = object,
    assay = "RNA",
    slot = "data",
    seeds = seed,
    ncores = ncores,
    min.cells = 3,
    min.feature = 0,
    custom = TRUE,
    geneset = geneset,
    msigdb = FALSE,
    species = "Homo sapiens",
    category = "H",
    subcategory = NULL,
    geneid = "symbol",
    method = c("singscore", "AUCell", "UCell"),
    aucell.MaxRank = NULL,
    ucell.MaxRank = NULL,
    kcdf = "Gaussian"
  )

  object[[paste0(prefix, "_singscore")]] <- as.numeric(object[["singscore"]]@data[1, colnames(object)])
  object[[paste0(prefix, "_AUCell")]] <- as.numeric(object[["AUCell"]]@data[1, colnames(object)])
  object[[paste0(prefix, "_UCell")]] <- as.numeric(object[["UCell"]]@data[1, colnames(object)])

  object[[paste0(prefix, "_Scoring")]] <-
    object[[paste0(prefix, "_singscore")]][, 1] +
    object[[paste0(prefix, "_AUCell")]][, 1] +
    object[[paste0(prefix, "_UCell")]][, 1] +
    object[[module_col]][, 1]

  # Keep the cell-level scores in metadata; remove temporary assays to avoid
  # collisions when additional gene sets are scored.
  for (assay_name in c("singscore", "AUCell", "UCell")) {
    if (assay_name %in% Assays(object)) object[[assay_name]] <- NULL
  }

  object
}

af <- score_gene_set(
  object = af,
  geneset_file = CHOLESTEROL_GENESET_FILE,
  prefix = "Cholesterol",
  ncores = N_CORES,
  seed = 2025
)

# Example for additional gene sets used in sensitivity analyses:
# af <- score_gene_set(af, file.path(PROJECT_DIR, "genesets", "lactate_genes.txt"), "Lactate", N_CORES)
# af <- score_gene_set(af, file.path(PROJECT_DIR, "genesets", "HALLMARK_GLYCOLYSIS.txt"), "Glycolysis", N_CORES)
# af <- score_gene_set(af, file.path(PROJECT_DIR, "genesets", "HALLMARK_HYPOXIA.txt"), "Hypoxia", N_CORES)

saveRDS(af, file.path(OUT_DIR, "01_global", "08_final_scRNA_object.rds"))

# --------------------------- 12. Reproducibility ------------------------------
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
message("Analysis completed. Results written to: ", OUT_DIR)

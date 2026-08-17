#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(irGSEA)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

# ------------------------------------------------------------------------------
# 1. Configuration
# ------------------------------------------------------------------------------

DATA_ROOT <- "data/spatial"
GENESET_DIR <- "gene_sets"
OUT_ROOT <- "results04"

SAMPLES <- c(
  "Colon1", "Colon2", "Colon3", "Colon4",
  "Liver1", "Liver2"
)

GENESET_FILES <- c(
  Lactate     = file.path(GENESET_DIR, "lactate_metabolism_genes.txt"),
  Cholesterol = file.path(GENESET_DIR, "cholesterol_metabolism_genes.txt"),
  Macrophage  = file.path(GENESET_DIR, "macrophage_signature_genes.txt")
)

MIN_FEATURES <- 200
SCALE_FACTOR <- 10000
HIGH_QUANTILE <- 0.75
N_CORES <- as.integer(Sys.getenv("N_CORES", "8"))

dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------------------------

read_gene_set <- function(file, features) {
  if (!file.exists(file)) {
    stop("Gene-set file not found: ", file)
  }

  genes <- read.table(
    file,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE
  )[[1]]

  genes <- unique(trimws(as.character(genes)))
  genes <- genes[nzchar(genes)]
  genes <- intersect(genes, features)

  if (length(genes) < 5) {
    stop("Too few genes from ", basename(file), " are present in the dataset.")
  }

  genes
}


extract_irGSEA_score <- function(object, assay_name) {
  mat <- object[[assay_name]]@data

  if (nrow(mat) != 1) {
    stop(
      "Expected one custom gene-set score in assay '",
      assay_name, "', but found ", nrow(mat), " rows."
    )
  }

  as.numeric(mat[1, colnames(object), drop = TRUE])
}



get_xy <- function(object) {
  coord <- GetTissueCoordinates(object)

  coord$barcode <- if ("barcode" %in% colnames(coord)) {
    coord$barcode
  } else if ("cell" %in% colnames(coord)) {
    coord$cell
  } else {
    rownames(coord)
  }

  if (all(c("x", "y") %in% colnames(coord))) {
    coord$x_plot <- coord$x
    coord$y_plot <- coord$y
  } else if (all(c("imagecol", "imagerow") %in% colnames(coord))) {
    coord$x_plot <- coord$imagecol
    coord$y_plot <- coord$imagerow
  } else if (all(c("col", "row") %in% colnames(coord))) {
    coord$x_plot <- coord$col
    coord$y_plot <- coord$row
  } else {
    stop(
      "Spatial coordinates were not recognized. Available columns: ",
      paste(colnames(coord), collapse = ", ")
    )
  }

  coord[, c("barcode", "x_plot", "y_plot")]
}


# ------------------------------------------------------------------------------
# 3. Gene-set scoring
# ------------------------------------------------------------------------------

score_one_gene_set <- function(object, genes, prefix) {

  # AddModuleScore is calculated from the normalized Spatial assay.
  object <- AddModuleScore(
    object = object,
    features = list(genes),
    name = paste0(prefix, "_AddModuleScore"),
    assay = "Spatial"
  )

  add_col <- paste0(prefix, "_AddModuleScore1")

  # irGSEA is run separately for each signature to avoid ambiguity when
  # extracting the resulting score matrices.
  object <- irGSEA.score(
    object = object,
    assay = "Spatial",
    slot = "data",
    seeds = 2025,
    ncores = N_CORES,
    min.cells = 3,
    min.feature = 0,
    custom = TRUE,
    geneset = list(genes),
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

  singscore_col <- paste0(prefix, "_singscore")
  aucell_col <- paste0(prefix, "_AUCell")
  ucell_col <- paste0(prefix, "_UCell")

  object[[singscore_col]] <- extract_irGSEA_score(object, "singscore")
  object[[aucell_col]] <- extract_irGSEA_score(object, "AUCell")
  object[[ucell_col]] <- extract_irGSEA_score(object, "UCell")

  # Primary composite score used in the study:
  # each algorithm-derived score enters with coefficient 1.
  composite_col <- paste0(prefix, "_Scoring")

  object[[composite_col]] <-
    object[[add_col]][, 1] +
    object[[singscore_col]][, 1] +
    object[[aucell_col]][, 1] +
    object[[ucell_col]][, 1]


  # Remove temporary irGSEA assays after scores have been safely copied to
  # metadata so the next gene set is processed independently.
  for (a in c("singscore", "AUCell", "UCell")) {
    if (a %in% Assays(object)) {
      object[[a]] <- NULL
    }
  }

  object
}


# ------------------------------------------------------------------------------
# 4. Continuous spatial-score plots
# ------------------------------------------------------------------------------

plot_continuous_scores <- function(object, sample_id, out_dir) {

  features <- c(
    "Lactate_Scoring",
    "Cholesterol_Scoring",
    "Macrophage_Scoring"
  )

  p <- SpatialFeaturePlot(
    object,
    features = features,
    assay = "Spatial",
    ncol = 3,
    pt.size.factor = 1.4
  )

  ggsave(
    file.path(out_dir, "continuous_signature_scores.pdf"),
    p,
    width = 12,
    height = 4
  )
}


# ------------------------------------------------------------------------------
# 5. Quantile-based overlap map
# ------------------------------------------------------------------------------

plot_quantile_overlap <- function(object, sample_id, out_dir) {

  score_cols <- c(
    Lactate = "Lactate_Scoring",
    Cholesterol = "Cholesterol_Scoring",
    Macrophage = "Macrophage_Scoring"
  )

  coord <- get_xy(object)

  meta <- object@meta.data %>%
    rownames_to_column("barcode")

  plot_df <- left_join(coord, meta, by = "barcode")

  thresholds <- sapply(
    score_cols,
    function(x) {
      quantile(
        plot_df[[x]],
        probs = HIGH_QUANTILE,
        na.rm = TRUE,
        names = FALSE
      )
    }
  )

  plot_df <- plot_df %>%
    mutate(
      Lactate_high =
        .data[[score_cols["Lactate"]]] >= thresholds["Lactate"],
      Cholesterol_high =
        .data[[score_cols["Cholesterol"]]] >= thresholds["Cholesterol"],
      Macrophage_high =
        .data[[score_cols["Macrophage"]]] >= thresholds["Macrophage"]
    )

  plot_df$score_group <- apply(
    plot_df[, c("Lactate_high", "Cholesterol_high", "Macrophage_high")],
    1,
    function(x) {
      hit <- c(
        if (x[1]) "Lactate",
        if (x[2]) "Cholesterol",
        if (x[3]) "Macrophage"
      )

      hit <- hit[!is.na(hit)]

      if (length(hit) == 0) "None" else paste(hit, collapse = "+")
    }
  )

  group_levels <- c(
    "None",
    "Lactate",
    "Cholesterol",
    "Macrophage",
    "Lactate+Cholesterol",
    "Lactate+Macrophage",
    "Cholesterol+Macrophage",
    "Lactate+Cholesterol+Macrophage"
  )

  plot_df$score_group <- factor(
    plot_df$score_group,
    levels = group_levels
  )

  overlap_colors <- c(
    "None" = "grey85",
    "Lactate" = "#E41A1C",
    "Cholesterol" = "#4DAF4A",
    "Macrophage" = "#377EB8",
    "Lactate+Cholesterol" = "#FFBF00",
    "Lactate+Macrophage" = "#984EA3",
    "Cholesterol+Macrophage" = "#00A9A5",
    "Lactate+Cholesterol+Macrophage" = "#222222"
  )

  p <- ggplot(
    plot_df,
    aes(x = x_plot, y = y_plot, color = score_group)
  ) +
    geom_point(size = 0.8, alpha = 0.9) +
    scale_color_manual(
      values = overlap_colors,
      drop = FALSE,
      name = NULL
    ) +
    scale_y_reverse() +
    coord_fixed() +
    labs(title = sample_id, x = NULL, y = NULL) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right"
    )

  ggsave(
    file.path(out_dir, "quantile_overlap_map.pdf"),
    p,
    width = 7,
    height = 6
  )

  write.csv(
    plot_df[, c(
      "barcode",
      "Lactate_Scoring",
      "Cholesterol_Scoring",
      "Macrophage_Scoring",
      "Lactate_high",
      "Cholesterol_high",
      "Macrophage_high",
      "score_group"
    )],
    file.path(out_dir, "spot_signature_classification.csv"),
    row.names = FALSE
  )

  data.frame(
    sample_id = sample_id,
    signature = names(thresholds),
    quantile = HIGH_QUANTILE,
    threshold = as.numeric(thresholds),
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------------------------
# 6. Per-sample analysis
# ------------------------------------------------------------------------------

analyze_spatial_sample <- function(sample_id) {

  message("Processing: ", sample_id)

  input_dir <- file.path(DATA_ROOT, sample_id)
  out_dir <- file.path(OUT_ROOT, sample_id)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (!dir.exists(input_dir)) {
    stop("Spatial input directory not found: ", input_dir)
  }

  # Standard Seurat loader. The public data directory should follow the
  # standard 10x Visium layout.
  object <- Load10X_Spatial(
    data.dir = input_dir,
    assay = "Spatial",
    slice = sample_id
  )

  # Transparent minimal QC criterion retained from the original workflow.
  object <- subset(
    object,
    subset = nFeature_Spatial > MIN_FEATURES
  )

  # Use conventional log-normalized Spatial expression for gene-set scoring.
  # No SCT assay is copied or relabeled as RNA.
  object <- NormalizeData(
    object,
    assay = "Spatial",
    normalization.method = "LogNormalize",
    scale.factor = SCALE_FACTOR,
    verbose = FALSE
  )

  genes_used <- list()

  for (prefix in names(GENESET_FILES)) {

    genes <- read_gene_set(
      GENESET_FILES[[prefix]],
      features = rownames(object)
    )

    genes_used[[prefix]] <- genes

    object <- score_one_gene_set(
      object,
      genes = genes,
      prefix = prefix
    )
  }

  plot_continuous_scores(
    object,
    sample_id = sample_id,
    out_dir = out_dir
  )

  threshold_df <- plot_quantile_overlap(
    object,
    sample_id = sample_id,
    out_dir = out_dir
  )

  # Pairwise Spearman correlations are exported numerically only.
  # They are not used to define high-score spots.
  score_df <- object@meta.data[, c(
    "Lactate_Scoring",
    "Cholesterol_Scoring",
    "Macrophage_Scoring"
  )]

  correlation_matrix <- cor(
    score_df,
    method = "spearman",
    use = "pairwise.complete.obs"
  )

  correlation_long <- as.data.frame(as.table(correlation_matrix))
  colnames(correlation_long) <- c(
    "signature_1", "signature_2", "spearman_rho"
  )
  correlation_long$sample_id <- sample_id

  genes_used_df <- bind_rows(
    lapply(
      names(genes_used),
      function(prefix) {
        data.frame(
          sample_id = sample_id,
          signature = prefix,
          gene = genes_used[[prefix]],
          stringsAsFactors = FALSE
        )
      }
    )
  )

  saveRDS(
    object,
    file.path(out_dir, paste0(sample_id, "_spatial_scored.rds"))
  )

  list(
    thresholds = threshold_df,
    correlations = correlation_long,
    genes_used = genes_used_df
  )
}


# ------------------------------------------------------------------------------
# 7. Run all spatial sections
# ------------------------------------------------------------------------------

results <- lapply(SAMPLES, analyze_spatial_sample)
names(results) <- SAMPLES

thresholds_all <- bind_rows(lapply(results, `[[`, "thresholds"))
correlations_all <- bind_rows(lapply(results, `[[`, "correlations"))
genes_used_all <- bind_rows(lapply(results, `[[`, "genes_used"))

write.csv(
  thresholds_all,
  file.path(OUT_ROOT, "spatial_score_thresholds.csv"),
  row.names = FALSE
)

write.csv(
  correlations_all,
  file.path(OUT_ROOT, "spatial_score_correlations.csv"),
  row.names = FALSE
)

write.csv(
  genes_used_all,
  file.path(OUT_ROOT, "genes_used_for_scoring.csv"),
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(OUT_ROOT, "sessionInfo.txt")
)

message("Spatial transcriptomics analysis completed.")

# ==============================================================================
# Dimensionality Reduction — bulkRNAseq Pipeline
# PCA and UMAP on VST-normalized bulk RNA-seq data
# ==============================================================================

#' Run PCA on VST-normalized expression matrix
#'
#' @param vst_obj   SummarizedExperiment from normalize_vst()
#' @param n_top     Number of most variable genes to use
#' @param n_pcs     Number of PCs to retain
#' @returns List with: pca (prcomp object), pca_df (scores + metadata),
#'   var_explained (numeric vector of % variance per PC)
run_pca <- function(vst_obj, n_top = N_TOP_GENES_PCA, n_pcs = N_PCS) {
  mat <- SummarizedExperiment::assay(vst_obj)
  rv  <- apply(mat, 1, var)
  top_genes <- order(rv, decreasing = TRUE)[seq_len(min(n_top, nrow(mat)))]
  mat_top   <- t(mat[top_genes, ])   # samples x genes for prcomp

  message(sprintf("[PCA] Running on %d genes x %d samples...", ncol(mat_top), nrow(mat_top)))
  pca <- prcomp(mat_top, center = TRUE, scale. = FALSE)

  n_pcs       <- min(n_pcs, ncol(pca$x))
  var_exp     <- (pca$sdev^2 / sum(pca$sdev^2) * 100)[seq_len(n_pcs)]
  names(var_exp) <- paste0("PC", seq_len(n_pcs))

  pca_df <- as.data.frame(pca$x[, seq_len(n_pcs), drop = FALSE])
  pca_df$sample <- rownames(pca_df)
  pca_df <- merge(pca_df, as.data.frame(SummarizedExperiment::colData(vst_obj)),
                  by.x = "sample", by.y = "row.names", sort = FALSE)

  message(sprintf("[PCA] PC1: %.1f%% | PC2: %.1f%% | PC3: %.1f%% variance",
                  var_exp[1], var_exp[2], var_exp[3]))
  list(pca = pca, pca_df = pca_df, var_explained = var_exp)
}

#' Scree plot: variance explained per principal component
#'
#' @param var_explained Named numeric vector from run_pca()$var_explained
#' @param n_show        Number of PCs to display
plot_scree <- function(var_explained, n_show = 20) {
  df <- data.frame(
    PC  = factor(names(var_explained)[seq_len(n_show)],
                 levels = names(var_explained)[seq_len(n_show)]),
    var = var_explained[seq_len(n_show)]
  )
  ggplot(df, aes(x = PC, y = var)) +
    geom_col(fill = "#377EB8", alpha = 0.8) +
    geom_line(aes(group = 1), colour = "black") +
    geom_point(size = 2) +
    labs(title = "PCA: variance explained", x = "Principal component",
         y = "% variance explained") +
    BASE_THEME +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#' PCA scatter plot coloured by a metadata variable
#'
#' @param pca_df    data.frame from run_pca()$pca_df
#' @param var_explained Named numeric vector from run_pca()$var_explained
#' @param colour_by Column in pca_df to colour points by
#' @param shape_by  Column in pca_df to map to point shape (optional)
#' @param pc_x      PC index for x-axis
#' @param pc_y      PC index for y-axis
#' @param palette   Named character vector of colours (optional)
plot_pca <- function(pca_df, var_explained,
                     colour_by = "genotype", shape_by = "tissue",
                     pc_x = 1, pc_y = 2, palette = NULL) {
  pc_x_name <- paste0("PC", pc_x)
  pc_y_name <- paste0("PC", pc_y)

  p <- ggplot(pca_df, aes(
    x      = .data[[pc_x_name]],
    y      = .data[[pc_y_name]],
    colour = .data[[colour_by]],
    shape  = if (!is.null(shape_by)) .data[[shape_by]] else NULL
  )) +
    geom_point(size = 3, alpha = 0.85) +
    labs(
      title  = sprintf("PCA — coloured by %s", colour_by),
      x      = sprintf("%s (%.1f%%)", pc_x_name, var_explained[pc_x]),
      y      = sprintf("%s (%.1f%%)", pc_y_name, var_explained[pc_y]),
      colour = colour_by,
      shape  = shape_by
    ) +
    BASE_THEME

  if (!is.null(palette)) p <- p + scale_colour_manual(values = palette)
  p
}

#' Run UMAP on the top N principal components
#'
#' @param pca_result  Output of run_pca()
#' @param n_pcs       Number of PCs to use as input
#' @param seed        Random seed for reproducibility
#' @returns data.frame with UMAP1, UMAP2, and all metadata columns
run_umap <- function(pca_result,
                     n_pcs       = N_PCS_UMAP,
                     n_neighbors = UMAP_N_NEIGHBORS,
                     min_dist    = UMAP_MIN_DIST,
                     seed        = UMAP_SEED) {
  pc_matrix <- as.matrix(pca_result$pca_df[, paste0("PC", seq_len(n_pcs))])

  message(sprintf("[UMAP] Running on %d PCs x %d samples...", n_pcs, nrow(pc_matrix)))
  set.seed(seed)
  umap_coords <- uwot::umap(
    pc_matrix,
    n_neighbors = n_neighbors,
    min_dist    = min_dist,
    metric      = "euclidean",
    verbose     = FALSE
  )

  umap_df <- pca_result$pca_df
  umap_df$UMAP1 <- umap_coords[, 1]
  umap_df$UMAP2 <- umap_coords[, 2]
  umap_df
}

#' UMAP scatter plot coloured by a metadata variable
#'
#' @param umap_df   data.frame from run_umap()
#' @param colour_by Column to colour by
#' @param shape_by  Column to map to shape (optional)
#' @param palette   Named colour vector (optional)
plot_umap <- function(umap_df, colour_by = "genotype", shape_by = "tissue",
                      palette = NULL) {
  p <- ggplot(umap_df, aes(
    x      = UMAP1,
    y      = UMAP2,
    colour = .data[[colour_by]],
    shape  = if (!is.null(shape_by)) .data[[shape_by]] else NULL
  )) +
    geom_point(size = 3, alpha = 0.85) +
    labs(
      title  = sprintf("UMAP — coloured by %s", colour_by),
      x      = "UMAP 1", y = "UMAP 2",
      colour = colour_by, shape = shape_by
    ) +
    BASE_THEME

  if (!is.null(palette)) p <- p + scale_colour_manual(values = palette)
  p
}

#' Build a grid of PCA/UMAP plots coloured by each key experimental variable
#'
#' @param embed_df   data.frame with embedding coordinates and metadata
#' @param x_col     Column name for x-axis ("PC1" or "UMAP1")
#' @param y_col     Column name for y-axis ("PC2" or "UMAP2")
#' @param var_explained Named numeric vector (for PCA axis labels); NULL for UMAP
plot_embedding_grid <- function(embed_df, x_col = "PC1", y_col = "PC2",
                                var_explained = NULL) {
  make_panel <- function(colour_by, palette = NULL) {
    x_lab <- if (!is.null(var_explained) && x_col %in% names(var_explained))
      sprintf("%s (%.1f%%)", x_col, var_explained[x_col]) else x_col
    y_lab <- if (!is.null(var_explained) && y_col %in% names(var_explained))
      sprintf("%s (%.1f%%)", y_col, var_explained[y_col]) else y_col

    p <- ggplot(embed_df, aes(x = .data[[x_col]], y = .data[[y_col]],
                               colour = .data[[colour_by]])) +
      geom_point(size = 2.5, alpha = 0.8) +
      labs(title = colour_by, x = x_lab, y = y_lab, colour = colour_by) +
      BASE_THEME +
      theme(legend.position = "bottom")
    if (!is.null(palette)) p <- p + scale_colour_manual(values = palette)
    p
  }

  p_genotype   <- make_panel("genotype",   GENOTYPE_COLORS)
  p_tissue     <- make_panel("tissue",     TISSUE_COLORS)
  p_sex        <- make_panel("sex",        SEX_COLORS)
  p_age        <- make_panel("age_months",
    setNames(AGE_COLORS, names(AGE_COLORS)))

  (p_genotype | p_tissue) / (p_sex | p_age)
}

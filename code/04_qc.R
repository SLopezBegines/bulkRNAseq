# QC Functions — bulkRNAseq Pipeline
# Library size, gene detection, sample correlation, and summary reporting

#' Compute per-sample QC metrics from a raw count matrix
#'
#' @param counts   Gene x sample count matrix
#' @param metadata Sample metadata data.frame
#' @returns data.frame with QC metrics joined to metadata
compute_qc_metrics <- function(counts, metadata) {
  qc <- data.frame(
    sample           = colnames(counts),
    library_size     = colSums(counts),
    genes_detected   = colSums(counts > 0),
    median_count     = apply(counts, 2, median),
    pct_top50_genes  = apply(counts, 2, function(x) {
      s  <- sort(x, decreasing = TRUE)
      sum(s[seq_len(min(50, length(s)))]) / sum(s) * 100
    }),
    stringsAsFactors = FALSE
  )
  # Flag samples below QC thresholds
  qc$low_depth    <- qc$library_size    < QC_MIN_LIBRARY_SIZE
  qc$low_coverage <- qc$genes_detected  < QC_MIN_GENES_DETECTED

  merge(qc, metadata, by = "sample", sort = FALSE)
}

#' Barplot of library sizes, coloured by genotype, flagging low-depth samples
#'
#' @param qc_df    Output of compute_qc_metrics()
#' @param min_size Horizontal line marking the minimum acceptable library size
plot_library_size <- function(qc_df, min_size = QC_MIN_LIBRARY_SIZE) {
  qc_df <- qc_df[order(qc_df$library_size), ]
  qc_df$sample_ord <- factor(qc_df$sample, levels = qc_df$sample)

  ggplot(qc_df, aes(x = sample_ord, y = library_size / 1e6,
                    fill = genotype, alpha = ifelse(low_depth, 0.4, 1))) +
    geom_col() +
    geom_hline(yintercept = min_size / 1e6, linetype = "dashed", colour = "red") +
    scale_fill_manual(values = GENOTYPE_COLORS) +
    scale_alpha_identity() +
    facet_wrap(~ tissue + age_months, scales = "free_x", ncol = 4) +
    labs(
      title = "Library size per sample",
      x     = NULL, y = "Total counts (millions)",
      fill  = "Genotype"
    ) +
    BASE_THEME +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

#' Boxplot of genes detected per group
plot_genes_detected <- function(qc_df) {
  ggplot(qc_df, aes(x = interaction(genotype, tissue), y = genes_detected,
                    fill = genotype, colour = genotype)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_jitter(aes(shape = sex), width = 0.15, size = 1.5) +
    geom_hline(yintercept = QC_MIN_GENES_DETECTED, linetype = "dashed", colour = "red") +
    scale_fill_manual(values   = GENOTYPE_COLORS) +
    scale_colour_manual(values = GENOTYPE_COLORS) +
    facet_wrap(~ age_months, nrow = 1, labeller = label_both) +
    labs(
      title = "Genes detected per sample",
      x     = "Genotype × Tissue", y = "Genes detected (count > 0)",
      fill  = "Genotype", colour = "Genotype", shape = "Sex"
    ) +
    BASE_THEME +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

#' Density plot of log10(library size) to spot outlier samples
plot_library_distribution <- function(qc_df) {
  ggplot(qc_df, aes(x = log10(library_size), colour = genotype,
                    linetype = tissue, group = interaction(genotype, tissue))) +
    geom_density(size = 0.8) +
    geom_vline(xintercept = log10(QC_MIN_LIBRARY_SIZE),
               linetype = "dashed", colour = "red") +
    scale_colour_manual(values = GENOTYPE_COLORS) +
    facet_wrap(~ age_months, nrow = 1, labeller = label_both) +
    labs(
      title   = "Distribution of library sizes",
      x       = expression(log[10]~"(total counts)"),
      y       = "Density",
      colour  = "Genotype",
      linetype = "Tissue"
    ) +
    BASE_THEME
}

#' Heatmap of Pearson correlation between samples on log2-normalized counts
#'
#' @param counts     Gene x sample count matrix (raw or normalized)
#' @param metadata   Sample metadata
#' @param top_n      Use only top_n most variable genes
#' @param log_transform Apply log2(x + 1) before computing correlation?
plot_sample_correlation <- function(counts, metadata,
                                    top_n = 1000, log_transform = TRUE) {
  mat <- counts
  if (log_transform) mat <- log2(mat + 1)
  # Select most variable genes
  rv  <- apply(mat, 1, var)
  mat <- mat[order(rv, decreasing = TRUE)[seq_len(min(top_n, nrow(mat)))], ]

  cor_mat <- cor(mat, method = "pearson")

  # Build annotation for the heatmap
  ann_df <- metadata[, c("genotype", "tissue", "age_months", "sex"), drop = FALSE]
  ann_df$age_months <- as.character(ann_df$age_months)
  ann_colors <- list(
    genotype   = GENOTYPE_COLORS,
    tissue     = TISSUE_COLORS,
    sex        = SEX_COLORS,
    age_months = AGE_COLORS
  )

  pheatmap::pheatmap(
    cor_mat,
    annotation_col  = ann_df,
    annotation_row  = ann_df,
    annotation_colors = ann_colors,
    show_rownames   = FALSE,
    show_colnames   = FALSE,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method        = "ward.D2",
    color  = colorRampPalette(c("#053061","white","#67001F"))(100),
    breaks = seq(0.7, 1.0, length.out = 101),
    main   = sprintf("Sample Pearson correlation (top %d variable genes)", top_n),
    border_color = NA,
    fontsize     = 6
  )
}

#' Print and return a text summary of QC flags
qc_summary <- function(qc_df) {
  n_total      <- nrow(qc_df)
  n_low_depth  <- sum(qc_df$low_depth)
  n_low_cov    <- sum(qc_df$low_coverage)

  cat(sprintf("=== QC Summary (%d samples) ===\n", n_total))
  cat(sprintf("  Low library size (< %.0fM): %d samples\n",
              QC_MIN_LIBRARY_SIZE / 1e6, n_low_depth))
  if (n_low_depth > 0) {
    cat("    Flagged:", paste(qc_df$sample[qc_df$low_depth], collapse = "\n             "), "\n")
  }
  cat(sprintf("  Low gene detection (< %d genes): %d samples\n",
              QC_MIN_GENES_DETECTED, n_low_cov))
  if (n_low_cov > 0) {
    cat("    Flagged:", paste(qc_df$sample[qc_df$low_coverage], collapse = "\n             "), "\n")
  }
  invisible(qc_df[qc_df$low_depth | qc_df$low_coverage, ])
}

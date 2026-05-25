# ==============================================================================
# Normalization & Filtering — bulkRNAseq Pipeline
# ==============================================================================

#' Pre-filter low-count genes before DE analysis
#'
#' Removes genes that do not meet a minimum expression threshold in a minimum
#' number of samples. This is the standard pre-filtering step recommended by
#' the DESeq2 vignette and edgeR user guide.
#'
#' @param counts          Gene x sample raw count matrix
#' @param min_cpm         Minimum CPM threshold (default MIN_COUNT_CPM)
#' @param min_samples     Minimum number of samples that must meet min_cpm
#' @param min_total       Minimum total counts across all samples
#' @returns Filtered count matrix
filter_low_counts <- function(counts,
                               min_cpm     = MIN_COUNT_CPM,
                               min_samples = MIN_SAMPLES_EXPRESSED,
                               min_total   = MIN_TOTAL_COUNTS) {
  lib_sizes  <- colSums(counts)
  cpm_matrix <- t(t(counts) / lib_sizes * 1e6)

  keep_cpm   <- rowSums(cpm_matrix >= min_cpm) >= min_samples
  keep_total <- rowSums(counts) >= min_total
  keep       <- keep_cpm & keep_total

  n_before <- nrow(counts)
  counts   <- counts[keep, , drop = FALSE]
  message(sprintf("[FILTER] Genes: %d -> %d (removed %d low-count genes)",
                  n_before, nrow(counts), n_before - nrow(counts)))
  counts
}

#' Apply edgeR's filterByExpr as an alternative gene-level filter
#'
#' @param counts    Gene x sample count matrix
#' @param metadata  Sample metadata with a 'group' column
#' @returns Filtered count matrix
filter_by_expr_edger <- function(counts, metadata) {
  group  <- metadata$group
  keep   <- edgeR::filterByExpr(counts, group = group)
  n_before <- nrow(counts)
  counts <- counts[keep, , drop = FALSE]
  message(sprintf("[FILTER edgeR] Genes: %d -> %d (removed %d)",
                  n_before, nrow(counts), n_before - nrow(counts)))
  counts
}

#' Normalize raw counts to CPM (Counts Per Million)
#' Simple normalization for visualization; not for DE analysis.
#'
#' @param counts Gene x sample count matrix
#' @returns CPM matrix (same dimensions)
counts_to_cpm <- function(counts) {
  lib_sizes <- colSums(counts)
  t(t(counts) / lib_sizes * 1e6)
}

#' Apply variance-stabilizing transformation (VST) from DESeq2
#' Recommended for PCA/UMAP visualization with n > 30 samples.
#'
#' @param dds     DESeqDataSet (after size factor estimation)
#' @param blind   Use blind = TRUE for QC; FALSE for downstream with known groups
#' @returns SummarizedExperiment with assay "vst" (genes x samples matrix)
normalize_vst <- function(dds, blind = TRUE) {
  message("[NORM] Computing VST...")
  tictoc::tic("VST")
  vst_obj <- DESeq2::vst(dds, blind = blind)
  tictoc::toc()
  message(sprintf("[NORM] VST complete: %d genes x %d samples", nrow(vst_obj), ncol(vst_obj)))
  vst_obj
}

#' Apply rlog transformation from DESeq2
#' Better than VST for small sample sizes (n < 30); slower.
#'
#' @param dds   DESeqDataSet
#' @param blind Blind to experimental design?
#' @returns SummarizedExperiment with rlog values
normalize_rlog <- function(dds, blind = TRUE) {
  message("[NORM] Computing rlog (slow for large n)...")
  tictoc::tic("rlog")
  rlog_obj <- DESeq2::rlog(dds, blind = blind)
  tictoc::toc()
  rlog_obj
}

#' Plot size factor estimates to check for outlier samples
#'
#' @param dds DESeqDataSet after estimateSizeFactors()
#' @param metadata Sample metadata
plot_size_factors <- function(dds, metadata) {
  sf_df <- data.frame(
    sample      = colnames(dds),
    size_factor = DESeq2::sizeFactors(dds),
    stringsAsFactors = FALSE
  )
  sf_df <- merge(sf_df, metadata, by = "sample", sort = FALSE)
  sf_df <- sf_df[order(sf_df$size_factor), ]
  sf_df$sample_ord <- factor(sf_df$sample, levels = sf_df$sample)

  ggplot(sf_df, aes(x = sample_ord, y = size_factor,
                    fill = genotype, shape = tissue)) +
    geom_col(alpha = 0.7) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_fill_manual(values = GENOTYPE_COLORS) +
    facet_wrap(~ age_months, scales = "free_x", nrow = 1, labeller = label_both) +
    labs(
      title = "DESeq2 size factors",
      x     = NULL, y = "Size factor",
      fill  = "Genotype", shape = "Tissue"
    ) +
    BASE_THEME +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

#' Compute and plot mean-variance relationship (before and after normalization)
#'
#' @param dds DESeqDataSet (before or after DESeq())
plot_mean_variance <- function(dds) {
  DESeq2::meanSdPlot(DESeq2::assay(DESeq2::vst(dds, blind = TRUE)))
}

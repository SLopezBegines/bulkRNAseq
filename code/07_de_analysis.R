# Differential Expression Analysis — bulkRNAseq Pipeline
# DESeq2 (primary) and edgeR (alternative), with cross-method comparison
# Adapted from snRNAseq_mouse/code/03_vulcano_plots.R

# 1. DESeq2 --------------------------------------------------------------------

#' Run full DESeq2 workflow and extract results for a contrast
#'
#' @param dds         DESeqDataSet (from build_dds(); genes already filtered)
#' @param contrast    Character vector of length 3: c("variable", "numerator", "denominator")
#'                    e.g. c("genotype", "5xFAD", "BL6")
#' @param shrink_type LFC shrinkage method: "apeglm" (recommended) | "ashr" | "normal"
#' @returns data.frame with gene, baseMean, log2FoldChange (shrunken), log2FoldChange_MLE,
#'          lfcSE, padj (|LFC|>FC_THRESH Wald test), svalue (apeglm only),
#'          direction, significant
run_deseq2 <- function(dds, contrast = c("genotype", "5xFAD", "BL6"),
                       shrink_type = "apeglm") {
  message(sprintf("[DESeq2] Running: %s vs %s (in %s)", contrast[2], contrast[3], contrast[1]))
  tictoc::tic("DESeq2")
  dds <- DESeq2::DESeq(dds, parallel = FALSE)
  tictoc::toc()

  message("[DESeq2] Extracting results...")
  # For apeglm shrinkage we need the coefficient name, not contrast
  if (shrink_type == "apeglm") {
    coef_name <- paste0(contrast[1], "_", contrast[2], "_vs_", contrast[3])
    coef_name <- gsub(";", ".", coef_name) # DESeq2 replaces ; with .
    # Test |LFC| > FC_THRESH directly; padj from this call encodes the lfcThreshold hypothesis
    res_raw <- DESeq2::results(dds, name = coef_name, alpha = P_VAL_THRESH,
      lfcThreshold = FC_THRESH, altHypothesis = "greaterAbs")
    # lfcShrink with lfcThreshold returns svalue — apeglm FDR for |LFC| > FC_THRESH
    res <- DESeq2::lfcShrink(dds,
      coef = coef_name, type = "apeglm",
      lfcThreshold = FC_THRESH, res = res_raw, quiet = TRUE
    )
  } else {
    res_raw <- DESeq2::results(dds, contrast = contrast, alpha = P_VAL_THRESH,
      lfcThreshold = FC_THRESH, altHypothesis = "greaterAbs")
    res <- DESeq2::lfcShrink(dds,
      contrast = contrast, type = shrink_type,
      res = res_raw, quiet = TRUE
    )
  }

  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  # apeglm lfcShrink with lfcThreshold replaces padj with svalue; restore padj for volcano
  if (!"padj" %in% names(res_df)) {
    res_df$padj <- res_raw$padj[match(rownames(res_df), rownames(res_raw))]
  }
  # Retain unshrunken MLE LFC for fair cross-method comparison (shrunken vs unshrunken is not a DESeq2/edgeR difference)
  res_df$log2FoldChange_MLE <- res_raw$log2FoldChange[match(rownames(res_df), rownames(res_raw))]

  # Significance: svalue for apeglm (resolves Issue 5); padj for other shrinkage types.
  # Both now test |LFC| > FC_THRESH — no redundant post-hoc |LFC| filter.
  if (shrink_type == "apeglm") {
    res_df$significant <- !is.na(res_df$svalue) & res_df$svalue < P_VAL_THRESH
  } else {
    res_df$significant <- !is.na(res_df$padj) & res_df$padj < P_VAL_THRESH
  }
  res_df$direction <- ifelse(
    res_df$significant & res_df$log2FoldChange > 0, "UP",
    ifelse(res_df$significant & res_df$log2FoldChange < 0, "DOWN", "NS")
  )

  n_up <- sum(res_df$direction == "UP", na.rm = TRUE)
  n_down <- sum(res_df$direction == "DOWN", na.rm = TRUE)
  message(sprintf(
    "[DESeq2] DE genes: %d UP, %d DOWN (svalue/padj<%.2f, |LFC|>%.1f threshold-tested)",
    n_up, n_down, P_VAL_THRESH, FC_THRESH
  ))

  res_df[order(res_df$padj, na.last = TRUE), ]
}

# 2. edgeR ---------------------------------------------------------------------

#' Run full edgeR quasi-likelihood workflow for a genotype contrast
#'
#' @param counts    Filtered gene x sample count matrix
#' @param metadata  Sample metadata data.frame
#' @param design    Model matrix formula as a string (e.g. "~ sex + age_months + genotype")
#'                  The last term should be the contrast variable.
#' @param contrast_col  Column in metadata for the main contrast
#' @param numerator     Level to compare (numerator)
#' @param denominator   Reference level (denominator)
#' @returns data.frame with gene, logFC, logCPM, F, PValue, FDR, direction, significant
run_edger <- function(counts, metadata,
                      design_formula = ~genotype,
                      contrast_col = "genotype",
                      numerator = "5xFAD",
                      denominator = "BL6") {
  message(sprintf("[edgeR] Running: %s vs %s", numerator, denominator))

  metadata[[contrast_col]] <- relevel(
    factor(metadata[[contrast_col]]),
    ref = denominator
  )
  # Match build_dds(): 5xFAD pathology is non-linear with age, so treat as factor in both methods
  if ("age_months" %in% all.vars(design_formula)) {
    metadata$age_months <- factor(metadata$age_months)
  }

  design <- model.matrix(design_formula, data = metadata)

  dge <- edgeR::DGEList(counts = counts, samples = metadata)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  dge <- edgeR::estimateDisp(dge, design, robust = TRUE)

  tictoc::tic("edgeR glmQLFit")
  fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
  tictoc::toc()

  # Coefficient for the contrast (last column of design matrix = numerator level)
  geno_cols <- grep(paste0("^", contrast_col), colnames(design), value = TRUE)
  stopifnot(length(geno_cols) == 1)
  coef_name <- geno_cols
  # glmTreat tests |LFC| > FC_THRESH directly (TREAT; FDR-controlled for effect-size)
  qlf <- edgeR::glmTreat(fit, coef = coef_name, lfc = FC_THRESH)

  res_df <- edgeR::topTags(qlf,
    n = Inf, adjust.method = "BH",
    sort.by = "PValue"
  )$table
  res_df$gene <- rownames(res_df)
  if (!"F" %in% names(res_df)) res_df$F <- NA_real_  # glmTreat uses LR, not F

  res_df$significant <- res_df$FDR < P_VAL_THRESH  # no redundant |LFC| post-filter
  res_df$direction <- ifelse(
    res_df$significant & res_df$logFC > 0, "UP",
    ifelse(res_df$significant & res_df$logFC < 0, "DOWN", "NS")
  )

  n_up <- sum(res_df$direction == "UP", na.rm = TRUE)
  n_down <- sum(res_df$direction == "DOWN", na.rm = TRUE)
  message(sprintf(
    "[edgeR] DE genes: %d UP, %d DOWN (FDR<%.2f, |LFC|>%.1f threshold-tested)",
    n_up, n_down, P_VAL_THRESH, FC_THRESH
  ))

  res_df[order(res_df$PValue, na.last = TRUE), ]
}

# 3. Visualization -------------------------------------------------------------

#' Volcano plot (adapted from snRNAseq_mouse/code/03_vulcano_plots.R)
#' Works with both DESeq2 and edgeR output.
#'
#' @param res_df      DE results data.frame
#' @param title       Plot title
#' @param lfc_col     Column with log2 fold-change values
#' @param pval_col    Column with (adjusted) p-values
#' @param label_top_n Number of top genes to label by |LFC|
plot_volcano <- function(res_df, title,
                         lfc_col = "log2FoldChange",
                         pval_col = "padj",
                         label_top_n = 20) {
  df <- res_df[!is.na(res_df[[pval_col]]) & res_df[[pval_col]] > 0, ]
  df$direction <- factor(df$direction, levels = c("UP", "DOWN", "NS"))

  # Top genes to label
  df_sig <- df[df$significant, ]
  top_ids <- df_sig$gene[order(abs(df_sig[[lfc_col]]), decreasing = TRUE)][
    seq_len(min(label_top_n, nrow(df_sig)))
  ]
  df$label <- ifelse(df$gene %in% top_ids, df$gene, NA_character_)

  colors <- c("UP" = "#E41A1C", "DOWN" = "#377EB8", "NS" = "grey70")

  ggplot(df, aes(
    x = .data[[lfc_col]],
    y = -log10(.data[[pval_col]]),
    colour = direction,
    label = label
  )) +
    geom_point(size = 1.2, alpha = 0.6) +
    geom_vline(
      xintercept = c(-FC_THRESH, FC_THRESH),
      linetype = "dashed", colour = "black", linewidth = 0.4
    ) +
    geom_hline(
      yintercept = -log10(P_VAL_THRESH),
      linetype = "dashed", colour = "black", linewidth = 0.4
    ) +
    ggrepel::geom_text_repel(
      na.rm = TRUE, size = 2.5, max.overlaps = 15,
      segment.colour = "grey50", segment.size = 0.3
    ) +
    scale_colour_manual(values = colors) +
    labs(
      title   = title,
      x       = expression(log[2] ~ "Fold Change"),
      y       = expression(-log[10] ~ "(adj. p-value)"),
      colour  = "Direction"
    ) +
    BASE_THEME
}

#' MA plot (log-ratio vs mean expression)
#'
#' @param res_df   DE results data.frame
#' @param title    Plot title
#' @param lfc_col  Column with log2FC
#' @param mean_col Column with mean expression (baseMean for DESeq2, logCPM for edgeR)
plot_ma <- function(res_df, title,
                    lfc_col = "log2FoldChange",
                    mean_col = "baseMean") {
  df <- res_df[!is.na(res_df[[lfc_col]]), ]
  x_var <- if (mean_col == "baseMean") log2(df[[mean_col]] + 1) else df[[mean_col]]

  df$x_plot <- x_var
  df$direction <- factor(df$direction, levels = c("UP", "DOWN", "NS"))
  colors <- c("UP" = "#E41A1C", "DOWN" = "#377EB8", "NS" = "grey70")

  ggplot(df, aes(x = x_plot, y = .data[[lfc_col]], colour = direction)) +
    geom_point(size = 1, alpha = 0.5) +
    geom_hline(yintercept = 0, colour = "black") +
    geom_hline(
      yintercept = c(-FC_THRESH, FC_THRESH),
      linetype = "dashed", colour = "black", linewidth = 0.4
    ) +
    scale_colour_manual(values = colors) +
    labs(
      title = title,
      x = if (mean_col == "baseMean") {
        expression(log[2] ~ "(mean expression + 1)")
      } else {
        "Average log CPM"
      },
      y = expression(log[2] ~ "Fold Change"),
      colour = "Direction"
    ) +
    BASE_THEME
}

#' Heatmap of top DE genes (VST expression across samples)
#'
#' @param res_df    DE results data.frame
#' @param vst_obj   SummarizedExperiment from normalize_vst()
#' @param metadata  Sample metadata
#' @param n_top     Number of top DE genes to show
#' @param title     Heatmap title
plot_de_heatmap <- function(res_df, vst_obj, metadata,
                            n_top = 50, title = "Top DE genes") {
  sig_genes <- res_df$gene[res_df$significant]
  if (length(sig_genes) == 0) {
    message("[HEATMAP] No significant genes to plot.")
    return(invisible(NULL))
  }
  top_genes <- sig_genes[seq_len(min(n_top, length(sig_genes)))]

  mat <- SummarizedExperiment::assay(vst_obj)
  mat <- mat[rownames(mat) %in% top_genes, , drop = FALSE]
  mat <- mat[match(top_genes[top_genes %in% rownames(mat)], rownames(mat)), ]
  mat <- t(scale(t(mat))) # z-score across samples

  ann_df <- metadata[colnames(mat), c("genotype", "tissue", "age_months", "sex"),
    drop = FALSE
  ]
  ann_df$age_months <- as.character(ann_df$age_months)

  pheatmap::pheatmap(
    mat,
    annotation_col = ann_df,
    annotation_colors = list(
      genotype   = GENOTYPE_COLORS,
      tissue     = TISSUE_COLORS,
      sex        = SEX_COLORS,
      age_months = AGE_COLORS
    ),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_method = "ward.D2",
    show_rownames = nrow(mat) <= 60,
    show_colnames = FALSE,
    color = colorRampPalette(c("#053061", "white", "#67001F"))(100),
    breaks = seq(-3, 3, length.out = 101),
    main = title,
    border_color = NA,
    fontsize_row = 6
  )
}

# 4. Cross-method comparison ---------------------------------------------------

#' Merge DESeq2 and edgeR results for direct comparison
#'
#' @param deseq2_res   Output of run_deseq2()
#' @param edger_res    Output of run_edger()
#' @returns Merged data.frame with columns from both methods, suffixed _deseq2 and _edger
merge_de_results <- function(deseq2_res, edger_res) {
  d2 <- deseq2_res[, c("gene", "log2FoldChange", "log2FoldChange_MLE", "padj", "significant", "direction")]
  names(d2)[2:6] <- paste0(names(d2)[2:6], "_deseq2")

  er <- edger_res[, c("gene", "logFC", "FDR", "significant", "direction")]
  names(er)[2:5] <- paste0(names(er)[2:5], "_edger")

  merged <- merge(d2, er, by = "gene", all = TRUE)
  merged$agreement <- with(
    merged,
    ifelse(significant_deseq2 & significant_edger & direction_deseq2 == direction_edger,
      "Concordant",
      ifelse(significant_deseq2 | significant_edger, "Method-specific", "NS")
    )
  )
  merged
}

#' Scatter plot comparing log2FC estimates from DESeq2 and edgeR
plot_lfc_comparison <- function(merged_df, title = "DESeq2 vs edgeR: unshrunken log2FC (MLE)") {
  df <- merged_df[!is.na(merged_df$log2FoldChange_MLE_deseq2) &
    !is.na(merged_df$logFC_edger), ]
  cor_val <- round(cor(df$log2FoldChange_MLE_deseq2, df$logFC_edger,
    use = "complete.obs"
  ), 3)
  colors <- c("Concordant" = "#4DAF4A", "Method-specific" = "#FF7F00", "NS" = "grey80")

  ggplot(df, aes(x = log2FoldChange_MLE_deseq2, y = logFC_edger, colour = agreement)) +
    geom_point(size = 1, alpha = 0.4) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_smooth(method = "lm", se = FALSE, colour = "black", linewidth = 0.6) +
    scale_colour_manual(values = colors) +
    annotate("text",
      x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
      label = sprintf("r = %s", cor_val), size = 4
    ) +
    labs(
      title  = title,
      x      = expression(log[2] ~ "FC (DESeq2, unshrunken MLE)"),
      y      = expression(log[2] ~ "FC (edgeR, unshrunken)"),
      colour = "Status"
    ) +
    BASE_THEME
}

#' Venn diagram of significant DE genes from both methods (text-based summary)
#' Returns a printable data.frame; use ggVennDiagram or VennDiagram for plots.
venn_de_overlap <- function(merged_df) {
  d2_sig <- merged_df$gene[!is.na(merged_df$significant_deseq2) & merged_df$significant_deseq2]
  er_sig <- merged_df$gene[!is.na(merged_df$significant_edger) & merged_df$significant_edger]
  shared <- intersect(d2_sig, er_sig)
  only_d2 <- setdiff(d2_sig, er_sig)
  only_er <- setdiff(er_sig, d2_sig)

  cat("=== DE overlap ===\n")
  cat(sprintf("  DESeq2 only:    %d genes\n", length(only_d2)))
  cat(sprintf("  edgeR only:     %d genes\n", length(only_er)))
  cat(sprintf("  Shared:         %d genes\n", length(shared)))
  cat(sprintf(
    "  Jaccard index:  %.3f\n",
    length(shared) / length(union(d2_sig, er_sig))
  ))

  invisible(list(deseq2_only = only_d2, edger_only = only_er, shared = shared))
}

# 5. Deferred: surrogate variable / hidden-batch analysis ----------------------
# Issue 6 — With 192 samples across many animals, technical structure not
# captured by sex/age/tissue/genotype is likely. A future module (code/09_batch.R)
# should estimate surrogate variables via sva::svaseq or RUVSeq and append them
# to the design formula. Deferred: implement as an optional, parameter-gated step
# outside the default path; document as exploratory until validated.

#' P-value rank comparison plot: DESeq2 padj vs edgeR FDR (ranked)
plot_pval_comparison <- function(merged_df, title = "Adjusted p-value ranks: DESeq2 vs edgeR") {
  df <- merged_df[!is.na(merged_df$padj_deseq2) & !is.na(merged_df$FDR_edger), ]
  df$rank_d2 <- rank(df$padj_deseq2)
  df$rank_er <- rank(df$FDR_edger)

  ggplot(df, aes(x = rank_d2, y = rank_er, colour = agreement)) +
    geom_point(size = 0.8, alpha = 0.4) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    scale_colour_manual(values = c(
      "Concordant" = "#4DAF4A",
      "Method-specific" = "#FF7F00",
      "NS" = "grey80"
    )) +
    labs(
      title = title,
      x = "Rank by DESeq2 padj", y = "Rank by edgeR FDR",
      colour = "Status"
    ) +
    BASE_THEME
}

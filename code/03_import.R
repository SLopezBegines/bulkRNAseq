# ==============================================================================
# Data Import — bulkRNAseq Pipeline
# GSE168137: 5xFAD mouse model bulk RNA-seq
#
# Column name convention:
#   5xFAD;BL6_{tissue}_{age}mon_{sex}_{id}   (transgenic)
#   BL6_{tissue}_{age}mon_{sex}_{id}          (wildtype)
# ==============================================================================

#' Parse sample metadata encoded in GSE168137 column names
#'
#' @param sample_names Character vector of column names from the count matrix
#' @returns data.frame with columns: sample, genotype, tissue, age_months,
#'   sex, animal_id, group (genotype:tissue:age:sex composite label)
parse_sample_metadata <- function(sample_names) {
  parse_one <- function(x) {
    genotype  <- ifelse(grepl("^5xFAD", x), "5xFAD", "BL6")
    stripped  <- sub("^5xFAD;BL6_", "", x)
    stripped  <- sub("^BL6_",       "", stripped)
    parts     <- strsplit(stripped, "_")[[1]]
    # parts: [1]=tissue [2]=Xmon [3]=sex [4]=animal_id
    list(
      sample     = x,
      genotype   = genotype,
      tissue     = parts[1],
      age_months = as.integer(sub("mon", "", parts[2])),
      sex        = parts[3],
      animal_id  = parts[4]
    )
  }
  meta <- do.call(rbind, lapply(sample_names, function(x) as.data.frame(parse_one(x))))
  meta$age_months <- as.integer(meta$age_months)
  meta$group      <- with(meta,
    paste(genotype, tissue, paste0(age_months, "mo"), sex, sep = "_"))
  rownames(meta)  <- meta$sample
  meta
}

#' Load the GSE168137 count matrix and paired metadata
#'
#' @param count_file  Path to GSE168137_countList.txt.gz
#' @param strip_version  Strip Ensembl gene version suffix (e.g. ".4")?
#' @returns List with:
#'   - counts:   integer matrix (genes x samples)
#'   - metadata: data.frame with sample annotations
load_counts <- function(count_file, strip_version = TRUE) {
  message("[IMPORT] Reading count matrix from: ", count_file)
  tictoc::tic("read count matrix")
  raw <- data.table::fread(count_file, sep = "\t", header = TRUE,
                           data.table = FALSE, check.names = FALSE)
  tictoc::toc()

  gene_ids <- raw[[1]]
  if (strip_version) gene_ids <- sub("\\.[0-9]+$", "", gene_ids)

  counts <- as.matrix(raw[, -1, drop = FALSE])
  # Round any fractional counts (one gene ENSMUSG00000000028 has 39.42 in one sample)
  counts <- round(counts)
  storage.mode(counts) <- "integer"
  rownames(counts) <- gene_ids

  # Remove duplicate gene IDs if any (keep first occurrence)
  dup <- duplicated(rownames(counts))
  if (any(dup)) {
    message(sprintf("[IMPORT] Removing %d duplicate gene IDs", sum(dup)))
    counts <- counts[!dup, ]
  }

  sample_names <- colnames(counts)
  metadata     <- parse_sample_metadata(sample_names)

  message(sprintf("[IMPORT] Loaded %d genes x %d samples", nrow(counts), ncol(counts)))
  message("[IMPORT] Genotypes: ", paste(unique(metadata$genotype),   collapse = ", "))
  message("[IMPORT] Tissues:   ", paste(unique(metadata$tissue),     collapse = ", "))
  message("[IMPORT] Ages (mo): ", paste(sort(unique(metadata$age_months)), collapse = ", "))
  message("[IMPORT] Sexes:     ", paste(unique(metadata$sex),        collapse = ", "))

  list(counts = counts, metadata = metadata)
}

#' Subset count matrix and metadata to a specific experimental group
#'
#' @param counts    Gene x sample count matrix
#' @param metadata  Sample metadata data.frame
#' @param tissue    "cortex" | "hippocampus" | NULL (keep both)
#' @param age       Integer vector of ages to keep, e.g. c(4, 8) | NULL (keep all)
#' @param sex       "Female" | "Male" | NULL (keep both)
#' @returns List with filtered counts and metadata
subset_experiment <- function(counts, metadata, tissue = NULL, age = NULL, sex = NULL) {
  keep <- rep(TRUE, nrow(metadata))
  if (!is.null(tissue)) keep <- keep & metadata$tissue     %in% tissue
  if (!is.null(age))    keep <- keep & metadata$age_months %in% age
  if (!is.null(sex))    keep <- keep & metadata$sex        %in% sex

  n_before <- ncol(counts)
  counts   <- counts[, keep, drop = FALSE]
  metadata <- metadata[keep, , drop = FALSE]
  message(sprintf("[SUBSET] %d -> %d samples after filtering", n_before, ncol(counts)))
  list(counts = counts, metadata = metadata)
}

#' Build a DESeqDataSet from a count matrix and metadata
#'
#' @param counts    Gene x sample integer matrix
#' @param metadata  Sample metadata data.frame (rownames must match colnames of counts)
#' @param design    Formula for the DESeq2 model, e.g. ~ sex + age_months + genotype
#' @returns DESeqDataSet object
build_dds <- function(counts, metadata, design = ~ genotype) {
  stopifnot(all(colnames(counts) == rownames(metadata)))
  # Relevel: BL6 (wildtype) as reference for genotype
  if ("genotype" %in% all.vars(design)) {
    metadata$genotype <- factor(metadata$genotype, levels = c("BL6", "5xFAD"))
  }
  if ("tissue" %in% all.vars(design)) {
    metadata$tissue <- factor(metadata$tissue)
  }
  if ("sex" %in% all.vars(design)) {
    metadata$sex <- factor(metadata$sex)
  }
  if ("age_months" %in% all.vars(design)) {
    metadata$age_months <- factor(metadata$age_months)
  }
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData   = metadata,
    design    = design
  )
  message(sprintf("[DDS] DESeqDataSet built: %d genes x %d samples", nrow(dds), ncol(dds)))
  dds
}

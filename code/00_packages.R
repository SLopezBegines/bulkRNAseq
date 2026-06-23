# ==============================================================================
# Package Loading — bulkRNAseq Pipeline
# Uses pak for all installations: parallel, handles Bioc + GitHub, fast.
# ==============================================================================

if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak", repos = sprintf(
    "https://r-lib.github.io/p/pak/stable/%s/%s/%s",
    .Platform$pkgType, R.Version()$os, R.Version()$arch
  ))
}
library(pak)

load_pkgs <- function(pkgs) {
  bare <- sub("^[^:]+::", "", pkgs)
  missing_idx <- !vapply(bare, requireNamespace, logical(1), quietly = TRUE)
  if (any(missing_idx)) {
    message("[INSTALL] Installing: ", paste(pkgs[missing_idx], collapse = ", "))
    pak::pkg_install(pkgs[missing_idx], ask = FALSE, upgrade = FALSE)
  }
  invisible(lapply(bare, function(p) {
    suppressPackageStartupMessages(suppressWarnings(
      library(p, character.only = TRUE, warn.conflicts = FALSE, quietly = TRUE)
    ))
  }))
}

# CRAN packages ----------------------------------------------------------------
message("[PACKAGES] Loading CRAN packages...")
cran_packages <- c(
  # Core data manipulation
  "dplyr", "tidyr", "tibble", "stringr", "purrr", "data.table",
  # Visualization
  "ggplot2", "patchwork", "ggrepel", "scales", "RColorBrewer",
  "pheatmap", "viridis",
  # Dimensionality reduction
  "uwot", # UMAP (pure R, no Python dependency)
  # I/O
  "readxl", "writexl",
  # Utilities
  "tictoc", "remotes",
  "kableExtra",
  "R.utils"
)
load_pkgs(cran_packages)

# Bioconductor packages --------------------------------------------------------
message("[PACKAGES] Loading Bioconductor packages...")
bioc_packages <- paste0("bioc::", c(
  "DESeq2",
  "edgeR",
  "limma",
  "EnhancedVolcano",
  "org.Mm.eg.db", # Mouse gene annotations
  "AnnotationDbi",
  "BiocParallel",
  "SummarizedExperiment",
  "ggbio",
  "apeglm",
  "vsn"
))
load_pkgs(bioc_packages)

message("\n[PACKAGES] All packages loaded.")
cat("R version:", R.version$version.string, "\n")
cat("DESeq2:", as.character(packageVersion("DESeq2")), "\n")
cat("edgeR:", as.character(packageVersion("edgeR")), "\n")
cat("apeglm:", as.character(packageVersion("apeglm")), "\n")

# renv::snapshot(prompt = FALSE, type = "explicit", lockfile = "renv.lock")
rm(load_pkgs, cran_packages, bioc_packages)

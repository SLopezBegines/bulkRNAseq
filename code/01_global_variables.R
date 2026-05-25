# ==============================================================================
# Global Variables — bulkRNAseq Pipeline
# Adapted from snRNAseq_mouse/code/global_variables.R
# ==============================================================================

# --- Output file extensions ---------------------------------------------------
tiff_extension <- ".tiff"
pdf_extension  <- ".pdf"

# --- Species / organism -------------------------------------------------------
species        <- 10090          # NCBI Taxonomy ID: Mus musculus
organism       <- "org.Mm.eg.db" # Bioconductor annotation package
kegg_organism  <- "mmu"          # KEGG organism code
keyType        <- "ENSEMBL"

# --- Statistical thresholds ---------------------------------------------------
P_VAL_THRESH   <- 0.05   # Adjusted p-value cutoff for DE significance
FC_THRESH      <- 1.0    # |log2 fold-change| cutoff for DE significance
P_VAL_STRICT   <- 0.01   # Stricter threshold for highlighting top hits

# --- Bulk RNA-seq filtering thresholds ----------------------------------------
# Applied during gene-level pre-filtering (before DESeq2/edgeR)
MIN_COUNT_CPM       <- 1     # Minimum CPM in at least MIN_SAMPLES_EXPRESSED samples
MIN_SAMPLES_EXPRESSED <- 3   # Minimum number of samples meeting MIN_COUNT_CPM
MIN_TOTAL_COUNTS    <- 10    # Minimum total counts across all samples per gene

# --- QC thresholds (bulk RNA-seq sample-level) --------------------------------
QC_MIN_LIBRARY_SIZE  <- 1e6    # Samples below this are flagged (low sequencing depth)
QC_MIN_GENES_DETECTED <- 5000  # Samples with fewer detected genes are flagged

# --- Normalization method ------------------------------------------------------
# "vst"  — variance-stabilizing transformation (DESeq2); fast, good for n>30
# "rlog" — regularized log (DESeq2); better for small n but slow
NORM_METHOD <- "vst"

# --- PCA / UMAP parameters ----------------------------------------------------
N_TOP_GENES_PCA  <- 500    # Top variable genes to use for PCA
N_PCS            <- 30     # Number of PCs to compute
N_PCS_UMAP       <- 15     # Number of PCs to pass to UMAP
UMAP_SEED        <- 42
UMAP_N_NEIGHBORS <- 15
UMAP_MIN_DIST    <- 0.3

# --- Plot aesthetics ----------------------------------------------------------
# Colour palettes for key experimental factors
GENOTYPE_COLORS <- c("5xFAD" = "#E41A1C", "BL6" = "#377EB8")
TISSUE_COLORS   <- c("cortex" = "#4DAF4A", "hippocampus" = "#984EA3")
SEX_COLORS      <- c("Female" = "#FF7F00", "Male" = "#A65628")
AGE_COLORS      <- c("4" = "#FED976", "8" = "#FD8D3C", "12" = "#E31A1C", "18" = "#800026")

BASE_THEME <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92"),
    legend.position  = "right"
  )

# --- Checkpoint naming --------------------------------------------------------
CHECKPOINT_PREFIX <- "checkpoint_"

# --- File counter (incremented by save_plot) ----------------------------------
image_number <- 1

# --- Directory setup ----------------------------------------------------------
create_directories <- function(base_path) {
  dirs <- c(
    "figures", "figures/QC", "figures/normalization",
    "figures/PCA", "figures/UMAP",
    "figures/DE", "figures/DE/DESeq2", "figures/DE/edgeR",
    "figures/DE/comparison",
    "tables", "tables/QC", "tables/DE",
    "reports",
    "RData"
  )
  for (d in dirs) {
    path <- file.path(base_path, d)
    if (!dir.exists(path)) {
      dir.create(path, recursive = TRUE, showWarnings = FALSE)
      message("  [DIR] Created: ", path)
    }
  }
  message("[OK] Directory structure ready at: ", base_path)
}

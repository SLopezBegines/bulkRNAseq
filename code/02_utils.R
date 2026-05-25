# ==============================================================================
# Utility Functions — bulkRNAseq Pipeline
# Ported and adapted from snRNAseq_mouse/code/01_aux_functions.R
# ==============================================================================

# 1. CHECKPOINT SYSTEM ---------------------------------------------------------
# Enables crash recovery: save progress -> reload on restart.
# Pattern: check_checkpoint() -> if not found, compute -> save_checkpoint()

#' Save a named R object as a checkpoint .rds file (atomic write)
save_checkpoint <- function(obj, name, base = output_path) {
  dir <- file.path(base, "RData")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  tmp   <- file.path(dir, paste0(CHECKPOINT_PREFIX, name, "_TMP.rds"))
  final <- file.path(dir, paste0(CHECKPOINT_PREFIX, name, ".rds"))
  saveRDS(obj, file = tmp)
  file.rename(tmp, final)
  message(sprintf("[CHECKPOINT] Saved: %s  (%s)", final, format(Sys.time(), "%H:%M:%S")))
}

#' Check whether a checkpoint file exists
check_checkpoint <- function(name, base = output_path) {
  file.exists(file.path(base, "RData", paste0(CHECKPOINT_PREFIX, name, ".rds")))
}

#' Load a checkpoint file and return the object
load_checkpoint <- function(name, base = output_path) {
  path <- file.path(base, "RData", paste0(CHECKPOINT_PREFIX, name, ".rds"))
  if (!file.exists(path)) stop("[CHECKPOINT] Not found: ", path)
  message(sprintf("[CHECKPOINT] Loaded: %s  (%s)", path, format(Sys.time(), "%H:%M:%S")))
  readRDS(path)
}

# 2. PLOT HELPERS --------------------------------------------------------------

#' Save a ggplot to TIFF and PDF; increment global image_number counter
#' @param plotname  Base name (no extension)
#' @param plot      A ggplot or patchwork object
#' @param width     Width in inches
#' @param height    Height in inches
#' @param subdir    Subdirectory under output_path/figures/ (e.g. "QC", "DE/DESeq2")
save_plot <- function(plotname, plot, width = 8, height = 6, subdir = "") {
  dir <- file.path(output_path, "figures", subdir)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  filename <- file.path(dir, sprintf("%03d_%s", image_number, plotname))
  tryCatch(
    {
      if (inherits(plot, c("ggplot", "patchwork"))) {
        ggsave(paste0(filename, tiff_extension), plot,
               width = width, height = height, units = "in", dpi = 300)
        ggsave(paste0(filename, pdf_extension),  plot,
               width = width, height = height, units = "in")
      } else if (is.function(plot)) {
        tiff(paste0(filename, tiff_extension), width = width, height = height, units = "in", res = 300)
        plot(); dev.off()
        pdf(paste0(filename, pdf_extension), width = width, height = height)
        plot(); dev.off()
      } else {
        warning(sprintf("[SAVE_PLOT] Unsupported plot type '%s', skipping.", class(plot)[1]))
        return(invisible(NULL))
      }
      image_number <<- image_number + 1
    },
    error = function(e) warning(sprintf("[SAVE_PLOT] Failed '%s': %s", filename, e$message))
  )
}

# 3. TABLE HELPERS -------------------------------------------------------------

#' Save a data frame as a TSV file with a timestamped filename
#' @param df     Data frame to save
#' @param name   Base name (no extension)
#' @param subdir Subdirectory under output_path/tables/
save_table <- function(df, name, subdir = "") {
  dir <- file.path(output_path, "tables", subdir)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  filename <- file.path(dir, paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", name, ".tsv"))
  tryCatch(
    {
      write.table(df, file = filename, sep = "\t", row.names = FALSE, quote = FALSE)
      message(sprintf("[SAVE_TABLE] Saved: %s", filename))
      invisible(filename)
    },
    error = function(e) {
      warning(sprintf("[SAVE_TABLE] Failed '%s': %s", filename, e$message))
      invisible(NULL)
    }
  )
}

# 4. MEMORY -------------------------------------------------------------------

#' Return current R session RAM usage in MB
ram_mb <- function() round(sum(gc(verbose = FALSE)[, 2]), 0)

# 5. LOGGING ------------------------------------------------------------------

#' Initialize a timestamped log file
#' @param script_name  Label for the log filename
#' @param log_dir      Directory for log files
#' @returns List with log_msg() and close_log()
setup_logging <- function(script_name, log_dir = "logs") {
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  log_file <- file.path(log_dir,
    paste0(script_name, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
  log_con <- file(log_file, open = "wt")

  log_msg <- function(...) {
    line <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(...))
    message(line)
    writeLines(line, log_con)
  }

  close_log <- function() {
    if (tryCatch(isOpen(log_con), error = function(e) FALSE)) {
      log_msg("[EXIT] Script finished.")
      close(log_con)
    }
  }

  list(log_msg = log_msg, close_log = close_log, log_file = log_file)
}

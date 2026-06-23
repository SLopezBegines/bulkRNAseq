# Task: Fix Statistical Errors in the bulkRNAseq Pipeline

**Audience:** Claude Code, working inside the `bulkRNAseq` repository.
**Goal:** Correct verified statistical-methodology errors in the differential-expression
workflow. This is a methodology fix, not a refactor. Do not restyle, rename, or
reorganise code that is not listed below.

---

## Ground rules (read first)

1. **Obey the repo `CLAUDE.md` safety rules.** Only edit the files named in each task.
   Use targeted edits, not full rewrites. Read the full function before changing it.
2. **Do not "fix" the things listed under [Do NOT change](#do-not-change-these-are-correct).**
   They were checked and are correct; changing them will introduce bugs.
3. **Each fix must keep the pipeline runnable.** After each change, confirm the
   function still sources without error (`source("code/07_de_analysis.R")`) and that
   the notebook chunks that call it still have matching argument names.
4. **Statistics over convenience.** Where a choice exists, prefer the option that
   controls the false discovery rate correctly, even if it changes the gene counts.
   Changed numbers are expected and acceptable.
5. **One commit per issue**, message prefixed `stats-fix:`. Do not push.

---

## Severity summary

| # | Issue | Severity | Files |
|---|-------|----------|-------|
| 1 | `tissue` confounder omitted from the default design while all tissues are analysed | **High** | `notebooks/analysis.qmd` |
| 2 | Significance uses a post-hoc \|log2FC\| ≥ 1 filter that does not control FDR for the effect-size hypothesis | **High** | `code/07_de_analysis.R` |
| 3 | DESeq2 vs edgeR are not the same model (age = factor in DESeq2, continuous in edgeR) | **High** | `code/07_de_analysis.R`, `code/03_import.R` |
| 4 | Cross-method comparison correlates *shrunken* DESeq2 LFC against *unshrunken* edgeR LFC | **Medium** | `code/07_de_analysis.R` |
| 5 | Significance call mixes a p-value (H0: LFC=0) with a threshold applied to the shrunken estimate | **Medium** | `code/07_de_analysis.R` |
| 6 | No hidden-batch / surrogate-variable handling | **Low (enhancement)** | new module |
| 7 | No correction for running many contrasts/strata | **Low (enhancement)** | `notebooks/analysis.qmd` |
| 8 | edgeR contrast coefficient selected by string concatenation (fragile) | **Low** | `code/07_de_analysis.R` |

---

## Issue 1 — `tissue` is a confounder but is omitted from the default design (High)

**Where:** `notebooks/analysis.qmd`, params chunk (~line 57 and ~line 66).

```r
TISSUE_FILTER <- NULL                       # keeps ALL tissues
DE_DESIGN     <- ~ sex + age_months + genotype   # no tissue term
```

**Why it is wrong:** With `TISSUE_FILTER = NULL` all 192 samples (cortex + hippocampus)
enter one model, but the design has no `tissue` term. PCA in this very pipeline shows
tissue is the dominant axis of variance — larger than the genotype effect. Omitting it
leaves cortex/hippocampus differences in the residual, inflating dispersion estimates
and making the genotype p-values anti-conservative. The inline comment even says
"if subsetting to a single tissue, drop tissue from the formula" — but the default does
the opposite: all tissues, no tissue term.

> The full, config-driven version of this fix (global variables + a `run_de()` wrapper
> that assembles the design and extracts contrasts) is specified in
> [Issue 9](#issue-9--de-design-one-fit-many-contrasts-driven-by-global-variables-high).
> Issue 1 below is the minimal guard; Issue 9 is the complete framework.

**Fix (preferred — region-stratified):** make the intended workflow explicit. Run the
notebook once per tissue with the tissue dropped from the formula, e.g.

```r
TISSUE_FILTER <- "cortex"        # then re-run with "hippocampus"
DE_DESIGN     <- ~ sex + age_months + genotype
```

Add a short guard so the two never disagree silently — if more than one tissue is
present in the analysis subset, `tissue` must be in the formula:

```r
n_tissues <- length(unique(meta_sub$tissue))
if (n_tissues > 1 && !grepl("tissue", deparse(DE_DESIGN))) {
  stop("Multiple tissues in the subset but 'tissue' is not in DE_DESIGN. ",
       "Either set TISSUE_FILTER to a single tissue, or add '+ tissue' to DE_DESIGN.")
}
```

**Acceptable alternative (combined, exploratory only):** if a pooled run is wanted,
use `DE_DESIGN <- ~ tissue + sex + age_months + genotype` and label the output clearly
as exploratory (a single shared dispersion across two very different tissues is still a
compromise; it cannot recover region-specific effects).

**Verify:** running with `TISSUE_FILTER = NULL` and a tissue-free formula now errors
with the guard message; per-tissue runs complete.

---

## Issue 2 — Post-hoc |log2FC| ≥ 1 filter does not control FDR for the effect-size hypothesis (High)

**Where:** `code/07_de_analysis.R`, `run_deseq2()` (lines ~44–45) and `run_edger()`
(line ~105).

```r
# DESeq2
res_df$significant <- !is.na(res_df$padj) &
  res_df$padj < P_VAL_THRESH & abs(res_df$log2FoldChange) >= FC_THRESH
# edgeR
res_df$significant <- res_df$FDR < P_VAL_THRESH & abs(res_df$logFC) >= FC_THRESH
```

**Why it is wrong:** `padj`/`FDR` here test H0: LFC = 0. Filtering those genes
afterwards by `|LFC| ≥ 1` does **not** give an FDR-controlled set for the hypothesis
"|LFC| ≥ 1". It is the classic "test against zero, then filter on fold change" mistake,
which inflates the effective false discovery rate among the genes you report. Both
DESeq2 and edgeR provide a correct test against a fold-change threshold.

**Fix — DESeq2 (`run_deseq2`):** test against the threshold directly instead of
filtering after the fact.

```r
# apeglm branch: ask lfcShrink to test |LFC| > FC_THRESH -> returns s-values
res <- DESeq2::lfcShrink(
  dds, coef = coef_name, type = "apeglm",
  lfcThreshold = FC_THRESH, res = res_raw, quiet = TRUE
)
# significance by s-value (svalue is the apeglm analogue of an FDR-controlled quantity)
res_df$significant <- !is.na(res_df$svalue) & res_df$svalue < P_VAL_THRESH
```

For the non-apeglm branch, use `results(..., lfcThreshold = FC_THRESH,
altHypothesis = "greaterAbs")` and keep `significant <- padj < P_VAL_THRESH`
(the padj now already encodes the threshold test, so do **not** re-filter on `|LFC|`).

**Fix — edgeR (`run_edger`):** replace `glmQLFTest` with `glmTreat`, which tests
against a log2-FC threshold (TREAT).

```r
qlf <- edgeR::glmTreat(fit, coef = coef_name, lfc = FC_THRESH)
# ... topTags as before ...
res_df$significant <- res_df$FDR < P_VAL_THRESH   # no extra |logFC| filter
```

**Important:** once you switch to threshold-based testing, **remove the redundant
`abs(...) >= FC_THRESH`** from the `significant` line in both functions, otherwise you
double-apply the threshold.

**Verify:** gene counts will drop relative to the old post-hoc filter — expected.
Confirm `direction` is still assigned from the sign of the (shrunken, DESeq2 /
TREAT, edgeR) LFC, and that downstream summaries/plots still run.

---

## Issue 3 — DESeq2 and edgeR are not fitting the same model (High)

**Where:** `code/03_import.R` `build_dds()` line ~121 vs `code/07_de_analysis.R`
`run_edger()` line ~85.

- `build_dds()` does `metadata$age_months <- factor(metadata$age_months)`, so **DESeq2
  treats age as a factor** (one term per timepoint).
- `run_edger()` does `model.matrix(design_formula, data = metadata)` on metadata where
  `age_months` is an **integer** (from `parse_sample_metadata`, line ~32), so **edgeR
  treats age as a continuous linear covariate**.

The notebook claims "we pass the same design formula to ensure comparability" — but the
two methods are adjusting for age differently, which biases the cross-method comparison
(Issue 4) and the concordance statistics reported.

**Decision required (pick one and apply consistently):**

- **Recommended — age as factor in both** (4/8/12/18 mo pathology is non-linear, so a
  linear age term is a poor adjustment). In `run_edger()`, factorise before building the
  model matrix:

  ```r
  if ("age_months" %in% all.vars(design_formula)) {
    metadata$age_months <- factor(metadata$age_months)
  }
  ```

- **Alternative — age continuous in both:** remove the `factor(age_months)` line in
  `build_dds()`. Only choose this if a specific monotonic-with-age hypothesis is intended.

Whichever is chosen, the two methods must end up with the same encoding for every
covariate. Add a one-line comment recording the decision and its rationale.

**Verify:** print `colnames(model.matrix(...))` for edgeR and `resultsNames(dds)` for
DESeq2; the covariate structure (number of age terms) must match.

---

## Issue 4 — Comparing shrunken vs unshrunken fold-changes (Medium)

**Where:** `code/07_de_analysis.R`, `merge_de_results()` (line ~269) and
`plot_lfc_comparison()` (line ~290).

**Why it is wrong:** `log2FoldChange_deseq2` is the **apeglm-shrunken** estimate, while
edgeR's `logFC` is **unshrunken**. Correlating them (the reported `r`) and labelling
"method-specific" genes partly measures the effect of shrinkage, not a real DESeq2-vs-
edgeR difference. A fair comparison uses comparable estimators.

**Fix:** keep the unshrunken DESeq2 MLE LFC for the comparison. In `run_deseq2()`,
retain `res_raw$log2FoldChange` as an extra column before returning, e.g.
`res_df$log2FoldChange_MLE <- res_raw$log2FoldChange[match(...)]`. Then in
`merge_de_results()`/`plot_lfc_comparison()`, use `log2FoldChange_MLE` (DESeq2) vs
`logFC` (edgeR) for the scatter and the correlation. Continue using the **shrunken**
LFC for the volcano/reporting.

**Verify:** the comparison scatter is now MLE-vs-MLE-like; the reported `r` may change.
Caption/notebook text should state both axes are unshrunken estimates.

---

## Issue 5 — Significance mixes a zero-null p-value with the shrunken estimate (Medium)

**Where:** `code/07_de_analysis.R`, `run_deseq2()` apeglm branch.

The current `significant` uses `padj` (from `res_raw`, testing LFC = 0) together with the
**shrunken** `log2FoldChange`. The p-value and the effect-size estimate refer to
different quantities.

**Fix:** this is resolved automatically by adopting Issue 2's `lfcThreshold` approach
(s-values from `lfcShrink(..., lfcThreshold=)` are internally consistent with the
shrunken estimate). No separate change needed if Issue 2 is implemented; just confirm
the apeglm branch no longer combines `res_raw$padj` with shrunken-LFC filtering.

---

## Issue 6 — No hidden-batch / surrogate-variable handling (Low, enhancement)

192 samples generated across many animals almost certainly carry technical structure not
captured by `sex`/`age`/`tissue`. Consider adding an optional step using `sva::svaseq`
(or `RUVSeq`) to estimate surrogate variables on the filtered counts and append them to
the design. Implement as a new, optional helper (e.g. `code/08_batch.R`) gated behind a
parameter; do not force it into the default path. Document the assumption that it is
exploratory.

---

## Issue 7 — Multiple testing across strata/contrasts (Low, enhancement)

Once region-stratified (Issue 1) and/or per-timepoint runs exist, several contrasts are
performed, each BH-corrected within ~18k genes. The number of contrasts is not accounted
for. Add a short note in the notebook (and, if a combined gene list is produced, apply a
second-level adjustment across contrasts). This is a documentation + reporting fix, not a
code-path rewrite.

---

## Issue 8 — Fragile edgeR coefficient selection (Low)

**Where:** `code/07_de_analysis.R`, `run_edger()` line ~96.

```r
coef_name <- paste0(contrast_col, numerator)   # e.g. "genotype5xFAD"
```

This depends on R's name-mangling and breaks if the contrast variable or level naming
changes. Safer: select the genotype coefficient from the actual design matrix, e.g.

```r
geno_cols <- grep(paste0("^", contrast_col), colnames(design), value = TRUE)
stopifnot(length(geno_cols) == 1)
coef_name <- geno_cols
```

**Verify:** `coef_name` matches a real column of `design`.

---

## Issue 9 — DE design: one fit, many contrasts, driven by global variables (High)

This section defines the analysis framework that makes Issues 1, 2, 3 and 5 fall out
correctly, and lets the user change *what is compared* by editing global variables only —
no edits to the modelling code.

### Principle: one fit per stratum, many contrasts — never pairwise re-fits

Do **not** subset to two groups and re-run DESeq2 for each comparison. Dispersion
(gene-wise variance) is estimated by sharing information across all samples in the fit;
re-fitting on small two-group subsets gives noisier dispersions, inconsistent size
factors, and a larger multiple-testing burden. Instead: **one `DESeq()` fit per
biological stratum, then extract each comparison as a contrast** from that fit.

The trade-off to respect:

- **Additive model** (`~ sex + age + tissue + genotype`) returns *one* effect per factor,
  averaged over the others. Clean and powerful, but assumes no interaction — which is
  biologically wrong for 5xFAD (the genotype effect is region- and age-dependent).
- **Interaction model** (`+ genotype:age`) lets the effect vary; more correct but costs
  power. Two-way interactions are feasible here (~6/cell); three- and four-way are not.
- **Stratification** (separate fit per tissue) is the honest default, because cortex and
  hippocampus have genuinely different variance structures and should not share one
  dispersion model. This is the recommended primary frame.

So: **stratify by tissue, fit `~ sex + age + genotype` within each tissue, and extract
genotype / age / sex as contrasts.** Add `genotype:age` within a tissue when the question
is age-dependence. Compare tissues either by contrasting the two stratified runs, or — if
a formal tissue test is wanted — one pooled model reported as exploratory.

### Global-variable configuration block

Add to `code/01_global_variables.R`. These supersede the ad-hoc `TISSUE_FILTER` /
`DE_DESIGN` params currently in `notebooks/analysis.qmd`; consolidate so the globals are
the single source of truth and the notebook only reads them.

```r
# --- DE: subset BEFORE modelling (NULL keeps all levels) ----------------------
SUBSET_GENOTYPE <- NULL            # "5xFAD" | "BL6" | NULL
SUBSET_TISSUE   <- NULL            # "cortex" | "hippocampus" | NULL
SUBSET_AGE      <- NULL            # e.g. c(4, 8) | NULL
SUBSET_SEX      <- NULL            # "Male" | "Female" | NULL

# --- DE: what to compare ------------------------------------------------------
DE_FACTOR       <- "genotype"      # factor whose differences you want: "genotype" | "tissue" | "sex" | "age_months"
DE_NUMERATOR    <- "5xFAD"         # Wald contrast numerator (level of DE_FACTOR)
DE_DENOMINATOR  <- "BL6"           # Wald contrast denominator (= reference level)
DE_COVARIATES   <- c("sex", "age_months")   # nuisance adjusters; must NOT contain DE_FACTOR
DE_INTERACTION  <- NULL            # e.g. "age_months" -> adds DE_FACTOR:age_months; NULL = additive
DE_TEST         <- "Wald"          # "Wald" (one 2-level contrast) | "LRT" (whole factor >2 levels, or interaction)
AGE_AS_FACTOR   <- TRUE            # treat age as factor (recommended; pathology is non-linear)
```

### Reference wrapper

Add as `code/08_de_run.R` (keep the optional batch helper from Issue 6 as
`code/09_batch.R` to avoid a filename clash). It assembles the design from the globals,
applies the Issue 1 guard, fits once, and returns a tidy results table using the
threshold-based testing from Issue 2.

```r
#' Config-driven DE: subset -> assemble design -> one fit -> extract contrast.
#' All behaviour is controlled by the DE_* / SUBSET_* globals.
run_de <- function(counts, metadata) {
  # 1. Subset (extend subset_experiment to also filter genotype)
  sub  <- subset_experiment(counts, metadata,
                            tissue = SUBSET_TISSUE, age = SUBSET_AGE, sex = SUBSET_SEX)
  cts  <- sub$counts; meta <- sub$metadata
  if (!is.null(SUBSET_GENOTYPE)) {
    keep <- meta$genotype %in% SUBSET_GENOTYPE
    cts  <- cts[, keep, drop = FALSE]; meta <- meta[keep, , drop = FALSE]
  }

  # 2. Guard against the Issue-1 confounder: if >1 tissue remains, tissue must be
  #    either the factor of interest or a covariate.
  if (length(unique(meta$tissue)) > 1 &&
      DE_FACTOR != "tissue" && !("tissue" %in% DE_COVARIATES)) {
    stop("Multiple tissues present but 'tissue' is neither DE_FACTOR nor a covariate. ",
         "Set SUBSET_TISSUE to one tissue, or add 'tissue' to DE_COVARIATES.")
  }

  # 3. Reference level of the factor of interest = denominator (so the apeglm coef
  #    name is well-defined). Drop covariates that collapsed to one level after subsetting.
  meta[[DE_FACTOR]] <- relevel(factor(meta[[DE_FACTOR]]), ref = DE_DENOMINATOR)
  covars <- setdiff(DE_COVARIATES, DE_FACTOR)
  covars <- covars[vapply(covars, function(v) length(unique(meta[[v]])) > 1, logical(1))]

  # 4. Assemble formula; factor of interest (and its interaction) go LAST.
  rhs <- paste(c(covars, DE_FACTOR), collapse = " + ")
  if (!is.null(DE_INTERACTION)) rhs <- paste0(rhs, " + ", DE_FACTOR, ":", DE_INTERACTION)
  design <- stats::as.formula(paste("~", rhs))
  message("[DE] design: ", deparse(design),
          "  |  test: ", DE_TEST, "  |  n = ", ncol(cts))

  # 5. Build dds (build_dds must: respect the already-set DE_FACTOR reference, and
  #    factorise age_months iff AGE_AS_FACTOR) and fit ONCE.
  dds <- build_dds(cts, meta, design = design)

  if (DE_TEST == "LRT") {
    drop_term   <- if (!is.null(DE_INTERACTION)) paste0(DE_FACTOR, ":", DE_INTERACTION) else DE_FACTOR
    reduced_rhs <- setdiff(strsplit(rhs, " \\+ ")[[1]], drop_term)
    reduced     <- stats::as.formula(paste("~", paste(reduced_rhs, collapse = " + ")))
    dds <- DESeq2::DESeq(dds, test = "LRT", reduced = reduced)
    res <- DESeq2::results(dds, alpha = P_VAL_THRESH)              # omnibus LRT padj
    res_df <- as.data.frame(res); res_df$gene <- rownames(res_df)
    res_df$significant <- !is.na(res_df$padj) & res_df$padj < P_VAL_THRESH
  } else {
    dds <- DESeq2::DESeq(dds)
    coef_name <- paste0(DE_FACTOR, "_", DE_NUMERATOR, "_vs_", DE_DENOMINATOR)
    stopifnot(coef_name %in% DESeq2::resultsNames(dds))           # fail loudly if mis-specified
    res_raw <- DESeq2::results(dds, name = coef_name,
                               lfcThreshold = FC_THRESH, altHypothesis = "greaterAbs",
                               alpha = P_VAL_THRESH)
    res <- DESeq2::lfcShrink(dds, coef = coef_name, type = "apeglm",
                             lfcThreshold = FC_THRESH, res = res_raw, quiet = TRUE)
    res_df <- as.data.frame(res); res_df$gene <- rownames(res_df)
    res_df$significant <- !is.na(res_df$svalue) & res_df$svalue < P_VAL_THRESH
  }
  res_df[order(res_df$padj, na.last = TRUE), ]
}
```

This requires two small supporting edits: extend `subset_experiment()` to accept a
`genotype` filter, and generalise `build_dds()` so it does **not** hard-override the
factor-of-interest levels (it currently forces genotype levels) — it should leave an
already-`factor()`-ed column as-is and only enforce `factor(age_months)` when
`AGE_AS_FACTOR`.

### Worked examples (change globals only)

**A. Cortex vs hippocampus (tissue main effect).** Both tissues must be in, so this is a
pooled contrast — report as exploratory (shared dispersion across tissues).

```r
SUBSET_TISSUE  <- NULL
DE_FACTOR      <- "tissue";  DE_NUMERATOR <- "cortex"; DE_DENOMINATOR <- "hippocampus"
DE_COVARIATES  <- c("sex", "age_months", "genotype")
DE_INTERACTION <- NULL;      DE_TEST <- "Wald"
```

**B. Hippocampus only — does the genotype effect change with age?** (genotype × age)

```r
SUBSET_TISSUE  <- "hippocampus"
DE_FACTOR      <- "genotype"; DE_NUMERATOR <- "5xFAD"; DE_DENOMINATOR <- "BL6"
DE_COVARIATES  <- c("sex", "age_months")
DE_INTERACTION <- "age_months"        # adds genotype:age_months
DE_TEST        <- "LRT"               # omnibus: is there ANY age-dependence of the 5xFAD effect?
```

To then get the 5xFAD effect *at a specific age*, follow up with a Wald run
(`DE_INTERACTION <- NULL`, `SUBSET_AGE <- 12`, `DE_TEST <- "Wald"`).

**C. Hippocampus only — overall 5xFAD vs BL6** (additive baseline).

```r
SUBSET_TISSUE  <- "hippocampus"
DE_FACTOR      <- "genotype"; DE_NUMERATOR <- "5xFAD"; DE_DENOMINATOR <- "BL6"
DE_COVARIATES  <- c("sex", "age_months"); DE_INTERACTION <- NULL; DE_TEST <- "Wald"
```

### Interpretation caveats (do not skip)

- **With an interaction in the model, the plain `DE_FACTOR` coefficient is the effect at
  the *reference* level of the interacting variable, not the average.** Extract per-level
  effects by combining the main term with the relevant interaction term.
- **LRT returns an omnibus p-value, not a single fold-change.** Use it to ask "does this
  factor / interaction matter at all?", then use Wald contrasts for the effect sizes.
- **Power:** keep to additive + at most one two-way interaction. Do not attempt
  `genotype:age:sex` or higher here.
- **edgeR mirror:** the same globals drive `run_edger()` — build the model matrix from the
  identical `design`, `glmQLFit` once, then `glmTreat(coef=, lfc=FC_THRESH)` for a Wald
  contrast or `glmQLFTest` over interaction columns for the LRT-equivalent. Age encoding
  must match DESeq2 (Issue 3).

---

## Do NOT change (these are correct)

- **apeglm uses `coef`, not `contrast`** in `run_deseq2()` (lines ~24–32). This is the
  correct apeglm API. Leave it.
- **BL6 is the reference level** — `build_dds()` line ~112 sets
  `factor(genotype, levels = c("BL6", "5xFAD"))`. Correct; do not relevel elsewhere.
- **VST `blind` usage** — `blind = TRUE` for QC/PCA, `blind = FALSE` for the DE heatmap.
  This matches DESeq2 guidance. Leave it.
- **edgeR normalisation** — TMM (`calcNormFactors`) + `estimateDisp(robust = TRUE)` +
  robust `glmQLFit`. Correct.
- **Gene pre-filtering** (`filter_low_counts`, `filterByExpr`). Standard. Leave it.
- **Plot aesthetics, colours, themes.** Out of scope.

---

## Definition of done

- [ ] Issue 1: guard added; per-tissue runs verified; pooled run requires `+ tissue`.
- [ ] Issue 2: DESeq2 uses `lfcThreshold` (s-values); edgeR uses `glmTreat`; redundant
      `|LFC|` post-filter removed in both.
- [ ] Issue 3: age encoding identical across DESeq2 and edgeR; decision documented.
- [ ] Issue 4: cross-method comparison uses unshrunken LFCs on both axes.
- [ ] Issue 5: confirmed resolved via Issue 2.
- [ ] Issues 6–8: implemented or explicitly deferred with a note.
- [ ] Issue 9: `DE_*`/`SUBSET_*` globals added; `run_de()` implemented; `subset_experiment()`
      and `build_dds()` adjusted; worked examples A–C run; notebook reads globals only.
- [ ] `source()` of every edited `.R` file succeeds.
- [ ] The notebook renders (or the DE chunks execute) end-to-end on at least one tissue.
- [ ] One `stats-fix:` commit per issue. Not pushed.

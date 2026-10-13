########################################
# Script Name: GLM_Script.R
# Description: Per-SNP modeling of allele frequency change in experimentally
#              evolved S. cerevisiae populations, within and across antifungal
#              treatments.
# Author: Megan Sandoval-Powers
# Last Updated: 10-2026
# R Version: 4.5.1
#
# INPUT   Input_files/*_only_filt.txt
#           One table per treatment, produced by the SNP filtering script.
#           Columns alternate af_<sample> and N_<sample>, with the ancestor
#           relabeled as timepoint 00 of that treatment. Early-fixed sites
#           have already been removed.
#
# OUTPUTS Output_files/model.results.Dose.txt      drug intensity comparisons
#         Output_files/model.results.Complex.txt   single vs combined drug
#
# Candidate classification, annotation, and enrichment are handled in
# SNP_annotation_and_analysis.R.
########################################


# =============== PACKAGES ====================================================

library(tidyverse)   # dplyr, tidyr, stringr, readr
library(emmeans)     # emtrends() and pairs() for per-SNP slope contrasts
library(pbapply)     # progress bars for the per-SNP loops


# =============== CONSTANTS ===================================================

# Chromosome sizes for the S288C reference genome (R64-1-1), used to build
# cumulative genome-wide positions for plotting.
CHR_SIZES <- c(I = 230218, II = 813184, III = 316620, IV = 1531933,
               V = 576874, VI = 270161, VII = 1090940, VIII = 562643,
               IX = 439888, X = 745751, XI = 666816, XII = 1078177,
               XIII = 924431, XIV = 784333, XV = 1091291, XVI = 948066)


# =============== READ TREATMENT TABLES =======================================

treatment_files <- c(
  "Low CAS"  = "Input_files/Low_CAS_only_filt.txt",
  "High CAS" = "Input_files/High_CAS_only_filt.txt",
  "Low CLO"  = "Input_files/Low_CLO_only_filt.txt",
  "High CLO" = "Input_files/High_CLO_only_filt.txt",
  "CASCLO"   = "Input_files/CASCLO_only_filt.txt",
  "Control"  = "Input_files/Control_only_filt.txt")

treatments <- map(treatment_files, ~ read.table(.x, header = TRUE))


# =============== BUILD MODELING INDICES ======================================

# Rather than melting each SNP table into long format, which would be large and
# slow, we build a small index describing every sample column once. The modeling
# functions then pull one SNP's allele frequencies and coverages at a time and
# bind them to this index.
build_index <- function(df) {
  af_cols  <- grep("^af_", colnames(df))
  cov_cols <- grep("^N_",  colnames(df))
  
  # af and coverage columns must pair up one to one and in the same order
  stopifnot(length(af_cols) == length(cov_cols))
  stopifnot(all(str_replace(colnames(df)[af_cols], "^af_", "") ==
                  str_replace(colnames(df)[cov_cols], "^N_", "")))
  
  id <- colnames(df)[af_cols]
  parts <- str_split_fixed(id, "_", 6)   # af_EEA_<treatment>_<dose>_<rep>_<time>
  
  out <- data.frame(
    idaf      = af_cols,
    idcov     = cov_cols,
    treatment = factor(parts[, 3]),
    dose      = factor(parts[, 4]),
    rep       = factor(str_remove(parts[, 5], "^rep")),
    time      = as.numeric(parts[, 6]))
  
  # Fail loudly if the sample naming scheme has changed
  stopifnot(!any(is.na(out$time)),
            all(nzchar(as.character(out$treatment))),
            all(nzchar(as.character(out$dose))))
  
  out
}


# =============== COMBINE TABLES FOR CROSS-TREATMENT MODELS ===================

# Control samples are relabeled with the focal drug name and a dose of "ND"
# (no drug) so that dose can be modeled as a single factor with three levels.
relabel_control <- function(df, drug) {
  colnames(df) <- str_replace(colnames(df), "_Control_ND_",
                              paste0("_", drug, "_ND_"))
  df
}

# CAS intensity: Low CAS, High CAS, Control
Low_High_CAS_comb <- treatments[["Low CAS"]] %>%
  inner_join(treatments[["High CAS"]], by = c("CHROM", "POS")) %>%
  inner_join(treatments[["Control"]],  by = c("CHROM", "POS")) %>%
  relabel_control("CAS")

# CLO intensity: Low CLO, High CLO, Control
Low_High_CLO_comb <- treatments[["Low CLO"]] %>%
  inner_join(treatments[["High CLO"]], by = c("CHROM", "POS")) %>%
  inner_join(treatments[["Control"]],  by = c("CHROM", "POS")) %>%
  relabel_control("CLO")

# Single vs combined drug: the high single-drug treatments are used here, not
# the low concentrations CASCLO contains, so that differences reflect the drug
# combination rather than a difference in the severity of selection.
CAS_CLO_CASCLO_comb <- treatments[["High CAS"]] %>%
  inner_join(treatments[["High CLO"]], by = c("CHROM", "POS")) %>%
  inner_join(treatments[["CASCLO"]],   by = c("CHROM", "POS")) %>%
  inner_join(treatments[["Control"]],  by = c("CHROM", "POS"))


# =============== PER-SNP MODEL FUNCTIONS =====================================

# Assemble the modeling frame for a single SNP: its allele frequencies and coverages, bound to the sample index.
build_snp_frame <- function(snp_row, index_df) {
  freq <- cbind(as.numeric(snp_row[index_df$idaf]),
                as.numeric(snp_row[index_df$idcov]),
                index_df)
  colnames(freq)[1:2] <- c("af", "cov")
  freq
}

# Single-treatment model: does allele frequency change over time within this
# environment? Returns the slope (beta) and its p-value for one SNP.
#
# Quasibinomial accounts for overdispersion in pooled sequencing data, and
# coverage is used as a weight so better-covered sites carry more information.
# Convergence warnings are captured rather than discarded so they can be
# summarized across all SNPs afterward.
fit_single_treatment <- function(snp_row, index_df) {
  freq <- build_snp_frame(snp_row, index_df)
  
  warn_msg <- NULL
  mod <- withCallingHandlers(
    glm(af ~ time, family = "quasibinomial", weights = cov, data = freq),
    warning = function(w) {
      warn_msg <<- c(warn_msg, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  
  label <- paste0(unique(index_df$treatment), "_", unique(index_df$dose))
  
  data.frame(
    beta    = coef(mod)["time"],
    pval    = summary(mod)$coefficients["time", "Pr(>|t|)"],
    warning = if (is.null(warn_msg)) NA_character_
    else paste(unique(warn_msg), collapse = " | ")) %>%
    rename_with(~ paste0(.x, "_", label), c("beta", "pval")) %>%
    rename_with(~ paste0("warning_", label), "warning")
}

# Cross-treatment model: does the rate of allele frequency change differ among
# treatments? Fits time x group, estimates the slope within each group, and
# returns pairwise contrasts between those slopes for one SNP.
#
# group_var is "dose" for the intensity comparisons and "treatment" for the
# single vs combined drug comparison.
#
# NOTE ON MULTIPLE TESTING: pairs(adjust = "fdr") corrects across the pairwise contrasts made within a single SNP. In the downstream
# classification script, these contrasts are applied only to SNPs that have
# already passed a genome-wide FDR threshold (the top 5% of BH-adjusted
# p-values from the single-treatment models).
fit_cross_treatment <- function(snp_row, index_df, group_var, prefix = "") {
  freq <- build_snp_frame(snp_row, index_df)
  
  mod <- glm(as.formula(paste("af ~ time *", group_var)),
             family = "quasibinomial", weights = cov, data = freq)
  
  trends <- emtrends(mod, specs = group_var, var = "time")
  pair   <- summary(pairs(trends, adjust = "fdr"))
  
  contrasts <- gsub(" ", "", pair$contrast)
  out <- c(pair$p.value, pair$estimate)
  names(out) <- c(paste0("pval_adj_", prefix, contrasts),
                  paste0("beta_",     prefix, contrasts))
  
  as.data.frame(t(out))
}


# =============== CHUNKED RUNNER ==============================================

# Apply a per-SNP model function across a full SNP table in chunks, which keeps
# peak memory use manageable. Works with either model function above.
run_chunked <- function(snp_df, index_df, fit_fun, chunk_size = 1000, ...) {
  starts  <- seq(1, nrow(snp_df), by = chunk_size)
  results <- vector("list", length(starts))
  
  for (i in seq_along(starts)) {
    end <- min(starts[i] + chunk_size - 1, nrow(snp_df))
    cat("Processing rows", starts[i], "to", end, "\n")
    
    chunk <- snp_df[starts[i]:end, , drop = FALSE]
    results[[i]] <- pbapply::pbapply(chunk, 1,
                                     function(x) fit_fun(x, index_df, ...))
    gc()
  }
  
  out <- as.data.frame(do.call(rbind, do.call(c, results)))
  cbind(snp_df[, c("CHROM", "POS")], out)
}


# =============== RUN SINGLE-TREATMENT MODELS =================================

single_results <- imap(treatments, function(df, name) {
  cat("\n=== Single-treatment model:", name, "===\n")
  run_chunked(df, build_index(df), fit_single_treatment)
})

# Summarize GLM convergence warnings before dropping the warning columns.
# Reviewers may ask whether the quasibinomial fits converged.
iwalk(single_results, function(res, name) {
  wcol <- grep("^warning_", colnames(res), value = TRUE)
  n_warn <- sum(!is.na(res[[wcol]]))
  cat(name, ":", n_warn, "of", nrow(res), "SNPs produced a GLM warning\n")
})

# FDR-adjust the single-treatment p-values across all SNPs, and drop the warning columns now that they have been summarized.
adjust_pvals <- function(df) {
  pval_cols <- grep("^pval_", colnames(df), value = TRUE)
  
  for (col in pval_cols) {
    df[[sub("^pval_", "pval_adj_", col)]] <- p.adjust(df[[col]], method = "BH")
    df[[col]] <- NULL
  }
  
  df %>% select(-starts_with("warning_"))
}

single_results <- map(single_results, adjust_pvals)


# =============== RUN CROSS-TREATMENT MODELS ==================================

cat("\n=== Cross-treatment model: CAS intensity ===\n")
res_CAS_dose <- run_chunked(Low_High_CAS_comb,
                            build_index(Low_High_CAS_comb),
                            fit_cross_treatment,
                            group_var = "dose", prefix = "CAS_")

cat("\n=== Cross-treatment model: CLO intensity ===\n")
res_CLO_dose <- run_chunked(Low_High_CLO_comb,
                            build_index(Low_High_CLO_comb),
                            fit_cross_treatment,
                            group_var = "dose", prefix = "CLO_")

cat("\n=== Cross-treatment model: single vs combined drug ===\n")
res_complex <- run_chunked(CAS_CLO_CASCLO_comb,
                           build_index(CAS_CLO_CASCLO_comb),
                           fit_cross_treatment,
                           group_var = "treatment")


# =============== ASSEMBLE MASTER TABLES ======================================

# Add a cumulative genome-wide position in Mb and a SNP identifier, both used
# by the downstream classification and plotting script.
add_position_and_id <- function(df) {
  chr_starts <- c(0, cumsum(CHR_SIZES[-length(CHR_SIZES)]))
  names(chr_starts) <- names(CHR_SIZES)
  
  df %>%
    mutate(MB     = (POS + chr_starts[CHROM]) / 1e6,
           SNP_ID = paste0(CHROM, "_", POS)) %>%
    relocate(MB, SNP_ID, .after = POS)
}

# Master table 1: drug intensity comparisons
model.res.Dose <- reduce(
  list(single_results[["Low CAS"]],
       single_results[["High CAS"]],
       single_results[["Low CLO"]],
       single_results[["High CLO"]],
       single_results[["Control"]],
       res_CAS_dose,
       res_CLO_dose),
  inner_join, by = c("CHROM", "POS")) %>%
  add_position_and_id()

# Master table 2: single vs combined drug
model.res.Complex <- reduce(
  list(single_results[["High CAS"]],
       single_results[["High CLO"]],
       single_results[["CASCLO"]],
       single_results[["Control"]],
       res_complex),
  inner_join, by = c("CHROM", "POS")) %>%
  add_position_and_id()

# write_tsv(model.res.Dose,    "Output_files/model.results.Dose.txt")
# write_tsv(model.res.Complex, "Output_files/model.results.Complex.txt")


# =============== SESSION INFO ================================================

sessionInfo()

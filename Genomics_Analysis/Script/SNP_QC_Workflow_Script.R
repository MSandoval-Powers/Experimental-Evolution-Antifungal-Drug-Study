# =============================================================================
# Genomic analysis of Saccharomyces cerevisiae populations evolved under
# antifungal stress differing in intensity and complexity
#
# INPUT
#   Input_files/filtered_snps.txt
#     126,215 raw SNPs from [upstream variant-calling pipeline]. Columns
#     alternate alt_<sample> (alternate allele count) and N_<sample> (total
#     coverage), and include the ancestor, the four founders, and
#     mitochondrial SNPs.
#
# OUTPUTS
#   Output_files/snps_filtered_poly.txt           64,764 filtered nuclear SNPs
#   Output_files/Coverage_Summary_Filtered.csv    Table S1
#   Output_files/fixed_snps_filt.txt              755 early-fixed sites
#   Figures/PCA/combined_PCA_plot_raw.svg         Fig. S2
#   Figures/SFS/SFS_multi_panel_*_t14_raw.svg     Fig. S3
#   Figures/Coverage/chrVIII_coverage_fig_raw.svg Fig. S6
#
# Downstream modeling (quasibinomial GLMs, candidate classification,
# annotation, GO enrichment) is handled in other scripts.
#
# R version and package versions are recorded by sessionInfo() at the end.
# =============================================================================


# ================= SETUP =====================================================

library(tidyverse)   # ggplot2, dplyr, tidyr, stringr, readr
library(cowplot)     # theme_cowplot(), plot_grid() for SFS multipanels
library(patchwork)   # / operator, plot_annotation() for multipanel figures
library(factoextra)  # fviz_screeplot()
library(svglite)     # backend for ggsave(".svg")
library(ggrastr)     # geom_point_rast() keeps exported coverage plots small

# Resolve masking once rather than prefixing dplyr:: throughout
select <- dplyr::select
filter <- dplyr::filter

theme_set(theme_cowplot())

# Treatment colors used throughout. Defined once so every figure matches.
TREATMENT_COLORS <- c(
  "Ancestor" = "black",
  "Low CAS"  = "gold",
  "High CAS" = "darkgoldenrod",
  "Low CLO"  = "hotpink",
  "High CLO" = "firebrick4",
  "CASCLO"   = "orangered",
  "Control"  = "grey50")

TREATMENT_LEVELS <- c("Low CAS", "High CAS", "Low CLO", "High CLO",
                      "CASCLO", "Control")

# Chromosome sizes for the S288C reference genome (R64-1-1), used to build
# cumulative genome-wide positions for coverage plots.
CHR_SIZES <- c(I = 230218, II = 813184, III = 316620, IV = 1531933,
               V = 576874, VI = 270161, VII = 1090940, VIII = 562643,
               IX = 439888, X = 745751, XI = 666816, XII = 1078177,
               XIII = 924431, XIV = 784333, XV = 1091291, XVI = 948066)


# ================= HELPER FUNCTIONS ==========================================

#Summarize per-sample coverage from any SNP table with N_ columns. Used before and after filtering.
summarize_coverage <- function(df) {
  cov <- df %>%
    select(starts_with("N_")) %>%
    mutate(across(everything(), as.numeric))
  
  tibble(
    Population         = colnames(cov),
    Mean_Coverage      = colMeans(cov, na.rm = TRUE),
    Standard_deviation = apply(cov, 2, sd, na.rm = TRUE))
}

# Add a cumulative genome-wide position (in Mb) to a table with CHROM and POS.
add_genome_position <- function(df) {
  chr_starts <- c(0, cumsum(CHR_SIZES[-length(CHR_SIZES)]))
  names(chr_starts) <- names(CHR_SIZES)
  
  df %>%
    mutate(POS = as.numeric(POS),
           MB  = (POS + chr_starts[CHROM]) / 1e6)
}

# Pull the single column matching a pattern, erroring if it is not unique.
# Used instead of hard-coded column indices so the script survives changes
# to the input file's column order.
get_one_col <- function(df, pattern) {
  hits <- grep(pattern, colnames(df), value = TRUE)
  stopifnot(length(hits) == 1)
  hits
}


# ================= READ AND FORMAT SNP TABLE =================================

snps <- read.table("Input_files/filtered_snps.txt", header = TRUE)  # 126,215 SNPs

# Missing genotypes are encoded as "." in the input, where a sample had no
# reads at a site. Setting both the alt count and the coverage to 0 means these
# sites fail the minimum-coverage filter below, rather than propagating NAs
# into the allele frequency calculations.
snps_edited <- as.data.frame(lapply(snps, function(x) gsub("\\.", "0", x)))

# Rename sample columns to the treatment labels used in the manuscript.
colnames(snps_edited) <- colnames(snps_edited) %>%
  # The control has no dose but carries the S01 suffix used for low-dose
  # treatments. Rename it first so the general S01 -> Low rule below does not
  # relabel it as "Low Control".
  str_replace("CTR_S01", "Control_ND") %>%
  str_replace_all("CTR", "Control") %>%
  str_replace_all("CAC", "CASCLO") %>%
  str_replace_all("CAP", "CAS") %>%
  str_replace_all("S01", "Low") %>%
  str_replace_all("S03", "High") %>%
  str_replace_all("D01", "T01") %>%
  str_replace_all("D07", "T07") %>%
  str_replace_all("D14", "T14")


# ================= COVERAGE BEFORE FILTERING =================================

# Experimental populations only (drop the ancestor and founder columns).
snps_exp_pops <- snps_edited %>%
  select(-starts_with(c("alt_EEA_anc", "N_EEA_anc",
                        "alt_EEA_hap", "N_EEA_hap")))

summary_cov <- summarize_coverage(snps_exp_pops)

range(summary_cov$Mean_Coverage)     # 27X to 279X
median(summary_cov$Mean_Coverage)    # 83X
sum(summary_cov$Mean_Coverage < 30)  # 1 sample
sum(summary_cov$Mean_Coverage < 40)  # 25 samples


# ================= FILTER FOR HIGH-QUALITY POLYMORPHIC SNPs ==================

# 1. Minimum coverage. Require at least 10X in every experimental population.
#    A 20X threshold would leave fewer than 40,000 SNPs.
exp_cov <- snps_exp_pops %>%
  select(starts_with("N_")) %>%
  mutate(across(everything(), as.numeric))

snps_cov_filt <- snps_edited %>%
  select(CHROM, POS, 6:ncol(.)) %>%
  mutate(mincov = apply(exp_cov, 1, min)) %>%
  filter(mincov > 10) %>%                       # 123,612 SNPs
  mutate(across(-c(CHROM, POS), as.numeric))

# 2. Convert alternate allele counts to frequencies (alt / N) in place.
alt_cols <- grep("^alt_", colnames(snps_cov_filt), value = TRUE)
n_cols   <- str_replace(alt_cols, "^alt_", "N_")
stopifnot(all(n_cols %in% colnames(snps_cov_filt)))

snps_cov_filt[alt_cols] <- snps_cov_filt[alt_cols] / snps_cov_filt[n_cols]
colnames(snps_cov_filt) <- str_replace(colnames(snps_cov_filt), "^alt_", "af_")

# 3. Retain only sites polymorphic in the ancestor. This restricts the analysis
#    to standing genetic variation rather than de novo mutations or errors.
anc_af <- get_one_col(snps_cov_filt, "^af_EEA_anc")

snps_poly <- snps_cov_filt %>%
  filter(.data[[anc_af]] > 0.02, .data[[anc_af]] < 0.98)   # 70,137 SNPs

# 4. Founder-based filtering (following Phillips et al. 2021).
#    a. The four founders are haploid and isogenic, so a site appearing
#       polymorphic within a founder reflects mapping or sequencing error.
founder_af <- grep("^af_EEA_hap", colnames(snps_poly), value = TRUE)

snps_poly <- snps_poly %>%
  filter(if_all(all_of(founder_af), ~ .x == 0 | .x == 1))  # 64,809 SNPs

#    b. Sites at which all four founders carry the same allele cannot have
#       contributed standing variation to the 4S population.
snps_poly <- snps_poly %>%
  filter(!if_all(all_of(founder_af[-1]),
                 ~ .x == .data[[founder_af[1]]]))          # 64,781 SNPs

# 5. Drop mitochondrial SNPs.
sum(snps_poly$CHROM == "mitochondrion")   # 17
snps_filtered <- snps_poly %>% filter(CHROM != "mitochondrion")  # 64,764 SNPs

# write_tsv(snps_filtered, "Output_files/snps_filtered_poly.txt")


# ================= COVERAGE AFTER FILTERING (TABLE S1) =======================

filter_summary_cov <- snps_filtered %>%
  select(-starts_with(c("af_EEA_anc", "N_EEA_anc",
                        "af_EEA_hap", "N_EEA_hap", "mincov"))) %>%
  summarize_coverage()

range(filter_summary_cov$Mean_Coverage)   # 27X to 290X
median(filter_summary_cov$Mean_Coverage)  # 86X

# Reformat sample names into the treatment labels used in the manuscript.
coverage_table <- filter_summary_cov %>%
  separate(Population,
           into = c("N", "Exp", "Population", "Dose", "rep", "Timepoint"),
           sep = "_") %>%
  mutate(
    Replicate      = as.numeric(gsub("rep", "", rep)),
    Timepoint      = as.numeric(gsub("T", "", Timepoint)),
    Treatment      = paste(Population, Dose, sep = "_"),
    Treatment      = recode(Treatment,
                            "CASCLO_Low" = "CASCLO",
                            "CAS_Low"    = "Low CAS",
                            "CAS_High"   = "High CAS",
                            "CLO_Low"    = "Low CLO",
                            "CLO_High"   = "High CLO",
                            "Control_ND" = "Control"),
    Treatment      = factor(Treatment, levels = TREATMENT_LEVELS)) %>%
  select(Treatment, Replicate, Timepoint,
         Avg_Genome_Cov = Mean_Coverage)

# write_csv(coverage_table, "Output_files/Coverage_Summary_Filtered.csv")


# ================= PCA OF FILTERED SNPs (FIG. S2) ============================

# Allele frequency matrix: experimental populations plus the ancestor,
# excluding founders.
af_matrix <- snps_filtered %>%
  select(starts_with("af_"), -starts_with("af_EEA_hap"))

# Metadata for every sample, parsed from the column names. Built once and
# reused by each per-treatment PCA below.
sample_meta <- tibble(Sample = colnames(af_matrix)) %>%
  mutate(
    tag = str_remove(Sample, "^af_EEA_"),
    Population = case_when(
      str_detect(tag, "^anc")          ~ "Ancestor",
      str_detect(tag, "^CASCLO")       ~ "CASCLO",
      str_detect(tag, "^Control")      ~ "Control",
      str_detect(tag, "^CAS_Low")      ~ "Low CAS",
      str_detect(tag, "^CAS_High")     ~ "High CAS",
      str_detect(tag, "^CLO_Low")      ~ "Low CLO",
      str_detect(tag, "^CLO_High")     ~ "High CLO"),
    Replicate = str_extract(tag, "rep\\d+"),
    Timepoint = if_else(Population == "Ancestor", "0",
                        str_remove(str_extract(tag, "T\\d+$"), "T")),
    Timepoint = factor(as.character(as.numeric(Timepoint)),
                       levels = c("0", "1", "7", "14"))) %>%
  select(-tag)

stopifnot(!any(is.na(sample_meta$Population)))  # catch unparsed sample names

# Run a PCA on a chosen subset of populations and return the plotting frame
# together with the percent variance explained by PC1 and PC2.
run_pca <- function(populations) {
  samples <- sample_meta %>%
    filter(Population %in% populations) %>%
    pull(Sample)
  
  pca <- prcomp(t(as.matrix(af_matrix[, samples])), scale. = TRUE)
  var_per <- round(pca$sdev^2 / sum(pca$sdev^2) * 100, 1)
  
  list(
    data = tibble(Sample = rownames(pca$x),
                  PC1 = pca$x[, 1],
                  PC2 = pca$x[, 2]) %>%
      left_join(sample_meta, by = "Sample"),
    var_per = var_per)
}

# Build one faceted PCA panel. Used for the three comparisons in Fig. S2.
plot_pca_panel <- function(pca_result, populations) {
  df <- pca_result$data %>%
    mutate(Population = factor(Population, levels = populations))
  
  ggplot(df, aes(PC1, PC2, color = Population, shape = Population)) +
    geom_point(size = 5) +
    scale_color_manual(values = TREATMENT_COLORS[populations]) +
    scale_shape_manual(values = c(16, 15, rep(17, length(populations) - 2))) +
    facet_wrap(~Timepoint, nrow = 1,
               labeller = labeller(Timepoint = function(x) paste0("T", x))) +
    labs(x = paste0("PC1 - ", pca_result$var_per[1], "%"),
         y = paste0("PC2 - ", pca_result$var_per[2], "%")) +
    theme(legend.position = "right",
          strip.text    = element_text(size = 16, face = "bold"),
          legend.text   = element_text(size = 14),
          legend.title  = element_text(size = 16),
          axis.text     = element_text(size = 14),
          axis.title    = element_text(size = 16),
          panel.spacing = unit(1, "lines"),
          panel.border  = element_rect(color = "black", fill = NA,
                                       linewidth = 0.8))
}

cas_pops    <- c("Ancestor", "Control", "Low CAS", "High CAS")
clo_pops    <- c("Ancestor", "Control", "Low CLO", "High CLO")
casclo_pops <- c("Ancestor", "Control", "CASCLO", "High CAS", "High CLO")

pca_cas    <- run_pca(cas_pops)
pca_clo    <- run_pca(clo_pops)
pca_casclo <- run_pca(casclo_pops)

combined_PCA_plot <-
  plot_pca_panel(pca_cas, cas_pops) /
  plot_pca_panel(pca_clo, clo_pops) /
  plot_pca_panel(pca_casclo, casclo_pops) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 30, face = "bold"))

# ggsave("Figures/PCA/combined_PCA_plot_raw.svg", combined_PCA_plot,
#        width = 14, height = 10, dpi = 300)


# ================= PER-TREATMENT SNP TABLES ==================================

# Strip the "T" from timepoint suffixes so they can be treated as integers in
# the downstream modeling script, and relabel the ancestor as timepoint 00 of
# each treatment so every treatment table carries its own T0.
poly_snps <- snps_filtered %>% select(-mincov)
colnames(poly_snps) <- colnames(poly_snps) %>%
  str_replace_all("_T", "_") %>%
  str_replace("af_EEA_anc_4SH_12", "af_EEA_ANC_S00_rep00_00") %>%
  str_replace("N_EEA_anc_4SH_12",  "N_EEA_ANC_S00_rep00_00")

# Subset one treatment, carrying the ancestor columns along and renaming them
# to match that treatment so the modeling script sees a complete time series.
subset_treatment <- function(pattern, anc_tag) {
  anc_cols <- grep("ANC_S00", colnames(poly_snps), value = TRUE)
  trt_cols <- grep(pattern, colnames(poly_snps), value = TRUE)
  
  df <- poly_snps %>% select(CHROM, POS, all_of(anc_cols), all_of(trt_cols))
  colnames(df) <- str_replace(colnames(df), "ANC_S00", anc_tag)
  df
}

Low_CAS_only  <- subset_treatment("_CAS_Low_",  "CAS_Low")
High_CAS_only <- subset_treatment("_CAS_High_", "CAS_High")
Low_CLO_only  <- subset_treatment("_CLO_Low_",  "CLO_Low")
High_CLO_only <- subset_treatment("_CLO_High_", "CLO_High")
CASCLO_only   <- subset_treatment("_CASCLO_",   "CASCLO_Low")
Control_only  <- subset_treatment("_Control_",  "Control_ND")

treatments <- list("Low CAS"  = Low_CAS_only,
                   "High CAS" = High_CAS_only,
                   "Low CLO"  = Low_CLO_only,
                   "High CLO" = High_CLO_only,
                   "CASCLO"   = CASCLO_only,
                   "Control"  = Control_only)


# ================= SITE FREQUENCY SPECTRA (FIG. S3) ==========================

# Histogram of allele frequencies in each replicate at a given timepoint.
plot_SFS <- function(df, treatment_name, fill_color,
                     timepoint = "14", binwidth = 0.02) {
  
  af_cols <- grep(paste0("^af_.*_rep[0-9]{2}_", timepoint, "$"),
                  colnames(df), value = TRUE)
  af_cols <- af_cols[!grepl("rep00_00", af_cols)]   # drop the ancestor
  
  sfs_df <- df %>%
    select(all_of(af_cols)) %>%
    pivot_longer(everything(), names_to = "Replicate", values_to = "AF") %>%
    mutate(Replicate = str_extract(Replicate, "rep[0-9]{2}"),
           Replicate = factor(Replicate, levels = sort(unique(Replicate))))
  
  ggplot(sfs_df, aes(AF)) +
    geom_histogram(binwidth = binwidth, fill = fill_color, color = NA) +
    facet_wrap(~Replicate, nrow = 2, ncol = 6, axes = "all_x") +
    scale_y_continuous(limits = c(0, 30000), expand = c(0, 0)) +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    labs(x = "SNP frequency", y = "Count",
         title = paste0(treatment_name, " (T", as.integer(timepoint), ")")) +
    theme(axis.text   = element_text(size = 14),
          axis.title  = element_text(size = 16),
          strip.text  = element_text(size = 14),
          plot.title  = element_text(size = 16, face = "bold"))
}

# One panel per treatment at T14.
sfs_t14 <- imap(treatments,
                ~ plot_SFS(.x, .y, TREATMENT_COLORS[[.y]], timepoint = "14"))

SFS_CAS    <- plot_grid(sfs_t14[["Low CAS"]],  sfs_t14[["High CAS"]],
                        ncol = 1, align = "v",
                        labels = c("A", "B"), label_size = 30)
SFS_CLO    <- plot_grid(sfs_t14[["Low CLO"]],  sfs_t14[["High CLO"]],
                        ncol = 1, align = "v",
                        labels = c("C", "D"), label_size = 30)
SFS_CASCLO <- plot_grid(sfs_t14[["CASCLO"]],   sfs_t14[["Control"]],
                        ncol = 1, align = "v",
                        labels = c("E", "F"), label_size = 30)

# ggsave("Figures/SFS/SFS_multi_panel_CAS_t14_raw.svg", SFS_CAS,
#        width = 12, height = 10, dpi = 300)
# ggsave("Figures/SFS/SFS_multi_panel_CLO_t14_raw.svg", SFS_CLO,
#        width = 12, height = 10, dpi = 300)
# ggsave("Figures/SFS/SFS_multi_panel_CASCLO_Control_t14_raw.svg", SFS_CASCLO,
#        width = 12, height = 10, dpi = 300)


# ================= EARLY-FIXED SITES =========================================

# A small subset of SNPs reached fixation or loss (AF of 1 or 0) in every
# replicate at every sampled timepoint, including T1. Because they show no
# temporal variance they cannot be modeled by the GLMs used downstream, so
# they are identified here and analyzed separately.

get_fixed_snps <- function(df) {
  af_cols  <- grep("^af_", colnames(df), value = TRUE)
  anc_col  <- grep("rep00_00$", af_cols, value = TRUE)[1]
  evo_cols <- setdiff(af_cols, anc_col)
  
  evo_mat <- as.matrix(df[, evo_cols, drop = FALSE])
  
  all_zero <- apply(evo_mat, 1, function(x) all(x == 0))
  all_one  <- apply(evo_mat, 1, function(x) all(x == 1))
  
  tibble(
    SNP_ID        = paste(df$CHROM, df$POS, sep = "_"),
    FixationState = case_when(all_zero ~ "loss", all_one ~ "gain"),
    Fixed         = if_else(all_zero | all_one, "fixed", "not_fixed"),
    Ancestor_AF   = df[[anc_col]],
    Evo_AF        = rowMeans(evo_mat, na.rm = TRUE),
    Delta_AF      = rowMeans(evo_mat, na.rm = TRUE) - df[[anc_col]])
}

fixed_snps <- imap_dfr(treatments,
                       ~ get_fixed_snps(.x) %>% mutate(Treatment = .y))

all_fixed_snps <- fixed_snps %>%
  filter(Fixed == "fixed") %>%
  pull(SNP_ID) %>%
  unique()

length(all_fixed_snps)   # 2,159 before filtering

# Counts by treatment and direction. Losses greatly outnumber gains.
fixed_snps %>%
  filter(Fixed == "fixed") %>%
  count(Treatment, FixationState)

# Remove early-fixed sites from each treatment table so the remainder can be
# modeled downstream.
treatments_filt <- map(treatments,
                       ~ .x %>% filter(!paste(CHROM, POS, sep = "_") %in%
                                         all_fixed_snps))

# walk2(treatments_filt, names(treatments_filt),
#       ~ write_tsv(.x, paste0("Output_files/",
#                              str_replace_all(.y, " ", "_"), "_filt.txt")))

# Restrict the early-fixed set to drug-associated changes by excluding sites
# that also fixed in Controls (laboratory adaptation) and sites that were
# already near fixation in the ancestor.
control_fixed <- fixed_snps %>%
  filter(Fixed == "fixed", Treatment == "Control") %>%
  pull(SNP_ID) %>%
  unique()

fixed_snps_filt <- fixed_snps %>%
  filter(Fixed == "fixed",
         !SNP_ID %in% control_fixed,
         abs(Delta_AF) > 0.10)

length(unique(fixed_snps_filt$SNP_ID))   # 755 (~1.2% of 64,764)

# Per-treatment totals reported in the manuscript:
# CASCLO 507, High CAS 551, High CLO 24, Low CAS 140, Low CLO 15
fixed_snps_filt %>% count(Treatment, name = "n_sites")

# write_tsv(fixed_snps_filt, "Output_files/fixed_snps_filt.txt")


# ================= CHROMOSOME VIII COVERAGE (FIG. S6) ========================

# Allele frequency data cannot detect changes in chromosome copy number, so we
# examine normalized read depth. Coverage at each SNP is divided by the sample's
# own genome-wide median, placing all populations on a common scale regardless
# of library depth.
cov_long <- snps_edited %>%
  select(CHROM, POS, starts_with("N_")) %>%
  add_genome_position() %>%
  filter(CHROM != "mitochondrion") %>%
  pivot_longer(starts_with("N_"), names_to = "Sample", values_to = "Coverage") %>%
  mutate(
    Coverage  = as.numeric(Coverage),
    Timepoint = str_extract(Sample, "T\\d+"),
    Replicate = str_extract(Sample, "rep\\d+"),
    Treatment = case_when(
      # CASCLO must be tested first: its sample names contain both "CAS" and
      # "CLO" and would otherwise be miscalled as a single-drug treatment.
      str_detect(Sample, "CASCLO")                           ~ "CASCLO",
      str_detect(Sample, "CAS") & str_detect(Sample, "Low")  ~ "Low CAS",
      str_detect(Sample, "CAS") & str_detect(Sample, "High") ~ "High CAS",
      str_detect(Sample, "CLO") & str_detect(Sample, "Low")  ~ "Low CLO",
      str_detect(Sample, "CLO") & str_detect(Sample, "High") ~ "High CLO",
      str_detect(Sample, "Control")                          ~ "Control")) %>%
  filter(!is.na(Treatment)) %>%   # drops the ancestor and founder columns
  group_by(Sample) %>%
  mutate(rel_cov = Coverage / median(Coverage, na.rm = TRUE)) %>%
  ungroup()

# Panel A: median chromosome VIII coverage per replicate over time.
# 1.0 is ancestral copy number; 1.5 is one additional copy in a diploid if
# present in all cells.
chromVIII_panelA <- cov_long %>%
  filter(CHROM == "VIII") %>%
  group_by(Treatment, Replicate, Timepoint) %>%
  summarise(rel_cov = median(rel_cov, na.rm = TRUE), .groups = "drop") %>%
  mutate(Treatment = factor(Treatment, levels = TREATMENT_LEVELS)) %>%
  ggplot(aes(Timepoint, rel_cov, group = Replicate, color = Treatment)) +
  geom_hline(yintercept = c(1, 1.5), linetype = "dashed", color = "grey60") +
  geom_line(alpha = 0.7) +
  geom_point(size = 1.5) +
  facet_wrap(~Treatment, nrow = 2) +
  scale_color_manual(values = TREATMENT_COLORS[TREATMENT_LEVELS]) +
  scale_y_continuous(limits = c(0.9, 1.6), breaks = c(1.0, 1.25, 1.5)) +
  labs(x = NULL, y = "Normalized coverage (chr VIII)") +
  theme(legend.position = "none",
        axis.title = element_text(size = 15),
        axis.text  = element_text(size = 13),
        strip.text = element_text(size = 13, face = "bold"))

# Panel B: coverage along chromosome VIII in 20 kb windows. Uniform elevation
# across windows indicates a whole-chromosome gain rather than a segmental
# duplication. Subtelomeric regions are trimmed because repetitive sequence
# produces coverage artifacts in all samples, including Controls.
clo_pops_cov <- c("Control", "Low CLO", "High CLO", "CASCLO")

chromVIII_panelB <- cov_long %>%
  filter(CHROM == "VIII",
         Treatment %in% clo_pops_cov,
         POS > 20000, POS < 540000) %>%
  mutate(window    = floor(POS / 20000) * 20000,
         Treatment = factor(Treatment, levels = clo_pops_cov)) %>%
  group_by(Treatment, Timepoint, window) %>%
  summarise(rel_cov = median(rel_cov, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(window / 1000, rel_cov, color = Treatment,
             alpha = Timepoint, linewidth = Timepoint)) +
  annotate("rect", xmin = 120.091, xmax = 121.683,
           ymin = -Inf, ymax = Inf, alpha = 0.25, fill = "blue") +
  annotate("text", x = 120.1, y = 1.63, label = "ERG11 region",
           size = 4, hjust = -0.15, color = "blue") +
  geom_line() +
  facet_wrap(~Treatment, nrow = 1) +
  scale_color_manual(values = TREATMENT_COLORS[clo_pops_cov]) +
  scale_alpha_manual(values = c("T01" = 0.35, "T07" = 0.65, "T14" = 1)) +
  scale_discrete_manual("linewidth",
                        values = c("T01" = 0.7, "T07" = 1.0, "T14" = 1.4)) +
  scale_y_continuous(limits = c(0.8, 1.7)) +
  labs(x = "Position on chromosome VIII (kb)", y = "Normalized coverage") +
  guides(color = "none", linewidth = "none",
         alpha = guide_legend(title = "Timepoint",
                              override.aes = list(linewidth = 1.2,
                                                  color = "black"))) +
  theme(legend.position = "right",
        axis.title   = element_text(size = 15),
        axis.text    = element_text(size = 13),
        strip.text   = element_text(size = 13, face = "bold"),
        legend.title = element_text(size = 13),
        legend.text  = element_text(size = 13))

ChrVIII_coverage <- chromVIII_panelA / chromVIII_panelB +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 30, face = "bold"))

# ggsave("Figures/Coverage/chrVIII_coverage_fig_raw.svg", ChrVIII_coverage,
#        width = 12, height = 8, dpi = 300)


          
# ================= SESSION INFO ==============================================

sessionInfo()         


          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
        
          
          

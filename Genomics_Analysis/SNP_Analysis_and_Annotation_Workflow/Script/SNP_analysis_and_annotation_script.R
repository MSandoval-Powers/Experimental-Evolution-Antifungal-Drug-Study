########################################
# Script Name: SNP_analysis_and_annotation_script.R
# Description: Classification, annotation, and visualization of candidate SNPs from the GLM outputs.
# Author: Megan Sandoval-Powers
# Last Updated: 10-2026
# R Version: 4.5.1
#
# INPUTS  Input_files/model.results.Dose.txt      GLM output, drug intensity
#         Input_files/model.results.Complex.txt   GLM output, single vs combined
#         Input_files/*_only_filt.txt             per-treatment allele frequencies
#         Input_files/snps_annotated.txt          SnpEff output
#         Input_files/fixed_snps_filt.txt         early-fixed sites
#         Input_files/FungAMR_db.csv              resistance gene database
#         Input_files/snps_filtered_poly.txt      filtered SNP table
#         Input_files/filtered_snps.txt           unfiltered SNP table (REF/ALT)
#         Input_files/saccharomyces_cerevisiae_R64-1-1_20110208.gff
#
# OUTPUTS Figures/Manhattan_and_Upset_Plots/*.svg          Figs. 3 and 4
#         Figures/Beta_Coefficient_Plot/*.svg              Fig. 2
#         Figures/FungAMR_and_GO_plots/*.svg               Fig. 5
#         Figures/FKS1_Gene_Track_Plot/*.svg               Fig. S4
#         Figures/AF_Trajectory_Plots/*.svg                Fig. S5
#         Output_files/*_gene_level_summary.xlsx           Tables S2-S5
#         Output_files/*_GO_enrichment_results.csv         Tables S6 and S7
########################################


# ======================= R PACKAGES ==========================================

library(tidyverse)      # ggplot2, dplyr, tidyr, stringr, purrr, readr
library(cowplot)        # theme_cowplot()
library(patchwork)      # plot composition with / and |
library(ggh4x)          # themed facet strips (strip_themed, facet_grid2)
library(ComplexUpset)   # UpSet plots
library(ggrastr)        # rasterize dense point layers to keep SVGs small
library(multcompView)   # compact letter display for pairwise tests
library(topGO)          # GO enrichment
library(rrvgo)          # semantic reduction of GO terms
library(org.Sc.sgd.db)  # S. cerevisiae gene annotations
library(rtracklayer)    # import() for the GFF
library(openxlsx)       # Excel export


# ======================= R ENVIRONMENT SETTINGS ==============================

theme_set(theme_cowplot())
options(scipen = 999)


# ======================= SHARED CONSTANTS ====================================

# Treatment colors used across every figure.
TREATMENT_COLORS <- c(
  "Low CAS"  = "gold",
  "High CAS" = "darkgoldenrod",
  "Low CLO"  = "hotpink",
  "High CLO" = "firebrick4",
  "CASCLO"   = "orangered",
  "Control"  = "grey50")

# The single vs combined drug analysis drops the dose suffix from its column
# names, so the underlying classes are "CAS" and "CLO". Those treatments are
# the HIGH concentrations, and figures label them accordingly.
COMPLEX_LABELS <- c("CAS" = "High CAS", "CLO" = "High CLO",
                    "CASCLO" = "CASCLO", "Control" = "Control",
                    "Shared" = "Shared")

CHROM_ORDER <- c("I", "II", "III", "IV", "V", "VI", "VII", "VIII",
                 "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI")


# ======================= DATAFILES ===========================================

# GLM results
model.results.Dose    <- read_tsv("Input_files/model.results.Dose.txt")
model.results.Complex <- read_tsv("Input_files/model.results.Complex.txt")

# Per-treatment allele frequency tables, used for trajectory plots
treatment_files <- c(
  "Low CAS"  = "Input_files/Low_CAS_only_filt.txt",
  "High CAS" = "Input_files/High_CAS_only_filt.txt",
  "Low CLO"  = "Input_files/Low_CLO_only_filt.txt",
  "High CLO" = "Input_files/High_CLO_only_filt.txt",
  "CASCLO"   = "Input_files/CASCLO_only_filt.txt",
  "Control"  = "Input_files/Control_only_filt.txt")

treatments <- map(treatment_files, read_tsv)

snps_ann           <- read_tsv("Input_files/snps_annotated.txt", comment = "##")
fixed_snps_filt    <- read_tsv("Input_files/fixed_snps_filt.txt")
fungamr_db         <- read.csv("Input_files/FungAMR_db.csv", fileEncoding = "latin1")
snps_filtered_poly <- read_tsv("Input_files/snps_filtered_poly.txt")
snps               <- read.table("Input_files/filtered_snps.txt", header = TRUE)


# ======================= CANDIDATE IDENTIFICATION ============================

# Candidates are the top 5% of SNPs by FDR-adjusted p-value within each
# treatment. These thresholds are reused later to draw the Manhattan
# significance lines.
top_pct <- 0.05

# Significance threshold for the pairwise slope contrasts. Note that these
# p-values were adjusted across the pairwise comparisons made within each SNP
# (see the GLM script).
pcutoff_pairwise <- 0.05

quant_Low_CAS  <- quantile(model.results.Dose$pval_adj_CAS_Low,     top_pct, na.rm = TRUE)
quant_High_CAS <- quantile(model.results.Dose$pval_adj_CAS_High,    top_pct, na.rm = TRUE)
quant_Low_CLO  <- quantile(model.results.Dose$pval_adj_CLO_Low,     top_pct, na.rm = TRUE)
quant_High_CLO <- quantile(model.results.Dose$pval_adj_CLO_High,    top_pct, na.rm = TRUE)
quant_Control  <- quantile(model.results.Dose$pval_adj_Control_ND,  top_pct, na.rm = TRUE)

model.results.Dose <- model.results.Dose %>%
  mutate(
    sig_Low_CAS  = !is.na(pval_adj_CAS_Low)    & pval_adj_CAS_Low    <= quant_Low_CAS,
    sig_High_CAS = !is.na(pval_adj_CAS_High)   & pval_adj_CAS_High   <= quant_High_CAS,
    sig_Low_CLO  = !is.na(pval_adj_CLO_Low)    & pval_adj_CLO_Low    <= quant_Low_CLO,
    sig_High_CLO = !is.na(pval_adj_CLO_High)   & pval_adj_CLO_High   <= quant_High_CLO,
    sig_Control  = !is.na(pval_adj_Control_ND) & pval_adj_Control_ND <= quant_Control,
    any_prioritized = sig_Low_CAS | sig_High_CAS | sig_Low_CLO |
      sig_High_CLO | sig_Control)


# Classify SNPs for one drug's intensity comparison.
#
# Treatment-specific: a candidate in that treatment only, with a rate of change
# significantly different from the other concentration AND from Controls.
# Shared: a candidate in both concentrations, with rates that do not differ
# significantly from each other but do differ from Controls.
#
# Requiring a significant difference from Controls means classified responses
# reflect drug exposure rather than laboratory adaptation. SNPs that are
# candidates somewhere but meet neither definition fall into "Other" and are
# excluded from the overlap comparisons.
classify_dose <- function(df, drug) {
  sig_low   <- paste0("sig_Low_",  drug)
  sig_high  <- paste0("sig_High_", drug)
  p_highlow <- paste0("pval_adj_", drug, "_High-Low")
  p_lowND   <- paste0("pval_adj_", drug, "_Low-ND")
  p_highND  <- paste0("pval_adj_", drug, "_High-ND")
  
  # TRUE when a contrast exists and is below the threshold
  lt <- function(col) !is.na(df[[col]]) & df[[col]] <  pcutoff_pairwise
  ge <- function(col) !is.na(df[[col]]) & df[[col]] >= pcutoff_pairwise
  
  case_when(
    !df$any_prioritized ~ "None",
    
    df[[sig_low]] & !df[[sig_high]] & !df$sig_Control &
      lt(p_highlow) & lt(p_lowND)                     ~ paste0("Low ",  drug, "-specific"),
    
    df[[sig_high]] & !df[[sig_low]] & !df$sig_Control &
      lt(p_highlow) & lt(p_highND)                    ~ paste0("High ", drug, "-specific"),
    
    df$sig_Control & !df[[sig_low]] & !df[[sig_high]] &
      lt(p_lowND) & lt(p_highND)                      ~ "Control-specific",
    
    df[[sig_low]] & df[[sig_high]] & !df$sig_Control &
      ge(p_highlow) & lt(p_lowND) & lt(p_highND)      ~ paste0(drug, " shared"),
    
    TRUE ~ "Other")
}

model.results.Dose <- model.results.Dose %>%
  mutate(sig_class_CAS_Dose = classify_dose(., "CAS"),
         sig_class_CLO_Dose = classify_dose(., "CLO"))

table(model.results.Dose$sig_class_CAS_Dose)
table(model.results.Dose$sig_class_CLO_Dose)


# Single vs combined drug. Dose suffixes are dropped from the column names
# because only the high single-drug treatments were used in this model.
model.results.Complex <- model.results.Complex %>%
  rename_with(~ gsub("_(High|Low|ND)$", "", .x))

quant_CAS       <- quantile(model.results.Complex$pval_adj_CAS,     top_pct, na.rm = TRUE)
quant_CLO       <- quantile(model.results.Complex$pval_adj_CLO,     top_pct, na.rm = TRUE)
quant_CASCLO    <- quantile(model.results.Complex$pval_adj_CASCLO,  top_pct, na.rm = TRUE)
quant_Control_2 <- quantile(model.results.Complex$pval_adj_Control, top_pct, na.rm = TRUE)

model.results.Complex <- model.results.Complex %>%
  mutate(
    sig_CAS     = pval_adj_CAS     <= quant_CAS,
    sig_CLO     = pval_adj_CLO     <= quant_CLO,
    sig_CASCLO  = pval_adj_CASCLO  <= quant_CASCLO,
    sig_Control = pval_adj_Control <= quant_Control_2,
    any_prioritized = sig_CAS | sig_CLO | sig_CASCLO | sig_Control)

# Same logic as classify_dose, extended to three drug treatments.
model.results.Complex <- model.results.Complex %>%
  mutate(sig_class_Complex = case_when(
    !any_prioritized ~ "None",
    
    sig_CAS & !sig_CLO & !sig_CASCLO & !sig_Control &
      `pval_adj_CAS-CLO`     < pcutoff_pairwise &
      `pval_adj_CAS-CASCLO`  < pcutoff_pairwise &
      `pval_adj_CAS-Control` < pcutoff_pairwise          ~ "CAS-specific",
    
    sig_CLO & !sig_CAS & !sig_CASCLO & !sig_Control &
      `pval_adj_CAS-CLO`      < pcutoff_pairwise &
      `pval_adj_CASCLO-CLO`   < pcutoff_pairwise &
      `pval_adj_CLO-Control`  < pcutoff_pairwise         ~ "CLO-specific",
    
    sig_CASCLO & !sig_CAS & !sig_CLO & !sig_Control &
      `pval_adj_CAS-CASCLO`     < pcutoff_pairwise &
      `pval_adj_CASCLO-CLO`     < pcutoff_pairwise &
      `pval_adj_CASCLO-Control` < pcutoff_pairwise       ~ "CASCLO-specific",
    
    sig_Control & !sig_CAS & !sig_CLO & !sig_CASCLO &
      `pval_adj_CAS-Control`    < pcutoff_pairwise &
      `pval_adj_CLO-Control`    < pcutoff_pairwise &
      `pval_adj_CASCLO-Control` < pcutoff_pairwise       ~ "Control-specific",
    
    sig_CAS & sig_CASCLO & !sig_CLO & !sig_Control &
      `pval_adj_CAS-CASCLO`     >= pcutoff_pairwise &
      `pval_adj_CAS-Control`     < pcutoff_pairwise &
      `pval_adj_CASCLO-Control`  < pcutoff_pairwise      ~ "CAS & CASCLO shared",
    
    sig_CLO & sig_CASCLO & !sig_CAS & !sig_Control &
      `pval_adj_CASCLO-CLO`      >= pcutoff_pairwise &
      `pval_adj_CLO-Control`      < pcutoff_pairwise &
      `pval_adj_CASCLO-Control`   < pcutoff_pairwise     ~ "CLO & CASCLO shared",
    
    sig_CAS & sig_CLO & !sig_CASCLO & !sig_Control &
      `pval_adj_CAS-CLO`       >= pcutoff_pairwise &
      `pval_adj_CAS-Control`    < pcutoff_pairwise &
      `pval_adj_CLO-Control`    < pcutoff_pairwise       ~ "CAS & CLO shared",
    
    sig_CAS & sig_CLO & sig_CASCLO & !sig_Control &
      `pval_adj_CAS-CLO`        >= pcutoff_pairwise &
      `pval_adj_CAS-CASCLO`     >= pcutoff_pairwise &
      `pval_adj_CASCLO-CLO`     >= pcutoff_pairwise &
      `pval_adj_CAS-Control`     < pcutoff_pairwise &
      `pval_adj_CLO-Control`     < pcutoff_pairwise &
      `pval_adj_CASCLO-Control`  < pcutoff_pairwise      ~ "CAS & CLO & CASCLO shared",
    
    TRUE ~ "Other"))

table(model.results.Complex$sig_class_Complex)


# ======================= MIRRORED MANHATTAN PLOTS ============================

# Each panel shows two treatments mirrored above and below a central gap that
# carries the chromosome labels. Grey points are the genomic background;
# colored points are classified candidates; larger white-outlined points mark
# shared classes so they stay visible.

YMAX <- 15           # maximum displayed -log10(padj)
GAP  <- 2.5          # half-height of the central chromosome-label gap

BAND_LIGHT <- "#F7F9FB"
BAND_DARK  <- "#EDF1F6"
SIG_COL    <- "firebrick2"
SIG_LWD    <- 0.9
MID_COL    <- "black"

PT_GREY <- setNames(rep(c("gray70", "gray35"), length.out = length(CHROM_ORDER)),
                    CHROM_ORDER)

model.results.Dose$CHROM    <- factor(model.results.Dose$CHROM,    levels = CHROM_ORDER)
model.results.Complex$CHROM <- factor(model.results.Complex$CHROM, levels = CHROM_ORDER)

# Approximate chromosome midpoints, used to place labels in the central gap.
get_axis_set <- function(df) {
  df %>%
    filter(!is.na(CHROM)) %>%
    group_by(CHROM) %>%
    summarize(center = mean(MB), .groups = "drop")
}

axis_set_dose    <- get_axis_set(model.results.Dose)
axis_set_complex <- get_axis_set(model.results.Complex)

# Candidate thresholds on the -log10 scale, for the significance lines
sigline_Low_CAS  <- -log10(quant_Low_CAS)
sigline_High_CAS <- -log10(quant_High_CAS)
sigline_Low_CLO  <- -log10(quant_Low_CLO)
sigline_High_CLO <- -log10(quant_High_CLO)
sigline_CAS      <- -log10(quant_CAS)
sigline_CLO      <- -log10(quant_CLO)
sigline_CASCLO   <- -log10(quant_CASCLO)


mirror_manhattan <- function(df, yup, ydn, class_col, axis_set,
                             up_cols, down_cols,
                             sig_up, sig_dn,
                             tag_up, tag_dn, tag_up_col, tag_dn_col,
                             title, emph_cols = NULL) {
  
  # Shift both halves away from zero to open the chromosome-label gap
  df[[yup]] <- df[[yup]] + GAP
  df[[ydn]] <- df[[ydn]] - GAP
  sig_up_p  <- sig_up + GAP
  sig_dn_p  <- sig_dn + GAP
  
  # Genomic span of each chromosome, used to draw alternating background bands
  spans <- df %>%
    filter(!is.na(CHROM)) %>%
    group_by(CHROM) %>%
    summarize(xmin = min(MB), xmax = max(MB), .groups = "drop") %>%
    mutate(band = as.integer(CHROM) %% 2 == 0)
  
  # Ordinary candidate points
  overlay <- function(cols, ycol) {
    lapply(names(cols), function(cat)
      geom_point(
        data = df %>% filter(.data[[class_col]] == cat, !is.na(.data[[ycol]])),
        aes(y = .data[[ycol]]), color = cols[[cat]], size = 1))
  }
  
  # Shared classes, drawn larger with a white outline so they stand out
  emph <- function(cols, ycol) {
    lapply(names(cols), function(cat)
      geom_point(
        data = df %>% filter(.data[[class_col]] == cat, !is.na(.data[[ycol]])),
        aes(y = .data[[ycol]]), shape = 21, fill = cols[[cat]],
        color = "white", size = 2.4, stroke = 0.4))
  }
  
  p <- ggplot(df, aes(x = MB)) +
    # Alternating chromosome bands, upper and lower halves
    geom_rect(data = filter(spans, !band), inherit.aes = FALSE,
              aes(xmin = xmin, xmax = xmax, ymin = GAP, ymax = YMAX + GAP),
              fill = BAND_LIGHT) +
    geom_rect(data = filter(spans, !band), inherit.aes = FALSE,
              aes(xmin = xmin, xmax = xmax, ymin = -(YMAX + GAP), ymax = -GAP),
              fill = BAND_LIGHT) +
    geom_rect(data = filter(spans, band), inherit.aes = FALSE,
              aes(xmin = xmin, xmax = xmax, ymin = GAP, ymax = YMAX + GAP),
              fill = BAND_DARK) +
    geom_rect(data = filter(spans, band), inherit.aes = FALSE,
              aes(xmin = xmin, xmax = xmax, ymin = -(YMAX + GAP), ymax = -GAP),
              fill = BAND_DARK) +
    # Genomic background
    geom_point(data = df %>% filter(!is.na(.data[[yup]])),
               aes(y = .data[[yup]], color = CHROM), size = 0.55, alpha = 0.55) +
    geom_point(data = df %>% filter(!is.na(.data[[ydn]])),
               aes(y = .data[[ydn]], color = CHROM), size = 0.55, alpha = 0.55) +
    scale_color_manual(values = PT_GREY, guide = "none") +
    overlay(up_cols, yup) +
    overlay(down_cols, ydn)
  
  # Shared classes drawn last so they sit on top
  if (!is.null(emph_cols)) {
    p <- p + emph(emph_cols, yup) + emph(emph_cols, ydn)
  }
  
  p +
    geom_hline(yintercept =  GAP, color = MID_COL, linewidth = 0.3) +
    geom_hline(yintercept = -GAP, color = MID_COL, linewidth = 0.3) +
    geom_text(data = axis_set, inherit.aes = FALSE,
              aes(x = center, y = 0, label = CHROM),
              size = 3.4, color = "grey20") +
    geom_hline(yintercept =  sig_up_p, linewidth = SIG_LWD, color = SIG_COL) +
    geom_hline(yintercept = -sig_dn_p, linewidth = SIG_LWD, color = SIG_COL) +
    # Threshold labels show the original -log10 value, not the shifted position
    annotate("text", x = Inf, y =  sig_up_p, label = sprintf("%.2f", sig_up),
             hjust = 1.05, vjust = -0.6, size = 3, color = SIG_COL) +
    annotate("text", x = Inf, y = -sig_dn_p, label = sprintf("%.2f", sig_dn),
             hjust = 1.05, vjust = 1.4, size = 3, color = SIG_COL) +
    annotate("text", x = -Inf, y = YMAX + GAP - 0.8, label = tag_up,
             hjust = -0.08, vjust = 1, fontface = "bold", size = 7,
             color = tag_up_col) +
    annotate("text", x = -Inf, y = -(YMAX + GAP) + 0.8, label = tag_dn,
             hjust = -0.08, vjust = 0, fontface = "bold", size = 7,
             color = tag_dn_col) +
    scale_x_continuous(name = "Chromosome", expand = c(0.01, 0)) +
    # Axis labels remove the GAP shift so readers see the true -log10 values
    scale_y_continuous(
      name   = bquote(-log[10](italic(p[adj]))),
      breaks = c(-(seq(0, YMAX, 5) + GAP), seq(0, YMAX, 5) + GAP),
      labels = c(seq(0, YMAX, 5), seq(0, YMAX, 5)),
      limits = c(-(YMAX + GAP), YMAX + GAP),
      expand = c(0, 0)) +
    ggtitle(title) +
    theme_classic() +
    theme(legend.position = "none",
          plot.title   = element_text(face = "bold", size = 18),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y  = element_text(size = 14),
          axis.title   = element_text(size = 16))
}


# Drug intensity panels. Low concentration plots upward, high downward.
dose_m <- model.results.Dose %>%
  filter(!is.na(CHROM)) %>%
  mutate(y_CAS_Low  =  -log10(pval_adj_CAS_Low),
         y_CAS_High = -(-log10(pval_adj_CAS_High)),
         y_CLO_Low  =  -log10(pval_adj_CLO_Low),
         y_CLO_High = -(-log10(pval_adj_CLO_High)))

p_CAS_dose <- mirror_manhattan(
  dose_m, "y_CAS_Low", "y_CAS_High", "sig_class_CAS_Dose", axis_set_dose,
  up_cols    = c("Low CAS-specific"  = TREATMENT_COLORS[["Low CAS"]]),
  down_cols  = c("High CAS-specific" = TREATMENT_COLORS[["High CAS"]]),
  sig_up     = sigline_Low_CAS,  sig_dn = sigline_High_CAS,
  tag_up     = "Low CAS",        tag_dn = "High CAS",
  tag_up_col = TREATMENT_COLORS[["Low CAS"]],
  tag_dn_col = TREATMENT_COLORS[["High CAS"]],
  title      = "Low CAS vs High CAS",
  emph_cols  = c("CAS shared" = "dodgerblue"))

p_CLO_dose <- mirror_manhattan(
  dose_m, "y_CLO_Low", "y_CLO_High", "sig_class_CLO_Dose", axis_set_dose,
  up_cols    = c("Low CLO-specific"  = TREATMENT_COLORS[["Low CLO"]]),
  down_cols  = c("High CLO-specific" = TREATMENT_COLORS[["High CLO"]]),
  sig_up     = sigline_Low_CLO,  sig_dn = sigline_High_CLO,
  tag_up     = "Low CLO",        tag_dn = "High CLO",
  tag_up_col = TREATMENT_COLORS[["Low CLO"]],
  tag_dn_col = TREATMENT_COLORS[["High CLO"]],
  title      = "Low CLO vs High CLO",
  emph_cols  = c("CLO shared" = "purple"))


# Single vs combined drug panels. CASCLO plots upward, each single drug
# downward. Labels use "High CAS" and "High CLO" because those are the
# treatments the model used.
complex_m <- model.results.Complex %>%
  filter(!is.na(CHROM)) %>%
  mutate(y_CASCLO =  -log10(pval_adj_CASCLO),
         y_CAS_dn = -(-log10(pval_adj_CAS)),
         y_CLO_dn = -(-log10(pval_adj_CLO)))

p_CAS_CASCLO <- mirror_manhattan(
  complex_m, "y_CASCLO", "y_CAS_dn", "sig_class_Complex", axis_set_complex,
  up_cols    = c("CASCLO-specific" = TREATMENT_COLORS[["CASCLO"]]),
  down_cols  = c("CAS-specific"     = TREATMENT_COLORS[["High CAS"]],
                 "CAS & CLO shared" = "mediumturquoise"),
  sig_up     = sigline_CASCLO, sig_dn = sigline_CAS,
  tag_up     = "CASCLO",       tag_dn = "High CAS",
  tag_up_col = TREATMENT_COLORS[["CASCLO"]],
  tag_dn_col = TREATMENT_COLORS[["High CAS"]],
  title      = "CASCLO vs High CAS",
  emph_cols  = c("CAS & CASCLO shared"       = "forestgreen",
                 "CAS & CLO & CASCLO shared" = "royalblue"))

p_CLO_CASCLO <- mirror_manhattan(
  complex_m, "y_CASCLO", "y_CLO_dn", "sig_class_Complex", axis_set_complex,
  up_cols    = c("CASCLO-specific" = TREATMENT_COLORS[["CASCLO"]]),
  down_cols  = c("CLO-specific"     = TREATMENT_COLORS[["High CLO"]],
                 "CAS & CLO shared" = "mediumturquoise"),
  sig_up     = sigline_CASCLO, sig_dn = sigline_CLO,
  tag_up     = "CASCLO",       tag_dn = "High CLO",
  tag_up_col = TREATMENT_COLORS[["CASCLO"]],
  tag_dn_col = TREATMENT_COLORS[["High CLO"]],
  title      = "CASCLO vs High CLO",
  emph_cols  = c("CLO & CASCLO shared"       = "magenta",
                 "CAS & CLO & CASCLO shared" = "royalblue"))


# ======================= UPSET PLOTS =========================================

# Color one UpSet intersection. Shared by both UpSet sections below.
qcol <- function(sets, col) {
  upset_query(intersect = sets, color = col, fill = col,
              only_components = c("intersections_matrix", "Intersection size"))
}

# Standard UpSet options so both plots look the same.
upset_opts <- list(
  sort_intersections_by = "cardinality",
  base_annotations = list(
    "Intersection size" = intersection_size(counts = TRUE, text = list(size = 3.2)) +
      theme(panel.grid = element_blank(),
            axis.line  = element_line(colour = "black", linewidth = 0.4))),
  matrix = intersection_matrix(geom = geom_point(size = 2.8),
                               segment = geom_segment(linewidth = 1.4)),
  set_sizes = FALSE,
  stripes = "grey96")


# Two-set UpSet comparing concentrations of one drug. Shared candidates belong
# to both sets, which is what produces the intersection bar.
dose_upset <- function(dat, class_col, low_lab, high_lab, shared_lab,
                       set_low, set_high, col_low, col_high, col_shared, title) {
  
  df <- dat %>%
    filter(.data[[class_col]] %in% c(low_lab, high_lab, shared_lab)) %>%
    transmute(!!set_low  := .data[[class_col]] %in% c(low_lab,  shared_lab),
              !!set_high := .data[[class_col]] %in% c(high_lab, shared_lab))
  
  queries <- list(qcol(set_low, col_low), qcol(set_high, col_high))
  
  # Only color the shared intersection if any shared candidates exist
  if (any(df[[set_low]] & df[[set_high]])) {
    queries <- c(queries, list(qcol(c(set_low, set_high), col_shared)))
  }
  
  do.call(upset, c(
    list(df, intersect = c(set_low, set_high), name = "concentration",
         queries = queries, width_ratio = 0.22, height_ratio = 0.55),
    upset_opts)) +
    patchwork::plot_annotation(title = title)
}

upset_CAS <- dose_upset(
  model.results.Dose, "sig_class_CAS_Dose",
  low_lab = "Low CAS-specific", high_lab = "High CAS-specific",
  shared_lab = "CAS shared",
  set_low = "Low CAS", set_high = "High CAS",
  col_low = "#F2B705", col_high = "#8A6D0B", col_shared = "dodgerblue",
  title = "CAS candidate overlap")

upset_CLO <- dose_upset(
  model.results.Dose, "sig_class_CLO_Dose",
  low_lab = "Low CLO-specific", high_lab = "High CLO-specific",
  shared_lab = "CLO shared",
  set_low = "Low CLO", set_high = "High CLO",
  col_low = "#EC5FA3", col_high = "#8E1616", col_shared = "purple",
  title = "CLO candidate overlap")


# Three-set UpSet for the single vs combined drug comparison. Set names use
# "High CAS" and "High CLO" to match the Manhattan panels.
class_to_sets <- tibble::tribble(
  ~sig_class_Complex,          ~`High CAS`, ~`High CLO`, ~CASCLO,
  "CAS-specific",               TRUE,        FALSE,       FALSE,
  "CLO-specific",               FALSE,       TRUE,        FALSE,
  "CASCLO-specific",            FALSE,       FALSE,       TRUE,
  "CAS & CLO shared",           TRUE,        TRUE,        FALSE,
  "CAS & CASCLO shared",        TRUE,        FALSE,       TRUE,
  "CLO & CASCLO shared",        FALSE,       TRUE,        TRUE,
  "CAS & CLO & CASCLO shared",  TRUE,        TRUE,        TRUE)

upset_df <- model.results.Complex %>%
  filter(sig_class_Complex %in% class_to_sets$sig_class_Complex) %>%
  left_join(class_to_sets, by = "sig_class_Complex") %>%
  dplyr::select(`High CAS`, `High CLO`, CASCLO)

# Intersection colors match the classes used in the Manhattan panels
complex_queries <- list(
  qcol("High CAS",                           TREATMENT_COLORS[["High CAS"]]),
  qcol("High CLO",                           TREATMENT_COLORS[["High CLO"]]),
  qcol("CASCLO",                             TREATMENT_COLORS[["CASCLO"]]),
  qcol(c("High CAS", "High CLO"),            "mediumturquoise"),
  qcol(c("High CAS", "CASCLO"),              "forestgreen"),
  qcol(c("High CLO", "CASCLO"),              "magenta"),
  qcol(c("High CAS", "High CLO", "CASCLO"),  "royalblue"))

upset_complex <- do.call(upset, c(
  list(upset_df, intersect = c("High CAS", "High CLO", "CASCLO"),
       name = "treatment", queries = complex_queries,
       width_ratio = 0.18, height_ratio = 0.55),
  upset_opts)) +
  patchwork::plot_annotation(title = "Complexity candidate overlap")


# ======================= MANHATTAN + UPSET COMBINED ==========================

# Rasterizing only the point layers keeps the exported SVGs small while axes,
# text, and UpSet elements stay vector.
man_CAS_r        <- rasterise(p_CAS_dose,   layers = "Point", dpi = 300)
man_CLO_r        <- rasterise(p_CLO_dose,   layers = "Point", dpi = 300)
man_CAS_CASCLO_r <- rasterise(p_CAS_CASCLO, layers = "Point", dpi = 300)
man_CLO_CASCLO_r <- rasterise(p_CLO_CASCLO, layers = "Point", dpi = 300)

# Drug intensity figure: one row per drug, Manhattan beside its UpSet
manplot_upset_fig_dose <-
  ((man_CAS_r | wrap_elements(upset_CAS)) + plot_layout(widths = c(3, 1))) /
  ((man_CLO_r | wrap_elements(upset_CLO)) + plot_layout(widths = c(3, 1)))

# Complexity figure: two Manhattans stacked on the left, one UpSet spanning
# both rows on the right. "#" marks intentionally empty cells.
manplot_upset_design <- "
A#
AC
BC
B#
"

manplot_upset_fig_complex <-
  man_CAS_CASCLO_r + man_CLO_CASCLO_r + wrap_elements(upset_complex) +
  plot_layout(design = manplot_upset_design, widths = c(3, 1))

# ggsave("Figures/Manhattan_and_Upset_Plots/dose_manhattan_upset_raw.svg",
#         manplot_upset_fig_dose, width = 12, height = 7.5, device = svglite::svglite)
# ggsave("Figures/Manhattan_and_Upset_Plots/complex_manhattan_upset_raw.svg",
#        manplot_upset_fig_complex, width = 12, height = 7.5, device = svglite::svglite)


# ======================= BETA COEFFICIENT PLOTS ==============================

# Compare the rate of allele frequency change (beta) across treatments for each
# candidate set. Three panels: CAS intensity, CLO intensity, and single vs
# combined drug.

USE_ASINH <- TRUE   # compress extreme beta values while keeping the sign
SIG       <- 0.1    # smaller values stretch the axis more around zero
Y_ZOOM    <- NULL   # NULL uses the full data range
Y_LAB     <- "Rate of allele frequency change (β)"

# Keep only SNPs measured in every treatment. The Friedman and paired Wilcoxon
# tests require complete blocks.
complete_blocks <- function(df, label = "") {
  n_t      <- n_distinct(df$Treatment)
  n_before <- n_distinct(df$SNP_ID)
  
  keep <- df %>%
    filter(!is.na(beta)) %>%
    count(SNP_ID) %>%
    filter(n == n_t) %>%
    pull(SNP_ID)
  
  out <- df %>% filter(SNP_ID %in% keep)
  
  cat(sprintf("%-22s %5d SNPs -> %5d complete blocks (%d dropped)\n",
              label, n_before, n_distinct(out$SNP_ID),
              n_before - n_distinct(out$SNP_ID)))
  out
}

# multcompLetters() needs pairwise p-values as a named vector, not a matrix.
pmat_to_vec <- function(pmat) {
  v <- c()
  for (r in rownames(pmat)) {
    for (cc in colnames(pmat)) {
      if (!is.na(pmat[r, cc])) v[paste(r, cc, sep = "-")] <- pmat[r, cc]
    }
  }
  v
}

# Rename significance letters so the leftmost treatment gets "a". Only the
# labels change; the statistical groupings are unaffected.
relabel_by_position <- function(lets, treat_lvls) {
  ord  <- treat_lvls[treat_lvls %in% names(lets)]
  seen <- character(0)
  
  for (l in lets[ord]) {
    for (ch in strsplit(l, "")[[1]]) {
      if (!ch %in% seen) seen <- c(seen, ch)
    }
  }
  
  map <- setNames(letters[seq_along(seen)], seen)
  vapply(lets, function(l) paste(map[strsplit(l, "")[[1]]], collapse = ""),
         character(1))
}

# Shared y-axis transformation so all three panels stay comparable
if (USE_ASINH) {
  tr   <- function(x) asinh(x / SIG)
  BRKS <- c(-3, -1, -0.5, 0, 0.5, 1, 3, 6)
} else {
  tr   <- function(x) x
  BRKS <- NULL
}

# Polarize beta to the allele that increased in the treatment where each SNP
# was identified. Afterwards a positive beta in any treatment means that same
# allele increased there.
polarize <- function(df, disc_map) {
  df %>%
    mutate(DiscTrt = disc_map[as.character(CandidateSet)]) %>%
    group_by(CandidateSet, SNP_ID) %>%
    mutate(flip = sign(beta[as.character(Treatment) == DiscTrt][1])) %>%
    ungroup() %>%
    filter(!is.na(flip), flip != 0) %>%
    mutate(beta_raw = beta,
           beta     = beta * flip) %>%
    dplyr::select(-DiscTrt, -flip)
}

# Friedman test plus paired Wilcoxon contrasts for one candidate set.
beta_stats <- function(df, set_name) {
  d <- df %>%
    filter(CandidateSet == set_name) %>%
    complete_blocks(label = set_name) %>%
    arrange(Treatment, SNP_ID)   # paired tests need a consistent order
  
  cat("\n===== ", set_name, " =====\n", sep = "")
  
  ft <- friedman.test(beta ~ Treatment | SNP_ID, data = d)
  print(ft)
  
  pw <- pairwise.wilcox.test(d$beta, d$Treatment, paired = TRUE,
                             p.adjust.method = "holm")
  print(pw)
  
  meds <- d %>%
    group_by(Treatment) %>%
    summarise(median_beta = median(beta),
              mean_beta   = mean(beta),
              q25         = quantile(beta, 0.25),
              q75         = quantile(beta, 0.75),
              n           = n_distinct(SNP_ID),
              .groups = "drop")
  print(as.data.frame(meds), digits = 4)
  
  # Treatments sharing a letter are not significantly different
  lets <- multcompLetters(pmat_to_vec(pw$p.value), threshold = 0.05)$Letters
  lets <- relabel_by_position(lets, levels(d$Treatment))
  
  list(friedman = ft, pairwise = pw, medians = meds,
       letters = tibble(CandidateSet = set_name,
                        Treatment = names(lets), label = unname(lets)),
       data = d)
}

# Grouped boxplot of beta for one comparison, reused for all three panels.
make_grouped_box <- function(res_list, pal, title, treat_lvls, y_rng) {
  
  df <- bind_rows(lapply(res_list, `[[`, "data"))
  
  letters_df <- bind_rows(lapply(res_list, `[[`, "letters")) %>%
    mutate(CandidateSet = factor(CandidateSet, levels = levels(df$CandidateSet)),
           Treatment    = factor(Treatment,    levels = treat_lvls),
           y = tr(y_rng[2]) - 0.03 * diff(tr(y_rng)))
  
  df$beta_t <- tr(df$beta)
  dodge <- position_dodge(width = 0.8, preserve = "single")
  
  p <- ggplot(df, aes(x = CandidateSet, y = beta_t, fill = Treatment)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_boxplot(position = dodge, width = 0.7, alpha = 0.85,
                 color = "black", linewidth = 0.3, coef = Inf) +
    geom_text(data = letters_df, inherit.aes = FALSE,
              aes(x = CandidateSet, y = y, label = label, group = Treatment),
              position = dodge, size = 5) +
    scale_fill_manual(values = pal, breaks = treat_lvls) +
    coord_cartesian(ylim = tr(y_rng)) +
    labs(x = NULL, y = Y_LAB, fill = NULL, title = title) +
    theme_classic(base_size = 11) +
    theme(legend.position = "bottom",
          plot.title   = element_text(size = 12, face = "bold"),
          axis.text.x  = element_text(size = 12),
          axis.title.y = element_text(size = 12))
  
  # Data are transformed, but the axis shows the original beta values
  if (USE_ASINH) p + scale_y_continuous(breaks = tr(BRKS), labels = BRKS) else p
}


# Build the long-format beta table for one comparison.
build_beta_table <- function(df, candidate_defs, beta_cols, treat_levels) {
  
  bind_rows(imap(candidate_defs, function(sig_col, set_name) {
    df %>% filter(.data[[sig_col]]) %>% mutate(CandidateSet = set_name)
  })) %>%
    dplyr::select(SNP_ID, CandidateSet, all_of(names(beta_cols))) %>%
    pivot_longer(all_of(names(beta_cols)),
                 names_to = "Treatment", values_to = "beta") %>%
    mutate(Treatment    = recode(Treatment, !!!beta_cols),
           Treatment    = factor(Treatment, levels = treat_levels),
           CandidateSet = factor(CandidateSet, levels = names(candidate_defs)))
}

# CAS intensity
cas_betas <- build_beta_table(
  model.results.Dose,
  candidate_defs = c("Low CAS candidates"  = "sig_Low_CAS",
                     "High CAS candidates" = "sig_High_CAS"),
  beta_cols = c(beta_Control_ND = "Control",
                beta_CAS_Low    = "Low CAS",
                beta_CAS_High   = "High CAS"),
  treat_levels = c("Control", "Low CAS", "High CAS")) %>%
  polarize(c("Low CAS candidates"  = "Low CAS",
             "High CAS candidates" = "High CAS"))

cas_res <- set_names(
  lapply(levels(cas_betas$CandidateSet), function(s) beta_stats(cas_betas, s)),
  levels(cas_betas$CandidateSet))

# CLO intensity
clo_betas <- build_beta_table(
  model.results.Dose,
  candidate_defs = c("Low CLO candidates"  = "sig_Low_CLO",
                     "High CLO candidates" = "sig_High_CLO"),
  beta_cols = c(beta_Control_ND = "Control",
                beta_CLO_Low    = "Low CLO",
                beta_CLO_High   = "High CLO"),
  treat_levels = c("Control", "Low CLO", "High CLO")) %>%
  polarize(c("Low CLO candidates"  = "Low CLO",
             "High CLO candidates" = "High CLO"))

clo_res <- set_names(
  lapply(levels(clo_betas$CandidateSet), function(s) beta_stats(clo_betas, s)),
  levels(clo_betas$CandidateSet))

# Single vs combined drug. Treatments are labeled High CAS and High CLO here
# because those are the concentrations the complexity model used.
complex_betas <- build_beta_table(
  model.results.Complex,
  candidate_defs = c("CASCLO candidates"   = "sig_CASCLO",
                     "High CAS candidates" = "sig_CAS",
                     "High CLO candidates" = "sig_CLO"),
  beta_cols = c(beta_Control = "Control",
                beta_CAS     = "High CAS",
                beta_CLO     = "High CLO",
                beta_CASCLO  = "CASCLO"),
  treat_levels = c("Control", "High CAS", "High CLO", "CASCLO")) %>%
  polarize(c("CASCLO candidates"   = "CASCLO",
             "High CAS candidates" = "High CAS",
             "High CLO candidates" = "High CLO"))

complex_res <- set_names(
  lapply(levels(complex_betas$CandidateSet),
         function(s) beta_stats(complex_betas, s)),
  levels(complex_betas$CandidateSet))


# One shared y-axis range so all three panels are directly comparable
all_beta <- c(bind_rows(lapply(cas_res,     `[[`, "data"))$beta,
              bind_rows(lapply(clo_res,     `[[`, "data"))$beta,
              bind_rows(lapply(complex_res, `[[`, "data"))$beta)

if (is.null(Y_ZOOM)) {
  global_rng    <- range(all_beta, na.rm = TRUE)
  global_rng[2] <- global_rng[2] * 1.12   # headroom for significance letters
} else {
  global_rng <- Y_ZOOM
}

cat("\ny range:", round(global_rng, 3), " asinh:", USE_ASINH, "\n")

beta_cas_plot <- make_grouped_box(
  cas_res, TREATMENT_COLORS[c("Control", "Low CAS", "High CAS")] %>%
    replace(1, "grey70"),
  "CAS", c("Control", "Low CAS", "High CAS"), global_rng)

beta_clo_plot <- make_grouped_box(
  clo_res, c("Control" = "grey70",
             "Low CLO"  = TREATMENT_COLORS[["Low CLO"]],
             "High CLO" = TREATMENT_COLORS[["High CLO"]]),
  "CLO", c("Control", "Low CLO", "High CLO"), global_rng)

beta_complex_plot <- make_grouped_box(
  complex_res, c("Control"  = "grey70",
                 "High CAS" = TREATMENT_COLORS[["High CAS"]],
                 "High CLO" = TREATMENT_COLORS[["High CLO"]],
                 "CASCLO"   = TREATMENT_COLORS[["CASCLO"]]),
  "Single vs combined drug",
  c("Control", "High CAS", "High CLO", "CASCLO"), global_rng)


# Descriptive check: is the CASCLO response larger than the sum of the two
# single-drug responses at the same loci? Reported in the Results.
cm <- setNames(complex_res[["CASCLO candidates"]]$medians$median_beta,
               complex_res[["CASCLO candidates"]]$medians$Treatment)

cat(sprintf("\nCASCLO median %.4f vs High CAS + High CLO = %.4f  (%.2fx)\n",
            cm["CASCLO"], cm["High CAS"] + cm["High CLO"],
            cm["CASCLO"] / (cm["High CAS"] + cm["High CLO"])))

# Only the left panel keeps the y-axis title so the label is not repeated
combined_beta_plots <-
  (beta_cas_plot + (beta_clo_plot + theme(axis.title.y = element_blank())) +
     plot_layout(widths = c(1, 1))) /
  beta_complex_plot +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(text = element_text(size = 14))

# ggsave("Figures/Beta_Coefficient_Plot/beta_comparison_plot_raw.svg",
#         combined_beta_plots, width = 9, height = 8, units = "in", dpi = 300)


# ======================= CANDIDATE SNP ANNOTATION ============================

# Build a VCF from the filtered SNP table for SnpEff. Chromosome names must
# match the GFF and genome files used to build the SnpEff database.
snps_for_snpeff <- snps_filtered_poly %>%
  left_join(snps %>% dplyr::select(CHROM, POS, REF, ALT), by = c("CHROM", "POS")) %>%
  mutate(CHROM = if_else(CHROM %in% CHROM_ORDER, paste0("chr", CHROM), CHROM),
         CHROM = gsub("mitochondrion", "chrMito", CHROM, ignore.case = TRUE),
         ID = ".", QUAL = ".", FILTER = ".", INFO = ".") %>%
  dplyr::select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO)

# write.table(snps_for_snpeff, "Output_files/snps.vcf",
#             sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


# Parse the SnpEff output. The ANN field holds one or more pipe-delimited
# annotations per variant, separated by commas.
ANN_FIELDS <- c("Allele", "Effect", "Impact", "Gene_Name", "Gene_ID",
                "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank",
                "HGVS_c", "HGVS_p", "cDNA_pos", "CDS_pos", "Protein_pos",
                "Distance", "Errors")

ann_long <- snps_ann %>%
  dplyr::select(CHROM = 1, POS = 2, REF = 4, ALT = 5, ANN = 8) %>%
  mutate(ANN = gsub("^ANN=", "", ANN),
         ANN = strsplit(ANN, ",")) %>%
  tidyr::unnest(ANN) %>%
  separate(ANN, into = ANN_FIELDS, sep = "\\|",
           fill = "right", extra = "merge") %>%
  filter(CHROM != "chrMito") %>%
  mutate(CHROM  = sub("^chr", "", CHROM),
         SNP_ID = paste(CHROM, POS, sep = "_")) %>%
  relocate(SNP_ID, .after = POS)


# SnpEff reports only systematic names, so add the standard gene name and
# description. Protein-coding genes resolve through ORF; snoRNAs, tRNAs, and
# other noncoding features resolve through ALIAS.
add_gene_annotations <- function(df, orf_col = "Gene_Name") {
  keys <- df[[orf_col]]
  
  lookup <- function(column) {
    out <- mapIds(org.Sc.sgd.db, keys = keys, column = column,
                  keytype = "ORF", multiVals = "first")
    missing <- is.na(out)
    if (any(missing)) {
      out[missing] <- mapIds(org.Sc.sgd.db, keys = keys[missing],
                             column = column, keytype = "ALIAS",
                             multiVals = "first")
    }
    out
  }
  
  df$Standard_Gene_Name <- lookup("GENENAME")
  df$Description        <- lookup("DESCRIPTION")
  df
}

ann_long <- add_gene_annotations(ann_long)
rm(snps_ann)


# Classify early-fixed sites the same way as the model-derived candidates,
# based on which treatments each site fixed in.
classify_fixed <- function(fixed_df, treatments_keep, labels) {
  label_vec <- unlist(labels)   # named character vector, NA for missing keys
  
  fixed_df %>%
    filter(Treatment %in% treatments_keep) %>%
    distinct(SNP_ID, Treatment) %>%
    group_by(SNP_ID) %>%
    summarise(key = paste(sort(Treatment), collapse = "|"), .groups = "drop") %>%
    mutate(Category = unname(label_vec[key])) %>%
    dplyr::select(SNP_ID, Category)
}

# CAS intensity
fixed_CAS <- classify_fixed(
  fixed_snps_filt,
  treatments_keep = c("Low CAS", "High CAS"),
  labels = list("High CAS"          = "High CAS-specific",
                "Low CAS"           = "Low CAS-specific",
                "High CAS|Low CAS"  = "CAS shared")) %>%
  dplyr::rename(Fixed_Category_CAS_Dose = Category)

# Single vs combined drug
fixed_Complex <- classify_fixed(
  fixed_snps_filt,
  treatments_keep = c("High CAS", "High CLO", "CASCLO"),
  labels = list("High CAS"                     = "CAS-specific",
                "High CLO"                     = "CLO-specific",
                "CASCLO"                       = "CASCLO-specific",
                "High CAS|High CLO"            = "CAS & CLO shared",
                "CASCLO|High CAS"              = "CAS & CASCLO shared",
                "CASCLO|High CLO"              = "CLO & CASCLO shared",
                "CASCLO|High CAS|High CLO"     = "CAS & CLO & CASCLO shared")) %>%
  dplyr::rename(Fixed_Category_Complex = Category)

# Attach every classification to the annotation table, then keep the
# highest-impact annotation per SNP. SnpEff already orders annotations by
# impact (High > Moderate > Low > Modifier).
ann_long_merged_filt <- ann_long %>%
  left_join(model.results.Dose %>%
              dplyr::select(SNP_ID, sig_class_CAS_Dose, sig_class_CLO_Dose),
            by = "SNP_ID") %>%
  left_join(model.results.Complex %>%
              dplyr::select(SNP_ID, sig_class_Complex), by = "SNP_ID") %>%
  left_join(fixed_CAS,     by = "SNP_ID") %>%
  left_join(fixed_Complex, by = "SNP_ID") %>%
  group_by(SNP_ID) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  arrange(CHROM, POS)


# Combine model-derived and early-fixed candidates into one table per analysis,
# recording which source each came from.
UNINFORMATIVE <- c("None", "Other", "Control-specific")

build_candidate_table <- function(model_col, fixed_col) {
  candidates <- bind_rows(
    ann_long_merged_filt %>%
      filter(!is.na(.data[[model_col]]), !.data[[model_col]] %in% UNINFORMATIVE) %>%
      transmute(SNP_ID, Group = .data[[model_col]], Source = "Model-derived"),
    ann_long_merged_filt %>%
      filter(!is.na(.data[[fixed_col]])) %>%
      transmute(SNP_ID, Group = .data[[fixed_col]], Source = "Fixed")) %>%
    distinct()
  
  ann_long_merged_filt %>% inner_join(candidates, by = "SNP_ID")
}

CAS_variants     <- build_candidate_table("sig_class_CAS_Dose", "Fixed_Category_CAS_Dose")
complex_variants <- build_candidate_table("sig_class_Complex",  "Fixed_Category_Complex")


# Summary tables of variant impact and effect per candidate group (Tables S2, S3)
variant_summary_funx <- function(df) {
  
  pct_table <- function(col, section) {
    df %>%
      group_by(Group, .data[[col]]) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(Group) %>%
      mutate(label = sprintf("%d (%.1f%%)", n, 100 * n / sum(n))) %>%
      ungroup() %>%
      dplyr::select(-n) %>%
      pivot_wider(names_from = Group, values_from = label,
                  values_fill = "0 (0.0%)") %>%
      mutate(Section = section, Feature = .data[[col]]) %>%
      dplyr::select(Section, Feature, everything(), -all_of(col))
  }
  
  count_table <- function(distinct_cols, name_from, feature) {
    df %>%
      distinct(across(all_of(distinct_cols))) %>%
      count(.data[[name_from]]) %>%
      pivot_wider(names_from = all_of(name_from), values_from = n, values_fill = 0) %>%
      mutate(Section = "Summary", Feature = feature) %>%
      mutate(across(-c(Section, Feature), as.character))
  }
  
  bind_rows(
    pct_table("Impact", "Effect by impact"),
    pct_table("Effect", "Effect by type"),
    count_table(c("SNP_ID", "Group"),    "Group",  "TOTAL SNPs"),
    count_table(c("Gene_Name", "Group"), "Group",  "UNIQUE GENES"),
    count_table(c("SNP_ID", "Source"),   "Source", "Candidate source"))
}

variant_summary_CAS     <- variant_summary_funx(CAS_variants)
variant_summary_complex <- variant_summary_funx(complex_variants)


# Collapse candidate SNPs to the gene level. A gene carrying variants from more
# than one category is assigned to Shared, since gene-level grouping is based
# on which categories its variants belong to rather than on the SNP-level class.
build_gene_summary <- function(model_col, fixed_col, category_rules) {
  
  ann_long_merged_filt %>%
    filter((!is.na(.data[[model_col]]) & !.data[[model_col]] %in% UNINFORMATIVE) |
             (!is.na(.data[[fixed_col]]) & !.data[[fixed_col]] %in% UNINFORMATIVE)) %>%
    group_by(Standard_Gene_Name) %>%
    summarise(
      Systematic_Gene_Name    = dplyr::first(Gene_Name),
      CHROM                   = dplyr::first(CHROM),
      Candidate_SNP_Positions = paste0(min(POS), "–", max(POS)),
      SNP_Effects             = paste(sort(unique(Effect)), collapse = "; "),
      SNP_Impacts             = paste(sort(unique(Impact)), collapse = "; "),
      `SNP-level_Categories`  = paste(sort(unique(na.omit(
        c(.data[[model_col]], .data[[fixed_col]])))), collapse = "; "),
      .groups = "drop") %>%
    distinct(Standard_Gene_Name, .keep_all = TRUE) %>%
    mutate(`Gene-level_Categories` = category_rules(`SNP-level_Categories`)) %>%
    arrange(Standard_Gene_Name)
}

# CAS intensity: Low, High, or Shared
CAS_Dose_gene_level_summary <- build_gene_summary(
  "sig_class_CAS_Dose", "Fixed_Category_CAS_Dose",
  function(cats) case_when(
    grepl("CAS shared", cats) |
      (grepl("High CAS-specific", cats) & grepl("Low CAS-specific", cats)) ~ "Shared",
    grepl("High CAS-specific", cats) &
      !grepl("Low CAS-specific", cats) & !grepl("CAS shared", cats)        ~ "High",
    grepl("Low CAS-specific", cats) &
      !grepl("High CAS-specific", cats) & !grepl("CAS shared", cats)       ~ "Low",
    TRUE ~ NA_character_))

# Single vs combined drug: CAS, CLO, CASCLO, or Shared
Complex_gene_level_summary <- build_gene_summary(
  "sig_class_Complex", "Fixed_Category_Complex",
  function(cats) case_when(
    grepl("shared", cats, ignore.case = TRUE) |
      (grepl("\\bCAS-specific\\b", cats) & grepl("\\bCLO-specific\\b", cats)) |
      (grepl("\\bCAS-specific\\b", cats) & grepl("\\bCASCLO-specific\\b", cats)) |
      (grepl("\\bCLO-specific\\b", cats) & grepl("\\bCASCLO-specific\\b", cats)) ~ "Shared",
    grepl("\\bCAS-specific\\b", cats)    ~ "CAS",
    grepl("\\bCLO-specific\\b", cats)    ~ "CLO",
    grepl("\\bCASCLO-specific\\b", cats) ~ "CASCLO",
    TRUE ~ NA_character_))


# ======================= FUNGAMR COMPARISON ==================================

# Keep S. cerevisiae records with a non-negative confidence score. Negative
# scores mark records that argue against a resistance link.
fungamr_genes <- fungamr_db %>%
  filter(species == "Saccharomyces cerevisiae", confidence.score >= 0) %>%
  separate_rows(gene.or.protein.name, sep = ",\\s*") %>%
  mutate(gene = toupper(trimws(gene.or.protein.name))) %>%
  distinct(gene) %>%
  pull(gene)

flag_fungamr <- function(gene_level_summary, reference = fungamr_genes) {
  gene_level_summary %>% mutate(FungAMR_Hit = Standard_Gene_Name %in% reference)
}

CAS_Dose_gene_level_summary <- flag_fungamr(CAS_Dose_gene_level_summary)
Complex_gene_level_summary  <- flag_fungamr(Complex_gene_level_summary)

CAS_Dose_gene_level_summary %>% count(FungAMR_Hit)
Complex_gene_level_summary  %>% count(FungAMR_Hit)


# Test whether any one candidate category is disproportionately made up of
# known resistance genes. Each category is tested against all other candidate
# categories pooled, so the odds ratio reads as enrichment relative to the rest
# of the candidate set rather than relative to the genome.
fungamr_enrichment <- function(gene_level_summary,
                               category_col = "Gene-level_Categories") {
  n_all  <- nrow(gene_level_summary)
  n_hits <- sum(gene_level_summary$FungAMR_Hit)
  
  gene_level_summary %>%
    group_by(.data[[category_col]]) %>%
    summarise(hits     = sum(FungAMR_Hit),
              non_hits = sum(!FungAMR_Hit),
              total    = n(),
              .groups  = "drop") %>%
    rowwise() %>%
    mutate(bg_hits     = n_hits - hits,
           bg_non_hits = (n_all - total) - (n_hits - hits),
           ft = list(fisher.test(matrix(c(hits, non_hits, bg_hits, bg_non_hits),
                                        nrow = 2))),
           fisher_p   = ft$p.value,
           odds_ratio = unname(ft$estimate)) %>%
    ungroup() %>%
    dplyr::select(-ft) %>%
    mutate(fisher_p_adj = p.adjust(fisher_p, method = "BH"))
}

CAS_dose_fungamr_fisher_test <- fungamr_enrichment(CAS_Dose_gene_level_summary)
Complex_fungamr_fisher_test  <- fungamr_enrichment(Complex_gene_level_summary)


# Antifungal class for each FungAMR hit. These assignments were made by hand
# from the FungAMR drug annotation together with the supporting literature,
# since several genes are linked to more than one drug class.
fungamr_cat_levels <- c("Echinocandin", "Azole", "Glutarimide", "Quinolone",
                        "Antimetabolite", "Serine protease inhibitor", "Other")

fungamr_cat_colors <- c("Echinocandin"              = "gold3",
                        "Azole"                     = "firebrick",
                        "Glutarimide"               = "lightseagreen",
                        "Quinolone"                 = "#984EA3",
                        "Antimetabolite"            = "#FF7F00",
                        "Serine protease inhibitor" = "#F781BF",
                        "Other"                     = "grey60")

fungamr_gene_categories <- tribble(
  ~Standard_Gene_Name, ~Resistance_Category,
  "FKS1",  "Echinocandin",                # caspofungin target
  "CSF1",  "Azole",
  "DHH1",  "Azole",
  "ERG7",  "Azole",
  "ERG25", "Azole",
  "HYM1",  "Azole",
  "PDR5",  "Azole",
  "SUR1",  "Azole",
  "TAO3",  "Azole",
  "ACT1",  "Glutarimide",
  "AFT1",  "Glutarimide",
  "BCK1",  "Glutarimide",
  "BCY1",  "Glutarimide",
  "BUL1",  "Glutarimide",
  "ERG20", "Glutarimide",
  "PMS1",  "Glutarimide",
  "TOP2",  "Glutarimide",
  "YRR1",  "Glutarimide",
  "FKS3",  "Quinolone",
  "PDR10", "Quinolone",
  "PDR11", "Quinolone",
  "PDR12", "Quinolone",
  "PDR15", "Quinolone",
  "GFA1",  "Antimetabolite",
  "PTK2",  "Serine protease inhibitor",
  "SKY1",  "Serine protease inhibitor",
  "POL2",  "Other",
  "TSC11", "Other",
  "UTP18", "Other",
  "VMA1",  "Other") %>%
  mutate(Resistance_Category = factor(Resistance_Category,
                                      levels = fungamr_cat_levels))


# ======================= GO ENRICHMENT =======================================

# topGO needs gene-to-GO mapping; the org.Sc.sgd.db object stores GO-to-gene.
ORF2GO <- inverseList(as.list(org.Sc.sgdGO2ORF))

# Gene universe: every gene carrying at least one SNP that maps to a valid ORF
ann_collapsed <- ann_long %>%
  filter(!is.na(Gene_Name),
         Gene_Name %in% keys(org.Sc.sgd.db, keytype = "ORF")) %>%
  group_by(SNP_ID) %>%
  slice_head(n = 1) %>%
  ungroup()

all_genes_uni_GO <- unique(ann_collapsed$Gene_Name)   # 5,940 genes

get_candidate_genes <- function(gene_level_summary, category_col, category_vals) {
  gene_level_summary %>%
    filter(!is.na(Systematic_Gene_Name),
           .data[[category_col]] %in% category_vals,
           Systematic_Gene_Name %in% keys(org.Sc.sgd.db, keytype = "ORF")) %>%
    pull(Systematic_Gene_Name)
}

candidate_gene_lists <- list(
  Low_CAS        = get_candidate_genes(CAS_Dose_gene_level_summary, "Gene-level_Categories", "Low"),
  High_CAS       = get_candidate_genes(CAS_Dose_gene_level_summary, "Gene-level_Categories", "High"),
  Shared_CAS     = get_candidate_genes(CAS_Dose_gene_level_summary, "Gene-level_Categories", "Shared"),
  CAS_complex    = get_candidate_genes(Complex_gene_level_summary,  "Gene-level_Categories", "CAS"),
  CLO_complex    = get_candidate_genes(Complex_gene_level_summary,  "Gene-level_Categories", "CLO"),
  CASCLO_complex = get_candidate_genes(Complex_gene_level_summary,  "Gene-level_Categories", "CASCLO"),
  Shared_complex = get_candidate_genes(Complex_gene_level_summary,  "Gene-level_Categories", "Shared"))

lengths(candidate_gene_lists)

# topGO expects a named 0/1 vector over the whole gene universe
gene_scores <- lapply(candidate_gene_lists, function(candidates)
  setNames(as.integer(all_genes_uni_GO %in% candidates), all_genes_uni_GO))

vapply(gene_scores, sum, integer(1))   # sanity check against lengths() above


# The weight01 algorithm accounts for the GO hierarchy by down-weighting
# redundant terms, so its p-values are nominal and not FDR-adjusted.
run_topGO <- function(gene_scores, ORF2GO, ontology = "BP",
                      algorithm = "weight01", stat = "fisher",
                      min_node_size = 10, p_cutoff = 0.05) {
  
  GOdata <- new("topGOdata",
                ontology = ontology,
                allGenes = gene_scores,
                geneSel  = function(x) x == 1,
                annot    = annFUN.gene2GO,
                gene2GO  = ORF2GO,
                nodeSize = min_node_size)
  
  resultFisher <- runTest(GOdata, algorithm = algorithm, statistic = stat)
  
  args <- list(GOdata, resultFisher, orderBy = "pvalue",
               topNodes = length(score(resultFisher)), numChar = 1000)
  names(args)[2] <- "pvalue"
  
  results_table <- do.call(GenTable, args) %>%
    as_tibble() %>%
    mutate(pvalue    = as.numeric(pvalue),
           GeneRatio = Significant / Annotated) %>%
    filter(pvalue < p_cutoff) %>%
    arrange(pvalue)
  
  list(GOdata = GOdata, results = results_table)
}

GO_results <- lapply(gene_scores, run_topGO, ORF2GO = ORF2GO)
vapply(GO_results, function(x) nrow(x$results), integer(1))


# Collapse semantically similar terms so each group is represented by its most
# significant member.
reduce_GO_terms <- function(go_results, orgdb = org.Sc.sgd.db,
                            ontology = "BP", threshold = 0.9, measure = "Wang") {
  
  sim_matrix <- calculateSimMatrix(go_results$GO.ID, orgdb = orgdb,
                                   ont = ontology, method = measure)
  
  reduced_terms <- reduceSimMatrix(
    sim_matrix,
    scores    = setNames(-log10(go_results$pvalue), go_results$GO.ID),
    threshold = threshold,
    orgdb     = orgdb)
  
  go_results %>%
    left_join(reduced_terms %>% dplyr::select(go, parentTerm),
              by = c("GO.ID" = "go")) %>%
    arrange(pvalue)
}

GO_reduced <- lapply(GO_results, function(x) reduce_GO_terms(x$results))
vapply(GO_reduced, function(x) n_distinct(x$parentTerm), integer(1))


# Tag each gene set with its display category and keep the best term per parent
collapse_GO <- function(reduced_list, mapping, level_order) {
  tagged <- Map(function(nm, cat) mutate(reduced_list[[nm]], Category = cat),
                names(mapping), unname(mapping))
  
  bind_rows(tagged) %>%
    mutate(Category = factor(Category, levels = level_order)) %>%
    filter(!is.na(parentTerm)) %>%
    group_by(Category, parentTerm) %>%
    slice_min(pvalue, n = 1, with_ties = TRUE) %>%
    ungroup()
}

CAS_dose_GO_collapsed <- collapse_GO(
  GO_reduced,
  c(Low_CAS = "Low", High_CAS = "High", Shared_CAS = "Shared"),
  c("Low", "High", "Shared"))

Complex_GO_collapsed <- collapse_GO(
  GO_reduced,
  c(CAS_complex = "CAS", CLO_complex = "CLO",
    CASCLO_complex = "CASCLO", Shared_complex = "Shared"),
  c("CAS", "CLO", "CASCLO", "Shared"))


# Which candidate genes sit under the drug-related GO terms discussed in the
# text. "sterol" also matches ergosterol and sterol import.
genes_in_GO_term <- function(go_obj, term_id, candidates) {
  AnnotationDbi::select(
    org.Sc.sgd.db,
    keys    = intersect(genesInTerm(go_obj$GOdata, term_id)[[1]], candidates),
    columns = c("ORF", "GENENAME"),
    keytype = "ORF")
}

drug_pattern <- paste("fungal-type cell wall organization", "glucan biosynthetic",
                      "ergosterol", "sterol", "xenobiotic", "toxic",
                      "phospholipid transport", sep = "|")

drug_term_genes <- imap_dfr(GO_results, function(res, set) {
  hits <- res$results %>% filter(grepl(drug_pattern, Term, ignore.case = TRUE))
  if (nrow(hits) == 0) return(NULL)
  
  map_dfr(seq_len(nrow(hits)), function(i) {
    id    <- hits$GO.ID[i]
    genes <- genes_in_GO_term(res, id, candidate_gene_lists[[set]])
    
    tibble(Gene_Set   = set,
           GO.ID      = id,
           Term       = hits$Term[i],
           pvalue     = hits$pvalue[i],
           parentTerm = GO_reduced[[set]]$parentTerm[GO_reduced[[set]]$GO.ID == id],
           n_genes    = nrow(genes),
           genes      = paste(na.omit(genes$GENENAME), collapse = ", "))
  })
})

drug_term_genes %>% dplyr::select(Gene_Set, Term, genes) %>% as.data.frame()


# ======================= FUNGAMR AND GO PANELS ===============================

# Both panels share one x-coordinate system so their columns line up when the
# figure is assembled. Complexity columns are labeled High CAS and High CLO
# because the complexity model used the high single-drug treatments.
STRENGTH_COLS <- c("Low CAS", "High CAS", "Shared")
COMPLEX_COLS  <- c("High CAS", "High CLO", "CASCLO", "Shared")
GROUP_GAP     <- 0.8

make_col_map <- function(strength_cols = STRENGTH_COLS,
                         complex_cols  = COMPLEX_COLS,
                         group_gap     = GROUP_GAP) {
  tibble(
    Analysis = c(rep("Drug strength", length(strength_cols)),
                 rep("Complexity",    length(complex_cols))),
    Column   = c(strength_cols, complex_cols),
    xpos     = c(seq_along(strength_cols),
                 length(strength_cols) + group_gap + seq_along(complex_cols)))
}

# Boundaries of each analysis group, for the grey header strips
make_strips <- function(col_map) {
  col_map %>%
    group_by(Analysis) %>%
    summarise(xmin = min(xpos) - 0.5,
              xmax = max(xpos) + 0.5, .groups = "drop") %>%
    mutate(xmid = (xmin + xmax) / 2)
}


# Panel A: FungAMR candidate genes across both analyses.
#   filled cell = candidate in that column, colored by drug class
#   white cell  = candidate in that analysis but not that column
#   grey cell   = not a candidate in that analysis
merged_fungamr_panel <- function(dose_summary, complex_summary,
                                 lookup     = fungamr_gene_categories,
                                 cat_colors = fungamr_cat_colors,
                                 col_map    = make_col_map(),
                                 lab_angle  = 40) {
  
  dose_hits <- dose_summary %>%
    filter(FungAMR_Hit) %>%
    transmute(Standard_Gene_Name,
              Analysis = "Drug strength",
              Column = case_when(`Gene-level_Categories` == "Low"  ~ "Low CAS",
                                 `Gene-level_Categories` == "High" ~ "High CAS",
                                 TRUE ~ `Gene-level_Categories`))
  
  complex_hits <- complex_summary %>%
    filter(FungAMR_Hit) %>%
    transmute(Standard_Gene_Name,
              Analysis = "Complexity",
              Column   = recode(`Gene-level_Categories`, !!!COMPLEX_LABELS))
  
  hits <- bind_rows(dose_hits, complex_hits) %>% mutate(present = TRUE)
  
  # Track which genes were candidates in each analysis so "not tested" can be
  # distinguished from "tested but absent from this column"
  in_dose    <- unique(dose_hits$Standard_Gene_Name)
  in_complex <- unique(complex_hits$Standard_Gene_Name)
  genes      <- union(in_dose, in_complex)
  
  missing <- setdiff(genes, lookup$Standard_Gene_Name)
  if (length(missing))
    warning("Not in fungamr_gene_categories: ", paste(missing, collapse = ", "))
  
  df <- expand_grid(Standard_Gene_Name = genes, col_map) %>%
    left_join(hits,   by = c("Standard_Gene_Name", "Analysis", "Column")) %>%
    left_join(lookup, by = "Standard_Gene_Name") %>%
    mutate(
      tested = (Analysis == "Drug strength" & Standard_Gene_Name %in% in_dose) |
        (Analysis == "Complexity"    & Standard_Gene_Name %in% in_complex),
      fill_col = case_when(
        !tested         ~ "not a candidate in this analysis",
        !is.na(present) ~ as.character(Resistance_Category),
        TRUE            ~ "absent"))
  
  # Order genes by drug class, then alphabetically within class
  gene_order <- lookup %>%
    filter(Standard_Gene_Name %in% genes) %>%
    arrange(Resistance_Category, Standard_Gene_Name) %>%
    pull(Standard_Gene_Name)
  n <- length(gene_order)
  
  # ggplot places the first discrete level at the bottom, so reverse the order
  df <- df %>%
    mutate(Standard_Gene_Name = factor(Standard_Gene_Name, levels = rev(gene_order)))
  
  # Horizontal separators between drug classes, converted to reversed y-coords
  cat_vec <- lookup$Resistance_Category[match(gene_order, lookup$Standard_Gene_Name)]
  changes <- which(cat_vec[-1] != cat_vec[-length(cat_vec)])
  seps    <- n - changes + 0.5
  
  strips  <- make_strips(col_map)
  xlo     <- min(col_map$xpos) - 0.5
  xhi     <- max(col_map$xpos) + 0.5
  y_lab   <- n + 0.75
  y_strip <- c(n + 2.4, n + 3.2)
  
  ggplot(df) +
    geom_tile(aes(x = xpos, y = Standard_Gene_Name, fill = fill_col),
              width = 1, height = 1, color = "grey55", linewidth = 0.3) +
    annotate("segment", x = xlo, xend = xhi, y = seps, yend = seps,
             color = "black", linewidth = 0.8) +
    annotate("text", x = col_map$xpos, y = y_lab, label = col_map$Column,
             angle = lab_angle, hjust = 0, vjust = 0, size = 4, fontface = "bold") +
    annotate("rect", xmin = strips$xmin, xmax = strips$xmax,
             ymin = y_strip[1], ymax = y_strip[2], fill = "grey88") +
    annotate("text", x = strips$xmid, y = mean(y_strip), label = strips$Analysis,
             size = 4.6, fontface = "bold") +
    scale_fill_manual(
      values = c(cat_colors, "absent" = "white",
                 "not a candidate in this analysis" = "grey93"),
      breaks = c(names(cat_colors), "not a candidate in this analysis"),
      name   = "FungAMR drug class") +
    scale_x_continuous(expand = expansion(add = 0.1)) +
    scale_y_discrete(expand = expansion(add = c(0.5, 4.2))) +
    theme_classic() +
    theme(axis.text.x     = element_blank(),
          axis.ticks      = element_blank(),
          axis.text.y     = element_text(size = 12, face = "italic"),
          axis.title      = element_blank(),
          axis.line       = element_blank(),
          legend.position = "bottom",
          legend.title    = element_text(size = 12, face = "bold"),
          legend.text     = element_text(size = 11)) +
    guides(fill = guide_legend(nrow = 3, title.position = "top"))
}

fungamr_panel <- merged_fungamr_panel(CAS_Dose_gene_level_summary,
                                      Complex_gene_level_summary)


# Panel C: enriched GO terms across both analyses as a bubble grid.
merged_GO_panel <- function(dose_GO, complex_GO,
                            col_map       = make_col_map(),
                            top_n_per_col = NULL,
                            keep_terms    = NULL,
                            lab_angle     = 40,
                            wrap_width    = NULL,
                            term_relabel  = NULL,
                            order_by      = c("pvalue", "alpha"),
                            term_size     = 11) {
  
  order_by <- match.arg(order_by)
  
  dose <- dose_GO %>%
    filter(!is.na(parentTerm)) %>%
    transmute(parentTerm, pvalue, GeneRatio,
              Analysis = "Drug strength",
              Column = case_when(Category == "Low"  ~ "Low CAS",
                                 Category == "High" ~ "High CAS",
                                 TRUE ~ as.character(Category)))
  
  cplx <- complex_GO %>%
    filter(!is.na(parentTerm)) %>%
    transmute(parentTerm, pvalue, GeneRatio,
              Analysis = "Complexity",
              Column   = recode(as.character(Category), !!!COMPLEX_LABELS))
  
  # A term can appear more than once per column; keep its best p-value
  df <- bind_rows(dose, cplx) %>%
    group_by(Analysis, Column, parentTerm) %>%
    slice_min(pvalue, n = 1, with_ties = TRUE) %>%
    ungroup()
  
  if (!is.null(top_n_per_col)) {
    keep <- df %>%
      group_by(Analysis, Column) %>%
      slice_min(pvalue, n = top_n_per_col, with_ties = TRUE) %>%
      ungroup() %>%
      distinct(parentTerm) %>%
      pull(parentTerm)
    
    keep <- union(keep, keep_terms)   # always retain manually chosen terms
    df <- df %>% filter(parentTerm %in% keep)
  }
  
  df <- df %>% left_join(col_map, by = c("Analysis", "Column"))
  
  # Shorten verbose rrvgo parent terms before ordering and wrapping
  if (!is.null(term_relabel)) {
    df <- df %>%
      mutate(parentTerm = ifelse(parentTerm %in% names(term_relabel),
                                 unname(term_relabel[parentTerm]), parentTerm))
  }
  
  # ggplot builds discrete axes bottom-up, so reverse to put the strongest
  # term at the top of the figure
  term_order <- if (order_by == "alpha") {
    sort(unique(df$parentTerm), decreasing = TRUE)
  } else {
    df %>%
      group_by(parentTerm) %>%
      summarise(best = min(pvalue), .groups = "drop") %>%
      arrange(desc(best)) %>%
      pull(parentTerm)
  }
  
  nt <- length(term_order)
  
  df <- if (!is.null(wrap_width)) {
    df %>% mutate(parentTerm = factor(str_wrap(parentTerm, wrap_width),
                                      levels = str_wrap(term_order, wrap_width)))
  } else {
    df %>% mutate(parentTerm = factor(parentTerm, levels = term_order))
  }
  
  strips  <- make_strips(col_map)
  y_lab   <- nt + 0.75
  y_strip <- c(nt + 2.4, nt + 3.2)
  
  # Alternating background rows make a long term list easier to read across
  bands <- do.call(rbind, lapply(seq(nt, 1, by = -2), function(y)
    data.frame(xmin = strips$xmin, xmax = strips$xmax,
               ymin = y - 0.5,     ymax = y + 0.5)))
  
  ggplot(df) +
    geom_rect(data = bands, inherit.aes = FALSE,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "grey95") +
    geom_point(aes(x = xpos, y = parentTerm,
                   size = GeneRatio, color = -log10(pvalue)), alpha = 0.9) +
    annotate("text", x = col_map$xpos, y = y_lab, label = col_map$Column,
             angle = lab_angle, hjust = 0, vjust = 0, size = 4, fontface = "bold") +
    annotate("rect", xmin = strips$xmin, xmax = strips$xmax,
             ymin = y_strip[1], ymax = y_strip[2], fill = "grey88") +
    annotate("text", x = strips$xmid, y = mean(y_strip), label = strips$Analysis,
             size = 4.6, fontface = "bold") +
    scale_color_gradient(low = "#A3C9E3", high = "#0B3D6B",
                         name = expression(-log[10](p))) +
    scale_size_continuous(name = "Gene Ratio", range = c(1, 8)) +
    scale_x_continuous(expand = expansion(add = 0.5)) +
    scale_y_discrete(position = "right", expand = expansion(add = c(0.6, 4.2))) +
    theme_classic() +
    theme(axis.text.x     = element_blank(),
          axis.ticks      = element_blank(),
          axis.text.y     = element_text(size = term_size, lineheight = 0.8),
          axis.title      = element_blank(),
          axis.line       = element_blank(),
          legend.position = "bottom",
          legend.box      = "horizontal",
          legend.text     = element_text(size = 11),
          legend.title    = element_text(size = 12, face = "bold")) +
    guides(color = guide_colorbar(order = 1, title.position = "top",
                                  title.hjust = 0.5, barwidth = 8),
           size  = guide_legend(order = 2, title.position = "top",
                                title.hjust = 0.5, nrow = 1))
}

# Shorter display labels for long GO terms. The full term is kept in the
# exported CSV.
go_short <- c(
  "positive regulation of RNA polymerase II transcription preinitiation complex assembly" = "pos. reg. of RNAPII PIC assembly",
  "negative regulation of growth of unicellular organism as a thread of attached cells"   = "neg. regulation of pseudohyphal growth",
  "positive regulation of filamentous growth of a population of unicellular organisms"    = "pos. regulation of filamentous growth",
  "proteasome-mediated ubiquitin-dependent protein catabolic process"                     = "proteasome-mediated protein catabolism",
  "positive regulation of transcription elongation by RNA polymerase I"                   = "pos. regulation of RNAPI transcription elongation",
  "positive regulation of small molecule metabolic process"                               = "pos. regulation of small molecule metabolism",
  "positive regulation of protein-containing complex assembly"                            = "pos. reg. of protein complex assembly",
  "negative regulation of DNA-templated DNA replication"                                  = "neg. reg. of DNA replication",
  "mitochondrial respiratory chain complex IV assembly"                                   = "mito. complex IV assembly",
  "(1->3)-beta-D-glucan biosynthetic process"                                             = "β-(1,3)-glucan biosynthesis",
  "mRNA cis splicing, via spliceosome"                                                    = "mRNA cis splicing",
  "vesicle docking involved in exocytosis"                                                = "vesicle docking (exocytosis)",
  "mitochondrial cytochrome c oxidase assembly"                                           = "mito. cytochrome c oxidase assembly",
  "negative regulation of gene expression"                                                = "neg. reg. of gene expression")

GO_panel <- merged_GO_panel(CAS_dose_GO_collapsed, Complex_GO_collapsed,
                            top_n_per_col = 3,
                            keep_terms    = NULL,
                            wrap_width    = NULL,
                            term_relabel  = go_short,
                            order_by      = "pvalue",
                            term_size     = 11)


# ======================= ALLELE FREQUENCY TRAJECTORIES =======================

TRAJ_GROUPS <- c("Low CAS", "High CAS", "High CLO", "CASCLO", "Control")
TRAJ_COLORS <- TREATMENT_COLORS[TRAJ_GROUPS]
TRAJ_TEXT   <- c("black", "white", "white", "white", "white")

all_pos <- paste(snps_filtered_poly$CHROM, snps_filtered_poly$POS)

# Missense variants in one gene, taken from either analysis ("union"), or from
# just one of them.
gene_missense_pos <- function(gene, which_set = c("union", "strength", "complex")) {
  which_set <- match.arg(which_set)
  
  strength_ok <- rlang::expr(
    sig_class_CAS_Dose %in% c("Low CAS-specific", "High CAS-specific", "CAS shared") |
      Fixed_Category_CAS_Dose %in% c("Low CAS-specific", "High CAS-specific", "CAS shared"))
  
  complex_ok <- rlang::expr(
    (!is.na(sig_class_Complex) & !sig_class_Complex %in% UNINFORMATIVE) |
      (!is.na(Fixed_Category_Complex) & !Fixed_Category_Complex %in% c("None", "Other")))
  
  keep <- switch(which_set,
                 strength = strength_ok,
                 complex  = complex_ok,
                 union    = rlang::expr(!!strength_ok | !!complex_ok))
  
  ann_long_merged_filt %>%
    filter(Standard_Gene_Name == gene, Effect == "missense_variant", !!keep) %>%
    distinct(CHROM, POS)
}

# Reshape allele frequencies for the selected variants into long format,
# expanding the ancestor across every treatment so T0 appears in each facet.
gene_traj_data <- function(pos_df, gene_label) {
  
  idx <- match(paste(pos_df$CHROM, pos_df$POS), all_pos)
  idx <- idx[!is.na(idx)]
  message(gene_label, ": ", length(idx), " missense SNPs matched")
  
  long <- snps_filtered_poly[idx, ] %>%
    dplyr::select(CHROM, POS, matches("^af_EEA_"), -matches("hap")) %>%
    mutate(SNP_ID = paste0(CHROM, "_", POS)) %>%
    pivot_longer(starts_with("af_EEA_"), names_to = "sample", values_to = "af") %>%
    separate(sample,
             into = c("Prefix", "Exp", "Treatment", "Dose", "Rep", "Timepoint"),
             sep = "_", fill = "right", extra = "merge") %>%
    mutate(Timepoint = as.numeric(gsub("T", "", Timepoint)),
           Timepoint = ifelse(Treatment == "anc", 0, Timepoint),
           Rep       = ifelse(Treatment == "anc", "rep00", Rep),
           Dose      = ifelse(Treatment == "anc", NA, Dose))
  
  anc_df    <- long %>% filter(Treatment == "anc")
  nonanc_df <- long %>% filter(Treatment != "anc")
  
  # Copy the ancestor into every treatment and replicate so each facet has a T0
  anc_expanded <- anc_df %>%
    dplyr::select(CHROM, POS, SNP_ID, af) %>%
    distinct() %>%
    crossing(nonanc_df %>% distinct(Treatment, Dose, Rep)) %>%
    mutate(Timepoint = 0)
  
  bind_rows(nonanc_df, anc_expanded) %>%
    filter((Treatment == "CAS"     & Dose == "Low")  |
             (Treatment == "CAS"     & Dose == "High") |
             (Treatment == "CLO"     & Dose == "High") |
             (Treatment == "CASCLO")                   |
             (Treatment == "Control" & Dose == "ND")) %>%
    mutate(Group = case_when(
      Treatment == "CAS"     & Dose == "Low"  ~ "Low CAS",
      Treatment == "CAS"     & Dose == "High" ~ "High CAS",
      Treatment == "CLO"     & Dose == "High" ~ "High CLO",
      Treatment == "CASCLO"                   ~ "CASCLO",
      Treatment == "Control"                  ~ "Control"),
      Group = factor(Group, levels = TRAJ_GROUPS),
      Gene  = gene_label)
}

# by_rep = FALSE averages across replicates, giving one line per SNP per
# treatment. by_rep = TRUE keeps replicates separate, which is useful when the
# question concerns independent evolutionary events.
traj_row_plot <- function(traj_df, by_rep = FALSE, strip_gene = FALSE) {
  
  d <- if (by_rep) {
    traj_df %>% mutate(line_id = paste(SNP_ID, Group, Rep, sep = "_"), y = af)
  } else {
    traj_df %>%
      group_by(Gene, SNP_ID, Group, Timepoint) %>%
      summarise(y = mean(af, na.rm = TRUE), .groups = "drop") %>%
      mutate(line_id = paste(SNP_ID, Group, sep = "_"))
  }
  
  strip_spec <- strip_themed(
    background_x = elem_list_rect(fill = unname(TRAJ_COLORS)),
    text_x       = elem_list_text(color = TRAJ_TEXT))
  
  p <- ggplot(d, aes(x = Timepoint, y = y, group = line_id, color = Group)) +
    geom_line(linewidth = 0.8, alpha = if (by_rep) 0.35 else 0.6) +
    geom_point(size = if (by_rep) 0.9 else 1.5, alpha = 0.7) +
    scale_color_manual(values = TRAJ_COLORS) +
    scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1),
                       name = "Mean ALT allele frequency") +
    scale_x_continuous(breaks = c(0, 1, 7, 14), name = "Timepoint")
  
  p <- if (strip_gene) {
    p + facet_grid2(Gene ~ Group, strip = strip_spec)
  } else {
    p + facet_wrap2(~ Group, ncol = 2, strip = strip_spec)
  }
  
  p + theme_cowplot() +
    theme(legend.position = "none",
          strip.text   = element_text(size = 12, face = "bold"),
          panel.border = element_rect(color = "grey70", fill = NA),
          axis.text    = element_text(size = 9),
          axis.title   = element_text(size = 11))
}

fks1_traj <- gene_traj_data(gene_missense_pos("FKS1", "union"), "FKS1")
n_fks1    <- n_distinct(fks1_traj$SNP_ID)

fks1_row <- traj_row_plot(fks1_traj) +
  ggtitle(bquote(italic("FKS1") ~ "missense variants (n =" ~ .(n_fks1) * ")")) +
  theme(plot.title = element_text(size = 12, face = "bold"))

hym1_sur1_traj <- bind_rows(
  gene_traj_data(gene_missense_pos("HYM1", "union"), "HYM1"),
  gene_traj_data(gene_missense_pos("SUR1", "union"), "SUR1"))

hym1_sur1_traj_plot <- traj_row_plot(hym1_sur1_traj, by_rep = FALSE,
                                     strip_gene = TRUE)

# ggsave("Figures/FungAMR_and_GO_plots/fungamr_GO_fig_panelA_raw.svg",
#         fungamr_panel, width = 6, height = 8, bg = "white")
#  ggsave("Figures/FungAMR_and_GO_plots/fungamr_GO_fig_panelB_raw.svg",
#         fks1_row, width = 4, height = 6, bg = "white")
#  ggsave("Figures/FungAMR_and_GO_plots/fungamr_GO_fig_panelC_raw.svg",
#         GO_panel, width = 8, height = 7, bg = "white")
# ggsave("Figures/AF_Trajectory_Plots/hym1_sur1_traj_fig_raw.svg",
#        hym1_sur1_traj_plot, width = 8, height = 4, bg = "white", dpi = 300)


# ======================= SUMMARY TABLES ======================================

# Attach the FungAMR drug class to each gene-level summary. Only hits receive a
# category; non-hits stay NA.
annotate_summary <- function(gene_level_summary,
                             lookup  = fungamr_gene_categories,
                             relabel = NULL) {
  out <- gene_level_summary %>%
    # Drop first so re-running does not create .x/.y duplicate columns
    dplyr::select(-any_of(c("Resistance_Category", "SNP_Source"))) %>%
    left_join(lookup, by = "Standard_Gene_Name") %>%
    mutate(Resistance_Category = ifelse(FungAMR_Hit,
                                        as.character(Resistance_Category),
                                        NA_character_))
  
  # Relabel gene-level categories for export, e.g. CAS -> High CAS
  if (!is.null(relabel)) {
    out <- out %>%
      mutate(`Gene-level_Categories` = recode(`Gene-level_Categories`, !!!relabel))
  }
  
  out
}

CAS_Dose_summary_out <- annotate_summary(CAS_Dose_gene_level_summary)
Complex_summary_out  <- annotate_summary(Complex_gene_level_summary,
                                         relabel = COMPLEX_LABELS)

Complex_GO_out <- Complex_GO_collapsed %>%
  mutate(Category = recode(as.character(Category), !!!COMPLEX_LABELS),
         Category = factor(Category,
                           levels = c("High CAS", "High CLO", "CASCLO", "Shared")))

# write.xlsx(CAS_Dose_summary_out, "Output_files/CAS_Dose_gene_level_summary.xlsx", rowNames = FALSE)
# write.xlsx(Complex_summary_out,  "Output_files/Complex_gene_level_summary.xlsx",  rowNames = FALSE)
# write.xlsx(CAS_dose_GO_collapsed, "Output_files/CAS_dose_GO_enrichment_results.xlsx", row.names = FALSE)
# write.xlsx(Complex_GO_out,  "Output_files/Complex_GO_enrichment_results.xlsx",  row.names = FALSE)


# ======================= FKS1 HOTSPOT AND GENE TRACK =========================

# Known echinocandin resistance hotspots, in amino acid coordinates
# (Perlin 2007; Johnson et al. 2011)
HS1_RANGE   <- c(641, 649)
HS2_RANGE   <- c(1345, 1365)
BIN_WIDTH   <- 25     # residues per bin in the synonymous histogram
ARROW_TAPER <- 25     # residues used for the arrowhead taper
HS_FILL     <- "firebrick"
CAS_COL     <- TREATMENT_COLORS[["High CAS"]]

gff <- import("Input_files/saccharomyces_cerevisiae_R64-1-1_20110208.gff")

# FKS1 / YLR342W is a single intronless CDS on chrXII (809,997-815,627), so
# amino acid position maps linearly onto genomic coordinate.
fks1_cds <- gff[which(gff$type == "CDS" & gff$gene == "FKS1")]
stopifnot(length(fks1_cds) == 1)

FKS1_LEN <- sum(width(fks1_cds)) / 3 - 1   # 1,876 residues, excluding the stop

# FKS1 is also annotated as CND1, CWH53, GSC1, and ETG1. Confirm no variants
# were filed under an alias and missed by the Standard_Gene_Name filter.
fks1_aliases <- unique(c("FKS1", unlist(gff[which(gff$gene == "FKS1")]$Alias)))
CAS_variants %>% filter(Standard_Gene_Name %in% fks1_aliases) %>% count(Standard_Gene_Name)

FKS1_variants <- CAS_variants %>%
  filter(Standard_Gene_Name == "FKS1") %>%
  mutate(
    # Protein_pos is stored as "645/1876"; take the residue number
    aa_pos   = as.numeric(gsub("/.*", "", Protein_pos)),
    prot_len = as.numeric(gsub(".*/", "", Protein_pos)),
    FKS1_hotspot = case_when(
      aa_pos >= HS1_RANGE[1] & aa_pos <= HS1_RANGE[2] ~ "HS1",
      aa_pos >= HS2_RANGE[1] & aa_pos <= HS2_RANGE[2] ~ "HS2",
      TRUE ~ "Outside"),
    Effect_grp = case_when(
      Effect == "synonymous_variant" ~ "Synonymous",
      Effect == "missense_variant"   ~ "Missense",
      TRUE                           ~ "Other")) %>%
  filter(!is.na(aa_pos))

# Agreement between the GFF and SnpEff confirms both used the same reference,
# so the HS1/HS2 residue numbering is directly comparable.
stopifnot(FKS1_LEN == max(FKS1_variants$prot_len, na.rm = TRUE))

# Counts reported in the text: 187 synonymous, 15 missense, 2 synonymous in
# HS1 and 1 in HS2
FKS1_variants %>% count(FKS1_hotspot, Effect_grp)

hs_df       <- tibble(hotspot = c("HS1", "HS2"),
                      xmin = c(HS1_RANGE[1], HS2_RANGE[1]),
                      xmax = c(HS1_RANGE[2], HS2_RANGE[2]))
missense_df <- FKS1_variants %>% filter(Effect_grp == "Missense")
synon_df    <- FKS1_variants %>% filter(Effect_grp == "Synonymous")

# One shared x-axis so both panels line up
x_scale <- scale_x_continuous(limits = c(1, FKS1_LEN),
                              breaks = c(1, 500, 1000, 1500, FKS1_LEN),
                              expand = expansion(mult = 0.01))

# Top panel: protein track with missense variants marked above it
p_track <- ggplot() +
  # Faint full-height bands tie the two panels together visually
  geom_rect(data = hs_df, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = HS_FILL, alpha = 0.15) +
  annotate("polygon",
           x = c(1, FKS1_LEN - ARROW_TAPER, FKS1_LEN, FKS1_LEN - ARROW_TAPER, 1),
           y = c(-0.15, -0.15, 0, 0.15, 0.15),
           fill = "grey88", colour = "grey40", linewidth = 0.3) +
  geom_rect(data = hs_df, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -0.15, ymax = 0.15),
            fill = HS_FILL, colour = NA) +
  annotate("text", x = FKS1_LEN / 2, y = 0, label = "FKS1",
           fontface = "bold.italic", colour = "grey20", size = 4.5) +
  geom_segment(data = missense_df, inherit.aes = FALSE,
               aes(x = aa_pos, xend = aa_pos, y = 0.15, yend = 0.50),
               colour = CAS_COL, linewidth = 0.5) +
  geom_point(data = missense_df, inherit.aes = FALSE,
             aes(x = aa_pos, y = 0.50), colour = CAS_COL, size = 2) +
  geom_segment(data = hs_df, inherit.aes = FALSE,
               aes(x = (xmin + xmax) / 2, xend = (xmin + xmax) / 2,
                   y = -0.15, yend = -0.38), colour = "grey40", linewidth = 0.3) +
  geom_text(data = hs_df, inherit.aes = FALSE,
            aes(x = (xmin + xmax) / 2, y = -0.46, label = hotspot),
            fontface = "bold", size = 3.8) +
  annotate("text", x = 1, y = 0.62,
           label = paste0("Missense (n = ", nrow(missense_df), ")"),
           hjust = 0, size = 3.5, colour = "grey25") +
  x_scale +
  scale_y_continuous(name = NULL, limits = c(-0.55, 0.68), breaks = NULL) +
  theme_cowplot() +
  theme(legend.position = "none",
        axis.line    = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin  = margin(5, 5, 0, 5))

# Bottom panel: synonymous variant density along the protein
p_hist <- ggplot(synon_df, aes(x = aa_pos)) +
  geom_rect(data = hs_df, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = HS_FILL, alpha = 0.18) +
  geom_histogram(binwidth = BIN_WIDTH, boundary = 0,
                 fill = "grey40", colour = "white", linewidth = 0.2) +
  annotate("text", x = FKS1_LEN * 0.02, y = Inf,
           label = paste0("Synonymous (n = ", nrow(synon_df), ")"),
           hjust = 0, vjust = 1.6, size = 3.5, colour = "grey30") +
  x_scale +
  scale_y_continuous(name = paste0("Variants per ", BIN_WIDTH, " aa"),
                     expand = expansion(mult = c(0, 0.15))) +
  labs(x = "FKS1 amino acid position") +
  theme_cowplot() +
  theme(plot.margin = margin(0, 5, 5, 5))

fks1_track_fig <- p_track / p_hist + plot_layout(heights = c(1, 1.1))

# ggsave("Figures/FKS1_Gene_Track_Plot/fks1_variant_track_fig_raw.svg",
#        fks1_track_fig, width = 10, height = 6, bg = "white", dpi = 300)


# ======================= SESSION INFO ========================================

sessionInfo()

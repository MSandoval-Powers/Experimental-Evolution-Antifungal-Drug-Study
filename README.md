# Experimental-Evolution-Antifungal-Drug-Study

Analysis code and data for **Sandoval-Powers M, Crestani G, Dugo H, Shearman E, Burke MK.** *Experimental evolution reveals distinct genomic trajectories under different intensities and combinations of antifungal selection in* Saccharomyces cerevisiae.

[![DOI](https://img.shields.io/badge/data-Figshare-blue)](https://doi.org/PLACEHOLDER)
[![Data](https://img.shields.io/badge/reads-NCBI%20SRA-green)](https://www.ncbi.nlm.nih.gov/bioproject/PLACEHOLDER)

---

## Overview

Replicate populations founded from a genetically diverse, outcrossed **4S** *Saccharomyces cerevisiae* population were experimentally evolved under six selective regimes that varied in both the **intensity** and the **complexity** of antifungal selection. Pooled whole-genome sequencing across three timepoints was used to track genome-wide allele frequency dynamics, and evolved populations were phenotyped for growth performance across drug environments.

This repository contains the scripts and small input/output files needed to reproduce the analyses and figures. Large input and output files related to the Genomics Analysis are archived on Figshare (see [Data availability](#data-availability)). 

Figures are included alongside the scripts that generate them. Where a figure was assembled or annotated for publication, both versions are provided: the .svg written directly by the script, and the final .png after editing in Inkscape. Scripts reproduce the .svg, not the edited version.

### Experimental design

| Treatment | Drug condition | Replicates | Transfers | Est. generations |
|---|---|---|---|---|
| Low CAS | 0.035 µM caspofungin | 12 | 13 | 56 |
| High CAS | 0.07 µM caspofungin | 12 | 12 | 54 |
| Low CLO | 0.35 µM clotrimazole | 12 | 13 | 64 |
| High CLO | 3.00 µM clotrimazole | 12 | 13 | 59 |
| CASCLO | 0.035 µM caspofungin + 0.35 µM clotrimazole | 12 | 12 | 61 |
| Control | 0.05% DMSO | 12 | 13 | 69 |

72 evolved populations, sampled at T1, T7, and T14 → **216 pooled sequencing samples**. After filtering: **64,764 biallelic SNPs** at a median genome-wide coverage of 86× (range 27×–290×).

---

## Repository structure

```
Experimental-Evolution-Antifungal-Drugs/
├── Drug_Dose_Determination_Analysis/      # Dose–response assays used to determine drug concentrations for experimental evolution
│   ├── Input_files/
│   ├── Scripts/
│   └── Figures/                           # Fig. 
├── Genomics_Analysis/
│   ├── SNP_QC_Workflow/                   # Coverage QC, filtering, PCA, SFS, identifying fixed SNPs
│   │   ├── Input_files/                   # See Figshare
│   │   ├── Output_files/                  # See Figshare, except Coverage_Summary_Filtered.csv
│   │   ├── Scripts/
│   │   └── Figures/
│   │       ├── Coverage/                  # Figs. 
│   │       ├── PCA/                       # Fig. 
│   │       └── SFS/                       # Fig. 
│   ├── GLM_AF_Change_Workflow/            # Quasibinomial GLMs, candidate SNP identification
│   │   ├── Input_files/                   # See Figshare
│   │   ├── Output_files/                  # See Figshare
│   │   └── Scripts/
│   └── SNP_Analysis_and_Annotation_Workflow/   # Determining significant SNPs, SnpEff annotation, FungAMR, GO enrichment
│       ├── Input_files/                   # partly Figshare
│       ├── Output_files/
│       ├── Scripts/
│       └── Figures/
│           ├── AF_Trajectory_Plots/       # Fig. 5B
│           ├── Beta_Coefficient_Plot/     # Fig. 2
│           ├── FKS1_Gene_Track_Plot/      # Fig. 
│           ├── FungAMR_and_GO_Plots/      # Fig. 5A, 5C
│           └── Manhattan/                 # Figs. 3, 4
├── Growth_Phenotyping_Analysis/           # Analysis of growth performance data of evolved populations and ancestor
│   ├── Input_files/
│   ├── Output_files/
│   ├── Scripts/
│   └── Figures/                           # Fig. 
└── README.md
```


### `Drug_Dose_Determination_Analysis/`

Dose–response assays of the ancestral 4S population, used to select the low and high concentrations of each drug and the combined-drug pairing. 

| Script | Purpose |
|---|---|
| `Drug_dose_assays_script.R` | Summarizes OD600 across drug doses and timepoints; generates Fig. S1 |


### `Genomics_Analysis/`

The core genomic workflow, run in three stages. Each stage consumes the previous stage's output. Analysis begins from the merged SNP table; read alignment and variant
calling were performed separately and are not included here.

**1 · `SNP_QC_Workflow/`**

Starts from the merged SNP table and filters to biallelic sites polymorphic in the ancestor (0.02 < AF < 0.98) with read depth ≥ 10 in every population, excluding low-confidence sites identified from the founder genotypes. Also covers coverage QC, exploratory PCA and site frequency spectra, identification of fixed sites, and
chromosome VIII coverage.

Input and output files for this stage are hosted on Figshare, except `Coverage_Summary_Filtered.csv`, which is in the repository.

| Script | Purpose |
|---|---|
| `SNP_QC_Workflow_Script.R` | Coverage QC, SNP filtering (126,215 → 64,764 SNPs), PCA, site frequency spectra, fixed-SNP identification, chromosome VIII coverage 


**2 · `GLM_AF_Change_Workflow/`**

Per-SNP quasibinomial GLMs weighted by read depth, fitted in two forms: 
single-treatment models (`AF ~ time`) testing for change within an environment, and cross-treatment models (`AF ~ time × treatment`) comparing rates of change across environments. Three cross-treatment contrasts: CAS selection strength, CLO selection strength, and single- vs. combined-drug. 

| Script | Purpose |
|---|---|
| `GLM_Script.R` | Fits per-SNP GLMs across all six treatments and three cross-treatment contrasts; writes two master results tables of β coefficients and FDR-adjusted p-values |

Produces model output only; figures from these results are generated in stage 3. Input and output files for this stage are hosted on Figshare.


**3 · `SNP_Analysis_and_Annotation_Workflow/`**

Applies the top 5% significance threshold to the GLM output, classifies candidate SNPs as treatment-specific or shared using the pairwise slope contrasts, and merges them with the fixed-SNP candidates from stage 1. Candidates are annotated with SnpEff, summarized at the gene level, cross-referenced against the *S. cerevisiae* subset of the FungAMR database (Bédard et al. 2025), and tested for GO biological process enrichment with `topGO` and redundancy reduction in `rrvgo`.

| Script | Purpose |
|---|---|
| `SNP_analysis_and_annotation_script.R` | Candidate classification, SnpEff annotation, FungAMR enrichment, GO enrichment, and all main-text figures (Figs. 2–5) |

Input files are hosted on Figshare; output files are in the repository.

### `Growth_Phenotyping_Analysis/`

48-hour growth assays of T14 evolved populations and the ancestor across control, single-drug, combined-drug, and amphotericin B media. 

| Script | Purpose |
|---|---|
| `Growth_phenotyping_script.R` | Fits logistic growth curves to OD600 data for six assays, models doubling time and carrying capacity, generates Figs. 6 and S6 |

The raw plate-reader exports are converted to formatted tables by a step near the top of the script that is commented out by default; the formatted tables it produces are included in `Output_files/`. Uncomment only if regenerating them from raw.

---

## Data availability

| Resource | Location |
|---|---|
| Raw sequencing reads, evolved populations | NCBI SRA, BioProject `[accession]` |
| Raw sequencing reads, 4S ancestor | NCBI SRA, [PRJNA678990](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA678990) (Burke et al. 2020) |
| Large input and output files | Figshare, [`[DOI]`](https://doi.org/PLACEHOLDER) |
| Analysis code (this repository, archived) | Zenodo, [`[DOI]`](https://doi.org/PLACEHOLDER) |
| Manuscript | `[DOI]` |

Files too large for version control are archived on Figshare with **folder names and relative paths matching this repository**. Download the deposit and place the folders into your clone; scripts will resolve their paths without editing.

---

## Requirements

**R** ≥ 4.5.1

```r
# Required — CRAN
install.packages(c(
  "tidyverse",      # ggplot2, dplyr, tidyr, stringr, readr, purrr, tibble
  "data.table", "reshape2", "glue", "scales",
  "cowplot", "patchwork", "ggh4x", "ggridges", "ghibli", "ggpattern",
  "ComplexUpset", "gggenes",
  "factoextra", "svglite", "ggrastr",
  "lme4", "lmerTest", "emmeans", "multcompView", "growthcurver",
  "pbapply", "openxlsx"
))

# Required — Bioconductor
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("topGO", "org.Sc.sgd.db", "rrvgo", "Rgraphviz", "rtracklayer"))

```


---

## Reproducing the analysis

1. Clone this repository.
2. Download the Figshare deposit and place its folders into the corresponding directories here.
3. Run `Drug_Dose_Determination_Analysis/` (independent of the others).
4. Run `Genomics_Analysis/` in order: `SNP_QC_Workflow/` → `GLM_AF_Change_Workflow/` → `SNP_Analysis_and_Annotation_Workflow/`. 
5. Run `Growth_Phenotyping_Analysis/` scripts (independent of the genomic workflow).

Write calls (`ggsave`, `write_tsv`, `write.xlsx`) are commented out by default — uncomment to regenerate outputs. The GLM stage takes a while; its results tables
are on Figshare if you want to skip refitting.

Sequence data processing from raw reads requires the SRA download and is not reproducible from this repository alone; start from the filtered SNP table on Figshare to reproduce downstream results.

---

## Citation

```
[Full citation once published]
```


---

## Contact

**Megan Sandoval-Powers** — sandovme@oregonstate.edu
**Molly K. Burke** — molly.burke@oregonstate.edu

Department of Integrative Biology, Oregon State University, Corvallis, OR 97331, USA

---


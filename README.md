# Experimental-Evolution-Antifungal-Drug-Study

Analysis code and data for **Sandoval-Powers M, Crestani G, Dugo H, Shearman E, Burke MK.** *Experimental evolution reveals distinct genomic trajectories under different intensities and combinations of antifungal selection in* Saccharomyces cerevisiae.

[![DOI](https://img.shields.io/badge/data-Figshare-blue)](https://doi.org/PLACEHOLDER)
[![Data](https://img.shields.io/badge/reads-NCBI%20SRA-green)](https://www.ncbi.nlm.nih.gov/bioproject/PLACEHOLDER)
[![License](https://img.shields.io/badge/license-MIT-lightgrey)](LICENSE)

---

## Overview

Replicate populations founded from the outcrossed **4S** *Saccharomyces cerevisiae* population were experimentally evolved under six selective regimes that varied in both the **intensity** and the **complexity** of antifungal selection. Pooled whole-genome sequencing across three timepoints was used to track genome-wide allele frequency dynamics, and evolved populations were phenotyped for growth performance across drug environments.

This repository contains the scripts and small input files needed to reproduce the analyses and figures. Large input and output files related to Genomics Analysis are archived on Figshare (see [Data availability](#data-availability)).

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
├── Drug_Dose_Determination_Analysis/   # Dose–response assays used to determine drug concentrations for experimental evolution
├── Genomics_Analysis/                  # SNP QC, allele frequency modeling, analysis and annotation
├── Growth_Phenotyping_Analysis/        # Doubling time and carrying capacity of evolved populations
└── README.md
```

### `Drug_Dose_Determination_Analysis/`

Dose–response assays of the ancestral 4S population, used to select the low and high concentrations of each drug and to identify a combined-drug pairing that permitted growth.

| Script | Purpose |
|---|---|
| `[script]` | `[processes raw plate reader output]` |
| `[script]` | `[fits dose–response curves; generates Fig. S1]` |

### `Genomics_Analysis/`

The core genomic workflow, run in three stages. Each stage consumes the previous stage's output.

**1 · SNP QC**

Reads aligned to the *S. cerevisiae* S288C reference (R64-1-1) with BWA; variants called with GATK v4.5.0.0 following best practices, with base quality scores recalibrated against known SNPs from 1,011 natural isolates (Peter et al. 2018). Allele frequencies computed per site as `AD/DP`. Sites were then filtered to retain biallelic variants polymorphic in the ancestor (0.02 < AF < 0.98) with read depth ≥ 10 in every population, excluding low-confidence sites identified from founder genotypes.

| Script | Purpose |
|---|---|
| `[script]` | `[extracts AD/DP from merged VCF]` |
| `[script]` | `[applies depth, ancestral polymorphism, and founder-genotype filters]` |

**2 · Allele frequency modeling**

Per-SNP quasibinomial GLMs weighted by read depth, fitted in two forms: single-treatment models (`AF ~ time`) testing for change within an environment, and cross-treatment models (`AF ~ time × treatment`) comparing rates of change across environments. Slopes extracted with `emtrends`, contrasts with `pairs`, p-values FDR-corrected. Candidate loci are the top 5% of SNPs by adjusted p-value within each treatment.

Three cross-treatment contrasts: CAS selection strength, CLO selection strength, and single- vs. combined-drug.

| Script | Purpose |
|---|---|
| `[script]` | `[fits single-treatment GLMs]` |
| `[script]` | `[fits cross-treatment GLMs; extracts slope contrasts]` |
| `[script]` | `[classifies candidate SNPs as treatment-specific or shared]` |
| `[script]` | `[Manhattan and UpSet plots; Figs. 3–4]` |
| `[script]` | `[β coefficient distributions; Fig. 2]` |
| `[script]` | `[chromosome VIII coverage analysis; Fig. S5]` |

**3 · Annotation and enrichment**

Candidate SNPs annotated with SnpEff and summarized at the gene level, cross-referenced against the *S. cerevisiae* subset of the FungAMR database (Bédard et al. 2025), and tested for GO biological process enrichment using `topGO` (weight01 algorithm, Fisher's exact test, minimum node size 10) with redundancy reduction in `rrvgo` (Wang similarity, threshold 0.9).

| Script | Purpose |
|---|---|
| `[script]` | `[SnpEff annotation and gene-level summary]` |
| `[script]` | `[FungAMR overlap and enrichment testing]` |
| `[script]` | `[GO enrichment and semantic reduction; Fig. 5]` |

### `Growth_Phenotyping_Analysis/`

48-hour growth assays of T14 evolved populations across control, single-drug, combined-drug, and amphotericin B media. Doubling time and carrying capacity modeled with linear mixed-effects models (`lme4`) with replicate as a random effect; pairwise contrasts from estimated marginal means (`emmeans`), Tukey-corrected.

| Script | Purpose |
|---|---|
| `[script]` | `[extracts growth curves; computes DT and K]` |
| `[script]` | `[fits mixed models; pairwise contrasts]` |
| `[script]` | `[generates Figs. 6 and S6]` |

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
install.packages(c(
  "tidyverse",      # ggplot2, dplyr, tidyr, stringr, readr, purrr, tibble
  "data.table", "reshape2", "glue", "scales",
  "cowplot", "patchwork", "ggpubr", "ggh4x", "ggridges", "ghibli",
  "ggvenn", "ggVennDiagram", "ComplexUpset", "gggenes",
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
3. Run `Drug_Dose_Determination_Analysis/` scripts (independent of the others).
4. Run `Genomics_Analysis/` scripts in order: SNP QC → allele frequency modeling → annotation.
5. Run `Growth_Phenotyping_Analysis/` scripts (independent of the genomic workflow).

Sequence data processing from raw reads requires the SRA download and is not reproducible from this repository alone; start from the filtered SNP table on Figshare to reproduce downstream results.

---

## Citation

```
[Full citation once published]
```

Please also cite the Figshare deposit if you reuse the data directly.

---

## Contact

**Megan Sandoval-Powers** — sandovme@oregonstate.edu
**Molly K. Burke** — molly.burke@oregonstate.edu

Department of Integrative Biology, Oregon State University, Corvallis, OR 97331, USA

---

## License

Code in this repository is released under the MIT License. Data archived on Figshare is released under CC-BY 4.0.

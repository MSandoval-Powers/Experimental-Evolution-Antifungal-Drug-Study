# Genomics_Analysis

Genome-wide analysis of allele frequency change in experimentally evolved
*Saccharomyces cerevisiae* populations. See the [top-level README](../README.md)
for study design, requirements, and data availability.

Three stages, run in order. Each consumes the previous stage's output. Analysis
begins from the merged SNP table; read alignment and variant calling were
performed separately and are not included here.

```
SNP_QC_Workflow  →  GLM_AF_Change_Workflow  →  SNP_Analysis_and_Annotation_Workflow
```

---

## 1 · SNP_QC_Workflow/

Filters the merged SNP table and characterizes genome-wide patterns before
modeling. SNPs are retained if polymorphic in the ancestor (0.02 < AF < 0.98),
covered at read depth ≥ 10 in every population, and not flagged as
low-confidence by the founder genotypes. Mitochondrial SNPs are excluded,
leaving 64,764 biallelic SNPs at a median coverage of 86×.

Also covers coverage QC, PCA, site frequency spectra, identification of
early-fixed sites, and chromosome VIII coverage.

**Script:** `Scripts/SNP_QC_Workflow_Script.R`

**Inputs**

| File | Location |
|---|---|
| `filtered_snps.txt` (merged SNP table) | Figshare |

**Outputs**

| File | Location | Manuscript |
|---|---|---|
| `snps_filtered_poly.txt` (filtered SNP table) | Figshare | — |
| `*_only_filt.txt` (per-treatment allele frequency tables) | Figshare | — |
| `fixed_snps_filt.txt` (early-fixed sites) | Figshare | — |
| `Coverage_Summary_Filtered.csv` | `Output_files/` | Table S1 |
| PCA, SFS, chromosome VIII figures | `Figures/` | Figs. S2, S3, S6 |

Early-fixed sites reached a frequency of 0 or 1 by T1 and stayed fixed through
T14 in every replicate. They have no temporal variance and cannot be modeled in
stage 2, so they are carried forward to stage 3 as a separate candidate set.

---

## 2 · GLM_AF_Change_Workflow/

Per-SNP quasibinomial GLMs of allele frequency over time, weighted by read
depth, fitted in two forms:

- **Single-treatment** (`AF ~ time`) — allele frequency change within one environment
- **Cross-treatment** (`AF ~ time × treatment`) — rates of change compared among environments

Three cross-treatment contrasts: CAS intensity (Low CAS, High CAS, Control),
CLO intensity (Low CLO, High CLO, Control), and single versus combined drug
(High CAS, High CLO, CASCLO, Control).

**Script:** `Scripts/GLM_Script.R`

**Inputs**

| File | Location |
|---|---|
| `*_only_filt.txt` (per-treatment allele frequency tables, stage 1) | Figshare |

**Outputs**

| File | Location |
|---|---|
| `model.results.Dose.txt` (drug intensity contrasts) | Figshare |
| `model.results.Complex.txt` (single vs combined drug) | Figshare |

Figures from these results are generated in stage 3.

This stage takes ~30 minutes. SNPs are processed in chunks to limit memory
use. The results tables are on Figshare if you want to skip refitting.

---

## 3 · SNP_Analysis_and_Annotation_Workflow/

Identifies candidate SNPs from the GLM output, then annotates and characterizes
them:

1. Apply the top 5% significance threshold within each treatment
2. Classify candidates as treatment-specific or shared using the pairwise slope contrasts
3. Merge with the early-fixed candidates from stage 1
4. Annotate with SnpEff and summarize at the gene level
5. Cross-reference against the *S. cerevisiae* subset of FungAMR (Bédard et al. 2025)
6. Test GO biological process enrichment with topGO, reducing redundancy with rrvgo

**Script:** `Scripts/SNP_analysis_and_annotation_script.R`

**Inputs**

| File | Location |
|---|---|
| `model.results.Dose.txt`, `model.results.Complex.txt` (stage 2) | Figshare |
| `*_only_filt.txt`, `snps_filtered_poly.txt`, `fixed_snps_filt.txt` (stage 1) | Figshare |
| `filtered_snps.txt` (merged SNP table) | Figshare |
| `snps_annotated.txt` (SnpEff output) | Figshare |
| `saccharomyces_cerevisiae_R64-1-1_20110208.gff` | Figshare |
| `FungAMR_db.csv` | `Input_files/` |

**Outputs**

| File | Location | Manuscript |
|---|---|---|
| Gene-level summaries | `Output_files/` | Tables S2–S5 |
| GO enrichment results | `Output_files/` | Tables S6, S7 |
| `snps.vcf` (written for SnpEff) | `Output_files/` | — |
| β coefficient plots | `Figures/Beta_Coefficient_Plot/` | Fig. 2 |
| Manhattan and UpSet plots | `Figures/Manhattan/` | Figs. 3, 4 |
| FungAMR and GO panels | `Figures/FungAMR_and_GO_Plots/` | Fig. 5 |
| *FKS1* gene track | `Figures/FKS1_Gene_Track_Plot/` | Fig. S4 |
| AF trajectories | `Figures/AF_Trajectory_Plots/` | Fig. S5 |

SnpEff runs outside this repository. The script writes `snps.vcf`, SnpEff is run
on it separately, and `snps_annotated.txt` is read back in.

---

## Naming

Sequencing files and NCBI records use the original scheme; the manuscript uses
readable treatment labels.

| Original | Manuscript |
|---|---|
| `EEA_CAP_S01_*` | Low CAS |
| `EEA_CAP_S03_*` | High CAS |
| `EEA_CLO_S01_*` | Low CLO |
| `EEA_CLO_S03_*` | High CLO |
| `EEA_CAC_S01_*` | CASCLO |
| `EEA_CTR_S01_*` | Control |

Timepoints `D01`, `D07`, `D14` correspond to T1, T7, T14. The ancestor appears
as `EEA_anc_4SH_12` and the four founders as `EEA_hap_*`.


In the single versus combined drug analysis, column names and classification
labels are shortened to `CAS` and `CLO`. These refer to the **high**
concentrations of each drug. Figures and exported tables use the full labels.

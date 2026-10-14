# Growth_Phenotyping_Analysis

48-hour growth assays of evolved populations from the final timepoint (T14) and
the 4S ancestor, across control, single-drug, combined-drug, and amphotericin B
media. See the [top-level README](../README.md) for study design, requirements,
and data availability.

Doubling time (DT) and carrying capacity (*K*) were estimated by fitting
logistic growth curves with growthcurver and used as proxies for population
fitness. Six assays were run: five compare an evolved population against the
ancestor (Low CAS, High CAS, Low CLO, High CLO, CASCLO), and the
cross-tolerance assay compares High CAS against CASCLO populations, including in
amphotericin B.

**Script:** `Scripts/Growth_phenotyping_script.R`

**Inputs**

| File | Location |
|---|---|
| Raw plate-reader exports, one per assay | `Input_files/` |
| Sample-name CSVs, one per assay | `Input_files/` |

**Outputs**

| File | Location | Manuscript |
|---|---|---|
| Formatted OD600 tables, one per assay | `Output_files/` | — |
| DT and *K* figures | `Figures/` | Figs. 6, S6 |

---

## Notes

Assay configuration lives in the `ASSAYS` list near the top of the script: file
paths, plate layout, treatment labels, and the outlier filters applied before
each model. To add or modify an assay, edit that list rather than the analysis
code.

DT and *K* are modeled separately and use different outlier filters. Rows are
dropped where the ancestor either failed to grow in drug media or never reached
stationary phase, since logistic growth parameters cannot be estimated
meaningfully in those wells. The CAS filters are deliberately asymmetric; see
the comments in `ASSAYS`.

`process_OD600_data()` converts raw plate-reader exports into the formatted
tables in `Output_files/`. It is commented out by default, since the formatted
tables are what the rest of the script reads. Uncomment only to regenerate them
from raw.

Models are fitted with lme4 and lmerTest, with replicate as a random effect and
pairwise contrasts from estimated marginal means. One exception: for CLO
carrying capacity, assay and population are perfectly confounded, so that model
uses a nested random effect (`Assay/Replicate`) and its three *p*-values come
from refitting with each assay environment as the reference level rather than
from emmeans.

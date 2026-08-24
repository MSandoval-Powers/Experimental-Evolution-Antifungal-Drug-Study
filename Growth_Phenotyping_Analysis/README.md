# Growth Phenotyping Analysis

48-hour growth assays of evolved populations from the final timepoint (T14) and the 4S
ancestor, across control, single-drug, combined-drug, and amphotericin B media. Doubling
time (DT) and carrying capacity (*K*) were estimated by fitting logistic growth curves with the `growthcurver` package
and used as proxies for population fitness.

Six assays were run. Five compare an evolved population against the ancestor (Low CAS,
High CAS, Low CLO, High CLO, CASCLO); the cross-tolerance assay compares High CAS against
CASCLO populations, including in amphotericin B, to test for broader resistance to a
mechanistically distinct antifungal class.

## Contents

| Folder | Contents |
|---|---|
| `Scripts/` | `Growth_phenotyping_script.R` — fits growth curves, models DT and *K*, generates Figs. 6 and S6 |
| `Input_files/` | Raw plate-reader exports and sample-name CSVs, one pair per assay |
| `Output_files/` | Formatted OD600 tables, one per assay |
| `Figures/` | Main and supplementary DT/*K* figures |

## Notes

Assay configuration lives in the `ASSAYS` list near the top of the script — file paths,
plate layout, treatment labels, and the outlier filters applied before each model. To add
or modify an assay, edit that list rather than the analysis code.

DT and *K* are modeled separately and use different outlier filters. Rows are dropped
where the ancestor either failed to grow in drug media or never reached stationary phase,
since logistic growth parameters can't be estimated meaningfully in those wells. The CAS
filters are deliberately asymmetric — see the comments in `ASSAYS`.

`process_OD600_data()` converts raw plate-reader exports into the formatted tables in
`Output_files/`. It is commented out by default, since the formatted tables are what the
rest of the script reads; uncomment only to regenerate them from raw.

Models are fitted with `lme4` and `lmerTest` (replicate as a random effect), with pairwise
contrasts from estimated marginal means. One exception: for CLO carrying capacity, assay
and population are perfectly confounded, so that model uses a nested random effect
(`Assay/Replicate`) and its three p-values come from refitting with each assay environment
as the reference level rather than from `emmeans`.

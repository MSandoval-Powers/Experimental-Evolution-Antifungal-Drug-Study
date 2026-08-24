# Drug Dose Determination Analysis

Dose–response assays of the ancestral 4S *Saccharomyces cerevisiae* population, used to determine the drug concentrations for the experimental evolution study.

Two concentrations per drug were chosen to impose distinct levels of selective pressure — Low CAS (0.035 µM), High CAS (0.07 µM), Low CLO (0.35 µM), High CLO (3.00 µM) — with low concentrations producing modest inhibition and high concentrations producing substantial but sublethal inhibition. 
For the combined treatment, the low-concentration pairing (0.035 µM CAS + 0.35 µM CLO) was selected because higher combined concentrations produced complete growth arrest.

## Contents

| Folder | Contents |
|---|---|
| `Scripts/` | `Drug_dose_assays_script.R` — summarizes OD600 across drug doses and timepoints; generates Fig. S1 |
| `Input_files/` | `OD600_sumdata_24well_drugdose.csv` — replicate OD600 values at 24 h and 48 h, relative to control |
| `Figures/` | `Drug_dose_response_assays_sumfig.png` - contains figure generated from script |

## Notes

The input CSV contains OD600 values already normalized to the no-drug control. Assays were run in 24-well plates under conditions mirroring the evolution experiment.

The `ggsave` call at the end of the script is commented out by default — uncomment to regenerate the figure.

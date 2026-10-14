# Drug_Dose_Determination_Analysis

Dose–response assays of the ancestral 4S *Saccharomyces cerevisiae* population,
used to determine the drug concentrations for the experimental evolution study.
See the [top-level README](../README.md) for study design, requirements, and
data availability.

Two concentrations per drug were chosen to impose distinct levels of selective
pressure, with low concentrations producing modest inhibition and high
concentrations producing substantial but sublethal inhibition: Low CAS
(0.035 µM), High CAS (0.07 µM), Low CLO (0.35 µM), and High CLO (3.00 µM). For
the combined treatment, the low-concentration pairing (0.035 µM CAS + 0.35 µM
CLO) was selected because higher combined concentrations produced complete
growth arrest.

**Script:** `Scripts/Drug_dose_assays_script.R`

**Inputs**

| File | Location |
|---|---|
| `OD600_sumdata_24well_drugdose.csv` (replicate OD600 at 24 and 48 hr, relative to control) | `Input_files/` |

**Outputs**

| File | Location | Manuscript |
|---|---|---|
| `Drug_dose_response_assays_sumfig.png` | `Figures/` | Fig. S1 |

---

## Notes

The input CSV contains OD600 values already normalized to the no-drug control.
Assays were run in 24-well plates under conditions mirroring the evolution
experiment.

The `ggsave` call at the end of the script is commented out by default.
Uncomment to regenerate the figure.

# Sample Design, Labeling & Preprocessing

## Overview

These scripts build the factorial experimental design for the pristine
(single-polymer) and colored (mixed-polymer) microplastic-contaminated soil
samples, label the raw Vis-NIR-SWIR reflectance spectra collected with an ASD
FieldSpec4 spectroradiometer against that design, and preprocess the
resulting spectral library (splice correction + Savitzky-Golay denoising) for
downstream analysis.

---

## Data availability

"/raw_spectra/FIELDSPEC4_data", "/raw_spectra/FIELDSPEC4_data2", and auxiliary data files are available in [INSERT DOI ONCE THE IEEE Data Descriptions  ARTICLE IS PUBLISHED]

<!> ATTENTION for discrepancies within the directories, between the published dataset and the paths used here <!>

---

## Sample Design

`pristine_sample_design.R` and `colored_sample_design.R` each build a full
factorial design and assign `SAMPLE_ID`/`SAMPLE_CODE` values.

| Design    | Script                     | Factors                                              | N   | `SAMPLE_ID` range |
|---|---|---|---|---|
| Pristine  | `pristine_sample_design.R` | 4 polymers (PP, PVC, PET, PE) × 3 sizes × 6 masses × 4 replicates | 288 | 1–288   |
| Colored   | `colored_sample_design.R`  | 1 mix ("ALL") × 3 sizes × 5 masses × 4 replicates    | 60  | 289–348 |

In the pristine design, combinations already prepared in an earlier round
(PP and PVC, the two coarsest size intervals, all masses) are flagged
`PREPARED = TRUE` and assigned `SAMPLE_ID`s first, ahead of the remaining
combinations — so ID order reflects preparation order, not just the
factorial order.

`SAMPLE_CODE` encodes polymer, size class, mass, and replicate, e.g.
`PVC_A1.2_1` (pristine) or `MIX_A3750_1` (colored), where the size letter
comes from:

```r
size_code_map = c("4.80–1.00" = "A", "1.00–0.60" = "B", "0.60–0.053" = "C")
```

Each script writes its design to a CSV under `raw_spectra/`.

---

## Data Labeling

`label_pristine_data.R` and `label_colored_data.R` merge the raw `.asd`
spectra with the corresponding sample design.

| Object     | Script                  | Contents                                             |
|---|---|---|
| `pristine` | `label_pristine_data.R` | Pristine design + raw reflectance spectra (`spc`)    |
| `colored`  | `label_colored_data.R`  | Colored design + raw reflectance spectra (`spc`)     |

### How it works

For each set, the script:

1. Reads every `.asd` file in the matching raw spectra folder
   (`FIELDSPEC4_data` for pristine, `FIELDSPEC4_data2` for colored).
2. Drops the first reading (spectralon blank reference).
3. Isolates the contaminated-soil readings — the first 288 (pristine) or 60
   (colored) rows, per `raw_spectra/notes [31.05.2025].txt` and
   `raw_spectra/notes [04.08.2025].txt` respectively.
4. Renames spectra `spec001`–`spec288` (pristine) or `spec289`–`spec348`
   (colored), matching the `SAMPLE_ID` ranges above.
5. Merges the spectra into the design table's `spc` column and saves the
   result.

`label_pristine_data.R` also defines `plot_spectra()`, a helper for
spot-checking spectra by `SAMPLE_CODE` prefix (e.g. `plot_spectra(pristine,
"PVC_A1.2")`).

---

## Spectral Preprocessing

`spectra_processing.R` loads both raw, labeled datasets and applies two
correction steps in sequence, with before/after `matplot()` checks at each
stage:

| Step               | Function                        | Key parameters                                         |
|---|---|---|
| Splice correction   | `prospectr::spliceCorrection()` | Splice wavelengths pulled from an ASD file's metadata (fields 30–31) |
| Denoising           | `prospectr::savitzkyGolay()`    | `m =` differentiation order,<br> `p =`  polynomial order (`poly_sg`),<br> `w =` window size (`swindow_sg`)   |

A signal-to-noise ratio check (mean / SD over the 1000–1100 nm region) is
also run before and after denoising to confirm the filter improves spectral
quality.

---

## Output Files

| File                                                    | Contents |
|---|---|
| `raw_spectra/pristine_sample_design.csv`                | Pristine factorial sample design |
| `raw_spectra/colored_sample_design.csv`                 | Colored factorial sample design |
| `raw_spectra/raw_pristine.rds`                          | Pristine design + raw reflectance spectra |
| `raw_spectra/raw_colored.rds`                           | Colored design + raw reflectance spectra |
| `preprocessed_data/pristine_splice_corrected.rds`       | Pristine spectra after splice correction |
| `preprocessed_data/colored_splice_corrected.rds`        | Colored spectra after splice correction |
| `preprocessed_data/pristine_denoised.rds`               | Pristine spectra after Savitzky-Golay denoising |
| `preprocessed_data/colored_denoised.rds`                | Colored spectra after Savitzky-Golay denoising |

---

## Dependencies

```r
library(dplyr)
library(tidyr)
library(ggplot2)
library(asdreader)
library(prospectr)
```

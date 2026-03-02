# Spectral Dissimilarity Analysis — Microplastic contaminated soil

## Overview

The `exploratory_analysis.R` script performs spectral dissimilarity analysis on Vis-NIR-SWIR reflectance spectra of pristine and colored microplastic samples. It computes pairwise
dissimilarities between mean group spectra using five complementary metrics,
organizes results by experimental context, and ranks polymer comparisons from
most to least dissimilar.

---

## Data

| Object     | Source file              | Contents                                          |
|---|---|---|
| `pristine`  | `datsoil.rds`            | Raw spectra — pristine polymers (used for PCA)    |
| `colored`   | `datsoil2.rds`           | Raw spectra — colored polymers (used for PCA)     |
| `pristine2` | `deNoised_data.rds`      | Denoised spectra — pristine (used for dissimilarity) |
| `colored2`  | `deNoised_newData.rds`   | Denoised spectra — colored (used for dissimilarity)  |

Spectra are converted to absorbance via `log(1/R)` before any analysis.

---

## Dissimilarity Metrics

Five metrics are computed on mean group spectra, grouped by `POLYMER × SIZE_CODE × MASS_mg`.
They are applied in priority order during ranking:

| Priority | Variable | Full name | Role |
|---|---|---|---|
| 1° | `mahD` | Mahalanobis Distance | Accounts for covariance structure — best for detecting outlier spectra |
| 2° | `EucD` | Euclidean Distance | Overall magnitude and shape difference between spectra |
| 3° | `samD` | Spectral Angle Mapper | Pure shape difference, invariant to intensity scaling |
| 4° | `cd1` | Correlation Dissimilarity | Global shape difference based on linear correlation |
| 5° | `mwcd` | Moving Window Correlation Dissimilarity | Localized shape differences across spectral segments |

---

## Contextual Analysis (Pristine vs. Pristine)

### CONFIG block — the only section you need to change

```r
CONTROL_VARS = c("MASS_mg", "SIZE_CODE")
VARY_LABEL   = function(row) paste0(row$POLYMER, "_", row$SIZE_CODE, row$MASS_mg)
```

**`CONTROL_VARS`** defines which two variables are *held fixed* to form a subset.
Within each subset, all pairwise comparisons are made across the remaining (varying) variable.

| Context       | `CONTROL_VARS`                  | Fixed           | Varying |
|---|---|---|---|
| `CONTEXT_i`   | `c("MASS_mg", "SIZE_CODE")`     | mass + size     | polymer |
| `CONTEXT_ii`  | `c("POLYMER", "SIZE_CODE")`     | polymer + size  | mass    |
| `CONTEXT_iii` | `c("MASS_mg", "POLYMER")`       | mass + polymer  | size    |

**`VARY_LABEL`** controls how comparison labels are printed and stored
(e.g. `PE_A1.2 vs PET_A1.2`). It should always reference the *varying* variable
so labels remain informative within a given context.

### How it works

`control_combos` enumerates every unique combination of the fixed variables.
For each combination, `compare_subset()` finds matching rows in `groups`,
generates all pairwise index combinations via `combn()`, and looks up the
corresponding value in each metric matrix. Results are printed to console and
stored in a nested list keyed as `CONTEXT_n_VAR1_val1_VAR2_val2`.

### Output

Results are tidied into `tidy_results` — a long-format tibble with columns
`Context`, `Metric`, `Comparison`, `Value` — then summarised and ranked.
To save, uncomment and rename the last line to match the active context:

```r
# |> write.csv(file = "CONTEXT_n.csv", row.names = FALSE)
```

---

## Cross Analysis (Pristine vs. Colored)

Compares each pristine polymer (`PE`, `PET`, `PP`, `PVC`) against its
corresponding colored mixture (`MIX_`) at matched size and mass. This is a
non-square cross-dissimilarity: each `MIX_` column is matched against its
four candidate polymer rows.

`_cross` metric variants are used throughout. For Mahalanobis distance
specifically, the colored samples are projected into the PCA space of the
pristine reference set before distances are computed.

Save via:

```r
# |> write.csv(file = "CONTEXT_iv.csv", row.names = FALSE)
```

---

## Output Files

| File             | Contents |
|---|---|
| `CONTEXT_i.csv`   | Pairwise polymer dissimilarities — mass + size fixed |
| `CONTEXT_ii.csv`  | Pairwise mass dissimilarities — polymer + size fixed |
| `CONTEXT_iii.csv` | Pairwise size dissimilarities — polymer + mass fixed |
| `CONTEXT_iv.csv`  | Cross dissimilarities — pristine polymers vs. colored mixtures |

---

## Dependencies

```r
library(purrr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(factoextra)
library(effectsize)
library(resemble)
library(tripack)
library(splancs)
library(gplots)
library(RColorBrewer)
```
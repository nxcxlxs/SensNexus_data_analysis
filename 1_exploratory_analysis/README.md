# Spectral Exploratory Analysis

## Overview

The `exploratory_analysis.R` script performs a two-part exploratory analysis
on Vis-NIR-SWIR reflectance spectra of pristine and colored microplastic
samples:

1. **Principal Component Analysis (PCA)** on raw pristine spectra, used to
   visualize sample structure, project colored samples into the pristine
   PCA space, and check whether colored samples fall inside the pristine
   predictive domain (convex hull).
2. **Spectral Dissimilarity Analysis** on denoised spectra, which computes
   pairwise dissimilarities between mean group spectra using five
   complementary metrics, organizes results by experimental context, and
   ranks polymer comparisons from most to least dissimilar.

---

## Data

| Object      | Source file                                    | Contents                                          |
|---|---|---|
| `pristine`  | `../raw_spectra/raw_pristine.rds`              | Raw spectra — pristine polymers (used for PCA)    |
| `colored`   | `../raw_spectra/raw_colored.rds`               | Raw spectra — colored polymers (used for PCA)     |
| `pristine2` | `../preprocessed_data/pristine_denoised.rds`   | Denoised spectra — pristine (used for dissimilarity) |
| `colored2`  | `../preprocessed_data/colored_denoised.rds`    | Denoised spectra — colored (used for dissimilarity)  |

Spectra are converted to absorbance via `log(1/R)` before any analysis.

---

## Part 1 — Principal Component Analysis

### PCA on raw pristine spectra

`prcomp()` is run on the absorbance-converted pristine spectra (`pristine$spcA`).
The script reports the variance explained by the first 10 components and
plots:

- **Scree plot** (`fviz_eig`) — eigenvalues vs. number of dimensions.
- **Scores plot** — PC1 vs. PC2 for all pristine samples.
- **Loadings plots** — first and second PC loadings vs. wavelength (nm),
  plus the last PC's loading (shown separately as it is mostly noise).

### Predictive domain check (pristine vs. colored)

This section checks whether the colored samples fall within the spectral
space defined by the pristine reference set:

1. A Delaunay triangulation (`tri.mesh`) is built from the pristine PC1/PC2
   scores, and its **convex hull** is extracted (`convex.hull`).
2. Colored sample spectra are centered and projected into the pristine PCA
   space using the pristine model's center and rotation matrix.
3. Point-in-polygon testing (`inpip`) determines which projected colored
   points fall inside vs. outside the pristine convex hull.
4. The script reports the percentage of colored samples that fall **outside**
   the predictive domain, and lists their `SAMPLE_ID`, `SAMPLE_CODE`,
   `SIZE_INTERVALS_mm`, `MASS_mg`, and `REPLICATE` for follow-up.
5. Both a base-R plot and a `ggplot2` version are produced, showing pristine
   samples, colored samples, points outside the domain (marked with an "x"),
   and the convex hull boundary.

### Structure / biplot check

- Scores and loadings are extracted into data frames.
- A biplot overlays sample scores (PC1 vs. PC2) with loading vectors for
  representative wavelengths across the Vis-NIR-SWIR range (350, 700, 1100,
  2500 nm), to relate the PCA structure back to specific regions of the
  spectrum.

### Grouped score plots

Three `fviz_pca_ind` panels color the pristine PCA scores by:

- **Polymer** (`p1`)
- **Size interval** (`p2`)
- **Mass (mg)** (`p3`)

plus a fourth panel (`p4`) repeating the mass grouping with 95% confidence
ellipses. These are composed into a single tagged panel figure with
`patchwork`.

### Beyond visual inspection — statistical checks

- **Correlation tests** (`cor.test`) between `MASS_mg` and PC1 / PC2.
- **ANOVA** (`aov`) of PC1 and PC2 by `POLYMER`, to test whether polymer
  identity explains score variation.
- **Effect sizes** (`eta_squared`) for `MASS_mg` and `POLYMER` on both PC1
  and PC2, to quantify how much variance each factor explains.

---

## Part 2 — Spectral Dissimilarity Analysis

### Dissimilarity Metrics

Five metrics are computed on mean group spectra, grouped by `POLYMER × SIZE_CODE × MASS_mg`.
They are applied in priority order during ranking:

| Priority | Variable | Full name | Role |
|---|---|---|---|
| 1° | `mahD` | Mahalanobis Distance | Accounts for covariance structure — best for detecting outlier spectra |
| 2° | `EucD` | Euclidean Distance | Overall magnitude and shape difference between spectra |
| 3° | `samD` | Spectral Angle Mapper | Pure shape difference, invariant to intensity scaling |
| 4° | `cd1` | Correlation Dissimilarity | Global shape difference based on linear correlation |
| 5° | `mwcd` | Moving Window Correlation Dissimilarity | Localized shape differences across spectral segments |

### Contextual Analysis (Pristine vs. Pristine)

#### CONFIG block — the only section you need to change

```r
CONTROL_VARS = c("MASS_mg", "SIZE_CODE")
VARY_LABEL = function(row) paste0(row$POLYMER, "_", row$SIZE_CODE, row$MASS_mg)
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

#### How it works

`control_combos` enumerates every unique combination of the fixed variables.
For each combination, `compare_subset()` finds matching rows in `groups`,
generates all pairwise index combinations via `combn()`, and looks up the
corresponding value in each metric matrix. Results are printed to console and
stored in a nested list keyed as `CONTEXT_n_VAR1_val1_VAR2_val2`.

#### Output

Results are tidied into `tidy_results` — a long-format tibble with columns
`Context`, `Metric`, `Comparison`, `Value` — then summarised and ranked.
To save, uncomment and rename the last line to match the active context:

```r
# |> write.csv(file = "CONTEXT_n.csv", row.names = FALSE)
```

#### Heatmaps

For each of the five metrics (`EucD`, `mahD`, `cd1`, `mwcd`, `samD`), a
symmetric `heatmap.2` matrix (rows/columns labeled `POLYMER_SIZE_MASS`) is
plotted with a fixed color key, giving a visual complement to the ranked
`tidy_results` table.

### Cross Analysis (Pristine vs. Colored)

Compares each pristine polymer (`PE`, `PET`, `PP`, `PVC`) against its
corresponding colored mixture (`MIX_`) at matched size and mass. This is a
non-square cross-dissimilarity: each `MIX_` column is matched against its
four candidate polymer rows.

`_cross` metric variants are used throughout. For Mahalanobis distance
specifically, the colored samples are projected into the PCA space of the
pristine reference set before distances are computed.

Results are tidied and ranked the same way as the contextual analysis. Save via:

```r
# |> write.csv(file = "CONTEXT_iv.csv", row.names = FALSE)
```

#### Heatmaps

The same five metrics are also plotted as non-square `heatmap.2` matrices
(pristine polymer groups as rows, colored mixture groups as columns) —
one heatmap per metric.

---

## Output Files

| File             | Contents |
|---|---|
| `CONTEXT_i.csv`   | Pairwise polymer dissimilarities — mass + size fixed |
| `CONTEXT_ii.csv`  | Pairwise mass dissimilarities — polymer + size fixed |
| `CONTEXT_iii.csv` | Pairwise size dissimilarities — polymer + mass fixed |
| `CONTEXT_iv.csv`  | Cross dissimilarities — pristine polymers vs. colored mixtures |

PCA outputs (scree plot, scores/loadings plots, biplot, grouped panels,
predictive-domain plots) and dissimilarity heatmaps are generated as R plots
within the script and are not written to disk unless explicitly saved.

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

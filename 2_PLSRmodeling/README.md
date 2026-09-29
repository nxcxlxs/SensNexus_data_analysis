# PLSR Mass Prediction from Vis-NIR-SWIR Spectra

## Overview

The `PLSRmodeling.R` script builds and evaluates Partial Least Squares
Regression (PLSR) models that predict microplastic **mass (`MASS_mg`)** from
Vis-NIR-SWIR reflectance spectra. Four spectral preprocessing treatments
(M1 to M4) are compared through four validation schemes of increasing
difficulty:

1. **75/25 hold-out** on pristine samples, with the optimal number of
   components chosen by cross-validation.
2. **Stratified k-fold cross-validation** on the full pristine dataset.
3. **External validation** — models trained on all pristine samples are used
   to predict colored samples.
4. **Tuned external validation** — the number of components is re-optimized on
   one half of the colored samples and evaluated on the other half.

---

## Data

| Object         | Source file                                  | Contents                                   |
|---|---|---|
| `pristine`     | `../preprocessed_data/pristine_denoised.rds` | Denoised spectra — pristine (M2, M3, M4)   |
| `pristine_raw` | `../raw_spectra/raw_pristine.rds`            | Raw spectra — pristine (M1)                |
| `colored`      | `../preprocessed_data/colored_denoised.rds`  | Denoised spectra — colored (M2, M3, M4)    |
| `colored_raw`  | `../raw_spectra/raw_colored.rds`             | Raw spectra — colored (M1)                 |

Each object must contain a reflectance matrix `spc` and the columns
`SAMPLE_ID`, `POLYMER` (pristine only), `MASS_mg` and `SIZE_CODE`.

Spectra are converted to absorbance with `log(1/R)` (natural logarithm) before
any modelling, and stored as `spcA`.

---

## Preprocessing Treatments

| Model | Description                                   | Spectral matrix used | Pipeline |
|---|---|---|---|
| M1    | Raw data                                      | `spcA` (raw objects)  | Absorbance only |
| M2    | Minimal preprocessing                         | `spcA` (denoised)     | Denoised spectra from `spectra_processing.R` |
| M3    | Intermediate preprocessing                    | `spcAmovav`           | M2 + SNV + moving average (`w = 11`) |
| M4    | Full preprocessing                            | `spcARmovav`          | M3 + 5 nm resampling |

The resampling wavelength axis (`newWavs`) is built from the pristine
wavelengths and reused for the colored samples, which therefore share the
same spectral range.

---

## Part 1 — 75/25 Hold-out

### Split

Samples are split with `sample_frac(0.75)` within strata defined by
`POLYMER × MASS_mg × SIZE_CODE` (`set.seed(1)`), giving a calibration set
(`datC`, `rawC`) and a validation set (`datV`, `rawV`). The distributions of
`POLYMER`, `MASS_mg` and `SIZE_CODE` in both sets are printed for checking.

### Model fitting and component selection

For each treatment, a `plsr()` model is fitted with `method = "oscorespls"`,
`ncomp = 50` and internal cross-validation (`validation = "CV"`). The optimal
number of components is read from the cross-validated RMSEP curve
(`RMSEP()`, `which.min()`).

Plots produced:

- **RMSEP validation plots** — one panel per model, with the optimal `ncomp`
  and minimum RMSEP marked.
- **Regression coefficients** — coefficients vs. wavelength (nm) at each
  model's optimal `ncomp`.

### Evaluation

Predictions are made on the calibration and validation sets at the optimal
`ncomp`, and ME, RMSE and R² are printed for each model (see
[Validation statistics](#validation-statistics)).

---

## Part 2 — K-fold Cross-Validation

`cv_plsr_model(data, spc_matrix, ncomp, nfolds = 4, seed = 999)` runs a
stratified k-fold CV on the full pristine dataset:

1. Fold labels are assigned within each `POLYMER` × `MASS_mg` × `SIZE_CODE`
   stratum, so every fold contains all strata.
2. For each fold, a PLSR model is trained on the remaining folds and used to
   predict the held-out fold.
3. Predictions from all folds are pooled and matched back through `SAMPLE_ID`.

The function returns `ME`, `RMSE`, `R2` and a `table` of observed and predicted
values. It is run once per treatment using the `ncomp` values from Part 1.

---

## Part 3 — External Validation (Pristine → Colored)

Models are refitted on the **full** pristine dataset and used, with the optimal
`ncomp` from Part 1, to predict all colored samples.
ME, RMSE and R² are reported for each treatment. No colored sample is used in
training or in the choice of `ncomp`.

---

## Part 4 — Tuned External Validation

The colored samples are split 50/50 within `MASS_mg × SIZE_CODE` strata
(`set.seed(42)`) into:

- **`colored1` / `colored_raw1`** — used to tune the number of components.
- **`colored2` / `colored_raw2`** — used for the final evaluation.

### Hyperparameter search

`evaluate_ncomp()` predicts `colored1` with every `ncomp` from 1 to 30 for
each pristine-trained model and returns ME, RMSE and R². The results are
stacked in `results`, and `best_nc` keeps the row with the lowest RMSE per
model.

### Final evaluation

Tuned models (`PLSR_mod_tuned` to `PLSR_mod_tuned4`) are predicted on
`colored2` / `colored_raw2` and ME, RMSE and R² are printed for each treatment.

### Residual plots

A two-panel `patchwork` figure summarizes the results by observed mass, with
one facet per model:

- **(a)** Absolute residuals (`Predicted − Observed`).
- **(b)** Relative residuals on the log scale (`log(Predicted) − log(Observed)`).

---

## Validation statistics

Custom functions are used throughout, all ignoring `NA` values:

| Statistic | Definition |
|---|---|
| `ME`   | Mean error, `mean(pred − obs)`. Positive values indicate over-prediction. |
| `RMSE` | Root mean squared error, in mg. |
| `R2`   | `1 − SSE / SST`, computed against the mean of the observed values. It can be negative when a model predicts worse than the mean. |

---

## Output

This script does not write files to disk. Results are printed to the console and
plotted in R:

| Object / output        | Contents |
|---|---|
| `calib_stats`, `valid_stats` | ME, RMSE and R² for the hold-out split |
| `cv_plsr_results`      | ME, RMSE, R² and prediction tables from k-fold CV |
| `external_stats`       | ME, RMSE and R² for the naive external validation |
| `results`, `best_nc`   | `ncomp` grid search results and the selected optimum per model |
| `res`, `res2`          | Absolute and relative residual boxplots |

---

## Dependencies

```r
library(prospectr)
library(dplyr)
library(pls)
library(ggplot2)
library(patchwork)
```

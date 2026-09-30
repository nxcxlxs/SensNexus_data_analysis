# Random Forest Mass Prediction from Vis-NIR-SWIR Spectra

## Overview

The `RFmodeling.R` script builds and evaluates Random Forest (RF) regression
models that predict microplastic **mass (`MASS_mg`)** from Vis-NIR-SWIR
reflectance spectra. It mirrors the PLSR workflow (`PLSRmodeling.R`): four
spectral preprocessing treatments (M1 to M4) are compared through four
validation schemes of increasing difficulty:

1. **75/25 hold-out** on pristine samples, with variable importance inspection.
2. **Stratified k-fold cross-validation** on the full pristine dataset.
3. **External validation** — models trained on all pristine samples are used
   to predict colored samples.
4. **Tuned external validation** — `ntree` and `mtry` are re-optimized on one
   half of the colored samples and evaluated on the other half.

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

| Model | Description                | Spectral matrix used | Pipeline |
|---|---|---|---|
| M1    | Raw data                   | `spcA` (raw objects) | Absorbance only |
| M2    | Minimal preprocessing      | `spcA` (denoised)    | Denoised spectra from `spectra_processing.R` |
| M3    | Intermediate preprocessing | `spcAmovav`          | M2 + SNV + moving average (`w = 11`) |
| M4    | Full preprocessing         | `spcARmovav`         | M3 + 5 nm resampling |

The resampling wavelength axis (`newWavs`) is built from the pristine
wavelengths and reused for the colored samples, which therefore share the
same spectral range.

---

## Random Forest Configuration

All models use `randomForest::randomForest()` with the formula
`plastMass ~ .`, where the predictors are the spectral bands (columns named
`spec.<wavelength>`).

| Parameter      | Default (Parts 1–3) | Part 4 (tuned)                          |
|---|---|---|
| `ntree`        | 150                 | 100, 150, 200, 250, 500 (grid)          |
| `mtry`         | 10                  | 5 to 100, 13 values (grid)              |

---

## Part 1 — 75/25 Hold-out

### Split

Samples are split with `sample_frac(0.75)` within strata defined by
`POLYMER × MASS_mg × SIZE_CODE` (`set.seed(1)`), giving a calibration set
(`datC`, `rawC`) and a validation set (`datV`, `rawV`). The distributions of
`POLYMER`, `MASS_mg` and `SIZE_CODE` in both sets are printed for checking.

### Model fitting

One RF model (`RF_mod_mass` to `RF_mod_mass4`) is fitted per treatment with the
default configuration.

Plots and outputs produced:

- **Variable importance plots** — `varImpPlot()` for each model (`%IncMSE` and
  `IncNodePurity`).
- **Top-30 bands** — the 30 most important wavelengths of Model 1, ranked by
  each importance metric.

### Evaluation

Predictions are made on the calibration and validation sets, and ME, RMSE and
R² are printed for each model (see [Validation statistics](#validation-statistics)).

---

## Part 2 — K-fold Cross-Validation

`cv_random_forest(data, spc_matrix, ntree = 150, mtry = 10, nfolds = 4, seed = 1)`
runs a stratified k-fold CV on the full pristine dataset:

1. Fold labels are assigned within each `POLYMER` × `MASS_mg` × `SIZE_CODE`
   stratum, so every fold contains all strata.
2. For each fold, an RF model is trained on the remaining folds and used to
   predict the held-out fold.
3. Predictions from all folds are pooled and matched back through `SAMPLE_ID`.

The function returns `ntree`, `mtry`, `ME`, `RMSE`, `R2` and a `table` of
observed and predicted values. It is run once per treatment with the default
`ntree` and `mtry`.

---

## Part 3 — External Validation (Pristine → Colored)

Models are refitted on the **full** pristine dataset (`rawFULL`, `dataFULL2`,
`dataFULL3`, `dataFULL4`; `set.seed(2)`) and used to predict all colored
samples. ME, RMSE and R² are reported for each treatment. No colored sample is
used in training or in the choice of hyperparameters.

---

## Part 4 — Tuned External Validation

The colored samples are split 50/50 within `MASS_mg × SIZE_CODE` strata
(`set.seed(23)`) into:

- **`colored1` / `colored_raw1`** — used to tune `ntree` and `mtry`.
- **`colored2` / `colored_raw2`** — used for the final evaluation.

### Hyperparameter search

For each treatment, a pristine-trained RF is fitted for every combination in
`param_grid_ext` (5 `ntree` × 13 `mtry` = 65 combinations, 260 fits in total)
and used to predict `colored1`. ME, RMSE and R² are stacked in `ext_results`.
The best combination per model is inspected with two criteria:

The best combination per model is the one with the **lowest RMSE** and is stored
in `best_param` (one row per model, ordered M1 to M4). The tuned models read
`ntree` and `mtry` directly from `best_param`.

The search loop takes time to run; the elapsed time is printed at the end.

> **Note on ties:** `best_param` is built with `filter(RMSE == min(RMSE))`,
> which keeps *every* combination that shares the minimum RMSE within a model.
> An exact tie is unlikely with continuous RMSE values, but it is possible
> (for example, when `mtry` values above the number of predictors are capped by
> `randomForest` and yield equivalent models). If a tie occurs, `best_param`
> will contain more than four rows, and because the tuned models read
> `best_param$ntree[i]` and `best_param$mtry[i]` by **position** (1 to 4), the
> parameters of the following models would be misaligned, silently. Check that
> `nrow(best_param) == 4` (or `table(best_param$Model)`) before fitting the
> tuned models, and break ties explicitly if needed (e.g.
> `slice_min(RMSE, n = 1, with_ties = FALSE)`).

### Final evaluation

Tuned models (`RF_mod_tuned` to `RF_mod_tuned4`) are predicted on
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

| Object / output              | Contents |
|---|---|
| `calib_stats`, `valid_stats` | ME, RMSE and R² for the hold-out split |
| `impPlot` to `impPlot4`      | Variable importance matrices of the hold-out models |
| `cv_results`                 | ME, RMSE, R² and prediction tables from k-fold CV |
| `external_stats`             | ME, RMSE and R² for the naive external validation |
| `ext_results`                | `ntree`/`mtry` grid search results for every model |
| `best_param`                 | Selected `ntree`/`mtry` (lowest RMSE) per model |
| `tuned_stats`                | ME, RMSE and R² of the tuned models on `colored2` |
| `res`, `res2`                | Absolute and relative residual boxplots |

---

## Dependencies

```r
library(prospectr)
library(dplyr)
library(randomForest)
library(ggplot2)
library(patchwork)
```

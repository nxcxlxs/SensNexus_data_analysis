# Cubist Mass Prediction from Vis-NIR-SWIR Spectra

## Overview

The `CUBISTmodeling.R` script builds and evaluates Cubist regression models
that predict microplastic **mass (`MASS_mg`)** from Vis-NIR-SWIR reflectance
spectra. It mirrors the PLSR and RF workflows (`PLSRmodeling.R`,
`RFmodeling.R`): four spectral preprocessing treatments (M1 to M4) are compared
through four validation schemes of increasing difficulty:

1. **75/25 hold-out** on pristine samples, with a model dimensionality check
   (number of wavelengths used by each model).
2. **Stratified k-fold cross-validation** on the full pristine dataset.
3. **External validation** — models trained on all pristine samples are used
   to predict colored samples.
4. **Tuned external validation** — `committees` and `neighbors` are
   re-optimized on one half of the colored samples and evaluated on the other
   half.

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

## Cubist Configuration

All models use `Cubist::cubist()` with the spectral matrix as predictors (`x`)
and `MASS_mg` as the response (`y`).

| Parameter    | Default (Parts 1–3) | Part 4 (tuned)                    |
|---|---|---|
| `committees` | 1                   | 1, 10, 25, 50, 75, 100 (grid)     |
| `neighbors`  | 0                   | 0, 1, 2, 3, 5, 7, 9 (grid)        |

`neighbors` is **only used at prediction time** in Cubist, not in training.
Passing it to `cubist()` (as done in the CV function and in the full-dataset
fits) has no effect, so Parts 2 and 3 are effectively run with `neighbors = 0`.
It is only truly applied in Part 4, where it is passed to `predict()`.

---

## Part 1 — 75/25 Hold-out

### Split

Samples are split with `sample_frac(0.75)` within strata defined by
`POLYMER × MASS_mg × SIZE_CODE` (`set.seed(1)`), giving a calibration set
(`datC`, `rawC`) and a validation set (`datV`, `rawV`). The distributions of
`POLYMER`, `MASS_mg` and `SIZE_CODE` in both sets are printed for checking.

### Model fitting

One Cubist model (`CUB_mod_mass` to `CUB_mod_mass4`) is fitted per treatment
with the default configuration.

As a dimensionality check, `sum(varImp(model) > 0)` reports how many
wavelengths each model actually uses.

### Evaluation

Predictions are made on the calibration and validation sets, and ME, RMSE and
R² are printed for each model (see [Validation statistics](#validation-statistics)).

---

## Part 2 — K-fold Cross-Validation

`cv_cubist_model(data, spc_matrix, committees = 1, neighbors = 0, nfolds = 4, seed = 999)`
runs a stratified k-fold CV on the full pristine dataset:

1. The minimum stratum size is checked, and `nfolds` is reduced to that size if
   it is smaller than the requested number of folds (a note is printed).
2. Fold labels are assigned within each `POLYMER` × `MASS_mg` × `SIZE_CODE`
   stratum, so every fold contains all strata.
3. For each fold, a Cubist model is trained on the remaining folds and used to
   predict the held-out fold. Spectra are converted to data frames with
   generic column names (`V1`, `V2`, ...) so training and validation columns
   match.
4. Predictions from all folds are pooled and matched back through `SAMPLE_ID`.

The function returns `committees`, `neighbors`, `nfolds_used`,
`min_stratum_size`, `ME`, `RMSE`, `R2` and a `table` of observed and predicted
values. It is run once per treatment with the default `committees` and
`neighbors`.

---

## Part 3 — External Validation (Pristine → Colored)

Models are refitted on the **full** pristine dataset (`CUB_mod_FULL` to
`CUB_mod_FULL4`) and used to predict all colored samples. ME, RMSE and R² are
reported for each treatment. No colored sample is used in training or in the
choice of hyperparameters.

---

## Part 4 — Tuned External Validation

The colored samples are split 50/50 within `MASS_mg × SIZE_CODE` strata
(`set.seed(22)`) into:

- **`colored1` / `colored_raw1`** — used to tune `committees` and `neighbors`.
- **`colored2` / `colored_raw2`** — used for the final evaluation.

### Hyperparameter search

For each treatment, a pristine-trained Cubist model is fitted **once per
`committees` value** (6 values, 24 fits in total). Each fitted model is then
used to predict `colored1` with every `neighbors` value (7 values), giving
6 × 7 = 42 combinations per model, 168 in total. ME, RMSE and R² are stacked in
`ext_results`.

The best combination per model is the one with the **lowest RMSE** and is stored
in `best_param` (one row per model, ordered M1 to M4). The tuned models read
`committees` and `neighbors` directly from `best_param`, so the tuned parameters
are not hard-coded.

The search loop takes time to run; the elapsed time is printed at the end.

> **Note on ties:** `best_param` is built with `filter(RMSE == min(RMSE))`,
> which keeps *every* combination that shares the minimum RMSE within a model.
> An exact tie is unlikely with continuous RMSE values, but it is possible here
> (for example, when different `neighbors` values yield identical predictions
> for a given model). If a tie occurs, `best_param` will contain more than four
> rows, and because the tuned models read `best_param$committees[i]` and
> `best_param$neighbors[i]` by **position** (1 to 4), the parameters of the
> following models would be misaligned, silently. Check that
> `nrow(best_param) == 4` (or `table(best_param$Model)`) before fitting the
> tuned models, and break ties explicitly if needed (e.g.
> `slice_min(RMSE, n = 1, with_ties = FALSE)`).

### Final evaluation

Tuned models (`CUB_ext_tuned` to `CUB_ext_tuned4`) are refitted on the full
pristine dataset with the selected `committees`, predicted on
`colored2` / `colored_raw2` with the selected `neighbors`, and ME, RMSE and R²
are printed for each treatment.

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

This script does not write files to disk. Results are printed to the console:

| Object / output              | Contents |
|---|---|
| `calib_stats`, `valid_stats` | ME, RMSE and R² for the hold-out split |
| `cv_cubist_results`          | ME, RMSE, R² and prediction tables from k-fold CV |
| `full_stats`                 | ME, RMSE and R² for the naive external validation |
| `ext_results`                | `committees`/`neighbors` grid search results for every model |
| `best_param`                 | Selected `committees`/`neighbors` (lowest RMSE) per model |
| `CUBpred_ext` to `CUBpred_ext4` | Tuned predictions on `colored2` / `colored_raw2` |

---

## Residual analysis

Unlike the PLSR and RF scripts, this script does not produce residual plots.
Since Cubist failed in our main experiment, the residual analysis was omitted.

---

## Dependencies

```r
library(prospectr)
library(dplyr)
library(Cubist)
```

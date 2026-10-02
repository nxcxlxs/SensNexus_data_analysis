<p align="center"><img src="https://github.com/user-attachments/assets/4912b73a-ec1a-4cc5-9268-3ca425f6ff0d" width="400"></p>

# Microplastic Mass Prediction in Soil from Vis-NIR-SWIR Spectra

## Overview

This repository contains the R code used to predict microplastic **mass
(`MASS_mg`)** in contaminated soil from Vis-NIR-SWIR reflectance spectra
(ASD/Malvern Panalytical FieldSpec4). Two sample sets are used:

| Set      | Description                          | Design                                                  | n   |
|---|---|---|---|
| Pristine | Single-polymer samples (PP, PVC, PET, PE) | 4 polymers × 3 sizes × 6 masses × 4 replicates     | 288 |
| Colored  | Mixed-polymer samples ("ALL")        | 1 mix × 3 sizes × 5 masses × 4 replicates               | 60  |

Models are trained on pristine samples and tested, in increasing order of
difficulty, on held-out pristine samples, through cross-validation, and on the
colored samples (external validation).

---

## Repository structure

Each folder contains its own README with the details of the scripts.

| Folder                    | Purpose |
|---|---|
| `0_data_preparation/`     | Sample design, labeling of raw `.asd` spectra, splice correction and Savitzky-Golay denoising |
| `1_exploratory_analysis/` | PCA, predictive domain check (pristine vs. colored) and spectral dissimilarity analysis |
| `2_PLSRmodeling/`         | Partial Least Squares Regression models |
| `3_RFmodeling/`           | Random Forest models |
| `4_CUBISTmodeling/`       | Cubist models |

Scripts are meant to be run in folder order, since each step uses the outputs
of the previous ones. Modeling scripts read their inputs from `../raw_spectra/`
and `../preprocessed_data/`, so these two data folders are expected at the same
level as the numbered folders.

---

## Data availability and associated publications

| Resource                                           | Reference |
|---|---|
| Dataset archive (IEEE DataPort)  | https://dx.doi.org/10.21227/vh42-0d98 |
| Data descriptor article (IEEE Data Descriptions)   | [DOI WILL BE INCLUDED AS SOON THE IEEE Data Descriptions ARTICLE IS PUBLISHED] |
| Main article (modeling and analysis)               | [DOI WILL BE INCLUDED AS SOON THE MAIN PAPER IS PUBLISHED] |

The dataset archive is associated with the data descriptor article. The code in
this repository belongs to the main article.

<!> ATTENTION for discrepancies within the directories, between the published dataset and the paths used here <!>

### How to cite

If you use this code or data, please cite:

- **Repository:** Ortiz, N., Silva, M., Andrade, C., & ten Caten, A. (2026). *SensNexus_data_analysis* [Computer software]. GitHub. https://github.com/nxcxlxs/SensNexus_data_analysis
- **Dataset:** Ortiz, N., Silva, M., Andrade, C., & ten Caten, A. (2026). *SensNexus Spectral Dataset (SensNexusDat)* [Data set]. IEEE Dataport. https://doi.org/10.21227/vh42-0d98
- **Data descriptor:** [DATA DESCRIPTOR CITATION]
- **Main article:** [MAIN PAPER CITATION]

---

## Shared modeling framework

The three modeling folders (PLSR, RF, Cubist) follow the same design so that
their results are comparable.

### Preprocessing treatments

Spectra are converted to absorbance with `log(1/R)` (natural logarithm), and
four treatments (M1 to M4) are compared:

| Model | Description                | Pipeline |
|---|---|---|
| M1    | Raw data                   | Absorbance only |
| M2    | Minimal preprocessing      | Splice correction + Savitzky-Golay denoising |
| M3    | Intermediate preprocessing | M2 + SNV + moving average (`w = 11`) |
| M4    | Full preprocessing         | M2 + 5 nm resampling + SNV + moving average (`w = 11`) |

### Validation schemes

1. **75/25 hold-out** on pristine samples (stratified by polymer × mass × size).
2. **Stratified k-fold cross-validation** on the full pristine dataset.
3. **External validation** — models trained on all pristine samples predict the
   colored samples.
4. **Tuned external validation** — hyperparameters are re-optimized on one half
   of the colored samples and evaluated on the other half.

Performance is reported with ME, RMSE (mg) and R². The hyperparameter tuned in
each algorithm is the number of components (PLSR), `ntree` and `mtry` (RF), and
`committees` and `neighbors` (Cubist), always selected by the lowest RMSE.

### Notes

- In the RF and Cubist scripts, the best parameters are selected with
  `filter(RMSE == min(RMSE))`, which can keep tied combinations. See the
  note on ties in the corresponding READMEs.
- Residual plots are produced for PLSR and RF. They were omitted for Cubist
  because it failed in our main experiment.
- <span style="color: darkgreen;">**100% artisanal, organic code🌱🌿♻️**</span><br>
  Written by hand, by someone who is learning on the go,
  so it is definitely *not* the most elegant or efficient.<br> *Please be patient.*

---

## Dependencies

R packages used across the project:

```r
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(patchwork)
library(asdreader)
library(prospectr)
library(pls)
library(randomForest)
library(Cubist)
library(factoextra)
library(effectsize)
library(resemble)
library(tripack)
library(splancs)
library(gplots)
library(RColorBrewer)
```

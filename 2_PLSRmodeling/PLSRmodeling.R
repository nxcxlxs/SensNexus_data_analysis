require(prospectr)
require(dplyr)
require(pls)
require(ggplot2)
require(patchwork)


# load processed data
pristine = readRDS("../preprocessed_data/pristine_denoised.rds")
pristine_raw = readRDS("../raw_spectra/raw_pristine.rds")
colored = readRDS("../preprocessed_data/colored_denoised.rds") 
colored_raw = readRDS("../raw_spectra/raw_colored.rds")

# convert spectra to absorbance
pristine$spcA = log(1/pristine$spc)
pristine_raw$spcA = log(1/pristine_raw$spc)
pristine_raw$spcA = as.matrix(pristine_raw$spcA)

colored$spcA = log(1/colored$spc)
colored_raw$spcA = log(1/colored_raw$spc)
colored_raw$spcA = as.matrix(colored_raw$spcA)

# plot spectra profiles
par(mfrow = c(1, 2))
matplot(colnames(pristine$spcA),
        t(pristine$spcA),
        main = "Absorbance spectra profile\n (pristine plastics)",
        xlab = "Wavelength (nm)",
        ylab = "Absorbance",
        type = "l",
        lty = 1,
        col = rgb(0.5, 0.5, 0.5,
                  alpha = 0.3))

matplot(colnames(colored$spcA),
        t(colored$spcA),
        main = "Absorbance spectra profile\n (colored plastics)",
        xlab = "Wavelength (nm)",
        ylab = "Absorbance",
        type = "l",
        lty = 1,
        col = rgb(0.1, 0.5, 0.1,
                  alpha = 0.3))


################################################################################
#                     PREPROCESSING TREATMENTS BOILERPLATE                     #
################################################################################

# M1 (raw data) = `raw_pristine.rds` & `raw_colored.rds`
# from `label_pristine_data.R` & `label_colored_data.R`

# M2 (minimal) = `pristine_denoised.rds` & `colored_denoised.rds`
# from `spectra_processing.R`

# M3 (in-between) = SGf + SNV + movav
pristine$spcAmovav = movav(standardNormalVariate(pristine$spcA), w = 11)
colored$spcAmovav = movav(standardNormalVariate(colored$spcA), w = 11)

# M4 (full preprocessing) = SGf + 5nm resample + SNV + movav
oldWavs = as.numeric(colnames(pristine$spcA))
newWavs = seq(min(oldWavs), max(oldWavs), by = 5) # same range for `colored`

pristine$spcAR = resample(pristine$spcA,
                          wav = oldWavs,
                          new.wav = newWavs,
                          interpol = "linear")
pristine$spcARmovav = movav(standardNormalVariate(pristine$spcAR), w = 11)

colored$spcAR = resample(colored$spcA,
                         wav = oldWavs,
                         new.wav = newWavs,
                         interpol = "linear")
colored$spcARmovav = movav(standardNormalVariate(colored$spcAR), w = 11)


## visualize diferences (e.g., pristine raw vs. full preprocessing)============#
par(mfrow = c(1, 2))

matplot(colnames(pristine_raw$spcA),
        t(pristine_raw$spcA),
        main = "Absorbance spectra profile\n (raw)",
        xlab = "Wavelength (nm)",
        ylab = "Absorbance",
        type = "l",
        lty = 1,
        col = rgb(0.5, 0.5, 0.5,
                  alpha = 0.3))

matplot(colnames(pristine$spcARmovav),
        t(pristine$spcARmovav),
        main = "Absorbance spectra profile\n (full preprocessing)",
        xlab = "Wavelength (nm)",
        ylab = "Absorbance",
        type = "l",
        lty = 1,
        col = rgb(0.5, 0.1, 0.1,
                  alpha = 0.3))


################################################################################
#                       SPLIT DATASET FOR 75/25 HOLD-OUT                       #
################################################################################
pristine$strata = interaction(pristine$POLYMER,
                              pristine$MASS_mg, # treatment` flaggin...
                              pristine$SIZE_CODE)

pristine_raw$strata = interaction(pristine_raw$POLYMER,
                                  pristine_raw$MASS_mg,
                                  pristine_raw$SIZE_CODE)

set.seed(1)

datC = pristine |>
    group_by(strata) |>
    sample_frac(0.75)

datV = pristine |>
    filter(!(SAMPLE_ID %in% datC$SAMPLE_ID))


rawC = pristine_raw |> 
    group_by(strata) |> 
    sample_frac(0.75)

rawV = pristine_raw |> 
    filter(!(SAMPLE_ID %in% rawC$SAMPLE_ID))


# check distribution in calibration and validation sets
cat("`POLYMER` distribution in:\n",
    "-- Calibration --",
    paste(capture.output(table(datC$POLYMER)), collapse = "\n"), "\n\n",
    "-- Validation --",
    paste(capture.output(table(datV$POLYMER)), collapse = "\n"))


cat("`MASS_mg` distribution in:\n",
    "-------- Calibration --------",
    paste(capture.output(table(datC$MASS_mg)), collapse = "\n"), "\n\n",
    "-------- Validation --------",
    paste(capture.output(table(datV$MASS_mg)), collapse = "\n"))


cat("`SIZE_CODE` distribution in:\n",
    "-- Calibration --",
    paste(capture.output(table(datC$SIZE_CODE)), collapse = "\n"), "\n\n",
    "-- Validation --",
    paste(capture.output(table(datV$SIZE_CODE)), collapse = "\n"))


# fitting 75/25 hold-out split models
PLSR_mod_mass = plsr(MASS_mg ~ spcA,
                      data = rawC,           # M1 (raw data)
                      method = "oscorespls",
                      ncomp = 50,
                      validation = "CV")

set.seed(21) # ensure publication reproducibility

PLSR_mod_mass2 = plsr(MASS_mg ~ spcA,
                     data = datC,
                     method = "oscorespls", # M2 (minimal preprocessing)
                     ncomp = 50,
                     validation = "CV")

PLSR_mod_mass3 = plsr(MASS_mg ~ spcAmovav,
                      data = datC,
                      method = "oscorespls", # M3 (intermediate preprocessing)
                      ncomp = 50,
                      validation = "CV")

PLSR_mod_mass4 = plsr(MASS_mg ~ spcARmovav,
                     data = datC,
                     method = "oscorespls", # M4 (full preprocessing)
                     ncomp = 50,
                     validation = "CV")


# check optima `ncomps` values and its respective cross-validated RMSEP
RMSEP(PLSR_mod_mass)  # 7 comps...
RMSEP(PLSR_mod_mass2) # 11 comps...
RMSEP(PLSR_mod_mass3) # 12 comps...
RMSEP(PLSR_mod_mass4) # 20 comps...

## visualize it
par(mfrow = c(2, 2), oma = c(3, 1, 0 , 0))

### M1
validationplot(PLSR_mod_mass,
               val.type = "RMSEP",
               main = "MODEL 1",
               cex.main = 2,
               lwd = 3,
               xlab = "",
               cex.lab = 2,
               cex.axis = 1.5,
               mgp = c(2.7, 1, 0))
grid()
abline(v = which.min(RMSEP(PLSR_mod_mass)$val["CV", , ][-1]),
       lty = 4, lwd = 2,
       col = adjustcolor("darkgreen", alpha.f = 0.6))
abline(h = min(RMSEP(PLSR_mod_mass)$val["CV", ,][-1]),
       lty = 4, lwd = 2,
       col = adjustcolor("darkgreen", alpha.f = 0.6))


### M2
validationplot(PLSR_mod_mass2,
               val.type = "RMSEP",
               main = "MODEL 2",
               cex.main = 2,
               lwd = 3,
               xlab = "",
               ylab = "",
               cex.lab = 2,
               cex.axis = 1.5)
grid()
abline(v = which.min(RMSEP(PLSR_mod_mass2)$val["CV", , ][-1]),
       lty = 4, lwd = 2,
       col = adjustcolor("darkgreen", alpha.f = 0.6))
abline(h = min(RMSEP(PLSR_mod_mass2)$val["CV", ,][-1]),
       lty = 4, lwd = 2,
       col = adjustcolor("darkgreen", alpha.f = 0.6))

### M3
validationplot(PLSR_mod_mass3,
               val.type = "RMSEP",
               main = "MODEL 3",
               cex.main = 2,
               lwd = 3,
               cex.lab = 2,
               cex.axis = 1.5,
               mgp = c(2.7, 1, 0))
grid()
abline(v = which.min(RMSEP(PLSR_mod_mass3)$val["CV", , ][-1]),
       lty = 4, lwd = 2,
       col = adjustcolor("darkgreen", alpha.f = 0.6))
abline(h = min(RMSEP(PLSR_mod_mass3)$val["CV", ,][-1]),
       lty = 4, lwd = 2,
       col = adjustcolor("darkgreen", alpha.f = 0.6))

### M4
validationplot(PLSR_mod_mass,
               val.type = "RMSEP",
               main = "MODEL 4",
               cex.main = 2,
               lwd = 3,
               ylab = NA,
               cex.lab = 2,
               cex.axis = 1.5,
               mgp = c(2.7, 1, 0))
grid()
abline(v = which.min(RMSEP(PLSR_mod_mass4)$val["CV", , ][-1]),
       lty = 4, lwd = 2,
       col = adjustcolor("darkgreen", alpha.f = 0.6))
abline(h = min(RMSEP(PLSR_mod_mass4)$val["CV", ,][-1]),
       lty = 4, lwd = 2,
       col = adjustcolor("darkgreen", alpha.f = 0.6))

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
legend("bottom",
       legend = c("CV", "adjusted CV", "Optimal ncomp"),
       lty = c(1, 2, 4),
       col = c("black", "#DF536B", "darkgreen"),
       lwd = c(3, 3, 2),
       bty = "n",
       horiz = TRUE,
       cex = 2,
       inset = -0.03)


# regression coefficients
par(mfrow = c(2, 2))

plot(as.numeric(colnames(pristine_raw$spcA)),
     PLSR_mod_mass$coefficients[,1,10],
     main = "Model 1",
     type = "l",
     xlab = "Wavelength (nm)",
     ylab = "Regression coefficient") +
    abline(h = 0, col = "red", lty = 2)
grid()

plot(as.numeric(colnames(pristine$spcA)),
     PLSR_mod_mass2$coefficients[,1,10],
     main = "Model 2",
     type = "l",
     xlab = "Wavelength (nm)",
     ylab = "Regression coefficient") +
    abline(h = 0, col = "red", lty = 2)
grid()

plot(as.numeric(colnames(pristine$spcAmovav)),
     PLSR_mod_mass3$coefficients[,1,10],
     main = "Model 3",
     type = "l",
     xlab = "Wavelength (nm)",
     ylab = "Regression coefficient") +
    abline(h = 0, col = "red", lty = 2)
grid()

plot(as.numeric(colnames(pristine$spcARmovav)),
     PLSR_mod_mass4$coefficients[,1,10],
     main = "Model 4",
     type = "l",
     xlab = "Wavelength (nm)",
     ylab = "Regression coefficient") +
    abline(h = 0, col = "red", lty = 2)
grid()


# predict on calibration dataset
PLSR_predC = predict(PLSR_mod_mass, ncomp = 7, newdata = rawC$spcA)
PLSR_predC2 = predict(PLSR_mod_mass2, ncomp = 11, newdata = datC$spcA)
PLSR_predC3 = predict(PLSR_mod_mass3, ncomp = 12, newdata = datC$spcAmovav)
PLSR_predC4 = predict(PLSR_mod_mass4, ncomp = 20, newdata = datC$spcARmovav)

# predict on validation dataset
PLSR_predV = predict(PLSR_mod_mass, ncomp = 7, newdata = rawV$spcA)
PLSR_predV2 = predict(PLSR_mod_mass2, ncomp = 11, newdata = datV$spcA)
PLSR_predV3 = predict(PLSR_mod_mass3, ncomp = 12, newdata = datV$spcAmovav)
PLSR_predV4 = predict(PLSR_mod_mass4, ncomp = 20, newdata = datV$spcARmovav)


# plots
plot(log(datC$MASS_mg), log(PLSR_predC),
     main = "Calibration", 
     xlab = "Observed",
     ylab = "Predicted",
     pch = 16,
     ylim = c(0, 10))
abline(0, 1)

plot(log(datV$MASS_mg), log(PLSR_predV),
     main ="Validation",
     xlab = "Observed",
     ylab = "Predicted",
     pch = 16,
     ylim = c(0, 10))
abline(0, 1)


# set validation statistics ===================================================#
ME = function(obs, pred){
    mean(pred - obs, na.rm = T)
}

RMSE = function(obs, pred){
    sqrt(mean((pred - obs)^2, na.rm = T))
}

R2 = function(obs, pred){
    SSE = sum((pred - obs)^2, na.rm = T) # squared sum error
    SST = sum((obs - mean(obs, na.rm = T))^2, na.rm = T) # squares sum total
    R2 = 1 - SSE / SST
    return(R2)
}
#==============================================================================#


# evaluate quality of predictions
## observed responses
calib_obs = c(rep(list(datC$MASS_mg), 3), list(rawC$MASS_mg))

valid_obs = c(rep(list(datV$MASS_mg), 3), list(rawV$MASS_mg))

## group calibration predictions
calib_preds = list(PLSR_predC, PLSR_predC2, PLSR_predC3, PLSR_predC4)
## group validation predictions
valid_preds = list(PLSR_predV, PLSR_predV2, PLSR_predV3, PLSR_predV4)

### compute calibration statistics
calib_stats = mapply(function(obs, pred) {
    c(
        ME   = ME(obs, pred),
        RMSE = RMSE(obs, pred),
        R2   = R2(obs, pred)
    )
}, calib_obs, calib_preds)

### compute validation statistics
valid_stats = mapply(function(obs, pred) {
    c(
        ME   = ME(obs, pred),
        RMSE = RMSE(obs, pred),
        R2   = R2(obs, pred)
    )
}, valid_obs, valid_preds)

### display results
for (i in seq_along(calib_preds)) {
    cat(paste0("\n======= MODEL ", i, " =======\n"))
    cat("##### CALIBRATION #####\n")
    cat(sprintf("ME   : %.7f\n", calib_stats["ME", i]))
    cat(sprintf("RMSE : %.7f\n", calib_stats["RMSE", i]))
    cat(sprintf("R²   : %.7f\n", calib_stats["R2", i]))
    
    cat("\n\n##### VALIDATION #####\n")
    cat(sprintf("ME   : %.7f\n", valid_stats["ME", i]))
    cat(sprintf("RMSE : %.7f\n", valid_stats["RMSE", i]))
    cat(sprintf("R²   : %.7f\n", valid_stats["R2", i]))
} 


################################################################################
#                            K-FOLD CROSS-VALIDATION                           #
################################################################################
cv_plsr_model = function(data, spc_matrix, ncomp, nfolds = 4, seed = 999) {
    set.seed(seed)                                # check `min(table(data$strata))`
    
    data$foldCV = NA
    strata = interaction(data$POLYMER, data$MASS_mg, data$SIZE_CODE)
    
    for (i in unique(strata)) {
        idx = which(strata == i)
        folds = sample(rep(1:nfolds, length.out = length(idx)))
        data$foldCV[idx] = folds
    }
    
    valTable = data.frame(SAMPLE_ID = data$SAMPLE_ID,
                          obs = data$MASS_mg,
                          pred = NA)
    
    for (i in 1:nfolds) {
        train = pristine[data$foldCV != i, ]
        valid = pristine[data$foldCV == i, ]
        
        train_spc = spc_matrix[data$foldCV != i, ]
        valid_spc = spc_matrix[data$foldCV == i, ]
        
        pls_model = plsr(train$MASS_mg ~ train_spc,
                         method = "oscorespls",
                         ncomp = max(ncomp),
                         validation = "none")
        
        preds = predict(pls_model, newdata = valid_spc, ncomp = ncomp)
        
        valTable$pred[valTable$SAMPLE_ID %in% valid$SAMPLE_ID] = drop(preds)
        
        cat("Fold number", i, "done for PLSR.\n")
    }
    
    list(
        ME = ME(valTable$obs, valTable$pred),
        RMSE = RMSE(valTable$obs, valTable$pred),
        R2 = R2(valTable$obs, valTable$pred),
        table = valTable
    )
}

cv_PLSRmod = cv_plsr_model(pristine_raw, pristine_raw$spcA, ncomp = 7)
cv_PLSRmod2 = cv_plsr_model(pristine, pristine$spcA, ncomp = 11)
cv_PLSRmod3 = cv_plsr_model(pristine, pristine$spcAmovav, ncomp = 12)
cv_PLSRmod4 = cv_plsr_model(pristine, pristine$spcARmovav, ncomp = 20)

cv_plsr_results = list(cv_PLSRmod, cv_PLSRmod2, cv_PLSRmod3, cv_PLSRmod4)

for (i in seq_along(cv_plsr_results)) {
    cat(paste0("\n======= PLSR MODEL ", i, " (k-fold CV) =======\n"))
    cat(sprintf("ME   : %.4f\n", cv_plsr_results[[i]]$ME))
    cat(sprintf("RMSE : %.4f\n", cv_plsr_results[[i]]$RMSE))
    cat(sprintf("R²   : %.4f\n", cv_plsr_results[[i]]$R2))
}


################################################################################
#             FULL INTERNAL DATASET TRANING/EXTERNAL VALIDATION                #
################################################################################
set.seed(666)


# "naive" models
PLSR_mod_mass = plsr(MASS_mg ~ spcA,
                     data = pristine_raw,
                     method = "oscorespls",
                     ncomp = 50,
                     validation = "CV")

PLSR_mod_mass2 = plsr(MASS_mg ~ spcA,
                      data = pristine,
                      method = "oscorespls",
                      ncomp = 50,
                      validation = "CV")

PLSR_mod_mass3 = plsr(MASS_mg ~ spcAmovav,
                      data = pristine,
                      method = "oscorespls",
                      ncomp = 50,
                      validation = "CV")

PLSR_mod_mass4 = plsr(MASS_mg ~ spcARmovav,
                      data = pristine,
                      method = "oscorespls",
                      ncomp = 50,
                      validation = "CV")

par(mfrow = c(2, 2))
validationplot(PLSR_mod_mass, val.type = "RMSEP", main = "MODEL 1")
validationplot(PLSR_mod_mass2, val.type = "RMSEP", main = "MODEL 2")
validationplot(PLSR_mod_mass3, val.type = "RMSEP", main = "MODEL 3")
validationplot(PLSR_mod_mass4, val.type = "RMSEP", main = "MODEL 4")

min(RMSEP(PLSR_mod_mass)$val["CV", ,][-1])
min(RMSEP(PLSR_mod_mass2)$val["CV", ,][-1])
min(RMSEP(PLSR_mod_mass3)$val["CV", ,][-1])
min(RMSEP(PLSR_mod_mass4)$val["CV", ,][-1])

which.min(RMSEP(PLSR_mod_mass)$val["CV", , ][-1])
which.min(RMSEP(PLSR_mod_mass2)$val["CV", , ][-1])
which.min(RMSEP(PLSR_mod_mass3)$val["CV", , ][-1])
which.min(RMSEP(PLSR_mod_mass4)$val["CV", , ][-1])


## predicting on new dataset
PLSR_predNew = predict(PLSR_mod_mass, ncomp = 7, newdata = colored_raw$spcA)
PLSR_predNew2 = predict(PLSR_mod_mass2, ncomp = 11, newdata = colored$spcA)
PLSR_predNew3 = predict(PLSR_mod_mass3, ncomp = 12, newdata = colored$spcAmovav)
PLSR_predNew4 = predict(PLSR_mod_mass4, ncomp = 20, newdata = colored$spcARmovav)

external_obs = c(list(raw_new$MASS_mg), rep(list(new$MASS_mg), 3))
external_preds = list(PLSR_predNew, PLSR_predNew2, PLSR_predNew3, PLSR_predNew4)

external_stats = mapply(function(obs, pred) {
    c(
        ME   = ME(obs, pred),
        RMSE = RMSE(obs, pred),
        R2   = R2(obs, pred)
    )
}, external_obs, external_preds)

for (i in seq_along(external_preds)) {
    cat(paste0("\n======= MODEL ", i, " =======\n"))
    cat(sprintf("ME   : %.7f\n", external_stats["ME", i]))
    cat(sprintf("RMSE : %.7f\n", external_stats["RMSE", i]))
    cat(sprintf("R²   : %.7f\n", external_stats["R2", i]))
}


################################################################################
#                            TUNED EXTERNAL VALIDATION                         #
################################################################################
set.seed(42)


colored$strata = interaction(colored$MASS_mg,
                             colored$SIZE_CODE)

colored1 = colored |> 
    group_by(strata) |> 
    sample_frac(0.5)

colored2 = colored |> 
    filter(!(SAMPLE_ID %in% colored1$SAMPLE_ID))

raw_colored$strata = interaction(raw_colored$MASS_mg,
                                 raw_colored$SIZE_CODE)


raw_colored1 = raw_colored |> 
    group_by(strata) |> 
    sample_frac(0.5)


raw_colored2 = raw_colored |> 
    filter(!(SAMPLE_ID %in% raw_colored1$SAMPLE_ID))


# hyperparameter optimization test ============================================#
evaluate_ncomp = function(model, newdata, obs, grid){
    
    out = lapply(grid, function(nc){
        
        pred = predict(model, ncomp = nc, newdata = newdata)
        pred = as.vector(pred)
        
        c(ncomp = nc,
          ME    = ME(obs, pred),
          RMSE  = RMSE(obs, pred),
          R2    = R2(obs, pred))
    })
    
    as.data.frame(do.call(rbind, out))
}

grid = 1:30

res_M1_ext = evaluate_ncomp(model = PLSR_mod_mass,
                            newdata = raw_colored1$spcA,
                            obs = raw_colored1$MASS_mg,
                            grid = grid)
res_M1_ext = res_M1_ext |> 
    mutate(Model = rep("M1", 30), .before = ncomp)


res_M2_ext = evaluate_ncomp(model = PLSR_mod_mass2,
                            newdata = colored1$spcA,
                            obs = colored1$MASS_mg,
                            grid = grid)
res_M2_ext = res_M2_ext |> 
    mutate(Model = rep("M2", 30), .before = ncomp)


res_M3_ext = evaluate_ncomp(model = PLSR_mod_mass3,
                            newdata = colored1$spcAmovav,
                            obs = colored1$MASS_mg,
                            grid = grid)
res_M3_ext = res_M3_ext |> 
    mutate(Model = rep("M3", 30), .before = ncomp)


res_M4_ext = evaluate_ncomp(model = PLSR_mod_mass4,
                            newdata = colored1$spcARmovav,
                            obs = colored1$MASS_mg,
                            grid = grid)
res_M4_ext = res_M4_ext |> 
    mutate(Model = rep("M4", 30), .before = ncomp)


results = rbind(res_M1_ext, res_M2_ext,
                res_M3_ext, res_M4_ext,
                make.row.names = F)


## check optima parameters
results |> 
    group_by(Model) |> 
    filter(RMSE == min(RMSE))
#==============================================================================#


## tuned models
PLSR_mod_tuned = plsr(MASS_mg ~ spcA,
                      data = pristine_raw,
                      method = "oscorespls",
                      ncomp = 10,
                      validation = "CV")

PLSR_mod_tuned2 = plsr(MASS_mg ~ spcA,
                       data = pristine,
                       method = "oscorespls",
                       ncomp = 12,
                       validation = "CV")

PLSR_mod_tuned3 = plsr(MASS_mg ~ spcAmovav,
                       data = pristine,
                       method = "oscorespls",
                       ncomp = 6,
                       validation = "CV")

PLSR_mod_tuned4 = plsr(MASS_mg ~ spcARmovav,
                       data = pristine,
                       method = "oscorespls",
                       ncomp = 6,
                       validation = "CV")


## final predictions
### M1                                       -> pay attention to the ncomp values
PLSR_pred_tuned = predict(PLSR_mod_tuned, ncomp = 10, newdata = raw_colored2$spcA)
cat(paste0("\n======= PLSR MODEL 1 (fine-tuned) =======\n",
           sprintf("ME   : %.2f\n", ME(raw_colored2$MASS_mg, PLSR_pred_tuned)),
           sprintf("RMSE : %.2f\n", RMSE(raw_colored2$MASS_mg, PLSR_pred_tuned)),
           sprintf("R²   : %.2f\n", R2(raw_colored2$MASS_mg, PLSR_pred_tuned))
))

### M2
PLSR_pred_tuned2 = predict(PLSR_mod_tuned2, ncomp = 12, newdata = colored2$spcA)
cat(paste0("\n======= PLSR MODEL 2 (fine-tuned) =======\n",
           sprintf("ME   : %.2f\n", ME(colored2$MASS_mg, PLSR_pred_tuned2)),
           sprintf("RMSE : %.2f\n", RMSE(colored2$MASS_mg, PLSR_pred_tuned2)),
           sprintf("R²   : %.2f\n", R2(colored2$MASS_mg, PLSR_pred_tuned2))
))

### M3
PLSR_pred_tuned3 = predict(PLSR_mod_tuned3, ncomp = 6, newdata = colored2$spcAmovav)
cat(paste0("\n======= PLSR MODEL 3 (fine-tuned) =======\n",
           sprintf("ME   : %.2f\n", ME(colored2$MASS_mg, PLSR_pred_tuned3)),
           sprintf("RMSE : %.2f\n", RMSE(colored2$MASS_mg, PLSR_pred_tuned3)),
           sprintf("R²   : %.2f\n", R2(colored2$MASS_mg, PLSR_pred_tuned3))
))

### M4
PLSR_pred_tuned = predict(PLSR_mod_tuned4, ncomp = 6, newdata = colored2$spcARmovav)
cat(paste0("\n======= PLSR MODEL 4 (fine-tuned) =======\n",
           sprintf("ME   : %.2f\n", ME(colored2$MASS_mg, PLSR_pred_tuned)),
           sprintf("RMSE : %.2f\n", RMSE(colored2$MASS_mg, PLSR_pred_tuned)),
           sprintf("R²   : %.2f\n", R2(colored2$MASS_mg, PLSR_pred_tuned))
))


# residual visualization
residTUNED = PLSR_pred_tuned - raw_new2$MASS_mg
residTUNED2 = PLSR_pred_tuned2 - new2$MASS_mg
residTUNED3 = PLSR_pred_tuned3 - new2$MASS_mg
residTUNED4 = PLSR_pred_tuned4 - new2$MASS_mg

residTUNEDlog = log(PLSR_pred_tuned) - log(raw_new2$MASS_mg)
residTUNEDlog2 = log(PLSR_pred_tuned2) - log(new2$MASS_mg)
residTUNEDlog3 = log(PLSR_pred_tuned3) - log(new2$MASS_mg)
residTUNEDlog4 = log(PLSR_pred_tuned4) - log(new2$MASS_mg)


df_m1 = data.frame(
    Model = "M1 (Raw data)",
    Observed_Mass = raw_new2$MASS_mg, # M1 uses raw_new
    Residual = as.vector(residTUNED4)
)

df_m2 = data.frame(
    Model = "M2 (Minimal preprocessing)",
    Observed_Mass = new2$MASS_mg, # M2 uses new
    Residual = as.vector(residTUNED2)
)

df_m3 = data.frame(
    Model = "M3 (Intermediate preprocessing)",
    Observed_Mass = new2$MASS_mg, # M3 uses new
    Residual = as.vector(residTUNED3)
)

df_m4 = data.frame(
    Model = "M4 (Full preprocessing)",
    Observed_Mass = new2$MASS_mg, # M4 uses new
    Residual = as.vector(residTUNED)
)


# combine data frames
plot_data = bind_rows(df_m1, df_m2, df_m3, df_m4)

# model factor in correct order of complexity
model_order = c("M1 (Raw data)", "M2 (Minimal preprocessing)",
                "M3 (Intermediate preprocessing)", "M4 (Full preprocessing)")
plot_data$Model = factor(plot_data$Model, levels = model_order)

# combined boxplot
res = ggplot(plot_data, aes(x = factor(Observed_Mass), y = Residual)) +
    geom_boxplot(outlier.size = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
    facet_wrap(~ Model, ncol = 2) +
    labs(
        x = NULL,
        y = "Residuals (Predicted - Observed Mass)",
        title = NULL) +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          strip.text = element_text(size = 11),
          axis.title.y = element_text(size = 13))


# data frames for relative residuals (log scale)
df_m1_log = data.frame(
    Model = "M1 (Raw data)",
    Observed_Mass = raw_new2$MASS_mg,
    Relative_Residual = as.vector(residTUNEDlog4)  # M1 = residTUNEDlog
)

df_m2_log = data.frame(
    Model = "M2 (Minimal preprocessing)", 
    Observed_Mass = new2$MASS_mg,
    Relative_Residual = as.vector(residTUNEDlog2)
)

df_m3_log = data.frame(
    Model = "M3 (Intermediate preprocessing)",
    Observed_Mass = new2$MASS_mg, 
    Relative_Residual = as.vector(residTUNEDlog3)
)

df_m4_log = data.frame(
    Model = "M4 (Full preprocessing)",
    Observed_Mass = new2$MASS_mg,
    Relative_Residual = as.vector(residTUNEDlog)  # M4 = residTUNEDlog4
)

# combine and process
plot_data_log = bind_rows(df_m1_log, df_m2_log, df_m3_log, df_m4_log)

# set factor order
model_order = c("M1 (Raw data)", "M2 (Minimal preprocessing)",
                "M3 (Intermediate preprocessing)", "M4 (Full preprocessing)")


# relative residuals plot
res2 = ggplot(plot_data_log, aes(x = factor(Observed_Mass), y = Relative_Residual)) +
    geom_boxplot(outlier.size = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
    facet_wrap(~ Model, ncol = 2) +
    labs(
        x = "Mass (mg)",
        y = "Relative Residuals [log(Predicted) - log(Observed)]",
        title = NULL) +
    theme(axis.text.x = element_text(size = 13),
          axis.text.y = element_text(size = 12),
          strip.text = element_text(size = 11),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 13))


res / res2 + plot_annotation(tag_levels = 'a',
                             tag_prefix = '(',
                             tag_suffix = ')') &
    theme(plot.tag = element_text(size = 22))

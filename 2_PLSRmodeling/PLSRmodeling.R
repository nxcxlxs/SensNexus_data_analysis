# load processed data
pristine = readRDS("C:/nico/Dissertação/SENSNEXUS_data/analysis_ready_data/deNoised_data.rds")
pristine_raw = readRDS("C:/nico/Dissertação/SENSNEXUS_data/preprocessed_data/datsoil.rds")
colored = readRDS("C:/nico/Dissertação/SENSNEXUS_data/analysis_ready_data/deNoised_newData.rds") 
colored_raw = readRDS("C:/nico/Dissertação/SENSNEXUS_data/preprocessed_data/datsoil2.rds")

# convert spectra to absorbance
pristine$spcA = log(1/pristine$spc)
pristine_raw$spcA = log(1/pristine_raw$spc)
pristine_raw$spcA = as.matrix(pristine_raw$spcA)

colored$spcA = log(1/colored$spc)
colored_raw$spcA = log(1/colored_raw$spc)
raw_new$spcA = as.matrix(raw_new$spcA)

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
#                           PREPROCESSING TREATMENTS                           #
################################################################################
require(prospectr)
# M1 (raw data) = `datsoil.rds` & `datsoil2.rds`
# from `compile_design_data.R`, both loaded as `pristine_raw` & `colored_raw`

# M2 (minimal) = `deNoised_data.rds` & `deNoised_newData.rds`
# from `noise_scatter_correction.R`, both loaded as `pristine` & `colored`

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

require(dplyr)
set.seed(1)
datC = pristine |>
    group_by(strata) |>
    sample_frac(0.75)


datV = pristine |>
    filter(!(SAMPLE_ID %in% datC$SAMPLE_ID))


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


# set basic validation statistics =============================================#
ME = function(obs, pred){
        mean(pred - obs, na.rm = T)
}

RMSE = function(obs, pred){
        sqrt(mean((pred - obs)^2, na.rm = T))
}

R2 = function(obs, pred){
        SSE = sum((pred - obs)^2, na.rm = T) # squared error sum
        SST = sum((obs - mean(obs, na.rm = T))^2, na.rm = T) # squares total sum
        R2 = 1 - SSE / SST
        return(R2)
}
#=============THIS SHOULD BE INCLUDED AFTER THE FIRST PREDICTIONS==============#

# fitting PLSR
set.seed(21)
require(pls)
PLSR_mod_mass = plsr(MASS_mg ~ spcARmovav,
                data = datC,
                method = "oscorespls",      # Wadoux's approach...
                ncomp = 50,
                validation = "CV")


PLSR_mod_mass2 = plsr(MASS_mg ~ spcA,
                     data = datC,
                     method = "oscorespls", # my approach...
                     ncomp = 50,
                     validation = "CV")

PLSR_mod_mass3 = plsr(MASS_mg ~ spcAmovav,
                      data = datC,
                      method = "oscorespls", # in-between Wadoux's and mine...
                      ncomp = 50,
                      validation = "CV")

raw$strata = interaction(raw$POLYMER,
                         raw$MASS_mg,
                         raw$SIZE_CODE)

set.seed(1)

rawC = raw |> 
    group_by(strata) |> 
    sample_frac(0.75)

rawV = raw |> 
    filter(!(SAMPLE_ID %in% rawC$SAMPLE_ID))

PLSR_mod_mass4 = plsr(MASS_mg ~ spcA,
                      data = rawC,           # raw spectra...
                      method = "oscorespls",
                      ncomp = 50,
                      validation = "CV")

par(mfrow = c(2, 2))
validationplot(PLSR_mod_mass, val.type = "RMSEP", main = "MODEL 1") # M4
validationplot(PLSR_mod_mass2, val.type = "RMSEP", main = "MODEL 2")
validationplot(PLSR_mod_mass3, val.type = "RMSEP", main = "MODEL 3")
validationplot(PLSR_mod_mass4, val.type = "RMSEP", main = "MODEL 4") # M1

min(RMSEP(PLSR_mod_mass)$val["CV", ,][-1])
min(RMSEP(PLSR_mod_mass2)$val["CV", ,][-1])
min(RMSEP(PLSR_mod_mass3)$val["CV", ,][-1])
min(RMSEP(PLSR_mod_mass4)$val["CV", ,][-1])

which.min(RMSEP(PLSR_mod_mass)$val["CV", , ][-1])
which.min(RMSEP(PLSR_mod_mass2)$val["CV", , ][-1])
which.min(RMSEP(PLSR_mod_mass3)$val["CV", , ][-1])
which.min(RMSEP(PLSR_mod_mass4)$val["CV", , ][-1])

cat("The smallest RMSEP is:", min(RMSEP(PLSR_mod_mass4)$val["CV", ,][-1]))
cat("The optimal number of components is:",
    which.min(RMSEP(PLSR_mod_mass4)$val["CV", , ][-1]))
# ...the 'elbow rule' should be good enough  ¯\_(ツ)_/¯

# par(mfrow = c(2, 2), oma = c(3, 1, 0 , 0))
# 
# validationplot(PLSR_mod_mass4,
#                val.type = "RMSEP",
#                main = "MODEL 1",
#                cex.main = 2,
#                lwd = 3,
#                xlab = "",
#                cex.lab = 2,
#                cex.axis = 1.5,
#                mgp = c(2.7, 1, 0))
# grid()
# abline(v = 7, lty = 4, lwd = 2, col = adjustcolor("darkgreen", alpha.f = 0.6))
# abline(h = 472.0287, lty = 4, lwd = 2, col = adjustcolor("darkgreen", alpha.f = 0.6))
# 
# 
# validationplot(PLSR_mod_mass2,
#                val.type = "RMSEP",
#                main = "MODEL 2",
#                cex.main = 2,
#                lwd = 3,
#                xlab = "",
#                ylab = "",
#                cex.lab = 2,
#                cex.axis = 1.5)
# grid()
# abline(v = 11, lty = 4, lwd = 2, col = adjustcolor("darkgreen", alpha.f = 0.6))
# abline(h = 478.4669, lty = 4, lwd = 2, col = adjustcolor("darkgreen", alpha.f = 0.6))
# 
# validationplot(PLSR_mod_mass3,
#                val.type = "RMSEP",
#                main = "MODEL 3",
#                cex.main = 2,
#                lwd = 3,
#                cex.lab = 2,
#                cex.axis = 1.5,
#                mgp = c(2.7, 1, 0))
# grid()
# abline(v = 12, lty = 4, lwd = 2, col = adjustcolor("darkgreen", alpha.f = 0.6))
# abline(h = 434.3338, lty = 4, lwd = 2, col = adjustcolor("darkgreen", alpha.f = 0.6))
# 
# 
# validationplot(PLSR_mod_mass,
#                val.type = "RMSEP",
#                main = "MODEL 4",
#                cex.main = 2,
#                lwd = 3,
#                ylab = NA,
#                cex.lab = 2,
#                cex.axis = 1.5,
#                mgp = c(2.7, 1, 0))
# grid()
# abline(v = 20, lty = 4, lwd = 2, col = adjustcolor("darkgreen", alpha.f = 0.6))
# abline(h = 425.7052, lty = 4, lwd = 2, col = adjustcolor("darkgreen", alpha.f = 0.6))
# 
# par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
# plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
# legend("bottom", 
#        legend = c("CV", "adjusted CV", "Optimal ncomp"),
#        lty = c(1, 2, 4),
#        col = c("black", "#DF536B", "darkgreen"),
#        lwd = c(3, 3, 2),
#        bty = "n",
#        horiz = TRUE,
#        cex = 2,
#        inset = -0.03)

RMSEP(PLSR_mod_mass)  # 20 comps...
RMSEP(PLSR_mod_mass2) # 11 comps...
RMSEP(PLSR_mod_mass3) # 12 comps...
RMSEP(PLSR_mod_mass4) # 7 comps...

par(mfrow = c(1, 2))
plot(PLSR_mod_mass4,
     ncomp = 10,
     main = "12 components",
     xlab = "Observed",
     ylab = "Predicted")

plot(PLSR_mod_mass4,
     ncomp = 7,
     main = "9 components",
     xlab = "Observed",
     ylab = "Predicted")

# three first loadings
plot(PLSR_mod_mass4,
     "loadings",
     comps = 1:3,
     xlab = "Index of the wavelength",
     ylab = "Loading value") +
        legend("top",
               legend = c("comp1", "comp2", "comp3"),
               col = 1:3,
               lty = 1,
               bty = 'n',
               cex = 0.8)

# regression coefficients
par(mfrow = c(2, 2))

plot(as.numeric(colnames(raw$spcA)), PLSR_mod_mass4$coefficients[,1,10],
     main = "Model 1", # M4
     type = "l",
     xlab = "Wavelength (nm)",
     ylab = "Regression coefficient") +
    abline(h = 0, col = "red", lty = 2)
grid()

plot(as.numeric(colnames(data$spcA)), PLSR_mod_mass2$coefficients[,1,10],
     main = "Model 2",
     type = "l",
     xlab = "Wavelength (nm)",
     ylab = "Regression coefficient") +
    abline(h = 0, col = "red", lty = 2)
grid()

plot(as.numeric(colnames(data$spcAmovav)), PLSR_mod_mass3$coefficients[,1,10],
     main = "Model 3",
     type = "l",
     xlab = "Wavelength (nm)",
     ylab = "Regression coefficient") +
    abline(h = 0, col = "red", lty = 2)
grid()

plot(as.numeric(colnames(data$spcARmovav)), PLSR_mod_mass$coefficients[,1,10],
     main = "Model 4", # M1
     type = "l",
     xlab = "Wavelength (nm)",
     ylab = "Regression coefficient") +
    abline(h = 0, col = "red", lty = 2)
grid()

# predict on calibration dataset
PLSR_predC = predict(PLSR_mod_mass, ncomp = 20, newdata = datC$spcARmovav)
PLSR_predC2 = predict(PLSR_mod_mass2, ncomp = 11, newdata = datC$spcA)
PLSR_predC3 = predict(PLSR_mod_mass3, ncomp = 12, newdata = datC$spcAmovav)
PLSR_predC4 = predict(PLSR_mod_mass4, ncomp = 7, newdata = rawC$spcA)

# predict on validation dataset
PLSR_predV = predict(PLSR_mod_mass, ncomp = 20, newdata = datV$spcARmovav)
PLSR_predV2 = predict(PLSR_mod_mass2, ncomp = 11, newdata = datV$spcA)
PLSR_predV3 = predict(PLSR_mod_mass3, ncomp = 12, newdata = datV$spcAmovav)
PLSR_predV4 = predict(PLSR_mod_mass4, ncomp = 7, newdata = rawV$spcA)

# predict on new dataset
PLSR_predNew = predict(PLSR_mod_mass, ncomp = 20, newdata = new$spcARmovav)
PLSR_predNew2 = predict(PLSR_mod_mass2, ncomp = 11, newdata = new$spcA)
PLSR_predNew3 = predict(PLSR_mod_mass3, ncomp = 12, newdata = new$spcAmovav)
PLSR_predNew4 = predict(PLSR_mod_mass4, ncomp = 7, newdata = raw_new$spcA)

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
    
    cat("##### VALIDATION #####\n")
    cat(sprintf("ME   : %.7f\n", valid_stats["ME", i]))
    cat(sprintf("RMSE : %.7f\n", valid_stats["RMSE", i]))
    cat(sprintf("R²   : %.7f\n", valid_stats["R2", i]))
} 


#### EXTERNAL VALIDATION (colored rigid household plastics - cRHP - dataset)
external_obs = c(rep(list(new$MASS_mg), 3), list(raw_new$MASS_mg))
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
# ## calibration
# ME(rawC$MASS_mg, PLSR_predC4)
# RMSE(rawC$MASS_mg, PLSR_predC4)
# R2(rawC$MASS_mg, PLSR_predC4)
# residC = PLSR_predC4 - rawC$MASS_mg # absolute residual = real error
# residClog = log(PLSR_predC4) - log(rawC$MASS_mg) # log residual = relative error (emphasizes all scales equally)
# boxplot(residClog ~ rawC$MASS_mg,
#         main = "Residuals by Mass (Calibration)") # absolute residuals overemphasize
#                                                   # large errors at higher response values
# 
# ## validation
# ME(rawV$MASS_mg, PLSR_predV4)
# RMSE(rawV$MASS_mg, PLSR_predV4)
# R2(rawV$MASS_mg, PLSR_predV4)
# residV = PLSR_predV4 - rawV$MASS_mg
# residVlog = log(PLSR_predV4) - log(rawV$MASS_mg)
# boxplot(residVlog ~ rawV$MASS_mg,
#         main = "Residuals by Mass (Validation)")

# residuals boxplots
## calibration
residC = PLSR_predC - datC$MASS_mg
residC2 = PLSR_predC2 - datC$MASS_mg
residC3 = PLSR_predC3 - datC$MASS_mg
residC4 = PLSR_predC4 - rawC$MASS_mg

# valid_idx = which(PLSR_predC > 0 & datC$MASS_mg > 0) if needed...

residClog = log(PLSR_predC) - log(datC$MASS_mg)
residClog2 = log(PLSR_predC2) - log(datC$MASS_mg)
residClog3 = log(PLSR_predC3) - log(datC$MASS_mg)
residClog4 = log(PLSR_predC4) - log(rawC$MASS_mg)

par(mfrow = c(1, 2))
boxplot(residC ~ datC$MASS_mg,
        main = "Residuals by Mass (Calibration)")

boxplot(residClog ~ datC$MASS_mg, # ATTENTION to the log transformation
        main = "Relative residuals by Mass (Calibration)")

## validation
residV = PLSR_predV - datV$MASS_mg
residV2 = PLSR_predV2 - datV$MASS_mg
residV3 = PLSR_predV3 - datV$MASS_mg
residV4 = PLSR_predV4 - rawV$MASS_mg

residVlog = log(PLSR_predV) - log(datV$MASS_mg)
residVlog2 = log(PLSR_predV2) - log(datV$MASS_mg)
residVlog3 = log(PLSR_predV3) - log(datV$MASS_mg)
residVlog4 = log(PLSR_predV4) - log(rawV$MASS_mg)

boxplot(residV4 ~ rawV$MASS_mg,
        main = "Residuals by Mass (Validation)")

boxplot(residVlog4 ~ rawV$MASS_mg,
        main = "Relative residuals by Mass (Validation)")

################################################################################

# there's the 10-fold cross-validation approach left...      # ALWAYS check `min(table(data$strata))`
cv_plsr_model = function(data, spc_matrix, ncomp, nfolds = 4, seed = 999) {
    set.seed(seed)
    
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
        train = data[data$foldCV != i, ]
        valid = data[data$foldCV == i, ]
        
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

cv_PLSRmod = cv_plsr_model(data, data$spcARmovav, ncomp = 20)
cv_PLSRmod2 = cv_plsr_model(data, data$spcA, ncomp = 11)
cv_PLSRmod3 = cv_plsr_model(data, data$spcAmovav, ncomp = 12)
cv_PLSRmod4 = cv_plsr_model(raw, raw$spcA, ncomp = 7)

cv_plsr_results = list(cv_PLSRmod, cv_PLSRmod2, cv_PLSRmod3, cv_PLSRmod4)

for (i in seq_along(cv_plsr_results)) {
    cat(paste0("\n======= PLSR MODEL ", i, " (10-fold CV) =======\n"))
    cat(sprintf("ME   : %.4f\n", cv_plsr_results[[i]]$ME))
    cat(sprintf("RMSE : %.4f\n", cv_plsr_results[[i]]$RMSE))
    cat(sprintf("R²   : %.4f\n", cv_plsr_results[[i]]$R2))
}


#==============================================================================#
#                     FULL INTERNAL DATASET TRANING!                           #
#==============================================================================#
set.seed(666)


## "naive" models
PLSR_mod_mass = plsr(MASS_mg ~ spcARmovav,
                     data = data,
                     method = "oscorespls",
                     ncomp = 20,
                     validation = "CV")

PLSR_mod_mass2 = plsr(MASS_mg ~ spcA,
                      data = data,
                      method = "oscorespls",
                      ncomp = 11,
                      validation = "CV")

PLSR_mod_mass3 = plsr(MASS_mg ~ spcAmovav,
                      data = data,
                      method = "oscorespls",
                      ncomp = 12,
                      validation = "CV")

PLSR_mod_mass4 = plsr(MASS_mg ~ spcA,
                      data = raw,
                      method = "oscorespls",
                      ncomp = 7,
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

# predicting on new dataset
PLSR_predNew = predict(PLSR_mod_mass, ncomp = 20, newdata = new$spcARmovav)
PLSR_predNew2 = predict(PLSR_mod_mass2, ncomp = 11, newdata = new$spcA)
PLSR_predNew3 = predict(PLSR_mod_mass3, ncomp = 12, newdata = new$spcAmovav)
PLSR_predNew4 = predict(PLSR_mod_mass4, ncomp = 7, newdata = raw_new$spcA)

external_obs = c(rep(list(new$MASS_mg), 3), list(raw_new$MASS_mg))
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

## "tuned-on-external" models
PLSR_mod_leak = plsr(MASS_mg ~ spcARmovav,
                     data = data,
                     method = "oscorespls",
                     ncomp = 6,
                     validation = "CV")

PLSR_mod_leak2 = plsr(MASS_mg ~ spcA,
                      data = data,
                      method = "oscorespls",
                      ncomp = 12,
                      validation = "CV")

PLSR_mod_leak3 = plsr(MASS_mg ~ spcAmovav,
                      data = data,
                      method = "oscorespls",
                      ncomp = 6,
                      validation = "CV")

PLSR_mod_leak4 = plsr(MASS_mg ~ spcA,
                      data = raw,
                      method = "oscorespls",
                      ncomp = 10,
                      validation = "CV")

# check prediction performance
PLSR_predLeak = predict(PLSR_mod_leak, ncomp = 6, newdata = new$spcARmovav)
ME(new$MASS_mg, PLSR_predLeak)
RMSE(new$MASS_mg, PLSR_predLeak)
R2(new$MASS_mg, PLSR_predLeak)

PLSR_predLeak2 = predict(PLSR_mod_leak2, ncomp = 12, newdata = new$spcA)
ME(new$MASS_mg, PLSR_predLeak2)
RMSE(new$MASS_mg, PLSR_predLeak2)
R2(new$MASS_mg, PLSR_predLeak2)

PLSR_predLeak3 = predict(PLSR_mod_leak3, ncomp = 6, newdata = new$spcAmovav)
ME(new$MASS_mg, PLSR_predLeak3)
RMSE(new$MASS_mg, PLSR_predLeak3)
R2(new$MASS_mg, PLSR_predLeak3)

PLSR_predLeak4 = predict(PLSR_mod_leak4, ncomp = 10, newdata = raw_new$spcA)
ME(new$MASS_mg, PLSR_predLeak4)
RMSE(new$MASS_mg, PLSR_predLeak4)
R2(new$MASS_mg, PLSR_predLeak4)

# residual visualization
residLEAK = PLSR_predLeak - new$MASS_mg
residLEAK2 = PLSR_predLeak2 - new$MASS_mg
residLEAK3 = PLSR_predLeak3 - new$MASS_mg
residLEAK4 = PLSR_predLeak4 - raw_new$MASS_mg

residLEAKlog = log(PLSR_predLeak) - log(new$MASS_mg)
residLEAKlog2 = log(PLSR_predLeak2) - log(new$MASS_mg)
residLEAKlog3 = log(PLSR_predLeak3) - log(new$MASS_mg)
residLEAKlog4 = log(PLSR_predLeak4) - log(raw_new$MASS_mg)

# ggplot(new, aes(x = factor(MASS_mg), y = residLEAK2)) +
#     geom_boxplot() +
#     labs(title = "MODEL 2",
#          x = "Mass (mg)",
#          y = "Residuals")
# 
# ggplot(new, aes(x = factor(MASS_mg), y = residLEAKlog4)) +
#     geom_boxplot() +
#     labs(title = "MODEL 2",
#          x = "Mass (mg)",
#          y = "Relative residuals")

# data frame for each model, ensuring the correct 'Observed_Mass' is used
df_m1 = data.frame(
    Model = "M1 (Raw data)",
    Observed_Mass = raw_new$MASS_mg, # M1 uses raw_new
    Residual = residLEAK4  # M1 = residLEAK4
)

df_m2 = data.frame(
    Model = "M2 (Minimal preprocessing)",
    Observed_Mass = new$MASS_mg, # M2 uses new
    Residual = residLEAK2
)

df_m3 = data.frame(
    Model = "M3 (Intermediate preprocessing)",
    Observed_Mass = new$MASS_mg, # M3 uses new
    Residual = residLEAK3
)

df_m4 = data.frame(
    Model = "M4 (Full preprocessing)",
    Observed_Mass = raw_new$MASS_mg, # M4 uses raw_new
    Residual = residLEAK  # M4 = residLEAK
)

# combine data frames
plot_data = bind_rows(df_m1, df_m2, df_m3, df_m4)

# model factor in correct order of complexity
model_order = c("M1 (Raw data)", "M2 (Minimal preprocessing)",
                "M3 (Intermediate preprocessing)", "M4 (Full preprocessing)")
plot_data$Model = factor(plot_data$Model, levels = model_order)

plot_data_long = plot_data  |> 
    pivot_longer(
        cols = c(MASS_mg.10.comps, MASS_mg.12.comps, MASS_mg.6.comps),
        names_to = "Model_Type", 
        values_to = "Residual"
    )  |> 
    filter(!is.na(Residual))

# combined boxplot
res = ggplot(plot_data_long, aes(x = factor(Observed_Mass), y = Residual)) +
    geom_boxplot(outlier.size = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
    facet_wrap(~ Model, ncol = 2) +
    labs(
        x = NULL,
        y = "Residuals (Predicted - Observed Mass)",
        title = NULL
    ) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())

#===============================================================================

# data frames for relative residuals (log scale)
df_m1_log = data.frame(
    Model = "M1 (Raw data)",
    Observed_Mass = raw_new$MASS_mg,
    Relative_Residual = residLEAKlog4  # M1 = residLEAKlog4
)

df_m2_log = data.frame(
    Model = "M2 (Minimal preprocessing)", 
    Observed_Mass = new$MASS_mg,
    Relative_Residual = residLEAKlog2
)

df_m3_log = data.frame(
    Model = "M3 (Intermediate preprocessing)",
    Observed_Mass = new$MASS_mg, 
    Relative_Residual = residLEAKlog3
)

df_m4_log = data.frame(
    Model = "M4 (Full preprocessing)",
    Observed_Mass = raw_new$MASS_mg,
    Relative_Residual = residLEAKlog  # M4 = residLEAKlog
)

# combine and process
plot_data_log = bind_rows(df_m1_log, df_m2_log, df_m3_log, df_m4_log)

# reshape from wide to long format
plot_data_log_long = plot_data_log |> 
    pivot_longer(
        cols = starts_with("MASS_mg"),  # Adjust if column names differ
        names_to = "Model_Type", 
        values_to = "Relative_Residual"
    ) |> 
    filter(!is.na(Relative_Residual))

# set factor order
model_order = c("M1 (Raw data)", "M2 (Minimal preprocessing)",
                "M3 (Intermediate preprocessing)", "M4 (Full preprocessing)")
plot_data_log_long$Model = factor(plot_data_log_long$Model, levels = model_order)

# relative residuals plot
res2 = ggplot(plot_data_log_long, aes(x = factor(Observed_Mass), y = Relative_Residual)) +
    geom_boxplot(outlier.size = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
    facet_wrap(~ Model, ncol = 2) +
    labs(
        x = "Mass (mg)",
        y = "Relative Residuals [log(Predicted) - log(Observed)]",
        title = NULL)

require(patchwork)

res / res2 + plot_annotation(tag_levels = 'a',
                             tag_prefix = '(',
                             tag_suffix = ')')

################################################################################
#                    MODELING FOR SINGLE POLYMERS QUANTIFICATION               #
################################################################################

PP_data = data |> 
    filter(POLYMER == "PP")

PP_raw = raw |> 
    filter(POLYMER == "PP")

PVC_data = data |> 
    filter(POLYMER == "PVC")

PVC_raw = raw |> 
    filter(POLYMER == "PVC")

PET_data = data |> 
    filter(POLYMER == "PET")

PET_raw = raw |> 
    filter(POLYMER == "PET")

PE_data = data |> 
    filter(POLYMER == "PE")

PE_raw = raw |> 
    filter(POLYMER == "PE")

new = new |> 
    mutate(PP_mg = MASS_mg/4,
           PVC_mg = MASS_mg/4,
           PET_mg = MASS_mg/4,
           PE_mg = MASS_mg/4)

raw_new = raw_new |> 
    mutate(PP_mg = MASS_mg/4,
           PVC_mg = MASS_mg/4,
           PET_mg = MASS_mg/4,
           PE_mg = MASS_mg/4)

set.seed(777)

PLSR_mod_PP = plsr(MASS_mg ~ spcARmovav,
                     data = PP_data,
                     method = "oscorespls",
                     ncomp = 50,
                     validation = "CV")

PLSR_mod_PP2 = plsr(MASS_mg ~ spcA,
                      data = PP_data,
                      method = "oscorespls",
                      ncomp = 50,
                      validation = "CV")

PLSR_mod_PP3 = plsr(MASS_mg ~ spcAmovav,
                      data = PP_data,
                      method = "oscorespls",
                      ncomp = 50,
                      validation = "CV")

PLSR_mod_PP4 = plsr(MASS_mg ~ spcA,
                      data = PP_raw,
                      method = "oscorespls",
                      ncomp = 50,
                      validation = "CV")

which.min(RMSEP(PLSR_mod_PP)$val["CV", , ][-1]) # 3
which.min(RMSEP(PLSR_mod_PP2)$val["CV", , ][-1]) # 2
which.min(RMSEP(PLSR_mod_PP3)$val["CV", , ][-1]) # 8
which.min(RMSEP(PLSR_mod_PP4)$val["CV", , ][-1]) # 2

## polymer-specific CV
cv_PLSR_PP = cv_plsr_model(PP_data, PP_data$spcARmovav, ncomp = 3)
cv_PLSR_PP2 = cv_plsr_model(PP_data, PP_data$spcA, ncomp = 2)
cv_PLSR_PP3 = cv_plsr_model(PP_data, PP_data$spcAmovav, ncomp = 8)
cv_PLSR_PP4 = cv_plsr_model(PP_raw, PP_raw$spcA, ncomp = 2)

cv_plsr_PPresults = list(cv_PLSR_PP, cv_PLSR_PP2, cv_PLSR_PP3, cv_PLSR_PP4)

for (i in seq_along(cv_plsr_PPresults)) {
    cat(paste0("\n======= PLSR MODEL ", i, " (10-fold CV) =======\n"))
    cat(sprintf("ME   : %.4f\n", cv_plsr_PPresults[[i]]$ME))
    cat(sprintf("RMSE : %.4f\n", cv_plsr_PPresults[[i]]$RMSE))
    cat(sprintf("R²   : %.4f\n", cv_plsr_PPresults[[i]]$R2))
}


## polymer-specific external validation
PLSR_predPP = predict(PLSR_mod_PP, ncomp = 3, newdata = new$spcARmovav)
PLSR_predPP2 = predict(PLSR_mod_PP2, ncomp = 2, newdata = new$spcA)
PLSR_predPP3 = predict(PLSR_mod_PP3, ncomp = 8, newdata = new$spcAmovav)
PLSR_predPP4 = predict(PLSR_mod_PP4, ncomp = 2, newdata = raw_new$spcA)

PP_obs = c(rep(list(new$PP_mg), 3), list(raw_new$PP_mg))
PP_preds = list(PLSR_predPP, PLSR_predPP2, PLSR_predPP3, PLSR_predPP4)

PP_stats = mapply(function(obs, pred) {
    c(
        ME   = ME(obs, pred),
        RMSE = RMSE(obs, pred),
        R2   = R2(obs, pred)
    )
}, PP_obs, PP_preds)

for (i in seq_along(PP_preds)) {
    cat(paste0("\n======= MODEL ", i, " =======\n"))
    cat(sprintf("ME   : %.7f\n", PP_stats["ME", i]))
    cat(sprintf("RMSE : %.7f\n", PP_stats["RMSE", i]))
    cat(sprintf("R²   : %.7f\n", PP_stats["R2", i]))
}

#==============================================================================#

PLSR_mod_PVC = plsr(MASS_mg ~ spcARmovav,
                   data = PVC_data,
                   method = "oscorespls",
                   ncomp = 50,
                   validation = "CV")

PLSR_mod_PVC2 = plsr(MASS_mg ~ spcA,
                    data = PVC_data,
                    method = "oscorespls",
                    ncomp = 50,
                    validation = "CV")

PLSR_mod_PVC3 = plsr(MASS_mg ~ spcAmovav,
                    data = PVC_data,
                    method = "oscorespls",
                    ncomp = 50,
                    validation = "CV")

PLSR_mod_PVC4 = plsr(MASS_mg ~ spcA,
                    data = PVC_raw,
                    method = "oscorespls",
                    ncomp = 50,
                    validation = "CV")

which.min(RMSEP(PLSR_mod_PVC)$val["CV", , ][-1]) # 36
which.min(RMSEP(PLSR_mod_PVC2)$val["CV", , ][-1]) # 15
which.min(RMSEP(PLSR_mod_PVC3)$val["CV", , ][-1]) # 10
which.min(RMSEP(PLSR_mod_PVC4)$val["CV", , ][-1]) # 4

## polymer-specific CV
cv_PLSR_PVC = cv_plsr_model(PVC_data, PVC_data$spcARmovav, ncomp = 3)
cv_PLSR_PVC2 = cv_plsr_model(PVC_data, PVC_data$spcA, ncomp = 2)
cv_PLSR_PVC3 = cv_plsr_model(PVC_data, PVC_data$spcAmovav, ncomp = 8)
cv_PLSR_PVC4 = cv_plsr_model(PVC_raw, PVC_raw$spcA, ncomp = 2)

cv_plsr_PVCresults = list(cv_PLSR_PVC, cv_PLSR_PVC2, cv_PLSR_PVC3, cv_PLSR_PVC4)

for (i in seq_along(cv_plsr_PVCresults)) {
    cat(paste0("\n======= PLSR MODEL ", i, " (10-fold CV) =======\n"))
    cat(sprintf("ME   : %.4f\n", cv_plsr_PVCresults[[i]]$ME))
    cat(sprintf("RMSE : %.4f\n", cv_plsr_PVCresults[[i]]$RMSE))
    cat(sprintf("R²   : %.4f\n", cv_plsr_PVCresults[[i]]$R2))
}


## polymer-specific external validation
PLSR_predPVC = predict(PLSR_mod_PVC, ncomp = 36, newdata = new$spcARmovav)
PLSR_predPVC2 = predict(PLSR_mod_PVC2, ncomp = 15, newdata = new$spcA)
PLSR_predPVC3 = predict(PLSR_mod_PVC3, ncomp = 10, newdata = new$spcAmovav)
PLSR_predPVC4 = predict(PLSR_mod_PVC4, ncomp = 4, newdata = raw_new$spcA)

PVC_obs = c(rep(list(new$PVC_mg), 3), list(raw_new$PVC_mg))
PVC_preds = list(PLSR_predPVC, PLSR_predPVC2, PLSR_predPVC3, PLSR_predPVC4)

PVC_stats = mapply(function(obs, pred) {
    c(
        ME   = ME(obs, pred),
        RMSE = RMSE(obs, pred),
        R2   = R2(obs, pred)
    )
}, PVC_obs, PVC_preds)

for (i in seq_along(PVC_preds)) {
    cat(paste0("\n======= MODEL ", i, " =======\n"))
    cat(sprintf("ME   : %.7f\n", PVC_stats["ME", i]))
    cat(sprintf("RMSE : %.7f\n", PVC_stats["RMSE", i]))
    cat(sprintf("R²   : %.7f\n", PVC_stats["R2", i]))
}

#==============================================================================#

PLSR_mod_PET = plsr(MASS_mg ~ spcARmovav,
                    data = PET_data,
                    method = "oscorespls",
                    ncomp = 50,
                    validation = "CV")

PLSR_mod_PET2 = plsr(MASS_mg ~ spcA,
                     data = PET_data,
                     method = "oscorespls",
                     ncomp = 50,
                     validation = "CV")

PLSR_mod_PET3 = plsr(MASS_mg ~ spcAmovav,
                     data = PET_data,
                     method = "oscorespls",
                     ncomp = 50,
                     validation = "CV")

PLSR_mod_PET4 = plsr(MASS_mg ~ spcA,
                     data = PET_raw,
                     method = "oscorespls",
                     ncomp = 50,
                     validation = "CV")

which.min(RMSEP(PLSR_mod_PET)$val["CV", , ][-1]) # 4
which.min(RMSEP(PLSR_mod_PET2)$val["CV", , ][-1]) # 7
which.min(RMSEP(PLSR_mod_PET3)$val["CV", , ][-1]) # 6
which.min(RMSEP(PLSR_mod_PET4)$val["CV", , ][-1]) # 8

## polymer-specific CV
cv_PLSR_PET = cv_plsr_model(PET_data, PET_data$spcARmovav, ncomp = 3)
cv_PLSR_PET2 = cv_plsr_model(PET_data, PET_data$spcA, ncomp = 2)
cv_PLSR_PET3 = cv_plsr_model(PET_data, PET_data$spcAmovav, ncomp = 8)
cv_PLSR_PET4 = cv_plsr_model(PET_raw, PET_raw$spcA, ncomp = 2)

cv_plsr_PETresults = list(cv_PLSR_PET, cv_PLSR_PET2, cv_PLSR_PET3, cv_PLSR_PET4)

for (i in seq_along(cv_plsr_PETresults)) {
    cat(paste0("\n======= PLSR MODEL ", i, " (10-fold CV) =======\n"))
    cat(sprintf("ME   : %.4f\n", cv_plsr_PETresults[[i]]$ME))
    cat(sprintf("RMSE : %.4f\n", cv_plsr_PETresults[[i]]$RMSE))
    cat(sprintf("R²   : %.4f\n", cv_plsr_PETresults[[i]]$R2))
}


## polymer-specific external validation
PLSR_predPET = predict(PLSR_mod_PET, ncomp = 4, newdata = new$spcARmovav)
PLSR_predPET2 = predict(PLSR_mod_PET2, ncomp = 7, newdata = new$spcA)
PLSR_predPET3 = predict(PLSR_mod_PET3, ncomp = 6, newdata = new$spcAmovav)
PLSR_predPET4 = predict(PLSR_mod_PET4, ncomp = 8, newdata = raw_new$spcA)

PET_obs = c(rep(list(new$PET_mg), 3), list(raw_new$PET_mg))
PET_preds = list(PLSR_predPET, PLSR_predPET2, PLSR_predPET3, PLSR_predPET4)

PET_stats = mapply(function(obs, pred) {
    c(
        ME   = ME(obs, pred),
        RMSE = RMSE(obs, pred),
        R2   = R2(obs, pred)
    )
}, PET_obs, PET_preds)

for (i in seq_along(PET_preds)) {
    cat(paste0("\n======= MODEL ", i, " =======\n"))
    cat(sprintf("ME   : %.7f\n", PET_stats["ME", i]))
    cat(sprintf("RMSE : %.7f\n", PET_stats["RMSE", i]))
    cat(sprintf("R²   : %.7f\n", PET_stats["R2", i]))
}

#==============================================================================#

PLSR_mod_PE = plsr(MASS_mg ~ spcARmovav,
                    data = PE_data,
                    method = "oscorespls",
                    ncomp = 50,
                    validation = "CV")

PLSR_mod_PE2 = plsr(MASS_mg ~ spcA,
                     data = PE_data,
                     method = "oscorespls",
                     ncomp = 50,
                     validation = "CV")

PLSR_mod_PE3 = plsr(MASS_mg ~ spcAmovav,
                     data = PE_data,
                     method = "oscorespls",
                     ncomp = 50,
                     validation = "CV")

PLSR_mod_PE4 = plsr(MASS_mg ~ spcA,
                     data = PE_raw,
                     method = "oscorespls",
                     ncomp = 50,
                     validation = "CV")

which.min(RMSEP(PLSR_mod_PE)$val["CV", , ][-1]) # 5
which.min(RMSEP(PLSR_mod_PE2)$val["CV", , ][-1]) # 1
which.min(RMSEP(PLSR_mod_PE3)$val["CV", , ][-1]) # 1
which.min(RMSEP(PLSR_mod_PE4)$val["CV", , ][-1]) # 1

## polymer-specific CV
cv_PLSR_PE = cv_plsr_model(PE_data, PE_data$spcARmovav, ncomp = 3)
cv_PLSR_PE2 = cv_plsr_model(PE_data, PE_data$spcA, ncomp = 2)
cv_PLSR_PE3 = cv_plsr_model(PE_data, PE_data$spcAmovav, ncomp = 8)
cv_PLSR_PE4 = cv_plsr_model(PE_raw, PE_raw$spcA, ncomp = 2)

cv_plsr_PEresults = list(cv_PLSR_PE, cv_PLSR_PE2, cv_PLSR_PE3, cv_PLSR_PE4)

for (i in seq_along(cv_plsr_PEresults)) {
    cat(paste0("\n======= PLSR MODEL ", i, " (10-fold CV) =======\n"))
    cat(sprintf("ME   : %.4f\n", cv_plsr_PEresults[[i]]$ME))
    cat(sprintf("RMSE : %.4f\n", cv_plsr_PEresults[[i]]$RMSE))
    cat(sprintf("R²   : %.4f\n", cv_plsr_PEresults[[i]]$R2))
}

## polymer-specific external validation
PLSR_predPE = predict(PLSR_mod_PE, ncomp = 5, newdata = new$spcARmovav)
PLSR_predPE2 = predict(PLSR_mod_PE2, ncomp = 1, newdata = new$spcA)
PLSR_predPE3 = predict(PLSR_mod_PE3, ncomp = 1, newdata = new$spcAmovav)
PLSR_predPE4 = predict(PLSR_mod_PE4, ncomp = 1, newdata = raw_new$spcA)

PE_obs = c(rep(list(new$PE_mg), 3), list(raw_new$PE_mg))
PE_preds = list(PLSR_predPE, PLSR_predPE2, PLSR_predPE3, PLSR_predPE4)

PE_stats = mapply(function(obs, pred) {
    c(
        ME   = ME(obs, pred),
        RMSE = RMSE(obs, pred),
        R2   = R2(obs, pred)
    )
}, PE_obs, PE_preds)

for (i in seq_along(PE_preds)) {
    cat(paste0("\n======= MODEL ", i, " =======\n"))
    cat(sprintf("ME   : %.7f\n", PE_stats["ME", i]))
    cat(sprintf("RMSE : %.7f\n", PE_stats["RMSE", i]))
    cat(sprintf("R²   : %.7f\n", PE_stats["R2", i]))
}


################################################################################
#                    NIR - SWIR WAVELENGTH RANGE ONLY                          #
################################################################################

VIS = which(as.numeric(colnames(data$spcA)) < 1000)
VIS2 = which(as.numeric(colnames(data$spcARmovav)) < 1000)


data$NIRspcARmovav = data$spcARmovav[, -VIS2]
data$NIRspcA = data$spcA[, -VIS]
data$NIRspcAmovav = data$spcAmovav[, -VIS]
raw$NIRspcA = raw$spcA[, -VIS]

new$NIRspcARmovav = new$spcARmovav[, -VIS2]
new$NIRspcA = new$spcA[, -VIS]
new$NIRspcAmovav = new$spcAmovav[, -VIS]
raw_new$NIRspcA = raw_new$spcA[,-VIS]


set.seed(42)

PLSR_mod_NIR = plsr(MASS_mg ~ NIRspcARmovav,
                     data = data,
                     method = "oscorespls",
                     ncomp = 50,
                     validation = "CV")

PLSR_mod_NIR2 = plsr(MASS_mg ~ NIRspcA,
                      data = data,
                      method = "oscorespls",
                      ncomp = 50,
                      validation = "CV")

PLSR_mod_NIR3 = plsr(MASS_mg ~ NIRspcAmovav,
                      data = data,
                      method = "oscorespls",
                      ncomp = 50,
                      validation = "CV")

PLSR_mod_NIR4 = plsr(MASS_mg ~ NIRspcA,
                      data = raw,
                      method = "oscorespls",
                      ncomp = 50,
                      validation = "CV")

par(mfrow = c(2, 2))
validationplot(PLSR_mod_NIR, val.type = "RMSEP", main = "MODEL 1")
validationplot(PLSR_mod_NIR2, val.type = "RMSEP", main = "MODEL 2")
validationplot(PLSR_mod_NIR3, val.type = "RMSEP", main = "MODEL 3")
validationplot(PLSR_mod_NIR4, val.type = "RMSEP", main = "MODEL 4")

min(RMSEP(PLSR_mod_NIR)$val["CV", ,][-1])
min(RMSEP(PLSR_mod_NIR2)$val["CV", ,][-1])
min(RMSEP(PLSR_mod_NIR3)$val["CV", ,][-1])
min(RMSEP(PLSR_mod_NIR4)$val["CV", ,][-1])

which.min(RMSEP(PLSR_mod_NIR)$val["CV", , ][-1]) # 19
which.min(RMSEP(PLSR_mod_NIR2)$val["CV", , ][-1]) # 15
which.min(RMSEP(PLSR_mod_NIR3)$val["CV", , ][-1]) # 18
which.min(RMSEP(PLSR_mod_NIR4)$val["CV", , ][-1]) # 15


# 10-fold cross-validation
cv_PLSR_NIR = cv_plsr_model(data, data$NIRspcARmovav, ncomp = 19)
cv_PLSR_NIR2 = cv_plsr_model(data, data$NIRspcA, ncomp = 15)
cv_PLSR_NIR3 = cv_plsr_model(data, data$NIRspcAmovav, ncomp = 18)
cv_PLSR_NIR4 = cv_plsr_model(raw, raw$NIRspcA, ncomp = 15)

cv_plsr_NIRresults = list(cv_PLSR_NIR, cv_PLSR_NIR2,
                          cv_PLSR_NIR3, cv_PLSR_NIR4)

for (i in seq_along(cv_plsr_NIRresults)) {
    cat(paste0("\n======= PLSR NIR_MODEL ", i, " (10-fold CV) =======\n"))
    cat(sprintf("ME   : %.4f\n", cv_plsr_NIRresults[[i]]$ME))
    cat(sprintf("RMSE : %.4f\n", cv_plsr_NIRresults[[i]]$RMSE))
    cat(sprintf("R²   : %.4f\n", cv_plsr_NIRresults[[i]]$R2))
}


# predicting on new dataset
PLSR_predNIR = predict(PLSR_mod_NIR, ncomp = 19, newdata = new$NIRspcARmovav)
PLSR_predNIR2 = predict(PLSR_mod_NIR2, ncomp = 15, newdata = new$NIRspcA)
PLSR_predNIR3 = predict(PLSR_mod_NIR3, ncomp = 18, newdata = new$NIRspcAmovav)
PLSR_predNIR4 = predict(PLSR_mod_NIR4, ncomp = 15, newdata = raw_new$NIRspcA)

NIR_obs = c(rep(list(new$MASS_mg), 3), list(raw_new$MASS_mg))
NIR_preds = list(PLSR_predNIR, PLSR_predNIR2, PLSR_predNIR3, PLSR_predNIR4)

NIR_stats = mapply(function(obs, pred) {
    c(
        ME   = ME(obs, pred),
        RMSE = RMSE(obs, pred),
        R2   = R2(obs, pred)
    )
}, NIR_obs, NIR_preds)

for (i in seq_along(NIR_preds)) {
    cat(paste0("\n======= MODEL ", i, " =======\n"))
    cat(sprintf("ME   : %.7f\n", NIR_stats["ME", i]))
    cat(sprintf("RMSE : %.7f\n", NIR_stats["RMSE", i]))
    cat(sprintf("R²   : %.7f\n", NIR_stats["R2", i]))
}

#==============================================================================#
#              MODELING FOR SINGLE POLYMERS QUANTIFICATION (NIR)               #
#==============================================================================#

PLSR_NIR_PP = plsr(MASS_mg ~ NIRspcARmovav,
                   data = PP_data,
                   method = "oscorespls",
                   ncomp = 50,
                   validation = "CV")

PLSR_NIR_PP2 = plsr(MASS_mg ~ NIRspcA,
                    data = PP_data,
                    method = "oscorespls",
                    ncomp = 50,
                    validation = "CV")

PLSR_NIR_PP3 = plsr(MASS_mg ~ NIRspcAmovav,
                    data = PP_data,
                    method = "oscorespls",
                    ncomp = 50,
                    validation = "CV")

PLSR_NIR_PP4 = plsr(MASS_mg ~ NIRspcA,
                    data = PP_raw,
                    method = "oscorespls",
                    ncomp = 50,
                    validation = "CV")

which.min(RMSEP(PLSR_NIR_PP)$val["CV", , ][-1]) # 1
which.min(RMSEP(PLSR_NIR_PP2)$val["CV", , ][-1]) # 10
which.min(RMSEP(PLSR_NIR_PP3)$val["CV", , ][-1]) # 1
which.min(RMSEP(PLSR_NIR_PP4)$val["CV", , ][-1]) # 7

## polymer-specific CV
cv_PLSR_PP_NIR = cv_plsr_model(PP_data, PP_data$NIRspcARmovav, ncomp = 1)
cv_PLSR_PP_NIR2 = cv_plsr_model(PP_data, PP_data$NIRspcA, ncomp = 10)
cv_PLSR_PP_NIR3 = cv_plsr_model(PP_data, PP_data$NIRspcAmovav, ncomp = 1)
cv_PLSR_PP_NIR4 = cv_plsr_model(PP_raw, PP_raw$NIRspcA, ncomp = 7)

cv_plsr_PP_NIRresults = list(cv_PLSR_PP_NIR, cv_PLSR_PP_NIR2,
                             cv_PLSR_PP_NIR3, cv_PLSR_PP_NIR4)

for (i in seq_along(cv_plsr_PP_NIRresults)) {
    cat(paste0("\n======= PLSR MODEL ", i, " (10-fold CV) =======\n"))
    cat(sprintf("ME   : %.4f\n", cv_plsr_PP_NIRresults[[i]]$ME))
    cat(sprintf("RMSE : %.4f\n", cv_plsr_PP_NIRresults[[i]]$RMSE))
    cat(sprintf("R²   : %.4f\n", cv_plsr_PP_NIRresults[[i]]$R2))
}


## external validation
PLSR_predPP_NIR = predict(PLSR_NIR_PP, ncomp = 1, newdata = new$NIRspcARmovav)
PLSR_predPP_NIR2 = predict(PLSR_NIR_PP2, ncomp = 10, newdata = new$NIRspcA)
PLSR_predPP_NIR3 = predict(PLSR_NIR_PP3, ncomp = 1, newdata = new$NIRspcAmovav)
PLSR_predPP_NIR4 = predict(PLSR_NIR_PP4, ncomp = 7, newdata = raw_new$NIRspcA)

NIR_PP_obs = c(rep(list(new$PP_mg), 3), list(raw_new$PP_mg))
NIR_PP_preds = list(PLSR_predPP_NIR, PLSR_predPP_NIR2,
                    PLSR_predPP_NIR3, PLSR_predPP_NIR4)

NIR_PP_stats = mapply(function(obs, pred) {
    c(
        ME   = ME(obs, pred),
        RMSE = RMSE(obs, pred),
        R2   = R2(obs, pred)
    )
}, NIR_PP_obs, NIR_PP_preds)

for (i in seq_along(NIR_PP_preds)) {
    cat(paste0("\n======= NIR_MODEL ", i, " =======\n"))
    cat(sprintf("ME   : %.7f\n", NIR_PP_stats["ME", i]))
    cat(sprintf("RMSE : %.7f\n", NIR_PP_stats["RMSE", i]))
    cat(sprintf("R²   : %.7f\n", NIR_PP_stats["R2", i]))
}

#==============================================================================#

PLSR_NIR_PVC = plsr(MASS_mg ~ NIRspcARmovav,
                    data = PVC_data,
                    method = "oscorespls",
                    ncomp = 50,
                    validation = "CV")

PLSR_NIR_PVC2 = plsr(MASS_mg ~ NIRspcA,
                     data = PVC_data,
                     method = "oscorespls",
                     ncomp = 50,
                     validation = "CV")

PLSR_NIR_PVC3 = plsr(MASS_mg ~ NIRspcAmovav,
                     data = PVC_data,
                     method = "oscorespls",
                     ncomp = 50,
                     validation = "CV")

PLSR_NIR_PVC4 = plsr(MASS_mg ~ NIRspcA,
                     data = PVC_raw,
                     method = "oscorespls",
                     ncomp = 50,
                     validation = "CV")

which.min(RMSEP(PLSR_NIR_PVC)$val["CV", , ][-1]) # 3
which.min(RMSEP(PLSR_NIR_PVC2)$val["CV", , ][-1]) # 14
which.min(RMSEP(PLSR_NIR_PVC3)$val["CV", , ][-1]) # 7
which.min(RMSEP(PLSR_NIR_PVC4)$val["CV", , ][-1]) # 16

## polymer-specific CV
cv_PLSR_PVC_NIR = cv_plsr_model(PVC_data, PVC_data$NIRspcARmovav, ncomp = 3)
cv_PLSR_PVC_NIR2 = cv_plsr_model(PVC_data, PVC_data$NIRspcA, ncomp = 14)
cv_PLSR_PVC_NIR3 = cv_plsr_model(PVC_data, PVC_data$NIRspcAmovav, ncomp = 7)
cv_PLSR_PVC_NIR4 = cv_plsr_model(PVC_raw, PVC_raw$NIRspcA, ncomp = 16)

cv_plsr_PVC_NIRresults = list(cv_PLSR_PVC_NIR, cv_PLSR_PVC_NIR2,
                              cv_PLSR_PVC_NIR3, cv_PLSR_PVC_NIR4)

for (i in seq_along(cv_plsr_PVC_NIRresults)) {
    cat(paste0("\n======= PLSR MODEL ", i, " (10-fold CV) =======\n"))
    cat(sprintf("ME   : %.4f\n", cv_plsr_PVC_NIRresults[[i]]$ME))
    cat(sprintf("RMSE : %.4f\n", cv_plsr_PVC_NIRresults[[i]]$RMSE))
    cat(sprintf("R²   : %.4f\n", cv_plsr_PVC_NIRresults[[i]]$R2))
}


PLSR_predPVC_NIR = predict(PLSR_NIR_PVC, ncomp = 3, newdata = new$NIRspcARmovav)
PLSR_predPVC_NIR2 = predict(PLSR_NIR_PVC2, ncomp = 14, newdata = new$NIRspcA)
PLSR_predPVC_NIR3 = predict(PLSR_NIR_PVC3, ncomp = 7, newdata = new$NIRspcAmovav)
PLSR_predPVC_NIR4 = predict(PLSR_NIR_PVC4, ncomp = 16, newdata = raw_new$NIRspcA)

NIR_PVC_obs = c(rep(list(new$PVC_mg), 3), list(raw_new$PVC_mg))
NIR_PVC_preds = list(PLSR_predPVC_NIR, PLSR_predPVC_NIR2,
                     PLSR_predPVC_NIR3, PLSR_predPVC_NIR4)

NIR_PVC_stats = mapply(function(obs, pred) {
    c(
        ME   = ME(obs, pred),
        RMSE = RMSE(obs, pred),
        R2   = R2(obs, pred)
    )
}, NIR_PVC_obs, NIR_PVC_preds)

for (i in seq_along(NIR_PVC_preds)) {
    cat(paste0("\n======= NIR_MODEL ", i, " =======\n"))
    cat(sprintf("ME   : %.7f\n", NIR_PVC_stats["ME", i]))
    cat(sprintf("RMSE : %.7f\n", NIR_PVC_stats["RMSE", i]))
    cat(sprintf("R²   : %.7f\n", NIR_PVC_stats["R2", i]))
}

#==============================================================================#

PLSR_NIR_PET = plsr(MASS_mg ~ NIRspcARmovav,
                    data = PET_data,
                    method = "oscorespls",
                    ncomp = 50,
                    validation = "CV")

PLSR_NIR_PET2 = plsr(MASS_mg ~ NIRspcA,
                     data = PET_data,
                     method = "oscorespls",
                     ncomp = 50,
                     validation = "CV")

PLSR_NIR_PET3 = plsr(MASS_mg ~ NIRspcAmovav,
                     data = PET_data,
                     method = "oscorespls",
                     ncomp = 50,
                     validation = "CV")

PLSR_NIR_PET4 = plsr(MASS_mg ~ NIRspcA,
                     data = PET_raw,
                     method = "oscorespls",
                     ncomp = 50,
                     validation = "CV")

which.min(RMSEP(PLSR_NIR_PET)$val["CV", , ][-1]) # 8
which.min(RMSEP(PLSR_NIR_PET2)$val["CV", , ][-1]) # 9
which.min(RMSEP(PLSR_NIR_PET3)$val["CV", , ][-1]) # 6
which.min(RMSEP(PLSR_NIR_PET4)$val["CV", , ][-1]) # 5

## polymer-specific CV
cv_PLSR_PET_NIR = cv_plsr_model(PET_data, PET_data$NIRspcARmovav, ncomp = 8)
cv_PLSR_PET_NIR2 = cv_plsr_model(PET_data, PET_data$NIRspcA, ncomp = 9)
cv_PLSR_PET_NIR3 = cv_plsr_model(PET_data, PET_data$NIRspcAmovav, ncomp = 6)
cv_PLSR_PET_NIR4 = cv_plsr_model(PET_raw, PET_raw$NIRspcA, ncomp = 5)

cv_plsr_PET_NIRresults = list(cv_PLSR_PET_NIR, cv_PLSR_PET_NIR2,
                              cv_PLSR_PET_NIR3, cv_PLSR_PET_NIR4)

for (i in seq_along(cv_plsr_PET_NIRresults)) {
    cat(paste0("\n======= PLSR MODEL ", i, " (10-fold CV) =======\n"))
    cat(sprintf("ME   : %.4f\n", cv_plsr_PET_NIRresults[[i]]$ME))
    cat(sprintf("RMSE : %.4f\n", cv_plsr_PET_NIRresults[[i]]$RMSE))
    cat(sprintf("R²   : %.4f\n", cv_plsr_PET_NIRresults[[i]]$R2))
}


PLSR_predPET_NIR = predict(PLSR_NIR_PET, ncomp = 8, newdata = new$NIRspcARmovav)
PLSR_predPET_NIR2 = predict(PLSR_NIR_PET2, ncomp = 9, newdata = new$NIRspcA)
PLSR_predPET_NIR3 = predict(PLSR_NIR_PET3, ncomp = 6, newdata = new$NIRspcAmovav)
PLSR_predPET_NIR4 = predict(PLSR_NIR_PET4, ncomp = 5, newdata = raw_new$NIRspcA)

NIR_PET_obs = c(rep(list(new$PET_mg), 3), list(raw_new$PET_mg))
NIR_PET_preds = list(PLSR_predPET_NIR, PLSR_predPET_NIR2,
                     PLSR_predPET_NIR3, PLSR_predPET_NIR4)

NIR_PET_stats = mapply(function(obs, pred) {
    c(
        ME   = ME(obs, pred),
        RMSE = RMSE(obs, pred),
        R2   = R2(obs, pred)
    )
}, NIR_PET_obs, NIR_PET_preds)

for (i in seq_along(NIR_PET_preds)) {
    cat(paste0("\n======= NIR_MODEL ", i, " =======\n"))
    cat(sprintf("ME   : %.7f\n", NIR_PET_stats["ME", i]))
    cat(sprintf("RMSE : %.7f\n", NIR_PET_stats["RMSE", i]))
    cat(sprintf("R²   : %.7f\n", NIR_PET_stats["R2", i]))
}

#==============================================================================#

PLSR_NIR_PE = plsr(MASS_mg ~ NIRspcARmovav,
                   data = PE_data,
                   method = "oscorespls",
                   ncomp = 50,
                   validation = "CV")

PLSR_NIR_PE2 = plsr(MASS_mg ~ NIRspcA,
                    data = PE_data,
                    method = "oscorespls",
                    ncomp = 50,
                    validation = "CV")

PLSR_NIR_PE3 = plsr(MASS_mg ~ NIRspcAmovav,
                    data = PE_data,
                    method = "oscorespls",
                    ncomp = 50,
                    validation = "CV")

PLSR_NIR_PE4 = plsr(MASS_mg ~ NIRspcA,
                    data = PE_raw,
                    method = "oscorespls",
                    ncomp = 50,
                    validation = "CV")

which.min(RMSEP(PLSR_NIR_PE)$val["CV", , ][-1]) # 6
which.min(RMSEP(PLSR_NIR_PE2)$val["CV", , ][-1]) # 5
which.min(RMSEP(PLSR_NIR_PE3)$val["CV", , ][-1]) # 6
which.min(RMSEP(PLSR_NIR_PE4)$val["CV", , ][-1]) # 2

## polymer-specific CV
cv_PLSR_PE_NIR = cv_plsr_model(PE_data, PE_data$NIRspcARmovav, ncomp = 6)
cv_PLSR_PE_NIR2 = cv_plsr_model(PE_data, PE_data$NIRspcA, ncomp = 5)
cv_PLSR_PE_NIR3 = cv_plsr_model(PE_data, PE_data$NIRspcAmovav, ncomp = 6)
cv_PLSR_PE_NIR4 = cv_plsr_model(PE_raw, PE_raw$NIRspcA, ncomp = 2)

cv_plsr_PE_NIRresults = list(cv_PLSR_PE_NIR, cv_PLSR_PE_NIR2,
                             cv_PLSR_PE_NIR3, cv_PLSR_PE_NIR4)

for (i in seq_along(cv_plsr_PE_NIRresults)) {
    cat(paste0("\n======= PLSR MODEL ", i, " (10-fold CV) =======\n"))
    cat(sprintf("ME   : %.4f\n", cv_plsr_PE_NIRresults[[i]]$ME))
    cat(sprintf("RMSE : %.4f\n", cv_plsr_PE_NIRresults[[i]]$RMSE))
    cat(sprintf("R²   : %.4f\n", cv_plsr_PE_NIRresults[[i]]$R2))
}


PLSR_predPE_NIR = predict(PLSR_NIR_PE, ncomp = 6, newdata = new$NIRspcARmovav)
PLSR_predPE_NIR2 = predict(PLSR_NIR_PE2, ncomp = 5, newdata = new$NIRspcA)
PLSR_predPE_NIR3 = predict(PLSR_NIR_PE3, ncomp = 6, newdata = new$NIRspcAmovav)
PLSR_predPE_NIR4 = predict(PLSR_NIR_PE4, ncomp = 2, newdata = raw_new$NIRspcA)

NIR_PE_obs = c(rep(list(new$PE_mg), 3), list(raw_new$PE_mg))
NIR_PE_preds = list(PLSR_predPE_NIR, PLSR_predPE_NIR2,
                    PLSR_predPE_NIR3, PLSR_predPE_NIR4)

NIR_PE_stats = mapply(function(obs, pred) {
    c(
        ME   = ME(obs, pred),
        RMSE = RMSE(obs, pred),
        R2   = R2(obs, pred)
    )
}, NIR_PE_obs, NIR_PE_preds)

for (i in seq_along(NIR_PE_preds)) {
    cat(paste0("\n======= NIR_MODEL ", i, " =======\n"))
    cat(sprintf("ME   : %.7f\n", NIR_PE_stats["ME", i]))
    cat(sprintf("RMSE : %.7f\n", NIR_PE_stats["RMSE", i]))
    cat(sprintf("R²   : %.7f\n", NIR_PE_stats["R2", i]))
}

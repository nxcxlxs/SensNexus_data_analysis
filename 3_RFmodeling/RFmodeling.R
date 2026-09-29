require(prospectr)
require(dplyr)
require(randomForest)
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


## visualize differences (e.g., pristine raw vs. full preprocessing)============#
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
                              pristine$MASS_mg,
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


# fitting RF
## prepare calibration data
rawCsub = data.frame(plastMass = rawC$MASS_mg, rawC$spcA)
colnames(rawCsub) = c("plastMass", paste0("spec.", colnames(rawC$spcA)))

datCsub2 = data.frame(plastMass = datC$MASS_mg, datC$spcA)
colnames(datCsub2) = c("plastMass", paste0("spec.", colnames(datC$spcA)))

datCsub3 = data.frame(plastMass = datC$MASS_mg, datC$spcAmovav)
colnames(datCsub3) = c("plastMass", paste0("spec.", colnames(datC$spcAmovav)))

datCsub4 = data.frame(plastMass = datC$MASS_mg, datC$spcARmovav)
colnames(datCsub4) = c("plastMass", paste0("spec.", colnames(datC$spcARmovav)))


# same tests as in PLSR...
set.seed(1)

RF_mod_mass = randomForest(plastMass ~ .,
                           data = rawCsub,
                           ntree = 150,
                           mtry = 10,
                           importance = T,
                           na.action = na.omit)

RF_mod_mass2 = randomForest(plastMass ~ .,
                           data = datCsub2,
                           ntree = 150,
                           mtry = 10,
                           importance = T,
                           na.action = na.omit)

RF_mod_mass3 = randomForest(plastMass ~ .,
                           data = datCsub3,
                           ntree = 150,
                           mtry = 10,
                           importance = T,
                           na.action = na.omit)

RF_mod_mass4 = randomForest(plastMass ~ .,
                            data = datCsub4,
                            ntree = 150,
                            mtry = 10,
                            importance = T,
                            na.action = na.omit)

## variables importance
### plots
impPlot = varImpPlot(RF_mod_mass, main = "Model 1")
impPlot2 = varImpPlot(RF_mod_mass2, main = "Model 2")
impPlot3 = varImpPlot(RF_mod_mass3, main = "Model 3")
impPlot4 = varImpPlot(RF_mod_mass4, main = "Model 4")

### values inspection
head(impPlot[order(impPlot[,"%IncMSE"], decreasing = TRUE), ], 30)
head(impPlot[order(impPlot[,"IncNodePurity"], decreasing = TRUE), ], 30)


# prepare validation data
rawVsub = data.frame(plastMass = rawV$MASS_mg, rawV$spcA)
colnames(rawVsub) = c("plastMass", paste0("spec.", colnames(rawV$spcA)))

datVsub2 = data.frame(plastMass = datV$MASS_mg, datV$spcA)
colnames(datVsub2) = c("plastMass", paste0("spec.", colnames(datV$spcA)))


datVsub3 = data.frame(plastMass = datV$MASS_mg, datV$spcAmovav)
colnames(datVsub3) = c("plastMass", paste0("spec.", colnames(datV$spcAmovav)))

datVsub4 = data.frame(plastMass = datV$MASS_mg, datV$spcARmovav)
colnames(datVsub4) = c("plastMass", paste0("spec.", colnames(datV$spcARmovav)))


# predictions
## calibration
RFpredC = predict(RF_mod_mass, rawCsub)
RFpredC2 = predict(RF_mod_mass2, datCsub2)
RFpredC3 = predict(RF_mod_mass3, datCsub3)
RFpredC4 = predict(RF_mod_mass4, datCsub4)

## validation
RFpredV = predict(RF_mod_mass, rawVsub)
RFpredV2 = predict(RF_mod_mass2, datVsub2)
RFpredV3 = predict(RF_mod_mass3, datVsub3)
RFpredV4 = predict(RF_mod_mass4, datVsub4)


# ## plots
# par(mfrow = c(1, 2))
# 
# plot(log(datC$MASS_mg), log(RFpredC),
#      main = "Calibration",
#      xlab = "log(Observed)",
#      ylab = "log(Predicted)",
#      ylim = c(0, 10),
#      xlim = c(0, 10),
#      pch = 16) +
#   abline(0, 1)
# 
# plot(log(datV$MASS_mg), log(RFpredV),
#      main = "Validation",
#      xlab = "log(Observed)",
#      ylab = "log(Predicted)",
#      ylim = c(0, 10),
#      pch = 16) +
#   abline(0, 1)


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
calib_obs = list(rawCsub$plastMass, datCsub2$plastMass,
                      datCsub3$plastMass, datCsub4$plastMass)
valid_obs = list(rawVsub$plastMass, datVsub2$plastMass,
                      datVsub3$plastMass, datVsub4$plastMass)

## group calibration predictions
calib_preds = list(RFpredC, RFpredC2, RFpredC3, RFpredC4)
## group validation predictions
valid_preds = list(RFpredV, RFpredV2, RFpredV3, RFpredV4)

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
  cat(sprintf("ME   : %.4f\n", calib_stats["ME", i]))
  cat(sprintf("RMSE : %.4f\n", calib_stats["RMSE", i]))
  cat(sprintf("R²   : %.4f\n", calib_stats["R2", i]))
  
  cat("##### VALIDATION #####\n")
  cat(sprintf("ME   : %.4f\n", valid_stats["ME", i]))
  cat(sprintf("RMSE : %.4f\n", valid_stats["RMSE", i]))
  cat(sprintf("R²   : %.4f\n", valid_stats["R2", i]))
} 


################################################################################
#                            K-FOLD CROSS-VALIDATION                           #
################################################################################
cv_random_forest = function(data, spc_matrix,  ntree = 150,
                            mtry = 10, nfolds = 4, seed = 1) {
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
    
    train_rf = data.frame(plastMass = train$MASS_mg, spc_matrix[data$foldCV != i, ])
    colnames(train_rf) = c("plastMass", paste0("spec.", colnames(spc_matrix)))
    
    valid_rf = spc_matrix[data$foldCV == i, ]
    colnames(valid_rf) = paste0("spec.", colnames(spc_matrix))
    
    rf_mod = randomForest(plastMass ~ ., data = train_rf, ntree = ntree, mtry = mtry)
    val_pred = predict(rf_mod, newdata = valid_rf)
    
    valTable$pred[valTable$SAMPLE_ID %in% valid$SAMPLE_ID] = val_pred
    
    cat("Fold number", i, "done for RF.\n")
  }
  
  list(
    ntree= ntree,
    mtry = mtry,
    ME = ME(valTable$obs, valTable$pred),
    RMSE = RMSE(valTable$obs, valTable$pred),
    R2 = R2(valTable$obs, valTable$pred),
    table = valTable
  )
}

cv_RFmod = cv_random_forest(pristine_raw, pristine_raw$spcA)
cv_RFmod2 = cv_random_forest(pristine, pristine$spcA)
cv_RFmod3 = cv_random_forest(pristine, pristine$spcAmovav)
cv_RFmod4 = cv_random_forest(pristine, pristine$spcARmovav)

cv_results = list(cv_RFmod, cv_RFmod2,
                  cv_RFmod3, cv_RFmod4)

for (i in seq_along(cv_results)) {
  cat(paste0("\n======= MODEL ", i, " (4-fold CV) =======\n"))
  cat(sprintf("ME   : %.4f\n", cv_results[[i]]$ME))
  cat(sprintf("RMSE : %.4f\n", cv_results[[i]]$RMSE))
  cat(sprintf("R²   : %.4f\n", cv_results[[i]]$R2))
}


################################################################################
#            FULL INTERNAL DATASET TRAINING/EXTERNAL VALIDATION                #
################################################################################

# prepare data
## prepare FULL calibration data
rawFULL = data.frame(plastMass = pristine_raw$MASS_mg, pristine_raw$spcA)
colnames(rawFULL) = c("plastMass", paste0("spec.", colnames(pristine_raw$spcA)))

dataFULL2 = data.frame(plastMass = pristine$MASS_mg, pristine$spcA)
colnames(dataFULL2) = c("plastMass", paste0("spec.", colnames(pristine$spcA)))

dataFULL3 = data.frame(plastMass = pristine$MASS_mg, pristine$spcAmovav)
colnames(dataFULL3) = c("plastMass", paste0("spec.", colnames(pristine$spcAmovav)))

dataFULL4 = data.frame(plastMass = pristine$MASS_mg, pristine$spcARmovav)
colnames(dataFULL4) = c("plastMass", paste0("spec.", colnames(pristine$spcARmovav)))


# "naive" models
set.seed(2)

RF_mod_full = randomForest(plastMass ~ .,
                           data = rawFULL,
                           ntree = 150,
                           mtry = 10,
                           importance = T,
                           na.action = na.omit)

RF_mod_full2 = randomForest(plastMass ~ .,
                            data = dataFULL2,
                            ntree = 150,
                            mtry = 10,
                            importance = T,
                            na.action = na.omit)

RF_mod_full3 = randomForest(plastMass ~ .,
                            data = dataFULL3,
                            ntree = 150,
                            mtry = 10,
                            importance = T,
                            na.action = na.omit)

RF_mod_full4 = randomForest(plastMass ~ .,
                            data = dataFULL4,
                            ntree = 150,
                            mtry = 10,
                            importance = T,
                            na.action = na.omit)


varImpPlot(RF_mod_full, main = "Model 1")
varImpPlot(RF_mod_full2, main = "Model 2")
varImpPlot(RF_mod_full3, main = "Model 3")
varImpPlot(RF_mod_full4, main = "Model 4")


# prepare new data (external validation)
new_sub = data.frame(plastMass = colored_raw$MASS_mg, colored_raw$spcA)
colnames(new_sub) = c("plastMass", paste0("spec.", colnames(colored_raw$spcA)))

new_sub2 = data.frame(plastMass = colored$MASS_mg, colored$spcA)
colnames(new_sub2) = c("plastMass", paste0("spec.", colnames(colored$spcA)))

new_sub3 = data.frame(plastMass = colored$MASS_mg, colored$spcAmovav)
colnames(new_sub3) = c("plastMass", paste0("spec.", colnames(colored$spcAmovav)))

new_sub4 = data.frame(plastMass = colored$MASS_mg, colored$spcARmovav)
colnames(new_sub4) = c("plastMass", paste0("spec.", colnames(colored$spcARmovav)))


## predicting on new dataset 
RFpred_full = predict(RF_mod_full, new_sub)
RFpred_full2 = predict(RF_mod_full2, new_sub2)
RFpred_full3 = predict(RF_mod_full3, new_sub3)
RFpred_full4 = predict(RF_mod_full4, new_sub4)

external_obs = list(new_sub$plastMass, new_sub2$plastMass,
                    new_sub3$plastMass, new_sub4$plastMass)
external_preds = list(RFpred_full, RFpred_full2, RFpred_full3, RFpred_full4)

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
set.seed(23)

colored$strata = interaction(colored$MASS_mg,
                             colored$SIZE_CODE)

colored1 = colored |> 
  group_by(strata) |> 
  sample_frac(0.5)

colored2 = colored |> 
  filter(!(SAMPLE_ID %in% colored1$SAMPLE_ID))


colored_raw$strata = interaction(colored_raw$MASS_mg,
                                 colored_raw$SIZE_CODE)

colored_raw1 = colored_raw |> 
  group_by(strata) |> 
  sample_frac(0.5)

colored_raw2 = colored_raw |> 
  filter(!(SAMPLE_ID %in% colored_raw1$SAMPLE_ID))


# prepare external data
raw_new1_M1 = data.frame(plastMass = colored_raw1$MASS_mg, colored_raw1$spcA)
colnames(raw_new1_M1) = c("plastMass", paste0("spec.", colnames(colored_raw1$spcA)))
raw_new2_M1 = data.frame(plastMass = colored_raw2$MASS_mg, colored_raw2$spcA)
colnames(raw_new2_M1) = c("plastMass", paste0("spec.", colnames(colored_raw2$spcA)))

new1_M2 = data.frame(plastMass = colored1$MASS_mg, colored1$spcA)
colnames(new1_M2) = c("plastMass", paste0("spec.", colnames(colored1$spcA)))
new2_M2 = data.frame(plastMass = colored2$MASS_mg, colored2$spcA)
colnames(new2_M2) = c("plastMass", paste0("spec.", colnames(colored2$spcA)))

new1_M3 = data.frame(plastMass = colored1$MASS_mg, colored1$spcAmovav)
colnames(new1_M3) = c("plastMass", paste0("spec.", colnames(colored1$spcAmovav)))
new2_M3 = data.frame(plastMass = colored2$MASS_mg, colored2$spcAmovav)
colnames(new2_M3) = c("plastMass", paste0("spec.", colnames(colored2$spcAmovav)))

new1_M4 = data.frame(plastMass = colored1$MASS_mg, colored1$spcARmovav)
colnames(new1_M4) = c("plastMass", paste0("spec.", colnames(colored1$spcARmovav)))
new2_M4 = data.frame(plastMass = colored2$MASS_mg, colored2$spcARmovav)
colnames(new2_M4) = c("plastMass", paste0("spec.", colnames(colored2$spcARmovav)))

param_grid_ext = expand.grid(
  ntree = c(100, 150, 200, 250, 500),
  mtry = c(5, 10, 20, 30, 40, 50, 55, 60, 70, 80, 90, 95, 100)
)

spectral_preproc_list_ext = list(
  M1 = list(train = rawFULL, new = raw_new1_M1),
  M2 = list(train = dataFULL2, new = new1_M2),
  M3 = list(train = dataFULL3, new = new1_M3),
  M4 = list(train = dataFULL4, new = new1_M4)
  )


# hyperparameter optimization test ============================================#
t0 = Sys.time()

ext_results = list()

counter = 1

for (model_name in names(spectral_preproc_list_ext)) {
  
  cat("\n======= External tuning for", model_name, "=======\n")
  
  train_data = spectral_preproc_list_ext[[model_name]]$train
  new_data   = spectral_preproc_list_ext[[model_name]]$new
  
  for (p in seq_len(nrow(param_grid_ext))) {
    cat("Combination", p, "of", nrow(param_grid_ext), "fitted\n")
    
    ntree_val = param_grid_ext$ntree[p]
    mtry_val  = param_grid_ext$mtry[p]
    
    rf_model = randomForest(
      plastMass ~ .,
      data  = train_data,
      ntree = ntree_val,
      mtry  = mtry_val
    )
    
    pred_ext = predict(rf_model, new_data)
    
    ext_results[[counter]] = data.frame(
      Model = model_name,
      ntree = ntree_val,
      mtry  = mtry_val,
      ME    = ME(new_data$plastMass, pred_ext),
      RMSE  = RMSE(new_data$plastMass, pred_ext),
      R2    = R2(new_data$plastMass, pred_ext)
    )
    
    counter = counter + 1
  }
}

ext_results = do.call(rbind, ext_results)

t1 = Sys.time()
cat("Training time:", round(t1 - t0, 2), "minutes\n")


## check optima parameters
ext_results |> 
  group_by(Model) |> 
  filter(RMSE == min(RMSE))

ext_results |>
  group_by(Model) |>
  filter(abs(ME) == min(abs(ME)))
#==============================================================================#


## tuned models
RF_mod_tuned = randomForest(plastMass ~ .,
                            data = rawFULL,
                            ntree = 200,
                            mtry = 70,
                            importance = T,
                            na.action = na.omit)

RF_mod_tuned2 = randomForest(plastMass ~ .,
                             data = dataFULL2,
                             ntree = 250,
                             mtry = 50,
                             importance = T,
                             na.action = na.omit)

RF_mod_tuned3 = randomForest(plastMass ~ .,
                             data = dataFULL3,
                             ntree = 150,
                             mtry = 80,
                             importance = T,
                             na.action = na.omit)

RF_mod_tuned4 = randomForest(plastMass ~ .,
                             data = dataFULL4,
                             ntree = 100,
                             mtry = 50,
                             importance = T,
                             na.action = na.omit)

RFpred_tuned = predict(RF_mod_tuned, raw_new2_M1)
RFpred_tuned2 = predict(RF_mod_tuned2, new2_M2)
RFpred_tuned3 = predict(RF_mod_tuned3, new2_M3)
RFpred_tuned4 = predict(RF_mod_tuned4, new2_M4)

tuned_obs = list(raw_new2_M1$plastMass, new2_M2$plastMass,
                 new2_M3$plastMass, new2_M4$plastMass)
tuned_preds = list(RFpred_tuned, RFpred_tuned2, RFpred_tuned3, RFpred_tuned4)

tuned_stats = mapply(function(obs, pred) {
  c(
    ME   = ME(obs, pred),
    RMSE = RMSE(obs, pred),
    R2   = R2(obs, pred)
  )
}, tuned_obs, tuned_preds)

for (i in seq_along(tuned_preds)) {
  cat(paste0("\n======= MODEL ", i, " =======\n"))
  cat(sprintf("ME   : %.7f\n", tuned_stats["ME", i]))
  cat(sprintf("RMSE : %.7f\n", tuned_stats["RMSE", i]))
  cat(sprintf("R²   : %.7f\n", tuned_stats["R2", i]))
}


# residual visualization
residTUNED = RFpred_tuned - colored_raw2$MASS_mg
residTUNED2 = RFpred_tuned2 - colored2$MASS_mg
residTUNED3 = RFpred_tuned3 - colored2$MASS_mg
residTUNED4 = RFpred_tuned4 - colored2$MASS_mg

residTUNEDlog = log(RFpred_tuned) - log(colored_raw2$MASS_mg)
residTUNEDlog2 = log(RFpred_tuned2) - log(colored2$MASS_mg)
residTUNEDlog3 = log(RFpred_tuned3) - log(colored2$MASS_mg)
residTUNEDlog4 = log(RFpred_tuned4) - log(colored2$MASS_mg)



df_m1 = data.frame(
  Model = "M1 (Raw data)",
  Observed_Mass = colored_raw2$MASS_mg,
  Residual = residTUNED
)

df_m2 = data.frame(
  Model = "M2 (Minimal preprocessing)",
  Observed_Mass = colored2$MASS_mg,
  Residual = residTUNED2
)

df_m3 = data.frame(
  Model = "M3 (Intermediate preprocessing)",
  Observed_Mass = colored2$MASS_mg,
  Residual = residTUNED3
)

df_m4 = data.frame(
  Model = "M4 (Full preprocessing)",
  Observed_Mass = colored2$MASS_mg,
  Residual = residTUNED4
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
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 11),
        axis.title.y = element_text(size = 13))
#==============================================================================#


# data frames for relative residuals (log scale)
df_m1_log = data.frame(
  Model = "M1 (Raw data)",
  Observed_Mass = colored_raw2$MASS_mg,
  Relative_Residual = residTUNEDlog
  )

df_m2_log = data.frame(
  Model = "M2 (Minimal preprocessing)", 
  Observed_Mass = colored2$MASS_mg,
  Relative_Residual = residTUNEDlog2
  )

df_m3_log = data.frame(
  Model = "M3 (Intermediate preprocessing)",
  Observed_Mass = colored2$MASS_mg, 
  Relative_Residual = residTUNEDlog3
  )

df_m4_log = data.frame(
  Model = "M4 (Full preprocessing)",
  Observed_Mass = colored2$MASS_mg,
  Relative_Residual = residTUNEDlog4
  )

# combine and process
plot_data_log = bind_rows(df_m1_log, df_m2_log, df_m3_log, df_m4_log)


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

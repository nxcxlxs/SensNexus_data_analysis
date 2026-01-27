# load processed data
data = readRDS("C:/nico/Dissertação/SENSNEXUS_data/analysis_ready_data/deNoised_data.rds")
raw = readRDS("C:/nico/Dissertação/SENSNEXUS_data/preprocessed_data/datsoil.rds")
new = readRDS("C:/nico/Dissertação/SENSNEXUS_data/analysis_ready_data/deNoised_newData.rds") 
raw_new = readRDS("C:/nico/Dissertação/SENSNEXUS_data/preprocessed_data/datsoil2.rds")

# convert spectra to absorbance
data$spcA = log(1/data$spc)
raw$spcA = log(1/raw$spc)

new$spcA = log(1/new$spc)
raw_new$spcA = log(1/raw_new$spc)

# plot spectra profile
matplot(colnames(raw_new$spcA), t(raw_new$spcA),
        type = "l",
        lty = 1,
        col = rgb(0.5, 0.5, 0.5, alpha = 0.3))

# smooth it out (for Wadoux's approach)
oldWavs = as.numeric(colnames(data$spcA))
newWavs = seq(min(oldWavs), max(oldWavs), by = 5)

data$spcAR = prospectr::resample(data$spcA,
                                 wav = oldWavs,
                                 new.wav = newWavs,
                                 interpol = "linear")

new$spcAR = prospectr::resample(new$spcA,
                                wav = oldWavs,
                                new.wav = newWavs,
                                interpol = "linear")

require(prospectr)
## SNV for baseline correction
data$spcARsnv = standardNormalVariate(data$spcAR)
data$spcAsnv = standardNormalVariate(data$spcA)

new$spcARsnv = standardNormalVariate(new$spcAR)
new$spcAsnv = standardNormalVariate(new$spcA)

## Moving Window Average to the SNV spectra
data$spcARmovav = movav(data$spcARsnv, w = 11)
data$spcAmovav = movav(data$spcAsnv, w = 11)

new$spcARmovav = movav(new$spcARsnv, w = 11)
new$spcAmovav = movav(new$spcAsnv, w = 11)

matplot(colnames(data$spcAmovav), t(data$spcAmovav),
        type = "l",
        lty = 1,
        col = rgb(0.5, 0.5, 0.5,
                  alpha = 0.3))

#==============================================================================#
# TRY `data$spcA`, `data$spcAmovav` and `data$spcARmovav` TO CHECK DIFFERENCES #
#==============================================================================#

# combined stratification key for calib./valid. data split
require(dplyr)
data$strata = interaction(data$POLYMER,
                          data$MASS_mg, # this is pretty much a `treatment` flaggin...
                          data$SIZE_CODE)

raw$strata = interaction(raw$POLYMER,
                         raw$MASS_mg,
                         raw$SIZE_CODE)

set.seed(1)
datC = data |> 
  group_by(strata) |> 
  sample_frac(0.75)

datV = data |> 
  filter(!(SAMPLE_ID %in% datC$SAMPLE_ID))

rawC = raw |> 
  group_by(strata) |> 
  sample_frac(0.75)

rawV = raw |> 
  filter(!(SAMPLE_ID %in% rawC$SAMPLE_ID))

vars = c("MASS_mg", "POLYMER", "SIZE_CODE")
par(mfrow = c(3, 2))

for (i in seq_along(vars)) {
  var = vars[i]
  
  barplot(table(datC[[var]]),
          main = "Calibration",
          xlab = var)
  
  barplot(table(datV[[var]]),
          main = "Validation",
          xlab = var)
  cat("\nCALIBRATION:\n", var)
  print(table(datC[[var]]))
  
  cat("\nVALIDATION:\n", var)
  print(table(datV[[var]]))
}

################################################################################

# validation metrics
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

################################################################################

# fitting RF
require(randomForest)

# prepare calibration data
datCsub = data.frame(plastMass = datC$MASS_mg, datC$spcARmovav)
colnames(datCsub) = c("plastMass", paste0("spec.", colnames(datC$spcARmovav)))

datCsub2 = data.frame(plastMass = datC$MASS_mg, datC$spcA)
colnames(datCsub2) = c("plastMass", paste0("spec.", colnames(datC$spcA)))


datCsub3 = data.frame(plastMass = datC$MASS_mg, datC$spcAmovav)
colnames(datCsub3) = c("plastMass", paste0("spec.", colnames(datC$spcAmovav)))

rawCsub = data.frame(plastMass = rawC$MASS_mg, rawC$spcA)
colnames(rawCsub) = c("plastMass", paste0("spec.", colnames(rawC$spcA)))

# same tests as in PLSR...
set.seed(1)

RF_mod_mass = randomForest(plastMass ~ .,
                           data = datCsub,
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
                            data = rawCsub,
                            ntree = 150,
                            mtry = 10,
                            importance = T,
                            na.action = na.omit)

# DOESN'T LOOK SO GOOD

# rownames(RF_mod_mass$importance) = paste0(gsub("spec\\.", "",
#                                                rownames(RF_mod_mass$importance)), " nm")
# rownames(RF_mod_mass2$importance) = paste0(gsub("spec\\.", "",
#                                                 rownames(RF_mod_mass2$importance)), " nm")
# rownames(RF_mod_mass3$importance) = paste0(gsub("spec\\.", "",
#                                                 rownames(RF_mod_mass3$importance)), " nm")
# rownames(RF_mod_mass4$importance) = paste0(gsub("spec\\.", "",
#                                                 rownames(RF_mod_mass4$importance)), " nm")


rownames(RF_mod_mass$importance) = paste0(gsub("nm", "",
                                               rownames(RF_mod_mass$importance)), "")

rownames(RF_mod_mass2$importance) = paste0(gsub("nm", "",
                                               rownames(RF_mod_mass2$importance)), "")

rownames(RF_mod_mass3$importance) = paste0(gsub("nm", "",
                                               rownames(RF_mod_mass3$importance)), "")

rownames(RF_mod_mass4$importance) = paste0(gsub("nm", "",
                                               rownames(RF_mod_mass4$importance)), "")

par(mfrow = c(2, 2), mar = c(4, 6, 5, 2))
par(cex.main = 1.8,      # title size
    cex.lab = 2,     # axis label size
    cex.axis = 1.5,
    mgp = c(2.7, 1, 0))

varImpPlot(RF_mod_mass4,
           type = 1,
           pch = 1,
           main = "Model 1")
title(ylab = "Wavelenghts (nm)", line = 3.5)


varImpPlot(RF_mod_mass2,
           type = 1,
           main = "Model 2")

varImpPlot(RF_mod_mass3,
           type = 1,
           main = "Model 3")
title(ylab = "Wavelenghts (nm)", line = 3.5)

varImpPlot(RF_mod_mass,
           type = 1,
           main = "Model 4")


varImpPlot(RF_mod_mass4, type = 2, main = "Model 1")
varImpPlot(RF_mod_mass2, type = 2, main = "Model 2")
varImpPlot(RF_mod_mass3, type = 2, main = "Model 3")
varImpPlot(RF_mod_mass,  type = 2, main = "Model 4")

impPlot4 = varImpPlot(RF_mod_mass, main = "Model 4")
impPlot2 = varImpPlot(RF_mod_mass2, main = "Model 2")
impPlot3 = varImpPlot(RF_mod_mass3, main = "Model 3")
impPlot1 = varImpPlot(RF_mod_mass4, main = "Model 1")

head(impPlot1[order(impPlot1[,"%IncMSE"], decreasing = TRUE), ], 30)
head(impPlot1[order(impPlot1[,"IncNodePurity"], decreasing = TRUE), ], 30)

# prepare validation data
datVsub = data.frame(plastMass = datV$MASS_mg, datV$spcARmovav)
colnames(datVsub) = c("plastMass", paste0("spec.", colnames(datV$spcARmovav)))

datVsub2 = data.frame(plastMass = datV$MASS_mg, datV$spcA)
colnames(datVsub2) = c("plastMass", paste0("spec.", colnames(datV$spcA)))


datVsub3 = data.frame(plastMass = datV$MASS_mg, datV$spcAmovav)
colnames(datVsub3) = c("plastMass", paste0("spec.", colnames(datV$spcAmovav)))

rawVsub = data.frame(plastMass = rawV$MASS_mg, rawV$spcA)
colnames(rawVsub) = c("plastMass", paste0("spec.", colnames(rawV$spcA)))

# prepare new data (external validation)
new_sub = data.frame(plastMass = new$MASS_mg, new$spcARmovav)
colnames(new_sub) = c("plastMass", paste0("spec.", colnames(new$spcARmovav)))

new_sub2 = data.frame(plastMass = new$MASS_mg, new$spcA)
colnames(new_sub2) = c("plastMass", paste0("spec.", colnames(new$spcA)))

new_sub3 = data.frame(plastMass = new$MASS_mg, new$spcAmovav)
colnames(new_sub3) = c("plastMass", paste0("spec.", colnames(new$spcAmovav)))

new_sub4 = data.frame(plastMass = raw_new$MASS_mg, raw_new$spcA)
colnames(new_sub4) = c("plastMass", paste0("spec.", colnames(raw_new$spcA)))

# predictions
## calibration
RFpredC = predict(RF_mod_mass, datCsub)
RFpredC2 = predict(RF_mod_mass2, datCsub2)
RFpredC3 = predict(RF_mod_mass3, datCsub3)
RFpredC4 = predict(RF_mod_mass4, rawCsub)

## validation
RFpredV = predict(RF_mod_mass, datVsub)
RFpredV2 = predict(RF_mod_mass2, datVsub2)
RFpredV3 = predict(RF_mod_mass3, datVsub3)
RFpredV4 = predict(RF_mod_mass4, rawVsub)

## external
RFpred_new = predict(RF_mod_mass, new_sub)
RFpred_new2 = predict(RF_mod_mass2, new_sub2)
RFpred_new3 = predict(RF_mod_mass3, new_sub3)
RFpred_new4 = predict(RF_mod_mass4, new_sub4)


### plot calibration
par(mfrow = c(1, 2))

plot(log(datC$MASS_mg), log(RFpredC),
     main = "Calibration",
     xlab = "log(Observed)",
     ylab = "log(Predicted)",
     ylim = c(0, 10),
     xlim = c(0, 10),
     pch = 16) +
  abline(0, 1)

### plot validation
plot(log(datV$MASS_mg), log(RFpredV),
     main = "Validation",
     xlab = "log(Observed)",
     ylab = "log(Predicted)",
     ylim = c(0, 10),
     pch = 16) +
  abline(0, 1)

# evaluate quality of predictions
## observed responses
calib_obs = list(datCsub$plastMass, datCsub2$plastMass,
                      datCsub3$plastMass, rawCsub$plastMass)
valid_obs = list(datVsub$plastMass, datVsub2$plastMass,
                      datVsub3$plastMass, rawVsub$plastMass)

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

#### EXTERNAL VALIDATION (colored rigid household plastics dataset)
external_obs = list(new_sub$plastMass, new_sub2$plastMass,
                    new_sub3$plastMass, new_sub4$plastMass)
external_preds = list(RFpred_new, RFpred_new2, RFpred_new3, RFpred_new4)

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
# ME(rawC$MASS_mg, RFpredC4)
# RMSE(rawC$MASS_mg, RFpredC4)
# R2(rawC$MASS_mg, RFpredC4)
# residC = RFpredC4 - rawC$MASS_mg # absolute residual = real error
# residClog = log(RFpredC4) - log(rawC$MASS_mg) # log residual = relative error (emphasizes all scales equally)
# boxplot(residClog ~ rawC$MASS_mg,
#         main = "Residuals by Mass (Calibration)") # absolute residuals overemphasize
#                                                   # large errors at higher response values
# 
# ## validation
# ME(rawV$MASS_mg, RFpredV4)
# RMSE(rawV$MASS_mg, RFpredV4)
# R2(rawV$MASS_mg, RFpredV4)
# residV = RFpredV4 - rawV$MASS_mg
# residVlog = log(RFpredV4) - log(rawV$MASS_mg)
# boxplot(residVlog ~ rawV$MASS_mg,
#         main = "Residuals by Mass (Validation)")


# residuals boxplots
## calibration
residC = RFpredC - datC$MASS_mg
residC2 = RFpredC2 - datC$MASS_mg
residC3 = RFpredC3 - datC$MASS_mg
residC4 = RFpredC4 - rawC$MASS_mg

residClog = log(RFpredC) - log(datC$MASS_mg)
residClog2 = log(RFpredC2) - log(datC$MASS_mg)
residClog3 = log(RFpredC3) - log(datC$MASS_mg)
residClog4 = log(RFpredC4) - log(rawC$MASS_mg)

par(mfrow = c(1, 2))
boxplot(residC ~ datC$MASS_mg,
        main = "Residuals by Mass (Calibration)")

boxplot(residClog ~ datC$MASS_mg, # ATTENTION to the log transformation
        main = "Relative residuals by Mass (Calibration)")

## validation
residV = RFpredV - datV$MASS_mg
residV2 = RFpredV2 - datV$MASS_mg
residV3 = RFpredV3 - datV$MASS_mg
residV4 = RFpredV4 - rawV$MASS_mg

residVlog = log(RFpredV) - log(datV$MASS_mg)
residVlog2 = log(RFpredV2) - log(datV$MASS_mg)
residVlog3 = log(RFpredV3) - log(datV$MASS_mg)
residVlog4 = log(RFpredV4) - log(rawV$MASS_mg)

boxplot(residV ~ datV$MASS_mg,
        main = "Residuals by Mass (Validation)")

boxplot(residVlog ~ datV$MASS_mg,
        main = "Relative residuals by Mass (Validation)")

################################################################################

# 10-fold cross-validation
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

cv_RFmod = cv_random_forest(data, data$spcARmovav)
cv_RFmod2 = cv_random_forest(data, data$spcA)
cv_RFmod3 = cv_random_forest(data, data$spcAmovav)
cv_RFmod4 = cv_random_forest(raw, raw$spcA)

cv_results = list(cv_RFmod, cv_RFmod2,
                  cv_RFmod3, cv_RFmod4)

for (i in seq_along(cv_results)) {
  cat(paste0("\n======= MODEL ", i, " (4-fold CV) =======\n"))
  cat(sprintf("ME   : %.4f\n", cv_results[[i]]$ME))
  cat(sprintf("RMSE : %.4f\n", cv_results[[i]]$RMSE))
  cat(sprintf("R²   : %.4f\n", cv_results[[i]]$R2))
}

################################################################################

## test parameterization
### spectral preprocessing datasets
spectral_preproc_list = list(
  M1 = data$spcARmovav,  # SGf + 5nm resample + SNV + moving average
  M2 = data$spcA,        # SGf only
  M3 = data$spcAmovav,  # SGf + SNV + moving average 
  M4 = raw$spcA          # raw data
)

param_grid = expand.grid(
  ntree = c(100, 150, 200),
  mtry = c(10, 50, 55, 100)
)

cv_results = list()

### loop over RF params:
for (p in seq_len(nrow(param_grid))) {
  ntree_val = param_grid$ntree[p]
  mtry_val = param_grid$mtry[p]
  
  cat("\n\n=== RF parameters: ntree =", ntree_val, ", mtry =", mtry_val, "===\n")
  
  #### loop over spectral preprocess models
  for (model_name in names(spectral_preproc_list)) {
    spc_mat = spectral_preproc_list[[model_name]]
    
    res = cv_random_forest(
      data = data,
      spc_matrix = spc_mat,
      ntree = ntree_val,
      mtry = mtry_val
    )
    
    ##### store results with descriptive key
    key = paste0("ntree", ntree_val, "_mtry", mtry_val, "_", model_name)
    cv_results[[key]] = res
    
    ##### display results
    cat("======= MODEL", model_name, "(10-fold CV) =======\n")
    cat(sprintf("ME   : %.4f\n", res$ME))
    cat(sprintf("RMSE : %.4f\n", res$RMSE))
    cat(sprintf("R²   : %.4f\n\n", res$R2))
  }
}

### optionally convert results to a summary dataframe
summary_df = do.call(rbind, lapply(names(cv_results), function(k) {
  res = cv_results[[k]]
  data.frame(
    Model = k,
    ntree = res$ntree,
    mtry = res$mtry,
    ME = res$ME,
    RMSE = res$RMSE,
    R2 = res$R2
  )
}))

print(summary_df)


#==============================================================================#
#                     FULL INTERNAL DATASET TRANING!                           #
#==============================================================================#

## prepare data
# prepare FULL calibration data
dataFULL = data.frame(plastMass = data$MASS_mg, data$spcARmovav)
colnames(dataFULL) = c("plastMass", paste0("spec.", colnames(data$spcARmovav)))

dataFULL2 = data.frame(plastMass = data$MASS_mg, data$spcA)
colnames(dataFULL2) = c("plastMass", paste0("spec.", colnames(data$spcA)))

dataFULL3 = data.frame(plastMass = data$MASS_mg, data$spcAmovav)
colnames(dataFULL3) = c("plastMass", paste0("spec.", colnames(data$spcAmovav)))

rawFULL = data.frame(plastMass = raw$MASS_mg, raw$spcA)
colnames(rawFULL) = c("plastMass", paste0("spec.", colnames(raw$spcA)))


## "naive" models
set.seed(2)

RF_mod_full = randomForest(plastMass ~ .,
                           data = dataFULL,
                           ntree = 100,
                           mtry = 10,
                           importance = T,
                           na.action = na.omit)

RF_mod_full2 = randomForest(plastMass ~ .,
                            data = dataFULL2,
                            ntree = 100,
                            mtry = 20,
                            importance = T,
                            na.action = na.omit)

RF_mod_full3 = randomForest(plastMass ~ .,
                            data = dataFULL3,
                            ntree = 100,
                            mtry = 10,
                            importance = T,
                            na.action = na.omit)

RF_mod_full4 = randomForest(plastMass ~ .,
                            data = rawFULL,
                            ntree = 100,
                            mtry = 20,
                            importance = T,
                            na.action = na.omit)

varImpPlot(RF_mod_full, main = "Model 1")
varImpPlot(RF_mod_full2, main = "Model 2")
varImpPlot(RF_mod_full3, main = "Model 3")
varImpPlot(RF_mod_full4, main = "Model 4")


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

## "tuned-on-external" models
set.seed(3)

RF_mod_leak = randomForest(plastMass ~ .,
                           data = dataFULL,
                           ntree = 100,
                           mtry = 50,
                           importance = T,
                           na.action = na.omit)

RF_mod_leak2 = randomForest(plastMass ~ .,
                            data = dataFULL2,
                            ntree = 100,
                            mtry = 16,
                            importance = T,
                            na.action = na.omit)

RF_mod_leak3 = randomForest(plastMass ~ .,
                            data = dataFULL3,
                            ntree = 100,
                            mtry = 50,
                            importance = T,
                            na.action = na.omit)

RF_mod_leak4 = randomForest(plastMass ~ .,
                            data = rawFULL,
                            ntree = 100,
                            mtry = 95, # check higher!!!
                            importance = T,
                            na.action = na.omit)

RFpred_leak = predict(RF_mod_leak, new_sub)
RFpred_leak2 = predict(RF_mod_leak2, new_sub2)
RFpred_leak3 = predict(RF_mod_leak3, new_sub3)
RFpred_leak4 = predict(RF_mod_leak4, new_sub4)

leak_obs = list(new_sub$plastMass, new_sub2$plastMass,
                    new_sub3$plastMass, new_sub4$plastMass)
leak_preds = list(RFpred_leak, RFpred_leak2, RFpred_leak3, RFpred_leak4)

leak_stats = mapply(function(obs, pred) {
  c(
    ME   = ME(obs, pred),
    RMSE = RMSE(obs, pred),
    R2   = R2(obs, pred)
  )
}, leak_obs, leak_preds)

for (i in seq_along(leak_preds)) {
  cat(paste0("\n======= MODEL ", i, " =======\n"))
  cat(sprintf("ME   : %.7f\n", leak_stats["ME", i]))
  cat(sprintf("RMSE : %.7f\n", leak_stats["RMSE", i]))
  cat(sprintf("R²   : %.7f\n", leak_stats["R2", i]))
}


# residual visualization
residLEAK = RFpred_leak - new$MASS_mg
residLEAK2 = RFpred_leak2 - new$MASS_mg
residLEAK3 = RFpred_leak3 - new$MASS_mg
residLEAK4 = RFpred_leak4 - raw_new$MASS_mg

residLEAKlog = log(RFpred_leak) - log(new$MASS_mg)
residLEAKlog2 = log(RFpred_leak2) - log(new$MASS_mg)
residLEAKlog3 = log(RFpred_leak3) - log(new$MASS_mg)
residLEAKlog4 = log(RFpred_leak4) - log(raw_new$MASS_mg)

ggplot(raw_new, aes(x = factor(MASS_mg), y = residLEAK)) +
    geom_boxplot() +
    labs(title = "MODEL 4",
         x = "Mass (mg)",
         y = "Residuals")

ggplot(raw_new, aes(x = factor(MASS_mg), y = residLEAKlog)) +
    geom_boxplot() +
    labs(title = "MODEL 4",
         x = "Mass (mg)",
         y = "Relative residuals")

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
require(dplyr)
plot_data = bind_rows(df_m1, df_m2, df_m3, df_m4)

# model factor in correct order of complexity
model_order = c("M1 (Raw data)", "M2 (Minimal preprocessing)",
                "M3 (Intermediate preprocessing)", "M4 (Full preprocessing)")
plot_data$Model = factor(plot_data$Model, levels = model_order)

# combined boxplot
require(ggplot2)
res = ggplot(plot_data, aes(x = factor(Observed_Mass), y = Residual)) +
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
    title = NULL)

require(patchwork)

res / res2 + plot_annotation(tag_levels = 'a',
                             tag_prefix = '(',
                             tag_suffix = ')')

################################################################################
#                   MODELLING FOR SINGLE POLYMERS QUANTIFICATION               #
################################################################################

# set polymer-specific subsets
polymers = c("PP", "PVC", "PET", "PE")

for(p in polymers) {
  assign(paste0(p, "_data"), filter(data, POLYMER == p))
  assign(paste0(p, "_raw"), filter(raw, POLYMER == p))
}

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

t0 = Sys.time()

## PP modeling
### prepare RF-ready datasets
PP_data1 = data.frame(plastMass = PP_data$MASS_mg, PP_data$spcARmovav)
colnames(PP_data1) = c("plastMass", paste0("spec.", colnames(PP_data$spcARmovav)))

PP_data2 = data.frame(plastMass = PP_data$MASS_mg, PP_data$spcA)
colnames(PP_data2) = c("plastMass", paste0("spec.", colnames(PP_data$spcA)))

PP_data3 = data.frame(plastMass = PP_data$MASS_mg, PP_data$spcAmovav)
colnames(PP_data3) = c("plastMass", paste0("spec.", colnames(PP_data$spcAmovav)))

PP_dataR = data.frame(plastMass = PP_raw$MASS_mg, PP_raw$spcA)
colnames(PP_dataR) = c("plastMass", paste0("spec.", colnames(raw$spcA)))


set.seed(4)

RF_mod_PP = randomForest(plastMass ~ .,
                           data = PP_data1,
                           ntree = 100,
                           mtry = 10,
                           importance = T,
                           na.action = na.omit)

RF_mod_PP2 = randomForest(plastMass ~ .,
                            data = PP_data2,
                            ntree = 100,
                            mtry = 20,
                            importance = T,
                            na.action = na.omit)

RF_mod_PP3 = randomForest(plastMass ~ .,
                            data = PP_data3,
                            ntree = 100,
                            mtry = 10,
                            importance = T,
                            na.action = na.omit)

RF_mod_PP4 = randomForest(plastMass ~ .,
                            data = PP_dataR,
                            ntree = 100,
                            mtry = 20,
                            importance = T,
                            na.action = na.omit)

RFpred_PP = predict(RF_mod_PP, new_sub)
RFpred_PP2 = predict(RF_mod_PP2, new_sub2)
RFpred_PP3 = predict(RF_mod_PP3, new_sub3)
RFpred_PP4 = predict(RF_mod_PP4, new_sub4)

PP_obs = list(new_sub$plastMass, new_sub2$plastMass,
             new_sub3$plastMass, new_sub4$plastMass)
PP_preds = list(RFpred_PP, RFpred_PP2, RFpred_PP3, RFpred_PP4)

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

### polymer-specific CV
cv_RF_PP = cv_random_forest(PP_data, PP_data$spcARmovav)
cv_RF_PP2 = cv_random_forest(PP_data, PP_data$spcA)
cv_RF_PP3 = cv_random_forest(PP_data, PP_data$spcAmovav)
cv_RF_PP4 = cv_random_forest(PP_raw, PP_raw$spcA)

cv_RF_PPresults = list(cv_RF_PP, cv_RF_PP2,
                  cv_RF_PP3, cv_RF_PP4)

for (i in seq_along(cv_RF_PPresults)) {
  cat(paste0("\n======= MODEL ", i, " (10-fold CV) =======\n"))
  cat(sprintf("ME   : %.4f\n", cv_RF_PPresults[[i]]$ME))
  cat(sprintf("RMSE : %.4f\n", cv_RF_PPresults[[i]]$RMSE))
  cat(sprintf("R²   : %.4f\n", cv_RF_PPresults[[i]]$R2))
}


#==============================================================================#

## PVC modeling
### prepare RF-ready datasets
PVC_data1 = data.frame(plastMass = PVC_data$MASS_mg, PVC_data$spcARmovav)
colnames(PVC_data1) = c("plastMass", paste0("spec.", colnames(PVC_data$spcARmovav)))

PVC_data2 = data.frame(plastMass = PVC_data$MASS_mg, PVC_data$spcA)
colnames(PVC_data2) = c("plastMass", paste0("spec.", colnames(PVC_data$spcA)))

PVC_data3 = data.frame(plastMass = PVC_data$MASS_mg, PVC_data$spcAmovav)
colnames(PVC_data3) = c("plastMass", paste0("spec.", colnames(PVC_data$spcAmovav)))

PVC_dataR = data.frame(plastMass = PVC_raw$MASS_mg, PVC_raw$spcA)
colnames(PVC_dataR) = c("plastMass", paste0("spec.", colnames(raw$spcA)))


RF_mod_PVC = randomForest(plastMass ~ .,
                          data = PVC_data1,
                          ntree = 100,
                          mtry = 10,
                          importance = T,
                          na.action = na.omit)

RF_mod_PVC2 = randomForest(plastMass ~ .,
                           data = PVC_data2,
                           ntree = 100,
                           mtry = 20,
                           importance = T,
                           na.action = na.omit)

RF_mod_PVC3 = randomForest(plastMass ~ .,
                           data = PVC_data3,
                           ntree = 100,
                           mtry = 10,
                           importance = T,
                           na.action = na.omit)

RF_mod_PVC4 = randomForest(plastMass ~ .,
                           data = PVC_dataR,
                           ntree = 100,
                           mtry = 20,
                           importance = T,
                           na.action = na.omit)

RFpred_PVC = predict(RF_mod_PVC, new_sub)
RFpred_PVC2 = predict(RF_mod_PVC2, new_sub2)
RFpred_PVC3 = predict(RF_mod_PVC3, new_sub3)
RFpred_PVC4 = predict(RF_mod_PVC4, new_sub4)

PVC_obs = list(new_sub$plastMass, new_sub2$plastMass,
               new_sub3$plastMass, new_sub4$plastMass)
PVC_preds = list(RFpred_PVC, RFpred_PVC2, RFpred_PVC3, RFpred_PVC4)

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

### polymer-specific CV
cv_RF_PVC = cv_random_forest(PVC_data, PVC_data$spcARmovav)
cv_RF_PVC2 = cv_random_forest(PVC_data, PVC_data$spcA)
cv_RF_PVC3 = cv_random_forest(PVC_data, PVC_data$spcAmovav)
cv_RF_PVC4 = cv_random_forest(PVC_raw, PVC_raw$spcA)

cv_RF_PVCresults = list(cv_RF_PVC, cv_RF_PVC2,
                  cv_RF_PVC3, cv_RF_PVC4)

for (i in seq_along(cv_RF_PVCresults)) {
  cat(paste0("\n======= MODEL ", i, " (10-fold CV) =======\n"))
  cat(sprintf("ME   : %.4f\n", cv_RF_PVCresults[[i]]$ME))
  cat(sprintf("RMSE : %.4f\n", cv_RF_PVCresults[[i]]$RMSE))
  cat(sprintf("R²   : %.4f\n", cv_RF_PVCresults[[i]]$R2))
}


#==============================================================================#

## PET modeling
### prepare RF-ready datasets
PET_data1 = data.frame(plastMass = PET_data$MASS_mg, PET_data$spcARmovav)
colnames(PET_data1) = c("plastMass", paste0("spec.", colnames(PET_data$spcARmovav)))

PET_data2 = data.frame(plastMass = PET_data$MASS_mg, PET_data$spcA)
colnames(PET_data2) = c("plastMass", paste0("spec.", colnames(PET_data$spcA)))

PET_data3 = data.frame(plastMass = PET_data$MASS_mg, PET_data$spcAmovav)
colnames(PET_data3) = c("plastMass", paste0("spec.", colnames(PET_data$spcAmovav)))

PET_dataR = data.frame(plastMass = PET_raw$MASS_mg, PET_raw$spcA)
colnames(PET_dataR) = c("plastMass", paste0("spec.", colnames(raw$spcA)))


RF_mod_PET = randomForest(plastMass ~ .,
                          data = PET_data1,
                          ntree = 100,
                          mtry = 10,
                          importance = T,
                          na.action = na.omit)

RF_mod_PET2 = randomForest(plastMass ~ .,
                           data = PET_data2,
                           ntree = 100,
                           mtry = 20,
                           importance = T,
                           na.action = na.omit)

RF_mod_PET3 = randomForest(plastMass ~ .,
                           data = PET_data3,
                           ntree = 100,
                           mtry = 10,
                           importance = T,
                           na.action = na.omit)

RF_mod_PET4 = randomForest(plastMass ~ .,
                           data = PET_dataR,
                           ntree = 100,
                           mtry = 20,
                           importance = T,
                           na.action = na.omit)

RFpred_PET = predict(RF_mod_PET, new_sub)
RFpred_PET2 = predict(RF_mod_PET2, new_sub2)
RFpred_PET3 = predict(RF_mod_PET3, new_sub3)
RFpred_PET4 = predict(RF_mod_PET4, new_sub4)

PET_obs = list(new_sub$plastMass, new_sub2$plastMass,
               new_sub3$plastMass, new_sub4$plastMass)
PET_preds = list(RFpred_PET, RFpred_PET2, RFpred_PET3, RFpred_PET4)

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

### polymer-specific CV
cv_RF_PET = cv_random_forest(PET_data, PET_data$spcARmovav)
cv_RF_PET2 = cv_random_forest(PET_data, PET_data$spcA)
cv_RF_PET3 = cv_random_forest(PET_data, PET_data$spcAmovav)
cv_RF_PET4 = cv_random_forest(PET_raw, PET_raw$spcA)

cv_RF_PETresults = list(cv_RF_PET, cv_RF_PET2,
                  cv_RF_PET3, cv_RF_PET4)

for (i in seq_along(cv_RF_PETresults)) {
  cat(paste0("\n======= MODEL ", i, " (10-fold CV) =======\n"))
  cat(sprintf("ME   : %.4f\n", cv_RF_PETresults[[i]]$ME))
  cat(sprintf("RMSE : %.4f\n", cv_RF_PETresults[[i]]$RMSE))
  cat(sprintf("R²   : %.4f\n", cv_RF_PETresults[[i]]$R2))
}


#==============================================================================#

## PE modeling
### prepare RF-ready datasets
PE_data1 = data.frame(plastMass = PE_data$MASS_mg, PE_data$spcARmovav)
colnames(PE_data1) = c("plastMass", paste0("spec.", colnames(PE_data$spcARmovav)))

PE_data2 = data.frame(plastMass = PE_data$MASS_mg, PE_data$spcA)
colnames(PE_data2) = c("plastMass", paste0("spec.", colnames(PE_data$spcA)))

PE_data3 = data.frame(plastMass = PE_data$MASS_mg, PE_data$spcAmovav)
colnames(PE_data3) = c("plastMass", paste0("spec.", colnames(PE_data$spcAmovav)))

PE_dataR = data.frame(plastMass = PE_raw$MASS_mg, PE_raw$spcA)
colnames(PE_dataR) = c("plastMass", paste0("spec.", colnames(raw$spcA)))


RF_mod_PE = randomForest(plastMass ~ .,
                         data = PE_data1,
                         ntree = 100,
                         mtry = 10,
                         importance = T,
                         na.action = na.omit)

RF_mod_PE2 = randomForest(plastMass ~ .,
                          data = PE_data2,
                          ntree = 100,
                          mtry = 20,
                          importance = T,
                          na.action = na.omit)

RF_mod_PE3 = randomForest(plastMass ~ .,
                          data = PE_data3,
                          ntree = 100,
                          mtry = 10,
                          importance = T,
                          na.action = na.omit)

RF_mod_PE4 = randomForest(plastMass ~ .,
                          data = PE_dataR,
                          ntree = 100,
                          mtry = 20,
                          importance = T,
                          na.action = na.omit)

RFpred_PE = predict(RF_mod_PE, new_sub)
RFpred_PE2 = predict(RF_mod_PE2, new_sub2)
RFpred_PE3 = predict(RF_mod_PE3, new_sub3)
RFpred_PE4 = predict(RF_mod_PE4, new_sub4)

PE_obs = list(new_sub$plastMass, new_sub2$plastMass,
              new_sub3$plastMass, new_sub4$plastMass)
PE_preds = list(RFpred_PE, RFpred_PE2, RFpred_PE3, RFpred_PE4)

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

### polymer-specific CV
cv_RF_PE = cv_random_forest(PE_data, PE_data$spcARmovav)
cv_RF_PE2 = cv_random_forest(PE_data, PE_data$spcA)
cv_RF_PE3 = cv_random_forest(PE_data, PE_data$spcAmovav)
cv_RF_PE4 = cv_random_forest(PE_raw, PE_raw$spcA)

cv_RF_PEresults = list(cv_RF_PE, cv_RF_PE2,
                  cv_RF_PE3, cv_RF_PE4)

for (i in seq_along(cv_RF_PEresults)) {
  cat(paste0("\n======= MODEL ", i, " (10-fold CV) =======\n"))
  cat(sprintf("ME   : %.4f\n", cv_RF_PEresults[[i]]$ME))
  cat(sprintf("RMSE : %.4f\n", cv_RF_PEresults[[i]]$RMSE))
  cat(sprintf("R²   : %.4f\n", cv_RF_PEresults[[i]]$R2))
}


t1 = Sys.time()

print(t1 - t0)


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

dataNIR = data.frame(plastMass = data$MASS_mg, data$NIRspcARmovav)
colnames(dataNIR) = c("plastMass", paste0("spec.", colnames(data$NIRspcARmovav)))

dataNIR2 = data.frame(plastMass = data$MASS_mg, data$NIRspcA)
colnames(dataNIR2) = c("plastMass", paste0("spec.", colnames(data$NIRspcA)))

dataNIR3 = data.frame(plastMass = data$MASS_mg, data$NIRspcAmovav)
colnames(dataNIR3) = c("plastMass", paste0("spec.", colnames(data$NIRspcAmovav)))

rawNIR = data.frame(plastMass = raw$MASS_mg, raw$NIRspcA)
colnames(rawNIR) = c("plastMass", paste0("spec.", colnames(raw$NIRspcA)))


new_NIR = data.frame(plastMass = new$MASS_mg, new$NIRspcARmovav)
colnames(new_NIR) = c("plastMass", paste0("spec.", colnames(new$NIRspcARmovav)))

new_NIR2 = data.frame(plastMass = new$MASS_mg, new$NIRspcA)
colnames(new_NIR2) = c("plastMass", paste0("spec.", colnames(new$NIRspcA)))

new_NIR3 = data.frame(plastMass = new$MASS_mg, new$NIRspcAmovav)
colnames(new_NIR3) = c("plastMass", paste0("spec.", colnames(new$NIRspcAmovav)))

new_NIR4 = data.frame(plastMass = raw_new$MASS_mg, raw_new$NIRspcA)
colnames(new_NIR4) = c("plastMass", paste0("spec.", colnames(raw_new$NIRspcA)))


## test parameterization for NIR-only wavelengths
NIR_preproc_list = list(
  M1_NIR = data$NIRspcARmovav,  # SGf + 5nm resample + SNV + moving average
  M2_NIR = data$NIRspcA,        # SGf only
  M3_NIR = data$NIRspcAmovav,  # SGf + SNV + moving average 
  M4_NIR = raw$NIRspcA          # raw data
)

param_grid = expand.grid(
  ntree = c(100, 150, 200),
  mtry = c(10, 50, 55, 100)
)

NIRcv_results = list()

t0 = Sys.time()

### loop over RF params:
for (p in seq_len(nrow(param_grid))) {
  ntree_val = param_grid$ntree[p]
  mtry_val = param_grid$mtry[p]
  
  cat("\n\n=== RF parameters: ntree =", ntree_val, ", mtry =", mtry_val, "===\n")
  
  #### loop over spectral preprocess models
  for (model_name in names(NIR_preproc_list)) {
    spc_mat = NIR_preproc_list[[model_name]]
    
    res = cv_random_forest(
      data = data,
      spc_matrix = spc_mat,
      ntree = ntree_val,
      mtry = mtry_val
    )
    
    ##### store results with descriptive key
    key = paste0("ntree", ntree_val, "_mtry", mtry_val, "_", model_name)
    NIRcv_results[[key]] = res
    
    ##### display results
    cat("======= MODEL", model_name, "(10-fold CV) =======\n")
    cat(sprintf("ME   : %.4f\n", res$ME))
    cat(sprintf("RMSE : %.4f\n", res$RMSE))
    cat(sprintf("R²   : %.4f\n\n", res$R2))
  }
}

t1 = Sys.time()
cat("Training time:", round(t1 - t0, 2), "minutes\n")

### convert results to a summary dataframe
NIRsummary_df = do.call(rbind, lapply(names(NIRcv_results), function(k) {
  res = NIRcv_results[[k]]
  data.frame(
    Model = k,
    ntree = res$ntree,
    mtry = res$mtry,
    ME = res$ME,
    RMSE = res$RMSE,
    R2 = res$R2
  )
}))

print(NIRsummary_df)


## "naive" NIR-only models

set.seed(42)

RF_mod_NIR = randomForest(plastMass ~ .,
                          data = dataNIR,
                          ntree = 150,
                          mtry = 50,
                          importance = T,
                          na.action = na.omit)

RF_mod_NIR2 = randomForest(plastMass ~ .,
                           data = dataNIR2,
                           ntree = 150,
                           mtry = 10,
                           importance = T,
                           na.action = na.omit)

RF_mod_NIR3 = randomForest(plastMass ~ .,
                           data = dataNIR3,
                           ntree = 200,
                           mtry = 50,
                           importance = T,
                           na.action = na.omit)

RF_mod_NIR4 = randomForest(plastMass ~ .,
                           data = rawNIR,
                           ntree = 100,
                           mtry = 10,
                           importance = T,
                           na.action = na.omit)

varImpPlot(RF_mod_NIR, main = "NIR_Model 1")
varImpPlot(RF_mod_NIR2, main = "NIR_Model 2")
varImpPlot(RF_mod_NIR3, main = "NIR_Model 3")
varImpPlot(RF_mod_NIR4, main = "NIR_Model 4")


RFpred_NIR = predict(RF_mod_NIR, new_NIR)
RFpred_NIR2 = predict(RF_mod_NIR2, new_NIR2)
RFpred_NIR3 = predict(RF_mod_NIR3, new_NIR3)
RFpred_NIR4 = predict(RF_mod_NIR4, new_NIR4)

NIR_obs = list(new_sub$plastMass, new_sub2$plastMass,
               new_sub3$plastMass, new_sub4$plastMass)
NIR_preds = list(RFpred_NIR, RFpred_NIR2, RFpred_NIR3, RFpred_NIR4)

NIR_stats = mapply(function(obs, pred) {
  c(
    ME   = ME(obs, pred),
    RMSE = RMSE(obs, pred),
    R2   = R2(obs, pred)
  )
}, NIR_obs, NIR_preds)

for (i in seq_along(NIR_preds)) {
  cat(paste0("\n======= NIR_MODEL ", i, " =======\n"))
  cat(sprintf("ME   : %.7f\n", NIR_stats["ME", i]))
  cat(sprintf("RMSE : %.7f\n", NIR_stats["RMSE", i]))
  cat(sprintf("R²   : %.7f\n", NIR_stats["R2", i]))
}

#==============================================================================#
#                MODELLING FOR SINGLE POLYMERS QUANTIFICATION (NIR)            #
#==============================================================================#

## PP modeling
### prepare RF-ready datasets
PP_data_NIR1 = data.frame(plastMass = PP_data$MASS_mg, PP_data$NIRspcARmovav)
colnames(PP_data_NIR1) = c("plastMass", paste0("spec.", colnames(PP_data$NIRspcARmovav)))

PP_data_NIR2 = data.frame(plastMass = PP_data$MASS_mg, PP_data$NIRspcA)
colnames(PP_data_NIR2) = c("plastMass", paste0("spec.", colnames(PP_data$NIRspcA)))

PP_data_NIR3 = data.frame(plastMass = PP_data$MASS_mg, PP_data$NIRspcAmovav)
colnames(PP_data_NIR3) = c("plastMass", paste0("spec.", colnames(PP_data$NIRspcAmovav)))

PP_data_NIR_R = data.frame(plastMass = PP_raw$MASS_mg, PP_raw$NIRspcA)
colnames(PP_data_NIR_R) = c("plastMass", paste0("spec.", colnames(raw$NIRspcA)))


set.seed(4)

RF_NIR_PP = randomForest(plastMass ~ .,
                         data = PP_data_NIR1,
                         ntree = 150,
                         mtry = 50,
                         importance = T,
                         na.action = na.omit)

RF_NIR_PP2 = randomForest(plastMass ~ .,
                          data = PP_data_NIR2,
                          ntree = 150,
                          mtry = 10,
                          importance = T,
                          na.action = na.omit)

RF_NIR_PP3 = randomForest(plastMass ~ .,
                          data = PP_data_NIR3,
                          ntree = 200,
                          mtry = 50,
                          importance = T,
                          na.action = na.omit)

RF_NIR_PP4 = randomForest(plastMass ~ .,
                          data = PP_data_NIR_R,
                          ntree = 100,
                          mtry = 10,
                          importance = T,
                          na.action = na.omit)

RFpred_NIR_PP = predict(RF_NIR_PP, new_sub)
RFpred_NIR_PP2 = predict(RF_NIR_PP2, new_sub2)
RFpred_NIR_PP3 = predict(RF_NIR_PP3, new_sub3)
RFpred_NIR_PP4 = predict(RF_NIR_PP4, new_sub4)

PP_NIR_obs = list(new_sub$plastMass, new_sub2$plastMass,
                  new_sub3$plastMass, new_sub4$plastMass)
PP_NIR_preds = list(RFpred_NIR_PP, RFpred_NIR_PP2,
                    RFpred_NIR_PP3, RFpred_NIR_PP4)

PP_NIR_stats = mapply(function(obs, pred) {
  c(
    ME   = ME(obs, pred),
    RMSE = RMSE(obs, pred),
    R2   = R2(obs, pred)
  )
}, PP_NIR_obs, PP_NIR_preds)

for (i in seq_along(PP_NIR_preds)) {
  cat(paste0("\n======= MODEL ", i, " =======\n"))
  cat(sprintf("ME   : %.7f\n", PP_NIR_stats["ME", i]))
  cat(sprintf("RMSE : %.7f\n", PP_NIR_stats["RMSE", i]))
  cat(sprintf("R²   : %.7f\n", PP_NIR_stats["R2", i]))
}

### polymer-specific CV
cv_RF_PP_NIR = cv_random_forest(PP_data, PP_data$NIRspcARmovav)
cv_RF_PP_NIR2 = cv_random_forest(PP_data, PP_data$NIRspcA)
cv_RF_PP_NIR3 = cv_random_forest(PP_data, PP_data$NIRspcAmovav)
cv_RF_PP_NIR4 = cv_random_forest(PP_raw, PP_raw$NIRspcA)

cv_RF_PP_NIRresults = list(cv_RF_PP_NIR, cv_RF_PP_NIR2,
                           cv_RF_PP_NIR3, cv_RF_PP_NIR4)

for (i in seq_along(cv_RF_PP_NIRresults)) {
  cat(paste0("\n======= MODEL ", i, " (10-fold CV) =======\n"))
  cat(sprintf("ME   : %.4f\n", cv_RF_PP_NIRresults[[i]]$ME))
  cat(sprintf("RMSE : %.4f\n", cv_RF_PP_NIRresults[[i]]$RMSE))
  cat(sprintf("R²   : %.4f\n", cv_RF_PP_NIRresults[[i]]$R2))
}

## PVC modeling
### prepare RF-ready datasets
PVC_data_NIR1 = data.frame(plastMass = PVC_data$MASS_mg, PVC_data$NIRspcARmovav)
colnames(PVC_data_NIR1) = c("plastMass", paste0("spec.", colnames(PVC_data$NIRspcARmovav)))

PVC_data_NIR2 = data.frame(plastMass = PVC_data$MASS_mg, PVC_data$NIRspcA)
colnames(PVC_data_NIR2) = c("plastMass", paste0("spec.", colnames(PVC_data$NIRspcA)))

PVC_data_NIR3 = data.frame(plastMass = PVC_data$MASS_mg, PVC_data$NIRspcAmovav)
colnames(PVC_data_NIR3) = c("plastMass", paste0("spec.", colnames(PVC_data$NIRspcAmovav)))

PVC_data_NIR_R = data.frame(plastMass = PVC_raw$MASS_mg, PVC_raw$NIRspcA)
colnames(PVC_data_NIR_R) = c("plastMass", paste0("spec.", colnames(raw$NIRspcA)))


set.seed(4)

RF_NIR_PVC = randomForest(plastMass ~ .,
                          data = PVC_data_NIR1,
                          ntree = 150,
                          mtry = 50,
                          importance = T,
                          na.action = na.omit)

RF_NIR_PVC2 = randomForest(plastMass ~ .,
                           data = PVC_data_NIR2,
                           ntree = 150,
                           mtry = 10,
                           importance = T,
                           na.action = na.omit)

RF_NIR_PVC3 = randomForest(plastMass ~ .,
                           data = PVC_data_NIR3,
                           ntree = 200,
                           mtry = 50,
                           importance = T,
                           na.action = na.omit)

RF_NIR_PVC4 = randomForest(plastMass ~ .,
                           data = PVC_data_NIR_R,
                           ntree = 100,
                           mtry = 10,
                           importance = T,
                           na.action = na.omit)

RFpred_NIR_PVC = predict(RF_NIR_PVC, new_sub)
RFpred_NIR_PVC2 = predict(RF_NIR_PVC2, new_sub2)
RFpred_NIR_PVC3 = predict(RF_NIR_PVC3, new_sub3)
RFpred_NIR_PVC4 = predict(RF_NIR_PVC4, new_sub4)

PVC_NIR_obs = list(new_sub$plastMass, new_sub2$plastMass,
                   new_sub3$plastMass, new_sub4$plastMass)
PVC_NIR_preds = list(RFpred_NIR_PVC, RFpred_NIR_PVC2,
                     RFpred_NIR_PVC3, RFpred_NIR_PVC4)

PVC_NIR_stats = mapply(function(obs, pred) {
  c(
    ME   = ME(obs, pred),
    RMSE = RMSE(obs, pred),
    R2   = R2(obs, pred)
  )
}, PVC_NIR_obs, PVC_NIR_preds)

for (i in seq_along(PVC_NIR_preds)) {
  cat(paste0("\n======= MODEL ", i, " =======\n"))
  cat(sprintf("ME   : %.7f\n", PVC_NIR_stats["ME", i]))
  cat(sprintf("RMSE : %.7f\n", PVC_NIR_stats["RMSE", i]))
  cat(sprintf("R²   : %.7f\n", PVC_NIR_stats["R2", i]))
}

### polymer-specific CV
cv_RF_PVC_NIR = cv_random_forest(PVC_data, PVC_data$NIRspcARmovav)
cv_RF_PVC_NIR2 = cv_random_forest(PVC_data, PVC_data$NIRspcA)
cv_RF_PVC_NIR3 = cv_random_forest(PVC_data, PVC_data$NIRspcAmovav)
cv_RF_PVC_NIR4 = cv_random_forest(PVC_raw, PVC_raw$NIRspcA)

cv_RF_PVC_NIRresults = list(cv_RF_PVC_NIR, cv_RF_PVC_NIR2,
                            cv_RF_PVC_NIR3, cv_RF_PVC_NIR4)

for (i in seq_along(cv_RF_PVC_NIRresults)) {
  cat(paste0("\n======= MODEL ", i, " (10-fold CV) =======\n"))
  cat(sprintf("ME   : %.4f\n", cv_RF_PVC_NIRresults[[i]]$ME))
  cat(sprintf("RMSE : %.4f\n", cv_RF_PVC_NIRresults[[i]]$RMSE))
  cat(sprintf("R²   : %.4f\n", cv_RF_PVC_NIRresults[[i]]$R2))
}

## PET modeling
### prepare RF-ready datasets
PET_data_NIR1 = data.frame(plastMass = PET_data$MASS_mg, PET_data$NIRspcARmovav)
colnames(PET_data_NIR1) = c("plastMass", paste0("spec.", colnames(PET_data$NIRspcARmovav)))

PET_data_NIR2 = data.frame(plastMass = PET_data$MASS_mg, PET_data$NIRspcA)
colnames(PET_data_NIR2) = c("plastMass", paste0("spec.", colnames(PET_data$NIRspcA)))

PET_data_NIR3 = data.frame(plastMass = PET_data$MASS_mg, PET_data$NIRspcAmovav)
colnames(PET_data_NIR3) = c("plastMass", paste0("spec.", colnames(PET_data$NIRspcAmovav)))

PET_data_NIR_R = data.frame(plastMass = PET_raw$MASS_mg, PET_raw$NIRspcA)
colnames(PET_data_NIR_R) = c("plastMass", paste0("spec.", colnames(raw$NIRspcA)))


set.seed(4)

RF_NIR_PET = randomForest(plastMass ~ .,
                          data = PET_data_NIR1,
                          ntree = 150,
                          mtry = 50,
                          importance = T,
                          na.action = na.omit)

RF_NIR_PET2 = randomForest(plastMass ~ .,
                           data = PET_data_NIR2,
                           ntree = 150,
                           mtry = 10,
                           importance = T,
                           na.action = na.omit)

RF_NIR_PET3 = randomForest(plastMass ~ .,
                           data = PET_data_NIR3,
                           ntree = 200,
                           mtry = 50,
                           importance = T,
                           na.action = na.omit)

RF_NIR_PET4 = randomForest(plastMass ~ .,
                           data = PET_data_NIR_R,
                           ntree = 100,
                           mtry = 10,
                           importance = T,
                           na.action = na.omit)

RFpred_NIR_PET = predict(RF_NIR_PET, new_sub)
RFpred_NIR_PET2 = predict(RF_NIR_PET2, new_sub2)
RFpred_NIR_PET3 = predict(RF_NIR_PET3, new_sub3)
RFpred_NIR_PET4 = predict(RF_NIR_PET4, new_sub4)

PET_NIR_obs = list(new_sub$plastMass, new_sub2$plastMass,
                   new_sub3$plastMass, new_sub4$plastMass)
PET_NIR_preds = list(RFpred_NIR_PET, RFpred_NIR_PET2,
                     RFpred_NIR_PET3, RFpred_NIR_PET4)

PET_NIR_stats = mapply(function(obs, pred) {
  c(
    ME   = ME(obs, pred),
    RMSE = RMSE(obs, pred),
    R2   = R2(obs, pred)
  )
}, PET_NIR_obs, PET_NIR_preds)

for (i in seq_along(PET_NIR_preds)) {
  cat(paste0("\n======= MODEL ", i, " =======\n"))
  cat(sprintf("ME   : %.7f\n", PET_NIR_stats["ME", i]))
  cat(sprintf("RMSE : %.7f\n", PET_NIR_stats["RMSE", i]))
  cat(sprintf("R²   : %.7f\n", PET_NIR_stats["R2", i]))
}

### polymer-specific CV
cv_RF_PET_NIR = cv_random_forest(PET_data, PET_data$NIRspcARmovav)
cv_RF_PET_NIR2 = cv_random_forest(PET_data, PET_data$NIRspcA)
cv_RF_PET_NIR3 = cv_random_forest(PET_data, PET_data$NIRspcAmovav)
cv_RF_PET_NIR4 = cv_random_forest(PET_raw, PET_raw$NIRspcA)

cv_RF_PET_NIRresults = list(cv_RF_PET_NIR, cv_RF_PET_NIR2,
                            cv_RF_PET_NIR3, cv_RF_PET_NIR4)

for (i in seq_along(cv_RF_PET_NIRresults)) {
  cat(paste0("\n======= MODEL ", i, " (10-fold CV) =======\n"))
  cat(sprintf("ME   : %.4f\n", cv_results[[i]]$ME))
  cat(sprintf("RMSE : %.4f\n", cv_results[[i]]$RMSE))
  cat(sprintf("R²   : %.4f\n", cv_results[[i]]$R2))
}

## PE modeling
### prepare RF-ready datasets
PE_data_NIR1 = data.frame(plastMass = PE_data$MASS_mg, PE_data$NIRspcARmovav)
colnames(PE_data_NIR1) = c("plastMass", paste0("spec.", colnames(PE_data$NIRspcARmovav)))

PE_data_NIR2 = data.frame(plastMass = PE_data$MASS_mg, PE_data$NIRspcA)
colnames(PE_data_NIR2) = c("plastMass", paste0("spec.", colnames(PE_data$NIRspcA)))

PE_data_NIR3 = data.frame(plastMass = PE_data$MASS_mg, PE_data$NIRspcAmovav)
colnames(PE_data_NIR3) = c("plastMass", paste0("spec.", colnames(PE_data$NIRspcAmovav)))

PE_data_NIR_R = data.frame(plastMass = PE_raw$MASS_mg, PE_raw$NIRspcA)
colnames(PE_data_NIR_R) = c("plastMass", paste0("spec.", colnames(raw$NIRspcA)))


set.seed(4)

RF_NIR_PE = randomForest(plastMass ~ .,
                         data = PE_data_NIR1,
                         ntree = 150,
                         mtry = 50,
                         importance = T,
                         na.action = na.omit)

RF_NIR_PE2 = randomForest(plastMass ~ .,
                          data = PE_data_NIR2,
                          ntree = 150,
                          mtry = 10,
                          importance = T,
                          na.action = na.omit)

RF_NIR_PE3 = randomForest(plastMass ~ .,
                          data = PE_data_NIR3,
                          ntree = 200,
                          mtry = 50,
                          importance = T,
                          na.action = na.omit)

RF_NIR_PE4 = randomForest(plastMass ~ .,
                          data = PE_data_NIR_R,
                          ntree = 100,
                          mtry = 10,
                          importance = T,
                          na.action = na.omit)

RFpred_NIR_PE = predict(RF_NIR_PE, new_sub)
RFpred_NIR_PE2 = predict(RF_NIR_PE2, new_sub2)
RFpred_NIR_PE3 = predict(RF_NIR_PE3, new_sub3)
RFpred_NIR_PE4 = predict(RF_NIR_PE4, new_sub4)

PE_NIR_obs = list(new_sub$plastMass, new_sub2$plastMass,
                  new_sub3$plastMass, new_sub4$plastMass)
PE_NIR_preds = list(RFpred_NIR_PE, RFpred_NIR_PE2,
                    RFpred_NIR_PE3, RFpred_NIR_PE4)

PE_NIR_stats = mapply(function(obs, pred) {
  c(
    ME   = ME(obs, pred),
    RMSE = RMSE(obs, pred),
    R2   = R2(obs, pred)
  )
}, PE_NIR_obs, PE_NIR_preds)

for (i in seq_along(PE_NIR_preds)) {
  cat(paste0("\n======= MODEL ", i, " =======\n"))
  cat(sprintf("ME   : %.7f\n", PE_NIR_stats["ME", i]))
  cat(sprintf("RMSE : %.7f\n", PE_NIR_stats["RMSE", i]))
  cat(sprintf("R²   : %.7f\n", PE_NIR_stats["R2", i]))
}

### polymer-specific CV
cv_RF_PE_NIR = cv_random_forest(PE_data, PE_data$NIRspcARmovav)
cv_RF_PE_NIR2 = cv_random_forest(PE_data, PE_data$NIRspcA)
cv_RF_PE_NIR3 = cv_random_forest(PE_data, PE_data$NIRspcAmovav)
cv_RF_PE_NIR4 = cv_random_forest(PE_raw, PE_raw$NIRspcA)

cv_RF_PE_NIRresults = list(cv_RF_PE_NIR, cv_RF_PE_NIR2,
                           cv_RF_PE_NIR3, cv_RF_PE_NIR4)

for (i in seq_along(cv_RF_PE_NIRresults)) {
  cat(paste0("\n======= MODEL ", i, " (10-fold CV) =======\n"))
  cat(sprintf("ME   : %.4f\n", cv_results[[i]]$ME))
  cat(sprintf("RMSE : %.4f\n", cv_results[[i]]$RMSE))
  cat(sprintf("R²   : %.4f\n", cv_results[[i]]$R2))
}

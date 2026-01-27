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
matplot(colnames(data$spcA), t(data$spcA),
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

## convert columns names (integer) to numeric
wavs = as.numeric(colnames(data$spcAmovav))

matplot(colnames(data$spcARmovav), t(data$spcARmovav),
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
  SST = sum((obs - mean(obs, na.rm = T))^2, na.rm = T) # squared total sum
  R2 = 1 - SSE / SST
  return(R2)
}

################################################################################

# DEFAULT PARAMETERIZATION (neighbors = 0; comittees = 1)
set.seed(1)
require(Cubist)
CUB_mod_mass = cubist(datC$spcARmovav, datC$MASS_mg)

CUB_mod_mass2 = cubist(datC$spcA, datC$MASS_mg)

CUB_mod_mass3 = cubist(datC$spcAmovav, datC$MASS_mg)

CUB_mod_mass4 = cubist(rawC$spcA, rawC$MASS_mg)

# predictions
## calibration
CUBpredC = predict(CUB_mod_mass, datC$spcARmovav)
CUBpredC2 = predict(CUB_mod_mass2, datC$spcA)
CUBpredC3 = predict(CUB_mod_mass3, datC$spcAmovav)
CUBpredC4 = predict(CUB_mod_mass4, rawC$spcA)

## validation
CUBpredV = predict(CUB_mod_mass, datV$spcARmovav)
CUBpredV2 = predict(CUB_mod_mass2, datV$spcA)
CUBpredV3 = predict(CUB_mod_mass3, datV$spcAmovav)
CUBpredV4 = predict(CUB_mod_mass4, rawV$spcA)

## external
CUBpred_new = predict(CUB_mod_mass, new$spcARmovav)
CUBpred_new2 = predict(CUB_mod_mass2, new$spcA) 
CUBpred_new3 = predict(CUB_mod_mass3, new$spcAmovav)
CUBpred_new4 = predict(CUB_mod_mass4, raw_new$spcA)

# ME(datC$MASS_mg, CUBpredC)
# RMSE(datC$MASS_mg, CUBpredC)
# R2(datC$MASS_mg, CUBpredC)
# 
# ME(datV$MASS_mg, CUBpredV)
# RMSE(datV$MASS_mg, CUBpredV)
# R2(datV$MASS_mg, CUBpredV)

# evaluate quality of predictions
## observed responses
calib_obs = c(rep(list(datC$MASS_mg), 3), list(rawC$MASS_mg))

valid_obs = c(rep(list(datV$MASS_mg), 3), list(rawV$MASS_mg))

## group calibration predictions
calib_preds = list(CUBpredC, CUBpredC2, CUBpredC3, CUBpredC4)
## group validation predictions
valid_preds = list(CUBpredV, CUBpredV2, CUBpredV3, CUBpredV4)

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
external_preds = list(CUBpred_new, CUBpred_new2, CUBpred_new3, CUBpred_new4)

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

# check models dimensionality (mosty important wavelengths)
require(caret)
sum(varImp(CUB_mod_mass) > 0) # 29...
sum(varImp(CUB_mod_mass2) > 0) # 30...
sum(varImp(CUB_mod_mass3) > 0) # 28...
sum(varImp(CUB_mod_mass4) > 0) # 22...

################################################################################

# TUNE PARAMETERIZATION
set.seed(12)

# createDataPartition() exists, but I already have splitted my sets...
# 10-fold cross-validation to tune the models over two parameters
grid = expand.grid(committees = c(1, 10, 50, 100),
                   neighbors = c(0, 1, 5, 9))

# boost shit up...
require(doParallel)

t0 = Sys.time()

cl = makeCluster(parallel::detectCores() - 1) # n. of cores - 1...
registerDoParallel(cl)

CUB_tuned = caret::train(
  x = datC$spcARmovav,
  y = datC$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T) # default n. folds for "cv" = 10...
)

ggplot(CUB_tuned) +
  labs(title = "SGf + 5nm resample + SNV + moving average Preprocessing")


# The object selected and fit the final model (with the best results).
## The test data are predicted using:
### predict(CUB_tuned, datC$spcARmovav) 


CUB_tuned2 = caret::train(
  x = datC$spcA,
  y = datC$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T)
)

ggplot(CUB_tuned2) +
  labs(title = "SGF only Preprocessing")


CUB_tuned3 = caret::train(
  x = datC$spcAmovav,
  y = datC$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T)
)

ggplot(CUB_tuned3) +
  labs(title = "SGf + SNV + moving average Preprocessing")


CUB_tuned4 = caret::train(
  x = rawC$spcA,
  y = rawC$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T)
)

ggplot(CUB_tuned4) +
  labs(title = "No Preprocessing (raw reflectance data)")


t1 = Sys.time()
cat("Training time:", round(t1 - t0, 2), "minutes\n")

stopCluster(cl)
registerDoSEQ()


#==============================================================================#
#                     FULL INTERNAL DATASET TRANING!                           #
#==============================================================================#

CUB_mod_FULL = cubist(data$spcARmovav, data$MASS_mg,
                      committees = 1, neighbors = 0)

CUB_mod_FULL2 = cubist(data$spcA, data$MASS_mg,
                       committees = 1, neighbors = 0)

CUB_mod_FULL3 = cubist(data$spcAmovav, data$MASS_mg,
                       committees = 100, neighbors = 5)

CUB_mod_FULL4 = cubist(raw$spcA, raw$MASS_mg,
                       committees = 1, neighbors = 0)

## external validation
CUBpred_FULL = predict(CUB_mod_FULL, new$spcARmovav)
CUBpred_FULL2 = predict(CUB_mod_FULL2, new$spcA) 
CUBpred_FULL3 = predict(CUB_mod_FULL3, new$spcAmovav)
CUBpred_FULL4 = predict(CUB_mod_FULL4, raw_new$spcA)

full_obs = c(rep(list(new$MASS_mg), 3), list(raw_new$MASS_mg))
full_preds = list(CUBpred_FULL, CUBpred_FULL2, CUBpred_FULL3, CUBpred_FULL4)

full_stats = mapply(function(obs, pred) {
  c(
    ME   = ME(obs, pred),
    RMSE = RMSE(obs, pred),
    R2   = R2(obs, pred)
  )
}, full_obs, full_preds)

for (i in seq_along(full_preds)) {
  cat(paste0("\n======= MODEL ", i, " =======\n"))
  cat(sprintf("ME   : %.7f\n", full_stats["ME", i]))
  cat(sprintf("RMSE : %.7f\n", full_stats["RMSE", i]))
  cat(sprintf("R²   : %.7f\n", full_stats["R2", i]))
}

################################################################################

# k-fold cross-validation
cv_cubist_model = function(data, spc_matrix, committees = 1, neighbors = 0, nfolds = 4, seed = 999) {
  set.seed(seed)
  
  # check minimum stratum size and adjust nfolds if needed
  strata = interaction(data$POLYMER, data$MASS_mg, data$SIZE_CODE)
  min_stratum_size = min(table(strata))
  
  if (nfolds > min_stratum_size) {
    cat("Note: Minimum stratum size is", min_stratum_size)
    nfolds = min_stratum_size
  }
  
  data$foldCV = NA
  
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
    
    # convert to proper format for Cubist
    train_spc = as.data.frame(as.matrix(spc_matrix)[data$foldCV != i, ])
    valid_spc = as.data.frame(as.matrix(spc_matrix)[data$foldCV == i, ])
    
    # ensure consistent column names
    colnames(train_spc) = colnames(valid_spc) = paste0("V", 1:ncol(train_spc))
    
    # train Cubist model
    cubist_model = cubist(x = train_spc, 
                          y = train$MASS_mg,
                          committees = committees,
                          neighbors = neighbors)
    
    preds = predict(cubist_model, newdata = valid_spc)
    
    valTable$pred[valTable$SAMPLE_ID %in% valid$SAMPLE_ID] = preds
    
    cat("Fold number", i, "done for Cubist. Validation samples:", nrow(valid), "\n")
  }
  
  list(
    committees = committees,
    neighbors = neighbors,
    nfolds_used = nfolds,
    min_stratum_size = min_stratum_size,
    ME = ME(valTable$obs, valTable$pred),
    RMSE = RMSE(valTable$obs, valTable$pred),
    R2 = R2(valTable$obs, valTable$pred),
    table = valTable
  )
}

cv_CUBISTmod = cv_cubist_model(data, data$spcARmovav)
cv_CUBISTmod2 = cv_cubist_model(data, data$spcA)
cv_CUBISTmod3 = cv_cubist_model(data, data$spcAmovav)
cv_CUBISTmod4 = cv_cubist_model(raw, raw$spcA)

cv_cubist_results = list(cv_CUBISTmod, cv_CUBISTmod2, cv_CUBISTmod3, cv_CUBISTmod4)

for (i in seq_along(cv_cubist_results)) {
  cat(paste0("\n======= CUBIST MODEL ", i, " (", cv_cubist_results[[i]]$nfolds_used, "-fold CV) =======\n"))
  cat(sprintf("ME         : %.4f\n", cv_cubist_results[[i]]$ME))
  cat(sprintf("RMSE       : %.4f\n", cv_cubist_results[[i]]$RMSE))
  cat(sprintf("R²         : %.4f\n", cv_cubist_results[[i]]$R2))
}

################################################################################


## retry parameterization tuning
t0 = Sys.time()

cl = makeCluster(parallel::detectCores() - 1)
registerDoParallel(cl)

CUB_FULLtuned = caret::train(
  x = data$spcARmovav,
  y = data$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T) # default n. folds for "cv" = 10...
)

ggplot(CUB_FULLtuned) +
  labs(title = "SGf + 5nm resample + SNV + moving average Preprocessing")
CUB_FULLtuned


CUB_FULLtuned2 = caret::train(
  x = data$spcA,
  y = data$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T)
)

ggplot(CUB_FULLtuned2) +
  labs(title = "SGF only Preprocessing")
CUB_FULLtuned2


CUB_FULLtuned3 = caret::train(
  x = data$spcAmovav,
  y = data$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T)
)

ggplot(CUB_FULLtuned3) +
  labs(title = "SGf + SNV + moving average Preprocessing")
CUB_FULLtuned3


CUB_FULLtuned4 = caret::train(
  x = raw$spcA,
  y = raw$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T)
)

ggplot(CUB_FULLtuned4) +
  labs(title = "No Preprocessing (raw reflectance data)")
CUB_FULLtuned4

t1 = Sys.time()
cat("Training time:", round(t1 - t0, 2), "minutes\n")

stopCluster(cl)
registerDoSEQ()

### retry full internal dataset training
set.seed(111) # before was `999``...`
CUB_mod_full = cubist(data$spcARmovav, data$MASS_mg,
                      committees = 100, neighbors = 0)

CUB_mod_full2 = cubist(data$spcA, data$MASS_mg,
                       committees = 50, neighbors = 0)

CUB_mod_full3 = cubist(data$spcAmovav, data$MASS_mg,
                       committees = 100, neighbors = 5)

CUB_mod_full4 = cubist(raw$spcA, raw$MASS_mg,
                       committees = 1, neighbors = 5)

#### external validation
CUBpred_full = predict(CUB_mod_full, new$spcARmovav)
CUBpred_full2 = predict(CUB_mod_full2, new$spcA) 
CUBpred_full3 = predict(CUB_mod_full3, new$spcAmovav)
CUBpred_full4 = predict(CUB_mod_full4, raw_new$spcA)

full_obs2 = c(rep(list(new$MASS_mg), 3), list(raw_new$MASS_mg))
full_preds2 = list(CUBpred_full, CUBpred_full2, CUBpred_full3, CUBpred_full4)

full_stats2 = mapply(function(obs, pred) {
  c(
    ME   = ME(obs, pred),
    RMSE = RMSE(obs, pred),
    R2   = R2(obs, pred)
  )
}, full_obs2, full_preds2)

for (i in seq_along(full_preds)) {
  cat(paste0("\n======= MODEL ", i, " =======\n"))
  cat(sprintf("ME   : %.7f\n", full_stats2["ME", i]))
  cat(sprintf("RMSE : %.7f\n", full_stats2["RMSE", i]))
  cat(sprintf("R²   : %.7f\n", full_stats2["R2", i]))
}


# residual visualization
residLEAK = CUBpred_full - new$MASS_mg
residLEAK2 = CUBpred_full2 - new$MASS_mg
residLEAK3 = CUBpred_full3 - new$MASS_mg
residLEAK4 = CUBpred_full4 - raw_new$MASS_mg

residLEAKlog = log(CUBpred_full) - log(new$MASS_mg)
residLEAKlog2 = log(CUBpred_full2) - log(new$MASS_mg)
residLEAKlog3 = log(CUBpred_full3) - log(new$MASS_mg)
residLEAKlog4 = log(CUBpred_full4) - log(raw_new$MASS_mg)

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


## CROSS-VALIDATION (PP)
t0 = Sys.time()

cl = makeCluster(parallel::detectCores() - 1)
registerDoParallel(cl)

CUB_PPtuned = caret::train(
  x = PP_data$spcARmovav,
  y = PP_data$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T) # default n. folds for "cv" = 10...
)

ggplot(CUB_PPtuned) +
  labs(title = "SGf + 5nm resample + SNV + moving average Preprocessing (PP)")
CUB_PPtuned


CUB_PPtuned2 = caret::train(
  x = PP_data$spcA,
  y = PP_data$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T)
)

ggplot(CUB_PPtuned2) +
  labs(title = "SGF only Preprocessing (PP)")
CUB_PPtuned2


CUB_PPtuned3 = caret::train(
  x = PP_data$spcAmovav,
  y = PP_data$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T)
)

ggplot(CUB_PPtuned3) +
  labs(title = "SGf + SNV + moving average Preprocessing (PP)")
CUB_PPtuned3


CUB_PPtuned4 = caret::train(
  x = PP_raw$spcA,
  y = PP_raw$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T)
)

ggplot(CUB_PPtuned4) +
  labs(title = "No Preprocessing - raw reflectance data (PP)")
CUB_PPtuned4

t1 = Sys.time()
cat("Training time:", round(t1 - t0, 2), "minutes\n")

stopCluster(cl)
registerDoSEQ()


## EXTERNAL VALIDATION (PP)
CUB_mod_PP = cubist(PP_data$spcARmovav, PP_data$MASS_mg,
                    committees = 50, neighbors = 1)

CUB_mod_PP2 = cubist(PP_data$spcA, PP_data$MASS_mg,
                     committees = 1, neighbors = 9)

CUB_mod_PP3 = cubist(PP_data$spcAmovav, PP_data$MASS_mg,
                     committees = 100, neighbors = 1)

CUB_mod_PP4 = cubist(PP_raw$spcA, PP_raw$MASS_mg,
                     committees = 1, neighbors = 0)

CUBpred_PP = predict(CUB_mod_PP, new$spcARmovav)
CUBpred_PP2 = predict(CUB_mod_PP2, new$spcA)
CUBpred_PP3 = predict(CUB_mod_PP3, new$spcAmovav)
CUBpred_PP4 = predict(CUB_mod_PP4, raw_new$spcA)

PP_obs = c(rep(list(new$PP_mg), 3), list(raw_new$PP_mg))
PP_preds = list(CUBpred_PP, CUBpred_PP2, CUBpred_PP3, CUBpred_PP4)

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
## CROSS-VALIDATION (PVC)
t0 = Sys.time()

cl = makeCluster(parallel::detectCores() - 1)
registerDoParallel(cl)

CUB_PVCtuned = caret::train(
  x = PVC_data$spcARmovav,
  y = PVC_data$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T) # default n. folds for "cv" = 10...
)

ggplot(CUB_PVCtuned) +
  labs(title = "SGf + 5nm resample + SNV + moving average Preprocessing (PVC)")
CUB_PVCtuned


CUB_PVCtuned2 = caret::train(
  x = PVC_data$spcA,
  y = PVC_data$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T)
)

ggplot(CUB_PVCtuned2) +
  labs(title = "SGF only Preprocessing (PVC)")
CUB_PVCtuned2


CUB_PVCtuned3 = caret::train(
  x = PVC_data$spcAmovav,
  y = PVC_data$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T)
)

ggplot(CUB_PVCtuned3) +
  labs(title = "SGf + SNV + moving average Preprocessing (PVC)")
CUB_PVCtuned3


CUB_PVCtuned4 = caret::train(
  x = PVC_raw$spcA,
  y = PVC_raw$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T)
)

ggplot(CUB_PVCtuned4) +
  labs(title = "No Preprocessing - raw reflectance data (PVC)")
CUB_PVCtuned4

t1 = Sys.time()
cat("Training time:", round(t1 - t0, 2), "minutes\n")

stopCluster(cl)
registerDoSEQ()


## EXTERNAL VALIDATION (PVC)
CUB_mod_PVC = cubist(PVC_data$spcARmovav, PVC_data$MASS_mg,
                     committees = 100, neighbors = 1)

CUB_mod_PVC2 = cubist(PVC_data$spcA, PVC_data$MASS_mg,
                      committees = 50, neighbors = 0)

CUB_mod_PVC3 = cubist(PVC_data$spcAmovav, PVC_data$MASS_mg,
                      committees = 100, neighbors = 1)

CUB_mod_PVC4 = cubist(PVC_raw$spcA, PVC_raw$MASS_mg,
                      committees = 100, neighbors = 0)

CUBpred_PVC = predict(CUB_mod_PVC, new$spcARmovav)
CUBpred_PVC2 = predict(CUB_mod_PVC2, new$spcA)
CUBpred_PVC3 = predict(CUB_mod_PVC3, new$spcAmovav)
CUBpred_PVC4 = predict(CUB_mod_PVC4, raw_new$spcA)

PVC_obs = c(rep(list(new$PVC_mg), 3), list(raw_new$PVC_mg))
PVC_preds = list(CUBpred_PVC, CUBpred_PVC2, CUBpred_PVC3, CUBpred_PVC4)

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
## CROSS-VALIDATION (PET)
t0 = Sys.time()

cl = makeCluster(parallel::detectCores() - 1)
registerDoParallel(cl)

CUB_PETtuned = caret::train(
  x = PET_data$spcARmovav,
  y = PET_data$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T) # default n. folds for "cv" = 10...
)

ggplot(CUB_PETtuned) +
  labs(title = "SGf + 5nm resample + SNV + moving average Preprocessing (PET)")
CUB_PETtuned


CUB_PETtuned2 = caret::train(
  x = PET_data$spcA,
  y = PET_data$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T)
)

ggplot(CUB_PETtuned2) +
  labs(title = "SGF only Preprocessing (PET)")
CUB_PETtuned2


CUB_PETtuned3 = caret::train(
  x = PET_data$spcAmovav,
  y = PET_data$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T)
)

ggplot(CUB_PETtuned3) +
  labs(title = "SGf + SNV + moving average Preprocessing (PET)")
CUB_PETtuned3


CUB_PETtuned4 = caret::train(
  x = PET_raw$spcA,
  y = PET_raw$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T)
)

ggplot(CUB_PETtuned4) +
  labs(title = "No Preprocessing - raw reflectance data (PET)")
CUB_PETtuned4

t1 = Sys.time()
cat("Training time:", round(t1 - t0, 2), "minutes\n")

stopCluster(cl)
registerDoSEQ()


## EXTERNAL VALIDATION (PET)
CUB_mod_PET = cubist(PET_data$spcARmovav, PET_data$MASS_mg,
                     committees = 10, neighbors = 9)

CUB_mod_PET2 = cubist(PET_data$spcA, PET_data$MASS_mg,
                      committees = 10, neighbors = 0)

CUB_mod_PET3 = cubist(PET_data$spcAmovav, PET_data$MASS_mg,
                      committees = 10, neighbors = 5)

CUB_mod_PET4 = cubist(PET_raw$spcA, PET_raw$MASS_mg,
                      committees = 1, neighbors = 0)

CUBpred_PET = predict(CUB_mod_PET, new$spcARmovav)
CUBpred_PET2 = predict(CUB_mod_PET2, new$spcA)
CUBpred_PET3 = predict(CUB_mod_PET3, new$spcAmovav)
CUBpred_PET4 = predict(CUB_mod_PET4, raw_new$spcA)

PET_obs = c(rep(list(new$PET_mg), 3), list(raw_new$PET_mg))
PET_preds = list(CUBpred_PET, CUBpred_PET2, CUBpred_PET3, CUBpred_PET4)

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
## CROSS-VALIDATION (PE)
t0 = Sys.time()

cl = makeCluster(parallel::detectCores() - 1)
registerDoParallel(cl)

CUB_PEtuned = caret::train(
  x = PE_data$spcARmovav,
  y = PE_data$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T) # default n. folds for "cv" = 10...
)

ggplot(CUB_PEtuned) +
  labs(title = "SGf + 5nm resample + SNV + moving average Preprocessing (PE)")
CUB_PEtuned


CUB_PEtuned2 = caret::train(
  x = PE_data$spcA,
  y = PE_data$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T)
)

ggplot(CUB_PEtuned2) +
  labs(title = "SGF only Preprocessing (PE)")
CUB_PEtuned2


CUB_PEtuned3 = caret::train(
  x = PE_data$spcAmovav,
  y = PE_data$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T)
)

ggplot(CUB_PEtuned3) +
  labs(title = "SGf + SNV + moving average Preprocessing (PE)")
CUB_PEtuned3


CUB_PEtuned4 = caret::train(
  x = PE_raw$spcA,
  y = PE_raw$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T)
)

ggplot(CUB_PEtuned4) +
  labs(title = "No Preprocessing - raw reflectance data (PE)")
CUB_PEtuned4

t1 = Sys.time()
cat("Training time:", round(t1 - t0, 2), "minutes\n")

stopCluster(cl)
registerDoSEQ()


## EXTERNAL VALIDATION (PE)
CUB_mod_PE = cubist(PE_data$spcARmovav, PE_data$MASS_mg,
                    committees = 50, neighbors = 0)

CUB_mod_PE2 = cubist(PE_data$spcA, PE_data$MASS_mg,
                     committees = 50, neighbors = 0)

CUB_mod_PE3 = cubist(PE_data$spcAmovav, PE_data$MASS_mg,
                     committees = 1, neighbors = 9)

CUB_mod_PE4 = cubist(PE_raw$spcA, PE_raw$MASS_mg,
                     committees = 50, neighbors = 5)

CUBpred_PE = predict(CUB_mod_PE, new$spcARmovav)
CUBpred_PE2 = predict(CUB_mod_PE2, new$spcA)
CUBpred_PE3 = predict(CUB_mod_PE3, new$spcAmovav)
CUBpred_PE4 = predict(CUB_mod_PE4, raw_new$spcA)

PE_obs = c(rep(list(new$PE_mg), 3), list(raw_new$PE_mg))
PE_preds = list(CUBpred_PE, CUBpred_PE2, CUBpred_PE3, CUBpred_PE4)

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

t0 = Sys.time()

cl = makeCluster(parallel::detectCores() - 1)
registerDoParallel(cl)

CUB_NIRtuned = caret::train(
  x = data$NIRspcARmovav,
  y = data$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T) # default n. folds for "cv" = 10...
)

ggplot(CUB_NIRtuned) +
  labs(title = "SGf + 5nm resample + SNV + moving average Preprocessing")
CUB_NIRtuned


CUB_NIRtuned2 = caret::train(
  x = data$NIRspcA,
  y = data$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T)
)

ggplot(CUB_NIRtuned2) +
  labs(title = "SGF only Preprocessing")
CUB_NIRtuned2


CUB_NIRtuned3 = caret::train(
  x = data$NIRspcAmovav,
  y = data$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T)
)

ggplot(CUB_NIRtuned3) +
  labs(title = "SGf + SNV + moving average Preprocessing")
CUB_NIRtuned3


CUB_NIRtuned4 = caret::train(
  x = raw$NIRspcA,
  y = raw$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T)
)

ggplot(CUB_NIRtuned4) +
  labs(title = "No Preprocessing (raw reflectance data)")
CUB_NIRtuned4

t1 = Sys.time()
cat("Training time:", round(t1 - t0, 2), "minutes\n")

stopCluster(cl)
registerDoSEQ()


### retry full internal dataset training
CUB_mod_NIR = cubist(data$NIRspcARmovav, data$MASS_mg,
                      committees = 100, neighbors = 0)

CUB_mod_NIR2 = cubist(data$NIRspcA, data$MASS_mg,
                       committees = 50, neighbors = 5)

CUB_mod_NIR3 = cubist(data$NIRspcAmovav, data$MASS_mg,
                       committees = 10, neighbors = 5)

CUB_mod_NIR4 = cubist(raw$NIRspcA, raw$MASS_mg,
                       committees = 10, neighbors = 1)

#### external validation
CUBpred_NIR = predict(CUB_mod_NIR, new$NIRspcARmovav)
CUBpred_NIR2 = predict(CUB_mod_NIR2, new$NIRspcA) 
CUBpred_NIR3 = predict(CUB_mod_NIR3, new$NIRspcAmovav)
CUBpred_NIR4 = predict(CUB_mod_NIR4, raw_new$NIRspcA)

NIR_obs = c(rep(list(new$MASS_mg), 3), list(raw_new$MASS_mg))
NIR_preds = list(CUBpred_NIR, CUBpred_NIR2, CUBpred_NIR3, CUBpred_NIR4)

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
#                MODELLING FOR SINGLE POLYMERS QUANTIFICATION (NIR)            #
#==============================================================================#
## CROSS-VALIDATION (PP)
t0 = Sys.time()

cl = makeCluster(parallel::detectCores() - 1)
registerDoParallel(cl)

CUB_PP_NIRtuned = caret::train(
  x = PP_data$NIRspcARmovav,
  y = PP_data$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T) # default n. folds for "cv" = 10...
)

ggplot(CUB_PP_NIRtuned) +
  labs(title = "SGf + 5nm resample + SNV + moving average Preprocessing (PP)")
CUB_PP_NIRtuned


CUB_PP_NIRtuned2 = caret::train(
  x = PP_data$NIRspcA,
  y = PP_data$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T)
)

ggplot(CUB_PP_NIRtuned2) +
  labs(title = "SGF only Preprocessing (PP)")
CUB_PP_NIRtuned2


CUB_PP_NIRtuned3 = caret::train(
  x = PP_data$NIRspcAmovav,
  y = PP_data$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T)
)

ggplot(CUB_PP_NIRtuned3) +
  labs(title = "SGf + SNV + moving average Preprocessing (PP)")
CUB_PP_NIRtuned3


CUB_PP_NIRtuned4 = caret::train(
  x = PP_raw$NIRspcA,
  y = PP_raw$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T)
)

ggplot(CUB_PP_NIRtuned4) +
  labs(title = "No Preprocessing - raw reflectance data (PP)")
CUB_PP_NIRtuned4

t1 = Sys.time()
cat("Training time:", round(t1 - t0, 2), "minutes\n")

stopCluster(cl)
registerDoSEQ()


## EXTERNAL VALIDATION (PP)
CUB_NIR_PP = cubist(PP_data$NIRspcARmovav, PP_data$MASS_mg,
                    committees = 1, neighbors = 1)

CUB_NIR_PP2 = cubist(PP_data$NIRspcA, PP_data$MASS_mg,
                     committees = 1, neighbors = 5)

CUB_NIR_PP3 = cubist(PP_data$NIRspcAmovav, PP_data$MASS_mg,
                     committees = 10, neighbors = 1)

CUB_NIR_PP4 = cubist(PP_raw$NIRspcA, PP_raw$MASS_mg,
                     committees = 100, neighbors = 1)

CUBpred_NIR_PP = predict(CUB_NIR_PP, new$NIRspcARmovav)
CUBpred_NIR_PP2 = predict(CUB_NIR_PP2, new$NIRspcA)
CUBpred_NIR_PP3 = predict(CUB_NIR_PP3, new$NIRspcAmovav)
CUBpred_NIR_PP4 = predict(CUB_NIR_PP4, raw_new$NIRspcA)

PP_NIR_obs = c(rep(list(new$PP_mg), 3), list(raw_new$PP_mg))
PP_NIR_preds = list(CUBpred_NIR_PP, CUBpred_NIR_PP2,
                    CUBpred_NIR_PP3, CUBpred_NIR_PP4)

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

#==============================================================================#
## CROSS-VALIDATION (PVC)
t0 = Sys.time()

cl = makeCluster(parallel::detectCores() - 1)
registerDoParallel(cl)

CUB_PVC_NIRtuned = caret::train(
  x = PVC_data$NIRspcARmovav,
  y = PVC_data$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T) # default n. folds for "cv" = 10...
)

ggplot(CUB_PVC_NIRtuned) +
  labs(title = "SGf + 5nm resample + SNV + moving average Preprocessing (PVC)")
CUB_PVC_NIRtuned


CUB_PVC_NIRtuned2 = caret::train(
  x = PVC_data$NIRspcA,
  y = PVC_data$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T)
)

ggplot(CUB_PVC_NIRtuned2) +
  labs(title = "SGF only Preprocessing (PVC)")
CUB_PVC_NIRtuned2


CUB_PVC_NIRtuned3 = caret::train(
  x = PVC_data$NIRspcAmovav,
  y = PVC_data$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T)
)

ggplot(CUB_PVC_NIRtuned3) +
  labs(title = "SGf + SNV + moving average Preprocessing (PVC)")
CUB_PVC_NIRtuned3


CUB_PVC_NIRtuned4 = caret::train(
  x = PVC_raw$NIRspcA,
  y = PVC_raw$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T)
)

ggplot(CUB_PVC_NIRtuned4) +
  labs(title = "No Preprocessing - raw reflectance data (PVC)")
CUB_PVC_NIRtuned4

t1 = Sys.time()
cat("Training time:", round(t1 - t0, 2), "minutes\n")

stopCluster(cl)
registerDoSEQ()


## EXTERNAL VALIDATION (PVC)
CUB_NIR_PVC = cubist(PVC_data$NIRspcARmovav, PVC_data$MASS_mg,
                     committees = 10, neighbors = 0)

CUB_NIR_PVC2 = cubist(PVC_data$NIRspcA, PVC_data$MASS_mg,
                      committees = 100, neighbors = 0)

CUB_NIR_PVC3 = cubist(PVC_data$NIRspcAmovav, PVC_data$MASS_mg,
                      committees = 100, neighbors = 5)

CUB_NIR_PVC4 = cubist(PVC_raw$NIRspcA, PVC_raw$MASS_mg,
                      committees = 50, neighbors = 0)

CUBpred_NIR_PVC = predict(CUB_NIR_PVC, new$NIRspcARmovav)
CUBpred_NIR_PVC2 = predict(CUB_NIR_PVC2, new$NIRspcA)
CUBpred_NIR_PVC3 = predict(CUB_NIR_PVC3, new$NIRspcAmovav)
CUBpred_NIR_PVC4 = predict(CUB_NIR_PVC4, raw_new$NIRspcA)

PVC_NIR_obs = c(rep(list(new$PVC_mg), 3), list(raw_new$PVC_mg))
PVC_NIR_preds = list(CUBpred_NIR_PVC, CUBpred_NIR_PVC2,
                     CUBpred_NIR_PVC3, CUBpred_NIR_PVC4)

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

#==============================================================================#
## CROSS-VALIDATION (PET)
t0 = Sys.time()

cl = makeCluster(parallel::detectCores() - 1)
registerDoParallel(cl)

CUB_PET_NIRtuned = caret::train(
  x = PET_data$NIRspcARmovav,
  y = PET_data$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T) # default n. folds for "cv" = 10...
)

ggplot(CUB_PET_NIRtuned) +
  labs(title = "SGf + 5nm resample + SNV + moving average Preprocessing (PET)")
CUB_PET_NIRtuned


CUB_PET_NIRtuned2 = caret::train(
  x = PET_data$NIRspcA,
  y = PET_data$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T)
)

ggplot(CUB_PET_NIRtuned2) +
  labs(title = "SGF only Preprocessing (PET)")
CUB_PET_NIRtuned2


CUB_PET_NIRtuned3 = caret::train(
  x = PET_data$NIRspcAmovav,
  y = PET_data$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T)
)

ggplot(CUB_PET_NIRtuned3) +
  labs(title = "SGf + SNV + moving average Preprocessing (PET)")
CUB_PET_NIRtuned3


CUB_PET_NIRtuned4 = caret::train(
  x = PET_raw$NIRspcA,
  y = PET_raw$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T)
)

ggplot(CUB_PET_NIRtuned4) +
  labs(title = "No Preprocessing - raw reflectance data (PET)")
CUB_PET_NIRtuned4

t1 = Sys.time()
cat("Training time:", round(t1 - t0, 2), "minutes\n")

stopCluster(cl)
registerDoSEQ()


## EXTERNAL VALIDATION (PET)
CUB_NIR_PET = cubist(PET_data$NIRspcARmovav, PET_data$MASS_mg,
                     committees = 50, neighbors = 0)

CUB_NIR_PET2 = cubist(PET_data$NIRspcA, PET_data$MASS_mg,
                      committees = 100, neighbors = 5)

CUB_NIR_PET3 = cubist(PET_data$NIRspcAmovav, PET_data$MASS_mg,
                      committees = 100, neighbors = 1)

CUB_NIR_PET4 = cubist(PET_raw$NIRspcA, PET_raw$MASS_mg,
                      committees = 1, neighbors = 0)

CUBpred_NIR_PET = predict(CUB_NIR_PET, new$NIRspcARmovav)
CUBpred_NIR_PET2 = predict(CUB_NIR_PET2, new$NIRspcA)
CUBpred_NIR_PET3 = predict(CUB_NIR_PET3, new$NIRspcAmovav)
CUBpred_NIR_PET4 = predict(CUB_NIR_PET4, raw_new$NIRspcA)

PET_NIR_obs = c(rep(list(new$PET_mg), 3), list(raw_new$PET_mg))
PET_NIR_preds = list(CUBpred_NIR_PET, CUBpred_NIR_PET2,
                     CUBpred_NIR_PET3, CUBpred_NIR_PET4)

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

#==============================================================================#
## CROSS-VALIDATION (PE)
t0 = Sys.time()

cl = makeCluster(parallel::detectCores() - 1)
registerDoParallel(cl)

CUB_PE_NIRtuned = caret::train(
  x = PE_data$NIRspcARmovav,
  y = PE_data$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T) # default n. folds for "cv" = 10...
)

ggplot(CUB_PE_NIRtuned) +
  labs(title = "SGf + 5nm resample + SNV + moving average Preprocessing (PE)")
CUB_PE_NIRtuned


CUB_PE_NIRtuned2 = caret::train(
  x = PE_data$NIRspcA,
  y = PE_data$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T)
)

ggplot(CUB_PE_NIRtuned2) +
  labs(title = "SGF only Preprocessing (PE)")
CUB_PE_NIRtuned2


CUB_PE_NIRtuned3 = caret::train(
  x = PE_data$NIRspcAmovav,
  y = PE_data$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T)
)

ggplot(CUB_PE_NIRtuned3) +
  labs(title = "SGf + SNV + moving average Preprocessing (PE)")
CUB_PE_NIRtuned3


CUB_PE_NIRtuned4 = caret::train(
  x = PE_raw$NIRspcA,
  y = PE_raw$MASS_mg,
  method = "cubist",
  tuneGrid = grid,
  trControl = trainControl(
    method = "cv",
    verboseIter = T,
    allowParallel = T)
)

ggplot(CUB_PE_NIRtuned4) +
  labs(title = "No Preprocessing - raw reflectance data (PE)")
CUB_PE_NIRtuned4

t1 = Sys.time()
cat("Training time:", round(t1 - t0, 2), "minutes\n")

stopCluster(cl)
registerDoSEQ()


## EXTERNAL VALIDATION (PE)
CUB_NIR_PE = cubist(PE_data$NIRspcARmovav, PE_data$MASS_mg,
                    committees = 100, neighbors = 0)

CUB_NIR_PE2 = cubist(PE_data$NIRspcA, PE_data$MASS_mg,
                     committees = 100, neighbors = 0)

CUB_NIR_PE3 = cubist(PE_data$NIRspcAmovav, PE_data$MASS_mg,
                     committees = 100, neighbors = 5)

CUB_NIR_PE4 = cubist(PE_raw$NIRspcA, PE_raw$MASS_mg,
                     committees = 10, neighbors = 0)

CUBpred_NIR_PE = predict(CUB_NIR_PE, new$NIRspcARmovav)
CUBpred_NIR_PE2 = predict(CUB_NIR_PE2, new$NIRspcA)
CUBpred_NIR_PE3 = predict(CUB_NIR_PE3, new$NIRspcAmovav)
CUBpred_NIR_PE4 = predict(CUB_NIR_PE4, raw_new$NIRspcA)

PE_NIR_obs = c(rep(list(new$PE_mg), 3), list(raw_new$PE_mg))
PE_NIR_preds = list(CUBpred_NIR_PE, CUBpred_NIR_PE2,
                    CUBpred_NIR_PE3, CUBpred_NIR_PE4)

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

require(prospectr)
require(dplyr)
require(Cubist)
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


# fitting Cubist
CUB_mod_mass = cubist(rawC$spcA, rawC$MASS_mg)

CUB_mod_mass2 = cubist(datC$spcA, datC$MASS_mg)

CUB_mod_mass3 = cubist(datC$spcAmovav, datC$MASS_mg)

CUB_mod_mass4 = cubist(datC$spcARmovav, datC$MASS_mg)

# check models dimensionality (mosty important wavelengths)
sum(varImp(CUB_mod_mass) > 0)
sum(varImp(CUB_mod_mass2) > 0)
sum(varImp(CUB_mod_mass3) > 0)
sum(varImp(CUB_mod_mass4) > 0)


# predictions
## calibration
CUBpredC = predict(CUB_mod_mass, rawC$spcA)
CUBpredC2 = predict(CUB_mod_mass2, datC$spcA)
CUBpredC3 = predict(CUB_mod_mass3, datC$spcAmovav)
CUBpredC4 = predict(CUB_mod_mass4, datC$spcARmovav)

## validation
CUBpredV = predict(CUB_mod_mass, rawV$spcA)
CUBpredV2 = predict(CUB_mod_mass2, datV$spcA)
CUBpredV3 = predict(CUB_mod_mass3, datV$spcAmovav)
CUBpredV4 = predict(CUB_mod_mass4, datV$spcARmovav)


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


################################################################################
#                            K-FOLD CROSS-VALIDATION                           #
################################################################################
cv_cubist_model = function(data, spc_matrix, committees = 1, neighbors = 0,
                           nfolds = 4, seed = 999) {
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

cv_CUBISTmod = cv_cubist_model(pristine_raw, pristine_raw$spcA)
cv_CUBISTmod2 = cv_cubist_model(pristine, pristine$spcA)
cv_CUBISTmod3 = cv_cubist_model(pristine, pristine$spcAmovav)
cv_CUBISTmod4 = cv_cubist_model(pristine, pristine$spcARmovav)

cv_cubist_results = list(cv_CUBISTmod, cv_CUBISTmod2, cv_CUBISTmod3, cv_CUBISTmod4)

for (i in seq_along(cv_cubist_results)) {
  cat(paste0("\n======= CUBIST MODEL ", i, " (", cv_cubist_results[[i]]$nfolds_used, "-fold CV) =======\n"))
  cat(sprintf("ME         : %.4f\n", cv_cubist_results[[i]]$ME))
  cat(sprintf("RMSE       : %.4f\n", cv_cubist_results[[i]]$RMSE))
  cat(sprintf("R²         : %.4f\n", cv_cubist_results[[i]]$R2))
}


################################################################################
#            FULL INTERNAL DATASET TRAINING/EXTERNAL VALIDATION                #
################################################################################
CUB_mod_FULL = cubist(pristine_raw$spcA, pristine_raw$MASS_mg,
                      committees = 1, neighbors = 0)

CUB_mod_FULL2 = cubist(pristine$spcA, pristine$MASS_mg,
                       committees = 1, neighbors = 0)

CUB_mod_FULL3 = cubist(pristine$spcAmovav, pristine$MASS_mg,
                       committees = 1, neighbors = 0)

CUB_mod_FULL4 = cubist(pristine$spcARmovav, pristine$MASS_mg,
                       committees = 1, neighbors = 0)

## external validation
CUBpred_FULL = predict(CUB_mod_FULL, colored_raw$spcA)
CUBpred_FULL2 = predict(CUB_mod_FULL2, colored$spcA) 
CUBpred_FULL3 = predict(CUB_mod_FULL3, colored$spcAmovav)
CUBpred_FULL4 = predict(CUB_mod_FULL4, colored$spcARmovav)

full_obs = c(list(colored_raw$MASS_mg), rep(list(colored$MASS_mg), 3))
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
#                            TUNED EXTERNAL VALIDATION                         #
################################################################################
set.seed(22)

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


committees = c(1, 10, 25, 50, 75, 100)
neighbors = c(0, 1, 2, 3, 5, 7, 9) # turns out this is used only in prediction
                                   # not in training...

spectral_preproc_list_ext = list(
  M1 = list(train_x = pristine_raw$spcA, train_y = pristine_raw$MASS_mg,
            new_x = colored_raw1$spcA, new_y = colored_raw1$MASS_mg),
  M2 = list(train_x = pristine$spcA, train_y = pristine$MASS_mg,
            new_x = colored_raw1$spcA, new_y = colored_raw1$MASS_mg),
  M3 = list(train_x = pristine$spcAmovav, train_y = pristine$MASS_mg,
            new_x = colored_raw1$spcAmovav, new_y = colored_raw1$MASS_mg),
  M4 = list(train_x = pristine$spcARmovav, train_y = pristine$MASS_mg,
            new_x = colored_raw1$spcARmovav, new_y = colored_raw1$MASS_mg)
)


# hyperparameter optimization test ============================================#
t0 = Sys.time()

ext_results = list()

counter = 1

for (model_name in names(spectral_preproc_list_ext)) {
  
  cat("\n======= External tuning for", model_name, "=======\n")
  
  train_x = spectral_preproc_list_ext[[model_name]]$train_x
  train_y = spectral_preproc_list_ext[[model_name]]$train_y
  new_x   = spectral_preproc_list_ext[[model_name]]$new_x
  new_y   = spectral_preproc_list_ext[[model_name]]$new_y
  
  for (committees_val in committees) {
    
    # train ONCE per committees value
    cub_model = cubist(train_x, train_y, committees = committees_val)
    cat("  Fitted committees =", committees_val, "\n")
    
    # then just predict 7 times with different neighbors
    for (neighbors_val in neighbors) {
      
      pred_ext = predict(cub_model, newdata = new_x, neighbors = neighbors_val)
      
      ext_results[[counter]] = data.frame(
        Model = model_name,
        committees = committees_val,
        neighbors = neighbors_val,
        ME = ME(new_y, pred_ext),
        RMSE = RMSE(new_y, pred_ext),
        R2 = R2(new_y, pred_ext)
      )
      counter = counter + 1
    }
  }
}
ext_results = do.call(rbind, ext_results)


t1 = Sys.time()
cat("Training time:", format(difftime(t1, t0, units = "mins")), "\n")


## check optima parameteres
best_param =
  ext_results |> 
  group_by(Model) |> 
  filter(RMSE == min(RMSE))
#==============================================================================#

# tuned models
# M1
CUB_ext_tuned = cubist(raw$spcA,
                       raw$MASS_mg,
                       committees = 1) # best_param$committes[1]

CUBpred_ext = predict(CUB_ext_tuned, newdata = raw_new2$spcA, neighbors = 5) # best_param$neighbors[1]
cat(paste0("\n======= CUBIST MODEL 1 (fine-tuned) =======\n",
           sprintf("ME   : %.2f\n", ME(raw_new2$MASS_mg, CUBpred_ext)),
           sprintf("RMSE : %.2f\n", RMSE(raw_new2$MASS_mg, CUBpred_ext)),
           sprintf("R²   : %.2f\n", R2(raw_new2$MASS_mg, CUBpred_ext))
           ))

# M2
CUB_ext_tuned2 = cubist(data$spcA,
                        data$MASS_mg,
                        committees = 100)

CUBpred_ext2 = predict(CUB_ext_tuned2, newdata = new2$spcA, neighbors = 1)
cat(paste0("\n======= CUBIST MODEL 2 (fine-tuned) =======\n",
           sprintf("ME   : %.2f\n", ME(new2$MASS_mg, CUBpred_ext2)),
           sprintf("RMSE : %.2f\n", RMSE(new2$MASS_mg, CUBpred_ext2)),
           sprintf("R²   : %.2f\n", R2(new2$MASS_mg, CUBpred_ext2))
           ))

# M3
CUB_ext_tuned3 = cubist(data$spcAmovav,
                        data$MASS_mg,
                        committees = 50)

CUBpred_ext3 = predict(CUB_ext_tuned3, newdata = new2$spcAmovav, neighbors = 9)
cat(paste0("\n======= CUBIST MODEL 3 (fine-tuned) =======\n",
           sprintf("ME   : %.2f\n", ME(new2$MASS_mg, CUBpred_ext3)),
           sprintf("RMSE : %.2f\n", RMSE(new2$MASS_mg, CUBpred_ext3)),
           sprintf("R²   : %.2f\n", R2(new2$MASS_mg, CUBpred_ext3))
           ))

# M4
CUB_ext_tuned4 = cubist(data$spcARmovav,
                        data$MASS_mg,
                        committees = 50)

CUBpred_ext4 = predict(CUB_ext_tuned4, newdata = new2$spcARmovav, neighbors = 7)
cat(paste0("\n======= CUBIST MODEL 4 (fine-tuned) =======\n",
           sprintf("ME   : %.2f\n", ME(new2$MASS_mg, CUBpred_ext4)),
           sprintf("RMSE : %.2f\n", RMSE(new2$MASS_mg, CUBpred_ext4)),
           sprintf("R²   : %.2f\n", R2(new2$MASS_mg, CUBpred_ext4))
           ))


#
#
# LAST EDITING 2026-10-02
#        17h21
#
#



# 
# # residual visualization
# residLEAK = CUBpred_full - new$MASS_mg
# residLEAK2 = CUBpred_full2 - new$MASS_mg
# residLEAK3 = CUBpred_full3 - new$MASS_mg
# residLEAK4 = CUBpred_full4 - raw_new$MASS_mg
# 
# residLEAKlog = log(CUBpred_full) - log(new$MASS_mg)
# residLEAKlog2 = log(CUBpred_full2) - log(new$MASS_mg)
# residLEAKlog3 = log(CUBpred_full3) - log(new$MASS_mg)
# residLEAKlog4 = log(CUBpred_full4) - log(raw_new$MASS_mg)
# 
# # ggplot(new, aes(x = factor(MASS_mg), y = residLEAK2)) +
# #     geom_boxplot() +
# #     labs(title = "MODEL 2",
# #          x = "Mass (mg)",
# #          y = "Residuals")
# # 
# # ggplot(new, aes(x = factor(MASS_mg), y = residLEAKlog4)) +
# #     geom_boxplot() +
# #     labs(title = "MODEL 2",
# #          x = "Mass (mg)",
# #          y = "Relative residuals")
# 
# # data frame for each model, ensuring the correct 'Observed_Mass' is used
# df_m1 = data.frame(
#   Model = "M1 (Raw data)",
#   Observed_Mass = raw_new$MASS_mg, # M1 uses raw_new
#   Residual = residLEAK4  # M1 = residLEAK4
# )
# 
# df_m2 = data.frame(
#   Model = "M2 (Minimal preprocessing)",
#   Observed_Mass = new$MASS_mg, # M2 uses new
#   Residual = residLEAK2
# )
# 
# df_m3 = data.frame(
#   Model = "M3 (Intermediate preprocessing)",
#   Observed_Mass = new$MASS_mg, # M3 uses new
#   Residual = residLEAK3
# )
# 
# df_m4 = data.frame(
#   Model = "M4 (Full preprocessing)",
#   Observed_Mass = raw_new$MASS_mg, # M4 uses raw_new
#   Residual = residLEAK  # M4 = residLEAK
# )
# 
# # combine data frames
# require(dplyr)
# plot_data = bind_rows(df_m1, df_m2, df_m3, df_m4)
# 
# # model factor in correct order of complexity
# model_order = c("M1 (Raw data)", "M2 (Minimal preprocessing)",
#                 "M3 (Intermediate preprocessing)", "M4 (Full preprocessing)")
# plot_data$Model = factor(plot_data$Model, levels = model_order)
# 
# 
# # combined boxplot
# require(ggplot2)
# res = ggplot(plot_data, aes(x = factor(Observed_Mass), y = Residual)) +
#   geom_boxplot(outlier.size = 0.8) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
#   facet_wrap(~ Model, ncol = 2) +
#   labs(
#     x = NULL,
#     y = "Residuals (Predicted - Observed Mass)",
#     title = NULL
#   ) +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank())
# 
# #===============================================================================
# 
# # data frames for relative residuals (log scale)
# df_m1_log = data.frame(
#   Model = "M1 (Raw data)",
#   Observed_Mass = raw_new$MASS_mg,
#   Relative_Residual = residLEAKlog4  # M1 = residLEAKlog4
# )
# 
# df_m2_log = data.frame(
#   Model = "M2 (Minimal preprocessing)", 
#   Observed_Mass = new$MASS_mg,
#   Relative_Residual = residLEAKlog2
# )
# 
# df_m3_log = data.frame(
#   Model = "M3 (Intermediate preprocessing)",
#   Observed_Mass = new$MASS_mg, 
#   Relative_Residual = residLEAKlog3
# )
# 
# df_m4_log = data.frame(
#   Model = "M4 (Full preprocessing)",
#   Observed_Mass = raw_new$MASS_mg,
#   Relative_Residual = residLEAKlog  # M4 = residLEAKlog
# )
# 
# # combine and process
# plot_data_log = bind_rows(df_m1_log, df_m2_log, df_m3_log, df_m4_log)
# 
# 
# # set factor order
# model_order = c("M1 (Raw data)", "M2 (Minimal preprocessing)",
#                 "M3 (Intermediate preprocessing)", "M4 (Full preprocessing)")
# 
# # relative residuals plot
# res2 = ggplot(plot_data_log, aes(x = factor(Observed_Mass), y = Relative_Residual)) +
#   geom_boxplot(outlier.size = 0.8) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
#   facet_wrap(~ Model, ncol = 2) +
#   labs(
#     x = "Mass (mg)",
#     y = "Relative Residuals [log(Predicted) - log(Observed)]",
#     title = NULL)
# 
# require(patchwork)
# 
# res / res2 + plot_annotation(tag_levels = 'a',
#                              tag_prefix = '(',
#                              tag_suffix = ')')

require(dplyr)
set.seed(23)

new$strata = interaction(new$MASS_mg,
                         new$SIZE_CODE)

new1 = new |> 
  group_by(strata) |> 
  sample_frac(0.5)

new2 = new |> 
  filter(!(SAMPLE_ID %in% new1$SAMPLE_ID))

raw_new$strata = interaction(raw_new$MASS_mg,
                             raw_new$SIZE_CODE)


raw_new1 = raw_new |> 
  group_by(strata) |> 
  sample_frac(0.5)


raw_new2 = raw_new |> 
  filter(!(SAMPLE_ID %in% raw_new1$SAMPLE_ID))

# prepare external data
raw_new1_M1 = data.frame(plastMass = raw_new1$MASS_mg, raw_new1$spcA)
colnames(raw_new1_M1) = c("plastMass", paste0("spec.", colnames(raw_new1$spcA)))
raw_new2_M1 = data.frame(plastMass = raw_new2$MASS_mg, raw_new2$spcA)
colnames(raw_new2_M1) = c("plastMass", paste0("spec.", colnames(raw_new2$spcA)))

new1_M2 = data.frame(plastMass = new1$MASS_mg, new1$spcA)
colnames(new1_M2) = c("plastMass", paste0("spec.", colnames(new1$spcA)))
new2_M2 = data.frame(plastMass = new2$MASS_mg, new2$spcA)
colnames(new2_M2) = c("plastMass", paste0("spec.", colnames(new2$spcA)))

new1_M3 = data.frame(plastMass = new1$MASS_mg, new1$spcAmovav)
colnames(new1_M3) = c("plastMass", paste0("spec.", colnames(new1$spcAmovav)))
new2_M3 = data.frame(plastMass = new2$MASS_mg, new2$spcAmovav)
colnames(new2_M3) = c("plastMass", paste0("spec.", colnames(new2$spcAmovav)))

new1_M4 = data.frame(plastMass = new1$MASS_mg, new1$spcARmovav)
colnames(new1_M4) = c("plastMass", paste0("spec.", colnames(new1$spcARmovav)))
new2_M4 = data.frame(plastMass = new2$MASS_mg, new2$spcARmovav)
colnames(new2_M4) = c("plastMass", paste0("spec.", colnames(new2$spcARmovav)))

param_grid_ext = expand.grid(
  ntree = c(100, 150, 200, 250, 500),
  mtry = c(5, 10, 20, 30, 40, 50, 55, 60, 70, 80, 90, 95, 100)
)

spectral_preproc_list_ext = list(
  M1 = list(train = rawFULL, new = raw_new1_M1),
  M2 = list(train = dataFULL2, new = new1_M2),
  M3 = list(train = dataFULL3, new = new1_M3),
  M4 = list(train = dataFULL, new = new1_M4)
)

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

# hyperparameter selection =====================================================
ext_results |> 
  group_by(Model) |> 
  filter(RMSE == min(RMSE))

ext_results |>
  group_by(Model) |>
  filter(abs(ME) == min(abs(ME)))
#===============================================================================

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
                            data = dataFULL,
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
residTUNED = RFpred_tuned - raw_new2$MASS_mg
residTUNED2 = RFpred_tuned2 - new2$MASS_mg
residTUNED3 = RFpred_tuned3 - new2$MASS_mg
residTUNED4 = RFpred_tuned4 - new2$MASS_mg

residTUNEDlog = log(RFpred_tuned) - log(raw_new2$MASS_mg)
residTUNEDlog2 = log(RFpred_tuned2) - log(new2$MASS_mg)
residTUNEDlog3 = log(RFpred_tuned3) - log(new2$MASS_mg)
residTUNEDlog4 = log(RFpred_tuned4) - log(new2$MASS_mg)



df_m1 = data.frame(
  Model = "M1 (Raw data)",
  Observed_Mass = raw_new2$MASS_mg, # M1 uses raw_new
  Residual = residTUNED
)

df_m2 = data.frame(
  Model = "M2 (Minimal preprocessing)",
  Observed_Mass = new2$MASS_mg, # M2 uses new
  Residual = residTUNED2
)

df_m3 = data.frame(
  Model = "M3 (Intermediate preprocessing)",
  Observed_Mass = new2$MASS_mg, # M3 uses new
  Residual = residTUNED3
)

df_m4 = data.frame(
  Model = "M4 (Full preprocessing)",
  Observed_Mass = new2$MASS_mg, # M4 uses new
  Residual = residTUNED4
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
    title = NULL) +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 11),
        axis.title.y = element_text(size = 13))

#===============================================================================

# data frames for relative residuals (log scale)
df_m1_log = data.frame(
  Model = "M1 (Raw data)",
  Observed_Mass = raw_new2$MASS_mg,
  Relative_Residual = residTUNEDlog  # M1 = residTUNEDlog
)

df_m2_log = data.frame(
  Model = "M2 (Minimal preprocessing)", 
  Observed_Mass = new2$MASS_mg,
  Relative_Residual = residTUNEDlog2
)

df_m3_log = data.frame(
  Model = "M3 (Intermediate preprocessing)",
  Observed_Mass = new2$MASS_mg, 
  Relative_Residual = residTUNEDlog3
)

df_m4_log = data.frame(
  Model = "M4 (Full preprocessing)",
  Observed_Mass = raw_new2$MASS_mg,
  Relative_Residual = residTUNEDlog4  # M4 = residTUNEDlog4
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

require(patchwork)

res / res2 + plot_annotation(tag_levels = 'a',
                             tag_prefix = '(',
                             tag_suffix = ')') &
  theme(plot.tag = element_text(size = 22))

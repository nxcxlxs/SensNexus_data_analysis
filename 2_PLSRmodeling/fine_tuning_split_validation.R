set.seed(42)

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



# NAIVE TRY ====================================================================
evaluate_ncomp = function(model, newdata, obs, grid){
  
  out = lapply(grid, function(nc){
    
    pred = predict(model, ncomp = nc, newdata = newdata)
    pred = as.vector(pred)
    
    c(
      ncomp = nc,
      ME    = ME(obs, pred),
      RMSE  = RMSE(obs, pred),
      R2    = R2(obs, pred)
    )
  })
  
  as.data.frame(do.call(rbind, out))
}

grid = 1:30 # gotta make sure the models trained on full dataset have ncomp = 50

res_M1_ext = evaluate_ncomp(model   = PLSR_mod_mass4,
                            newdata = raw_new1$spcA,
                            obs     = raw_new1$MASS_mg,
                            grid    = grid)
res_M1_ext = res_M1_ext |> 
  mutate(Model = rep("M1", 30), .before = ncomp)



res_M2_ext = evaluate_ncomp(model   = PLSR_mod_mass2,
                            newdata = new1$spcA,
                            obs     = new1$MASS_mg,
                            grid    = grid)
res_M2_ext = res_M2_ext |> 
  mutate(Model = rep("M2", 30), .before = ncomp)

res_M3_ext = evaluate_ncomp(model   = PLSR_mod_mass3,
                            newdata = new1$spcAmovav,
                            obs     = new1$MASS_mg,
                            grid    = grid)
res_M3_ext = res_M3_ext |> 
  mutate(Model = rep("M3", 30), .before = ncomp)
    
res_M4_ext = evaluate_ncomp(model   = PLSR_mod_mass,
                            newdata = new1$spcARmovav,
                            obs     = new1$MASS_mg,
                            grid    = grid)
res_M4_ext = res_M4_ext |> 
  mutate(Model = rep("M4", 30), .before = ncomp)


results = rbind(res_M1_ext, res_M2_ext,
                res_M3_ext, res_M4_ext,
                make.row.names = F)


results = results |> 
  mutate(R2 = round(R2, 8))


results |> 
  group_by(Model) |> 
  filter(RMSE == min(RMSE))


# FINE TUNING ==================================================================
# M1                                       -> pay attention to the ncomp values
PLSR_pred_tuned4 = predict(PLSR_mod_leak4, ncomp = 10, newdata = raw_new2$spcA)
cat(paste0("\n======= PLSR MODEL 1 (fine-tuned) =======\n",
           sprintf("ME   : %.2f\n", ME(raw_new2$MASS_mg, PLSR_pred_tuned4)),
           sprintf("RMSE : %.2f\n", RMSE(raw_new2$MASS_mg, PLSR_pred_tuned4)),
           sprintf("R²   : %.2f\n", R2(raw_new2$MASS_mg, PLSR_pred_tuned4))
))

# M2
PLSR_pred_tuned2 = predict(PLSR_mod_leak2, ncomp = 12, newdata = new2$spcA)
cat(paste0("\n======= PLSR MODEL 2 (fine-tuned) =======\n",
           sprintf("ME   : %.2f\n", ME(new2$MASS_mg, PLSR_pred_tuned2)),
           sprintf("RMSE : %.2f\n", RMSE(new2$MASS_mg, PLSR_pred_tuned2)),
           sprintf("R²   : %.2f\n", R2(new2$MASS_mg, PLSR_pred_tuned2))
))

# M3
PLSR_pred_tuned3 = predict(PLSR_mod_leak3, ncomp = 6, newdata = new2$spcAmovav)
cat(paste0("\n======= PLSR MODEL 3 (fine-tuned) =======\n",
           sprintf("ME   : %.2f\n", ME(new2$MASS_mg, PLSR_pred_tuned3)),
           sprintf("RMSE : %.2f\n", RMSE(new2$MASS_mg, PLSR_pred_tuned3)),
           sprintf("R²   : %.2f\n", R2(new2$MASS_mg, PLSR_pred_tuned3))
))
                   
# M4
PLSR_pred_tuned = predict(PLSR_mod_leak, ncomp = 6, newdata = new2$spcARmovav)
cat(paste0("\n======= PLSR MODEL 4 (fine-tuned) =======\n",
           sprintf("ME   : %.2f\n", ME(new2$MASS_mg, PLSR_pred_tuned)),
           sprintf("RMSE : %.2f\n", RMSE(new2$MASS_mg, PLSR_pred_tuned)),
           sprintf("R²   : %.2f\n", R2(new2$MASS_mg, PLSR_pred_tuned))
))


## this pretty much just inflate the values...
## although, the M4 coefficient of determination looks better, its ME is still 
## worse than M2's; so the ranking of the models stays the same...

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
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 11),
        axis.title.y = element_text(size = 13))

#===============================================================================

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

require(patchwork)

res / res2 + plot_annotation(tag_levels = 'a',
                             tag_prefix = '(',
                             tag_suffix = ')') &
  theme(plot.tag = element_text(size = 22))

        
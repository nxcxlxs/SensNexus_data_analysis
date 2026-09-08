set.seed(22)

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

#===============================================================================

spectral_preproc_list_ext = list(
  M1 = list(train_x = raw$spcA, train_y = raw$MASS_mg,
            new_x = raw_new1$spcA, new_y = raw_new1$MASS_mg),
  M2 = list(train_x = data$spcA, train_y = data$MASS_mg,
            new_x = new1$spcA, new_y = new1$MASS_mg),
  M3 = list(train_x = data$spcAmovav, train_y = data$MASS_mg,
            new_x = new1$spcAmovav, new_y = new1$MASS_mg),
  M4 = list(train_x = data$spcARmovav, train_y = data$MASS_mg,
            new_x = new1$spcARmovav, new_y = new1$MASS_mg)
)

# hyperparameters
committees = c(1, 10, 25, 50, 75, 100)
neighbors = c(0, 1, 2, 3, 5, 7, 9) # turns out this is used only in prediction
                                   # not in training...

#===============================================================================
require(Cubist)

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
cat("Training time:", round(t1 - t0, 2), "minutes\n")


#===============================================================================
require(dplyr)

ext_results |> 
  group_by(Model) |> 
  filter(RMSE == min(RMSE))

ext_results |>
  group_by(Model) |>
  filter(abs(ME) == min(abs(ME)))

#===============================================================================
# M1
CUB_ext_tuned = cubist(raw$spcA,
                       raw$MASS_mg,
                       committees = 1)

CUBpred_ext = predict(CUB_ext_tuned, newdata = raw_new2$spcA, neighbors = 5)
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

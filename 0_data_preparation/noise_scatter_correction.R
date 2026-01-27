# load spliced corrected data
data = readRDS("./preprocessed_data/splice_corrected.rds")


# average the replicates (?)
require(dplyr)

avg_data = bind_cols(
        data[, c("POLYMER",
                 "SIZE_INTERVALS_mm",
                 "SIZE_CODE",
                 "MASS_mg")],
        data$spc
)

avg_data = avg_data |> 
        group_by(POLYMER, SIZE_INTERVALS_mm, SIZE_CODE, MASS_mg) |> 
        summarise(across(everything(), mean,
                         .names = "{.col}"),
                         .groups = "drop")


spc = avg_data[, 5:ncol(avg_data)]
avg_data = avg_data[, 1:4]
avg_data$spc = as.data.frame(spc)

# visualize side by side
par(mfrow = c(1, 2))

# vizualize spectra
matplot(x = colnames(avg_data$spc), y = t(avg_data$spc),
        main = "Raw averaged sample spectra",
        xlab = "Wavelength (nm)",
        ylab = "Reflectance",
        type = "l",
        lty = 1,
        col = rgb(red = 1, green = 0.5, blue = 0.5, alpha = 0.3))

################################################################################
#                             MOVING AVERAGE WINDOW                            #
################################################################################
# noise removal
require(prospectr)

swindow_ma = 21 # bands

#  wavelengths lost (begining and end) (window size - 1) / 2
avg_data$spc_ma = movav(avg_data$spc, swindow_ma)

# plot de-noised spectra
matplot(x = colnames(avg_data$spc_ma), y = t(avg_data$spc_ma),
        main = "De-noised averaged sample spectra (moving average)",
        xlab = "Wavelength (nm)",
        ylab = "Reflectance",
        type = "l",
        lty = 1,
        col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.3))

################################################################################
#                             SAVITSKY-GOLAY FILTER                            #
################################################################################
swindow_sg = 21

poly_sg = 2

data$spc_sg = savitzkyGolay(data$spc,
                            m = 0,
                            p = poly_sg,
                            w = swindow_sg)

# plot de-noised spectra 
matplot(x = colnames(data$spc_sg), y = t(data$spc_sg),
        main = "De-noised sample spectra (Savitzky-Golay filter_2)",
        xlab = "Wavelength (nm)",
        ylab = "Reflectance",
        type = "l",
        lty = 1,
        col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.3))

wav2 = as.numeric(colnames(data$spc_sg)) # use reflectance!

snr_region2 = which(wav2 > 1000 & wav2 < 1100)

snr2 = rowMeans(data$spc_sg[, snr_region2]) / apply(data$spc_sg[, snr_region2], 1, sd)
hist(snr2, breaks = 30, main = "SNR (1000–1100 nm) after")

identical(snr, snr2) # different?

sum(snr2>snr) # how many improved?

mean((snr2 - snr) / snr) * 100  # mean % improvement across spectra

hist(snr2 - snr, breaks = 30, main = "SNR Improvement (SG Filter_2)", xlab = "ΔSNR")

################################################################################
#                      MULTIPLICATICE SCATTER CORRECTION                       #
################################################################################
require(pls)

color = rgb(0, 0, 0, alpha = 0.3 )

avg_data$spc_msc = msc(avg_data$spc)

matplot(x = colnames(avg_data$spc_msc), y = t(avg_data$spc_msc),
        main = "Scatter corrected averaged sample spectra (MSC)",
        xlab = "Wavelength (nm)",
        ylab = "msc(Reflectance)",
        type = "l",
        lty = 1,
        col = color)

################################################################################
#                              STANDARD NORMAL VARIATE                         #
################################################################################
avg_data$spc_snv = standardNormalVariate(avg_data$spc)

matplot(x = colnames(avg_data$spc_snv), y = t(avg_data$spc_snv),
        main = "Scatter corrected averaged sample spectra (SNV)",
        xlab = "Wavelength (nm)",
        ylab = "snv(Reflectance)",
        type = "l",
        lty = 1,
        col = color)

################################################################################
#                               CENTERING AND SCALING                          #
################################################################################
centering = TRUE
scaling = FALSE

# centering just the average
avg_data$spc_cnt = scale(avg_data$spc, center = centering, scale = scaling)

matplot(x = colnames(avg_data$spc_cnt), y = t(avg_data$spc_cnt),
        main = "Centered averaged sample spectra",
        xlab = "Wavelength (nm)",
        ylab = "Reflectance",
        type = "l",
        lty = 1,
        col = color)

# centering the SNV
avg_data$spc_snv_cnt = scale(avg_data$spc_snv, center = centering, scale = scaling) 

matplot(x = colnames(avg_data$spc_snv_cnt), y = t(avg_data$spc_snv_cnt),
        main = "Centered scatter corrected averaged sample spectra (SNV)",
        xlab = "Wavelength (nm)",
        ylab = "snv(Reflectance)",
        type = "l",
        lty = 1,
        col = color)


################################################################################
#                                 CONTINUUM REMOVAL                            #
################################################################################
tp = "R" # reflectance, or "A" for absorbance

avg_data$spc_cr = continuumRemoval(avg_data$spc,
                                   wav = as.numeric(colnames(avg_data$spc)),
                                   type = tp)

matplot(x = colnames(avg_data$spc_cr), y = t(avg_data$spc_cr),
        main = "Scaled  averaged sample spectra",
        xlab = "Wavelength (nm)",
        ylab = "Reflectance",
        type = "l",
        lty = 1,
        col = color)


#==============================================================================#
# TURNS OUT THAT I DON'T REALLY THINK THE DATA I HAVE TO WORK NEED SO MUCH PRE-#
#-PROCESSING. OF COURSE I JUST EVALUATED IT VISUALLY OVER THE SPECTRA AVERAGE. #
# ANYWAY, I'LL JUST APPLY SAVITSKY-GOLAY FILTER AND FOLLOW ON TO MODELLING     #
#==============================================================================#

deNoised_data = data

swindow_sg = 21
poly_sg = 2

deNoised_data$spc = savitzkyGolay(data$spc,
                                m = 0,
                                p = poly_sg,
                                w = swindow_sg)

# plot de-noised spectra 
matplot(x = colnames(deNoised_data$spc), y = t(deNoised_data$spc),
        main = "De-noised full sample spectra (Savitzky-Golay filter)",
        xlab = "Wavelength (nm)",
        ylab = "Reflectance",
        type = "l",
        lty = 1,
        col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.3))


saveRDS(deNoised_data, "./analysis_ready_data/deNoised_data.rds")


paste("Lost", ncol(data$spc) - ncol(deNoised_data$spc),
      "nanometers of wavelength.",
      "Still have", ncol(deNoised_data$spc) * nrow(deNoised_data),
      "data points to work with, though.")


ref = readRDS("./preprocessed_data/ref_corrected.rds")

deNoised_ref = lapply(ref, function(x){
  spec_matrix = as.matrix(x)
  smoothed = savitzkyGolay(spec_matrix,
                           m = 0,
                           p = 2,
                           w = 21)
  as.data.frame(smoothed)
})

saveRDS(deNoised_ref, "./analysis_ready_data/deNoised_ref.rds")


#==============================================================================#
#         PREPROCESS NEW DATA (COLORED RIGID HOUSEHOLD PLASTIC - cRHP)         #
#==============================================================================#

new_data = readRDS("./preprocessed_data/splice_corrected2.rds")

matplot(x = colnames(new_data$spc), y = t(new_data$spc),
        main = "cRHP spectra",
        xlab = "Wavelength (nm)",
        ylab = "Reflectance",
        type = "l",
        lty = 1,
        col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.3))

deNoised_newData = new_data

deNoised_newData$spc = savitzkyGolay(new_data$spc,
                                  m = 0,
                                  p = poly_sg,
                                  w = swindow_sg)

par(mfrow = c(1, 2))

matplot(x = colnames(deNoised_newData$spc), y = t(deNoised_newData$spc),
        main = "De-noised cRHP spectra",
        xlab = "Wavelength (nm)",
        ylab = "Reflectance",
        type = "l",
        lty = 1,
        col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.3))

wav = as.numeric(colnames(new_data$spc))

snr_region = which(wav > 1000 & wav < 1100)

snr = rowMeans(new_data$spc[, snr_region]) / apply(new_data$spc[, snr_region], 1, sd) # raw data
snr2 = rowMeans(deNoised_newData$spc[, snr_region]) / apply(deNoised_newData$spc[, snr_region], 1, sd) # post correction
 
hist(snr, breaks = 30, main = "SNR (1000–1100 nm) after")
hist(snr2, breaks = 30, main = "SNR (1000–1100 nm) after")

identical(snr, snr2) # different?

sum(snr2>snr) # how many improved?

mean((snr2 - snr) / snr) * 100  # mean % improvement across spectra

par(mfrow = c(1, 1))

hist(snr2 - snr, breaks = 30, main = "SNR Improvement (SG Filter_2)", xlab = "ΔSNR")

saveRDS(deNoised_newData, "./analysis_ready_data/deNoised_newData.rds")

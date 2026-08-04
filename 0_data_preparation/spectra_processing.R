require(asdreader)
require(prospectr)


################################################################################
#                               SPLICE CORRECTION                              #
################################################################################

# load soil dataset
pristine = readRDS("../raw_spectra/raw_pristine.rds")
colored = readRDS("../raw_spectra/raw_colored.rds")

# visual inspection
matplot(colnames(pristine$spc),
        t(pristine$spc),
        main = "Raw",
        xlab = 'Wavelength',
        ylab = 'Reflectance',
        type = 'l',
        lwd = 2,
        lty = 1,
        col = rgb(red = 0.5,
                  green = 0.5,
                  blue = 0.5,
                  alpha = 0.3))


# get spectral files metadata (same for `pristine` and `colored`)
raw = "../raw_spectra/FIELDSPEC4_data/spec00000.asd"

str(get_metadata(raw))


# check splice artifacts
abline(v = c(980, 1020, 1780, 1820),
       col = adjustcolor("red", alpha.f = 0.5),
       lty = 2,
       lwd = 2)


# correct splices
correct_pristine = spliceCorrection(pristine$spc,
                                    as.numeric(colnames(pristine$spc)),
                                    splice = c(as.numeric(get_metadata(raw)[30]),
                                               as.numeric(get_metadata(raw)[31])))

correct_colored = spliceCorrection(colored$spc,
                                   as.numeric(colnames(colored$spc)),
                                   splice = c(as.numeric(get_metadata(raw)[30]),
                                              as.numeric(get_metadata(raw)[31])))

# side-by-side check
par(mfrow = c(1, 2))

matplot(colnames(pristine$spc),
        t(pristine$spc),
        main = "Raw",
        xlab = 'Wavelength',
        ylab = 'Reflectance',
        type = 'l',
        lwd = 2,
        lty = 1,
        col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.3))

matplot(colnames(correct_pristine$spc),
        t(correct_pristine$spc),
        main = "Corrected",
        xlab = 'Wavelength',
        ylab = 'Reflectance',
        type = 'l',
        lwd = 2,
        lty = 1,
        col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.3))


# save it
pristine$spc = as.data.frame(correct_pristine)
colored$spc = as.data.frame(correct_colored)

saveRDS(pristine, "../preprocessed_data/pristine_splice_corrected.rds")
saveRDS(colored, "../preprocessed_data/colored_splice_corrected.rds")


################################################################################
#                     DENOISE DATA (SAVITSKY-GOLAY FILTER)                     #
################################################################################

# load splice corrected data (for sanity sake)
pristine_corrected = readRDS("../preprocessed_data/pristine_splice_corrected.rds")
colored_corrected = readRDS("../preprocessed_data/colored_splice_corrected.rds")

# set parameters
swindow_sg = 21
poly_sg = 2


pristine_denoised = pristine_corrected
pristine_denoised$spc = savitzkyGolay(pristine_corrected$spc,
                                      m = 0,
                                      p = poly_sg,
                                      w = swindow_sg)

colored_denoised = colored_corrected
colored_denoised$spc = savitzkyGolay(colored_corrected$spc,
                                     m = 0,
                                     p = poly_sg,
                                     w = swindow_sg)

# side-by-side check
par(mfrow = c(1, 2))

matplot(colnames(pristine_corrected$spc),
        t(pristine_corrected$spc),
        main = "Corrected",
        xlab = 'Wavelength',
        ylab = 'Reflectance',
        type = 'l',
        lwd = 2,
        lty = 1,
        col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.3))

matplot(colnames(pristine_denoised$spc),
        t(pristine_denoised$spc),
        main = "Denoised",
        xlab = 'Wavelength',
        ylab = 'Reflectance',
        type = 'l',
        lwd = 2,
        lty = 1,
        col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.3))

paste("Lost", ncol(pristine_corrected$spc) - ncol(pristine_denoised$spc),
      "wavelengths.",
      "Still have", ncol(pristine_denoised$spc) * nrow(pristine_denoised),
      "data points to work with, though.")


# quick signal-to-noise ratio check
wav = as.numeric(colnames(pristine_corrected$spc))

snr_region = which(wav > 1000 & wav < 1100) # relatively stable region

snr = rowMeans(pristine_corrected$spc[, snr_region]) / apply(pristine_corrected$spc[, snr_region], 1, sd) # raw data
snr2 = rowMeans(pristine_denoised$spc[, snr_region]) / apply(pristine_denoised$spc[, snr_region], 1, sd) # post denoising

par(mfrow = c(1, 2))

hist(snr, breaks = 30, main = "SNR (1000–1100 nm) before")
hist(snr2, breaks = 30, main = "SNR (1000–1100 nm) after")

identical(snr, snr2) # different?

sum(snr2>snr) # how many improved?

mean((snr2 - snr) / snr) * 100  # mean % improvement across spectra

## visualize it
par(mfrow = c(1, 2))

matplot(colnames(pristine_corrected$spc[, snr_region]),
        t(pristine_corrected$spc[, snr_region]),
        type = "l",
        lwd = 2,
        lty = 1,
        col = rgb(0.5, 0.5, 0.5, alpha = 0.5))

matplot(colnames(pristine_denoised$spc[, snr_region]),
        t(pristine_denoised$spc[, snr_region]),
        type = "l",
        lwd = 2,
        lty = 1,
        col = rgb(0.5, 0.5, 0.5, alpha = 0.5))


# save it
saveRDS(pristine_denoised, "../preprocessed_data/pristine_denoised.rds")
saveRDS(colored_denoised, "../preprocessed_data/colored_denoised.rds")

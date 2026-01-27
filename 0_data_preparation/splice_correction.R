# load soil dataset
data = readRDS("./preprocessed_data/datsoil.rds")

# splice removal
require(asdreader)
require(prospectr)

## check raw spectra metadatata get splice location
raw = "./raw_spectra/spec00000.asd"

correct_data = spliceCorrection(data$spc, as.numeric(colnames(data$spc)),
                                splice = c(as.numeric(get_metadata(raw)[30]),
                                           as.numeric(get_metadata(raw)[31])))

# visualize data
par(mfrow = c(1, 2))

matplot(colnames(data$spc),
        t(data$spc),
        main = "Raw",
        xlab = 'Wavelength',
        ylab = 'Reflectance',
        type = 'l',
        lwd = 2,
        lty = 1,
        col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.3))

matplot(colnames(correct_data),
        t(correct_data),,
        main = "Corrected",
        xlab = 'Wavelength',
        ylab = 'Reflectance',
        type = 'l',
        lwd = 2,
        lty = 1,
        col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.3))


data$spc = as.data.frame(correct_data)

saveRDS(data, "./preprocessed_data/splice_corrected.rds")


# load references dataset
ref = readRDS("./preprocessed_data/ref_spec.rds")

ref_corrected = lapply(ref, function(data){
  corrected = spliceCorrection(data,
                               as.numeric(colnames(data)),
                               splice = c(as.numeric(get_metadata(raw)[30]),
                                          as.numeric(get_metadata(raw)[31])))
  as.data.frame(corrected)
})

par(mfrow = c(1, 2))

matplot(colnames(ref$soil),
        t(ref$soil),
        main = "Raw",
        xlab = 'Wavelength',
        ylab = 'Reflectance',
        type = 'l',
        lwd = 2,
        lty = 1,
        col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.3))

matplot(colnames(ref_corrected$soil),
        t(ref_corrected$soil),
        main = "Corrected",
        xlab = 'Wavelength',
        ylab = 'Reflectance',
        type = 'l',
        lwd = 2,
        lty = 1,
        col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.3))

saveRDS(ref_corrected, "./preprocessed_data/ref_corrected.rds")

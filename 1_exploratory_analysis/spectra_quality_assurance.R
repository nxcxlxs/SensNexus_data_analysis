require(resemble)
require(scatterplot3d)

# load dataset
data = readRDS("./preprocessed_data/splice_corrected.rds")

#==============================================================================#
#                 1. OUTLIER DETECTION BASED ON SPECTRAL SHAPE                .#
#==============================================================================#
# convert reflectance to absorbance
data$spcA = log(1/data$spc) # it stabilizes variance
                            # it avoid compressing informative low-reflectance regions
                            #  it emphasizes shape over absolute magnitude

# reduce dimensionality (approach for `p >> n` problem)
max_expl_var = 0.99

pcspectra = pc_projection(Xr = as.matrix(data$spcA),
                          pc_selection = list("cumvar", max_expl_var),
                          method = "pca",
                          center = TRUE,
                          scale = FALSE)

# average of PC scores
pcspectraCentre = colMeans(pcspectra$scores)
pcspectraCentre = t(as.matrix(pcspectraCentre))


# Mahalanobis Distance between scores centre and spectra scores
mahD = f_diss(Xr = pcspectra$scores,
               Xu = pcspectraCentre,
               diss_method = "mahalanobis",
               center = FALSE, scale = FALSE)

## index of the spectra against the Mahalanobis distance
plot(mahD,
     pch = 16,
     col = rgb(red = 0, green = 0.4, blue = 0.8, alpha = 0.5),
     ylab = "Mahalanobis distance")

## threshold dinamically adjusted based on the dataset (99th percentile)
threshold = quantile(mahD, 0.99) 
abline(h = threshold, col = "red")

### obtain the indices of the outliers
indxOutM = which(mahD > threshold)

### how many potential outliers?
length(indxOutM)

### visualize outliers
sct3d = scatterplot3d(pcspectra$scores[,1],
                      pcspectra$scores[,2],
                      pcspectra$scores[,3],
                      xlab = "PC 1",
                      ylab = "PC 2",
                      zlab = "PC 3",
                      color=rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.5),
                      pch = 16, grid=TRUE, angle=50)

### location of the outliers in red
sct3d$points3d(pcspectra$scores[indxOutM,1],
               pcspectra$scores[indxOutM,2],
               pcspectra$scores[indxOutM,3],
               pch = "X",
               col = "red")

### visualize outilers spectra
matplot(x = colnames(data$spcA), y = t(data$spcA),
        xlab = "Wavelength (nm)",
        ylab = "Absorbance",
        type = "l",
        lty = 1,
        col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.3))
  
  matlines(x = colnames(data$spcA), y = t(data$spcA)[,indxOutM],
           type = "l",
           lty = 1,
           col = 'red')

### name it...
data$SAMPLE_CODE[indxOutM]

#==============================================================================#
#                           2.SIGNAL-TO-NOISE RATIO                            #
#==============================================================================#
# simple SNR estimate: mean / sd across a stable region
wav = as.numeric(colnames(data$spc)) # use reflectance!

snr_region = which(wav > 1000 & wav < 1100)

matplot(wav[snr_region], t(data$spc[ ,snr_region]),
        type = "l",
        lty = 1,
        col = rgb(0.5, 0.5, 0.5, alpha = 0.3))

snr = rowMeans(data$spc[, snr_region]) / apply(data$spc[, snr_region], 1, sd)
hist(snr, breaks = 30, main = "SNR (1000–1100 nm)")
which(snr < 50)  #  depends on instrument

#==============================================================================#
#                     3.DETECT NEGATIVE VALUES OR SPIKES                       #
#==============================================================================#
neg_reflectance = which(data$spc < 0, arr.ind = TRUE) # use reflectance!
row_with_neg = unique(neg_reflectance[, 1])

first_deriv = t(apply(as.matrix(data$spc), 1, diff))
spike_flags = apply(abs(first_deriv), 1, max) > 0.1  # arbitrary threshold
sum(spike_flags)
names(spike_flags)[spike_flags]

#==============================================================================#
#                     4.STANDARD NORMAL VARIATE (SNV) PROFILE                  #
#==============================================================================#
# use absorbance!
data_snv = t(apply(data$spcA, 1, function(x) (x - mean(x)) / sd(x)))
snv_sd = apply(data_snv, 1, sd)
hist(snv_sd, main = "SD after SNV")
flag_snv = which(snv_sd > 1.5 | snv_sd < 0.5)  # abnormal scaling

#==============================================================================#
#                        5.COSINE SIMILARITY/SECTRAL ANGLE                     #
#==============================================================================#
cos_sim = function(x, y) sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2)))

ref_spectrum = colMeans(data$spcA) # use absorbance!
similarity = apply(data$spcA, 1, cos_sim, y = ref_spectrum)
hist(similarity, main = "Cosine Similarity to Mean Spectrum")
low_sim = which(similarity < 0.98)


#==============================================================================#
#                                   6.FOLLOW UP                                #
#==============================================================================#

# compile (possible) dissimilar/outliers/bad spectra 
suspect_indices = unique(c(indxOutM, low_sim))

is_bad_spectra = data$spc[sort(suspect_indices),]

# plot against the mean to visually inspect it
matplot(x = as.numeric(colnames(is_bad_spectra)),
        y = t(is_bad_spectra),
        xlab = "Wavelength (nm)",
        ylab = "Reflectance",
        type = "l",
        lty = 1,
        col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.3))

matlines(x = as.numeric(colnames(data$spc)),
         y = colMeans(data$spc),
         type = "l",
         lty = 1,
         col = 'red')

# probably not outliers, but actually really different spectra...
data$SAMPLE_CODE[suspect_indices]
data$SAMPLE_CODE[data$MASS_mg == 3750]


paste0(
  (length(data$SAMPLE_CODE[suspect_indices]) * 100) / length(data$MASS_mg[data$MASS_mg == 3750]),
  "% of samples with highest plastic masses plus ", data$SAMPLE_CODE[suspect_indices][29])


# high plastic masses and large sizes that are not flagged as 'bad spectra':
data$SAMPLE_CODE[data$MASS_mg == 3750 & !(data$SAMPLE_CODE %in% data$SAMPLE_CODE[suspect_indices])]


# ... since it's more likely to be informative deviations:
high_plastic_spectra = is_bad_spectra
rm(is_bad_spectra)


#==============================================================================#
#                     7.CHECK SIMILARITY USING SIDE INFO                       #
#==============================================================================#

# take plastic masses and associated spectra
specPlast = data[, c("MASS_mg", "spcA")]

# vector with wavelength
wavs = as.numeric(colnames(specPlast$spcA))

# two disjoint subsets
ns = nrow(specPlast)

# percentage selected
pd = 0.25

# compute number of samples
nsamples = round(ns * pd, digits = 0)

# vector of samples indices
origInd = 1:ns

# randomly select the samples
set.seed(666)
selectSmp = sample(origInd, size = nsamples) # just the indexes

# split samples
xuData = specPlast[selectSmp, ]

xrData = specPlast[-selectSmp, ]

# compute full set PCs
combX = rbind(xrData$spcA, xuData$spcA)

PCspectra = pc_projection(as.matrix(combX), pc_selection = list('cumvar', 0.999),
                          method = 'pca', center = TRUE, scale = FALSE)
names(PCspectra)

# plot fisrt two PCs scores of xrData
plot(x = PCspectra$scores[1:nrow(xrData), 1],
     y = PCspectra$scores[1:nrow(xrData), 2],
     xlab = 'PC 1',
     ylab = 'PC 2',
     type = 'p',
     pch = 16,
     col = rgb(red = 0, green = 0.4, blue = 0.8, alpha = 0.5))
grid()

# add first two PCs scores of xuData
points(x = PCspectra$scores[-c(1:nrow(xrData)), 1],
       y = PCspectra$scores[-c(1:nrow(xrData)),2],
       xlab = 'PC 1',
       ylab = 'PC 2',
       pch = 16,
       col = rgb(red = 0.8, green = 0.4, blue = 0, alpha = 0.5))

# compute Mahalanobis pairwise distance from xrData
mdXr = f_diss(PCspectra$scores,
              PCspectra$scores,
              diss_method = "mahalanobis",
              center = FALSE, scale = FALSE)

# mdXr is the dissimilarity matrix of the spectra computed in its PC space...
#...select for each spectrum in xuData its most spectrally similar in xrData 

nearN = NULL

for(i in (nrow(xrData)+1):nrow(mdXr)){
  nn_i = order(mdXr[1:nrow(xrData), ][, i])[2]
  
  nearN = c(nearN, nn_i)
}
nearN[1:10]

plot(xuData[, 'MASS_mg'],
     xrData[nearN, 'MASS_mg'],
     xlab = 'Polymer mass reference spectra (mg)',
     ylab = 'Polymer mass from nearest neighbour (mg)',
     pch = 16,
     col = rgb(red = 1, green = 0.2, blue = 0.2, alpha = 0.5))
grid()

# linear tendency shows similarity between the reference and nearest neighbour 
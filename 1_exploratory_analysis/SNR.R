raw_snr = raw$spc
minimal_snr = data$spc # SGf only...
between_snr = movav(standardNormalVariate(data$spc), w = 11)

oldWavs = as.numeric(colnames(data$spcA))
newWavs = seq(min(oldWavs), max(oldWavs), by = 5)
full_snr = movav(standardNormalVariate(resample(data$spc,
                    wav = oldWavs,
                    new.wav = newWavs,
                    interpol = "linear")), w = 11)


wav = as.numeric(colnames(data$spc)) # use reflectance!
wav2 = as.numeric(colnames(minimal_snr))
wav3 = as.numeric(colnames(between_snr))
wav4 = as.numeric(colnames(full_snr))


snr_region = which(wav > 1000 & wav < 1100)
snr_region2 = which(wav2 > 1000 & wav2 < 1100)
snr_region3 = which(wav3 > 1000 & wav3 < 1100)
snr_region4 = which(wav4 > 1000 & wav4 < 1100)

snr = rowMeans(raw_snr[, snr_region]) / apply(raw_snr[, snr_region], 1, sd)
snr2 = rowMeans(minimal_snr[, snr_region2]) / apply(minimal_snr[, snr_region2], 1, sd)
snr3 = rowMeans(between_snr[, snr_region3]) / apply(between_snr[, snr_region3], 1, sd)
snr4 = rowMeans(full_snr[, snr_region4]) / apply(full_snr[, snr_region4], 1, sd)

round(mean(snr), 2)
round(mean(snr2), 2)
round(mean(snr3), 2)
round(mean(snr4), 2)

mean((snr2 - snr) / snr) * 100  # mean % improvement across spectra

mean((snr3 - snr) / snr) * 100

mean((snr4 - snr3) / snr3) * 100

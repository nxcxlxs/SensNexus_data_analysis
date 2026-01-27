require(asdreader)
require(dplyr)
require(tidyr)
require(ggplot2)

# import every spectra file 
setwd("C:/nico/Dissertação/SENSNEXUS_data/raw_spectra")
list_asd = list.files("./FIELDSPEC4_data",
                      pattern = "\\.asd$",
                      full.names = TRUE)

spec = get_spectra(list_asd, type = "reflectance")

# remove spectralon reading (blank reference)
spec = spec[2:nrow(spec), ]

# tidy up the spectra names and convert it to data.frame
rownames(spec) = sprintf("spec%03d", 1:nrow(spec)) # prefix `spec`+ 3 seq digits
spec = as.data.frame(spec)

# isolate references
soil_spec = spec[289:292, ]
PP_spec = spec[293:294, ] 
PVC_spec = spec[295:296, ]
PET_spec = spec[297:298, ]
PE_spec = spec[299:300, ]

# set contaminated soil data
spec = spec[1:288, ]

# load the sample design spreadsheet
design = read.csv2("C:/nico/Dissertação/SENSNEXUS_data/sample_design.csv")
design = design[, -c(1, 9)] # remove row names and prep check boxes

design = design |>
  mutate(MASS_mg = as.numeric(MASS_mg)) |> 
  mutate_if(~ !is.numeric(.x), as.factor) |> 
  mutate(SAMPLE_CODE = as.character(SAMPLE_CODE))

# merge datasets
design$spc = spec
datsoil = design

# save compiled data
saveRDS(datsoil, "./preprocessed_data/datsoil.rds")

datsoil = readRDS("./preprocessed_data/datsoil.rds")

# save references
ref_spec = list(
  soil = soil_spec,
  PP   = PP_spec,
  PVC  = PVC_spec,
  PET  = PET_spec,
  PE   = PE_spec
)

saveRDS(ref_spec, "./preprocessed_data/ref_spec.rds")

ref_spec = readRDS("./preprocessed_data/ref_spec.rds")

# visualization
## plot function for references
plot_references = function(df, title) {
  df_long = df |> 
    mutate(Sample = rownames(df)) |> 
    pivot_longer(
      cols = -Sample,
      names_to = "Wavelength",
      values_to = "Reflectance") |> 
    mutate(Wavelength = as.numeric(Wavelength))
  
  ggplot(df_long, aes(x = Wavelength, y = Reflectance, color = factor(Sample))) +
    geom_line(linewidth = 1) +
    # scale_x_continuous(limits = c(1800, 1850)) + # to check shifts or so...
    # scale_x_continuous(breaks = seq(1750, 1800, by = 10)) +
    # theme(axis.text.x = element_text(angle = 90)) +
    labs(
      x = "Wavelength (nm)",
      y = "Reflectance",
      color = "Sample",
      title = title)
}

p = plot_references(soil_spec, "Soils Spectra Data")
# note the “shifts” of these soil spectra (1000 and 1830 nm) <- Wadoux
## it can be corrected by prospectr::spliceCorretion()

wav = as.numeric(colnames(soil_spec))

# try to improve with the splice correction
require(prospectr)
require(patchwork)
require(asdreader)

raw = "./raw_spectra/spec00000.asd"

soil_spec2 = spliceCorrection(soil_spec, wav,
                              splice = c(as.numeric(get_metadata(raw)[30]), # 1000
                                         as.numeric(get_metadata(raw)[31]))) # 1800

p2 = plot_references(as.data.frame(soil_spec2), "Soils Spectra Data Corrected")

plot_references(PP_spec, "Polypropylene Reference Spectra Data")
plot_references(PVC_spec, "Polyvinyl Chloride Reference Spectra Data")
plot_references(PET_spec, "Polyethylene Terephthalate Reference Spectra Data")
plot_references(PE_spec, "High-density Polyethylene Reference Spectra Data")


## plot function for samples
plot_spectra = function(data, code_prefix) {
  rows = which(startsWith(data$SAMPLE_CODE, code_prefix))
  df = data$spc[rows, ]
  
  df_long = df |>
    mutate(Sample = data$SAMPLE_ID[rows]) |>
    pivot_longer(
      cols = -Sample,
      names_to = "Wavelength",
      values_to = "Reflectance") |>
    mutate(Wavelength = as.numeric(Wavelength))
  
  ggplot(df_long, aes(x = Wavelength, y = Reflectance, color = factor(Sample))) +
    geom_line(linewidth = 1) +
    labs(
      x = "Wavelength (nm)",
      y = "Reflectance",
      color = "SAMPLE_ID",
      title = code_prefix)
}

plot_spectra(datsoil, "PVC_A1.2")

require(asdreader)
require(dplyr)
require(tidyr)
require(ggplot2)


# import every spectra file)
list_asd = list.files("../raw_spectra/FIELDSPEC4_data",
                      pattern = "\\.asd$",
                      full.names = TRUE)

spec = get_spectra(list_asd, type = "reflectance")

# remove spectralon reading (blank reference)
spec = spec[2:nrow(spec), ]

# tidy up the spectra names and convert it to data.frame
rownames(spec) = sprintf("spec%03d", 1:nrow(spec)) # prefix `spec`+ 3 seq digits
spec = as.data.frame(spec)

# isolate contaminated soil data (see `raw_spectra/notes [31.05.2025].txt`)
spec = spec[1:288, ]

# load the sample design file
design = read.csv2("../raw_spectra/pristine_sample_design.csv")
design = design[, -8] # remove row names and prep check boxes

design = design |>
  mutate(MASS_mg = as.numeric(MASS_mg)) |> 
  mutate_if(~ !is.numeric(.x), as.factor) |> 
  mutate(SAMPLE_CODE = as.character(SAMPLE_CODE))

# merge datasets
design$spc = spec
pristine = design

# save compiled data
saveRDS(pristine, "../raw_spectra/raw_pristine.rds")


## visualization
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

plot_spectra(datsoil, "PVC_A1.2") # e.g.

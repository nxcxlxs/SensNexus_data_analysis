require(tidyverse)
require(ggplot2)
require(tidyr)

# load reference spectra
ref_spec = readRDS("./preprocessed_data/ref_corrected.rds")

# function to plot easily
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
    scale_color_manual(values = soil$POLYMER) +
    # scale_color_manual(values = c("black", "grey40")) +
    labs(
      x = "Wavelength (nm)",
      y = "Reflectance",
      color = "Sample",
      title = title)
}

# HDPE
plot_references(ref_spec$PE, "PE characteristic NIR wavelengths") +
  geom_vline(xintercept = c(1210, 1720), linetype = "dashed",
             color = "red", linewidth = 0.8)

# PET
plot_references(ref_spec$PET, "PET characteristic NIR wavelengths") +
  geom_vline(xintercept = c(1200, 1420, 1660, 1730, 1910),
             linetype = "dashed", color = "red", linewidth = 0.8)

# PVC
plot_references(ref_spec$PVC, "PVC characteristic NIR wavelengths") +
  geom_vline(xintercept = c(930, 1040, 1210, 1418, 1420, 1715, 1730, 1740),
             linetype = "dashed", color = "red", linewidth = 0.8)

# PP
range_list = list(seq(914, 934, 1),
             seq(1020, 1052, 1),
             seq(1158, 1214, 1),
             seq(1707, 1763, 1),
             seq(1824, 1983, 1),
             seq(2067, 2466, 1))

range = unlist(range_list)

range_means = sapply(range_list, mean)

plot_references(ref_spec$PP, "PP characteristic NIR wavelengths") +
  geom_vline(xintercept = c(1630, range),
             linetype = "dashed", color = "darkorange", alpha = 0.3) +
  geom_vline(xintercept = c(1630, range_means),
             linetype = "dashed", color = "red",
             linewidth = 0.8) 

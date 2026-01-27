# load samples dataset
data = readRDS("C:/nico/Dissertação/SENSNEXUS_data/preprocessed_data/splice_corrected.rds")

## plot function for samples
require(dplyr)
require(tidyr)
require(ggplot2)
require(patchwork)

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

plot_spectra(data, "PP_A3750")

# extract tratments
treatments = unique(sub("(_\\d+)$", "", data$SAMPLE_CODE))

# save plots
for (t in treatments) {
  p = plot_spectra(data, t)
  ggsave(
    filename = paste0("plots/", t, ".png"),
    plot = p,
    dpi = 300)
}

# load references dataset
ref = readRDS("./preprocessed_data/ref_corrected2.rds")

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

plot_references(ref$soil, "Soil Spectra Data")

for (name in names(ref)) {
  data = ref[[name]]
  p = plot_references(data, title = name)
  ggsave(
    filename = paste0("plots/reference_spectra/", "colored_", name, ".png"),
    plot = p,
    dpi = 300)
}

################################################################################
#                            COMPILED TREATMENTS                               #
################################################################################
library(RColorBrewer)

mass_levels = levels(as.factor(data$MASS_mg))
# size_levels = levels(as.factor(data$SIZE_CODE))
n_mass = length(mass_levels)
# n_size = length(size_levels)
cols = brewer.pal(n = n_mass, name = "RdBu")
# cols = brewer.pal(n = n_size, name = "RdBu")
names(cols) = mass_levels
# names(cols) = size_levels


# loop over each polymer
for (poly in levels(factor(data$POLYMER))) {
  
  sub_data = data[data$POLYMER == poly, ]
  color_vector = cols[as.character(sub_data$MASS_mg)] # SIZE_CODE
  
  # set up layout with extra space for legend
  par(mar = c(5, 4, 4, 6), xpd = TRUE)  # Extra space on the right
  
  matplot(x = as.numeric(colnames(sub_data$spc)),
          y = t(sub_data$spc),
          type = "l",
          lty = 1,
          col = color_vector,
          xlab = "Wavelength /nm",
          ylab = "Reflectance",
          main = "MIX") # paste(poly)
  
  # add legend outside plot area
  legend("topright",
         inset = c(-0.18, 0),
         legend = names(cols),
         col = cols,
         pch = 19,
         title = "MASS (mg)", # SIZE_CODE
         cex = 0.8)
}

# check specific spectra (via SAMPLE_CODE)
matplot(x = colnames(data$spc), y = t(data$spc),
        xlab = "Wavelength (nm)",
        ylab = "Reflectance",
        type = "l",
        lty = 1,
        col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.3)) +
  
  matlines(x = colnames(data$spc), y = t(data$spc)[, data$SAMPLE_CODE == "MIX_C6_4"],
           col = 'red')

################################################################################

raw = readRDS("C:/nico/Dissertação/SENSNEXUS_data/preprocessed_data/datsoil.rds")
raw_new = readRDS("C:/nico/Dissertação/SENSNEXUS_data/preprocessed_data/datsoil2.rds")

raw_combined = rbind(raw, raw_new)

raw_long = bind_cols(
  raw_combined |> select(-spc),
  raw_combined$spc
) |>
  pivot_longer(
    cols = as.character(350:2500),
    names_to = "Wavelength",
    values_to = "Reflectance"
  ) |>
  mutate(Wavelength = as.numeric(Wavelength)) |> 
  mutate(POLYMER = recode(POLYMER,
                          "PE" = "HDPE",
                          "ALL" = "MIX"))

p1 = ggplot(raw_long, aes(x = Wavelength, y = Reflectance,
                     group = SAMPLE_ID, color = POLYMER)) +
  geom_line(alpha = 0.15) +
  scale_color_grey(start = 0, end = 0.7) +
  # scale_color_manual(values = c("darkgreen", "black", "red", "dodgerblue2", "goldenrod2")) +
  labs(
    x = "Wavelength (nm)",
    y = "Reflectance") +
  theme(legend.title = element_blank(),
        legend.position = "none") +
  guides(color = guide_legend(override.aes = list(alpha = 1)))


r = "C:/nico/Dissertação/SENSNEXUS_data/raw_spectra/FIELDSPEC4_data/spec00000.asd"

require(prospectr)
require(asdreader)
correct_raw = spliceCorrection(raw_combined$spc,
                               as.numeric(colnames(raw_combined$spc)),
                                splice = c(as.numeric(get_metadata(r)[30]),
                                           as.numeric(get_metadata(r)[31])))

raw_combined2 = raw_combined
raw_combined2$spc = as.data.frame(correct_raw)

corr_long = bind_cols(
  raw_combined2 |> select(-spc),
  raw_combined2$spc
) |>
  pivot_longer(
    cols = as.character(350:2500),
    names_to = "Wavelength",
    values_to = "Reflectance"
  ) |>
  mutate(Wavelength = as.numeric(Wavelength)) |> 
  mutate(POLYMER = recode(POLYMER,
                          "ALL" = "MIX"))

p2 = ggplot(corr_long, aes(x = Wavelength, y = Reflectance,
                     group = SAMPLE_ID, color = POLYMER)) +
  geom_line(alpha = 0.15) +
  scale_color_grey(start = 0, end = 0.7) +
  scale_color_manual(values = c("darkgreen", "black", "red", "dodgerblue2", "goldenrod2")) +
  labs(
    x = "Wavelength (nm)",
    y = "Reflectance") +
  theme(legend.title = element_blank(),
        legend.position = "top") +
  guides(color = guide_legend(override.aes = list(alpha = 1)))
 
p1/p2

# just pristine polymer soil contaminated spectra
polymer_type = ggplot(corr_long |> filter(!POLYMER == "MIX"),
                      aes(Wavelength, Reflectance,
                          group = SAMPLE_ID, color = factor(POLYMER))) +
                geom_line(aes(alpha = 0.15), linewidth = 0.45) +
                # scale_color_brewer(palette = "Spectral", direction = -1) +
                scale_color_manual(values = c("#3E88B3", "#54BF34",
                                              "#FF7F32", "darkred"),
                                   name = "Polymer") +
                scale_alpha_identity(guide = "none") +
                theme(axis.title.y = element_text(size = 15),
                      axis.text.y = element_text(size = 15),
                      axis.title.x = element_blank(),
                      axis.text.x = element_text(size = 15),
                      legend.position = "bottom",
                      legend.title = element_text(size = 15),
                      legend.text = element_text(size = 15),
                      legend.key.size = unit(1, 'cm')) +
                guides(color = guide_legend(nrow = 1))

mass_classes = ggplot(corr_long |> filter(!POLYMER == "MIX"),
                      aes(Wavelength, Reflectance,
                          group = SAMPLE_ID, color = factor(MASS_mg))) +
                geom_line(aes(alpha = 0.30), linewidth = 0.45) +
                scale_color_brewer(palette = "Spectral", direction = -1,
                                   name = "Mass (mg)") +
                scale_alpha_identity(guide = "none") +
                theme(axis.title.y = element_text(size = 15),
                      axis.text.y = element_text(size = 15),
                      axis.title.x = element_blank(),
                      axis.text.x = element_text(size = 15),
                      legend.position = "bottom",
                      legend.title = element_text(size = 15),
                      legend.text = element_text(size = 15),
                      legend.key.size = unit(1, 'cm')) +
                      guides(color = guide_legend(nrow = 1))

size_fractions = ggplot(corr_long |> filter(!POLYMER == "MIX"),
                        aes(Wavelength, Reflectance,
                            group = SAMPLE_ID, color = factor(SIZE_INTERVALS_mm))) +
                  geom_line(aes(alpha = 0.50), linewidth = 0.45) +
                  scale_color_brewer(palette = "Spectral", direction = 1,
                                     name = "Size intervals (mm)") +
                  scale_alpha_identity(guide = "none") +
                  theme(axis.title.y = element_blank(),
                        axis.text.y = element_text(size = 15),
                        axis.title.x = element_blank(),
                        axis.text.x = element_text(size = 15),
                        legend.position = "bottom",
                        legend.title = element_text(size = 15),
                        legend.text = element_text(size = 15),
                        legend.key.size = unit(1, 'cm'))

  # just colored household polymer soil contaminated spectra
mass_classes2 = ggplot(corr_long |> filter(POLYMER == "MIX"),
                      aes(Wavelength, Reflectance,
                          group = SAMPLE_ID, color = factor(MASS_mg))) +
                geom_line(aes(alpha = 0.30), linewidth = 0.45) +
                scale_color_brewer(palette = "Spectral", direction = -1,
                                   name = "Mass (mg)") +
                scale_alpha_identity(guide = "none") +
                theme(axis.title.y = element_text(size = 15),
                      axis.text.y = element_text(size = 15),
                      axis.title.x = element_text(size = 15),
                      axis.text.x = element_text(size = 15),
                      legend.position = "bottom",
                      legend.title = element_text(size = 15),
                      legend.text = element_text(size = 15),
                      legend.key.size = unit(1, 'cm'))

size_fractions2 = ggplot(corr_long |> filter(POLYMER == "MIX"),
                        aes(Wavelength, Reflectance,
                            group = SAMPLE_ID, color = factor(SIZE_INTERVALS_mm))) +
                  geom_line(aes(alpha = 0.50), linewidth = 0.45) +
                  scale_color_brewer(palette = "Spectral", direction = 1,
                                     name = "Size intervals (mm)") +
                  scale_alpha_identity(guide = "none") +
                  theme(axis.title.y = element_blank(),
                        axis.text.y = element_text(size = 15),
                        axis.title.x = element_text(size = 15),
                        axis.text.x = element_text(size = 15),
                        legend.position = "bottom",
                        legend.title = element_text(size = 15),
                        legend.text = element_text(size = 15),
                        legend.key.size = unit(1, 'cm')) 

prist_color = ggplot(
                corr_long,
                aes(
                  Wavelength,
                  Reflectance,
                  group = SAMPLE_ID,
                  color = ifelse(POLYMER == "MIX",
                                 "Pristine",
                                 "Colored household")
                )
              ) +
                geom_line(alpha = 0.25, linewidth = 0.45) +
                scale_color_manual(
                  name   = NULL,
                  values = c("Pristine" = "gray16",
                             "Colored household" = "darkred")) +
                guides(color = guide_legend(override.aes = list(linewidth = 1,
                                                                alpha = 1))) +
                  theme(axis.title.y = element_blank(),
                        axis.text.y = element_text(size = 15),
                        axis.title.x = element_blank(),
                        axis.text.x = element_text(size = 15),
                        legend.position = "bottom",
                        legend.title = element_text(size = 15),
                        legend.text = element_text(size = 15),
                        legend.key.size = unit(1, 'cm')) +
                  scale_alpha_identity(guide = "none")

require(patchwork)
(polymer_type + prist_color) / (mass_classes + size_fractions) / (mass_classes2 + size_fractions2) +
  plot_annotation(tag_levels = 'a',
                  # tag_levels = list(c('b', 'c', 'd', 'e')),
                  tag_prefix = '(', tag_suffix = ')')

# polymer_type + 
#   plot_annotation(tag_levels = 'a',
#                   tag_prefix = '(', tag_suffix = ')')

################################################################################

# just pristine polymer soil contaminated spectra
polymer_type = ggplot(corr_long |> filter(!POLYMER == "MIX"),
                      aes(Wavelength, Reflectance,
                          group = SAMPLE_ID, color = factor(POLYMER))) +
  geom_line(aes(alpha = 0.30), linewidth = 0.45) +
  # scale_color_manual(values = c("darkgreen", "black",
                                # "darkred", "steelblue"),
                     # name = "Polymer type") +
  # scale_color_viridis_d(option = "viridis",
  #                       name = "Polymer type") +
  scale_color_brewer(palette = "Spectral", direction = -1,
                     name = "Polymer") +
  scale_alpha_identity(guide = "none") +
  theme(axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 15),
        legend.position = "bottom",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1, 'cm')) +
  guides(color = guide_legend(nrow = 1)) +
  guides(color = guide_legend(override.aes = list(linewidth = 1, alpha = 1),
                              nrow = 1))

mass_classes = ggplot(corr_long |> filter(!POLYMER == "MIX"),
                      aes(Wavelength, Reflectance,
                          group = SAMPLE_ID, color = factor(MASS_mg))) +
  geom_line(aes(alpha = 0.30), linewidth = 0.45) +
  scale_color_brewer(palette = "Spectral", direction = -1,
                     name = "Mass (mg)") +
  # scale_color_viridis_d(option = "viridis",
  #                       name = "Mass (mg)") +
  scale_alpha_identity(guide = "none") +
  theme(axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 15),
        legend.position = "bottom",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1, 'cm')) +
  guides(color = guide_legend(nrow = 1)) +
  guides(color = guide_legend(override.aes = list(linewidth = 1, alpha = 1),
                              nrow = 1))

spectral_modified3 = c("#2B83BA", "#ABDDA4", "#D7191C")

size_fractions = ggplot(corr_long |> filter(!POLYMER == "MIX"),
                        aes(Wavelength, Reflectance,
                            group = SAMPLE_ID, color = factor(SIZE_INTERVALS_mm))) +
  geom_line(aes(alpha = 0.50), linewidth = 0.45) +
  scale_color_manual(values = spectral_modified3,
                     name = "Size intervals (mm)") +
  # scale_color_brewer(palette = "Spectral", direction = 1,
  #                    name = "Size intervals (mm)") +
  # scale_color_viridis_d(option = "viridis",
  #                       name = "Size intervals (mm)") +
  scale_alpha_identity(guide = "none") +
  theme(axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        legend.position = "bottom",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1, 'cm')) +
  guides(color = guide_legend(override.aes = list(linewidth = 1, alpha = 1),
                              nrow = 1))

# just colored household polymer soil contaminated spectra
mass_classes2 = ggplot(corr_long |> filter(POLYMER == "MIX"),
                       aes(Wavelength, Reflectance,
                           group = SAMPLE_ID, color = factor(MASS_mg))) +
  geom_line(aes(alpha = 0.50), linewidth = 0.45) +
  scale_color_brewer(palette = "Spectral", direction = -1,
                     name = "Mass (mg)") +
  # scale_color_viridis_d(option = "viridis",
  #                       name = "Mass (mg)") +
  scale_alpha_identity(guide = "none") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 15),
        legend.position = "bottom",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1, 'cm')) +
  guides(color = guide_legend(override.aes = list(linewidth = 1, alpha = 1),
                              nrow = 1))

size_fractions2 = ggplot(corr_long |> filter(POLYMER == "MIX"),
                         aes(Wavelength, Reflectance,
                             group = SAMPLE_ID, color = factor(SIZE_INTERVALS_mm))) +
  geom_line(aes(alpha = 0.50), linewidth = 0.45) +
  scale_color_manual(values = spectral_modified3,
                     name = "Size intervals (mm)") +
  # scale_color_brewer(palette = "Spectral", direction = 1,
  #                    name = "Size intervals (mm)") +
  # scale_color_viridis_d(option = "viridis",
  #                       name = "Size intervals (mm)") +
  scale_alpha_identity(guide = "none") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        legend.position = "bottom",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1, 'cm')) +
  guides(color = "none")

prist_color = ggplot(
  corr_long,
  aes(
    Wavelength,
    Reflectance,
    group = SAMPLE_ID,
    color = ifelse(!POLYMER == "MIX",
                   "Pristine",
                   "Colored")
  )
) +
  geom_line(alpha = 0.25, linewidth = 0.45) +
  scale_color_manual(
    name   = NULL,
    values = c("Pristine" = "black",
               "Colored" = "darkred")) +
  scale_alpha_identity(guide = "none") +
  guides(color = guide_legend(override.aes = list(linewidth = 1, alpha = 1),
                              nrow = 1)) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 15),
        legend.position = "bottom",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1, 'cm'))

# Final plot arrangement
(polymer_type + prist_color) / (mass_classes + mass_classes2) / (size_fractions + size_fractions2) +
  plot_annotation(tag_levels = 'a',
                  tag_prefix = '(', tag_suffix = ')') &
  theme(plot.tag = element_text(face = 'bold', size = 20))

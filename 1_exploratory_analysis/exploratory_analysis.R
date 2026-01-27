# load processed data
data = readRDS("C:/nico/Dissertação/SENSNEXUS_data/preprocessed_data/datsoil.rds")
data2 = readRDS("C:/nico/Dissertação/SENSNEXUS_data/analysis_ready_data/deNoised_data.rds")
new = readRDS("C:/nico/Dissertação/SENSNEXUS_data/preprocessed_data/datsoil2.rds")
new2 = readRDS("C:/nico/Dissertação/SENSNEXUS_data/analysis_ready_data/deNoised_newData.rds") 


# convert spectra to absorbance
data$spcA = log(1/data$spc)
data2$spcA = log(1/data2$spc)
new$spcA = log(1/new$spc)
new2$spcA = log(1/new2$spc)


# preprocess approaches
## 5 nm resample
oldWavs = as.numeric(colnames(data2$spcA))
newWavs = seq(min(oldWavs), max(oldWavs), by = 5)

data2$spcAR = prospectr::resample(data2$spcA,
                                 wav = oldWavs,
                                 new.wav = newWavs,
                                 interpol = "linear")

new2$spcAR = prospectr::resample(new2$spcA,
                                wav = oldWavs,
                                new.wav = newWavs,
                                interpol = "linear")

require(prospectr)
## SNV for baseline correction
data2$spcARsnv = standardNormalVariate(data2$spcAR)
data2$spcAsnv = standardNormalVariate(data2$spcA)

new2$spcARsnv = standardNormalVariate(new2$spcAR)
new2$spcAsnv = standardNormalVariate(new2$spcA)

## Moving Window Average to the SNV spectra
data2$spcARmovav = movav(data2$spcARsnv, w = 11)
data2$spcAmovav = movav(data2$spcAsnv, w = 11)

new2$spcARmovav = movav(new2$spcARsnv, w = 11)
new2$spcAmovav = movav(new2$spcAsnv, w = 11)

# try and remove "noisy" edges...

# edges = which(as.numeric(colnames(data$spcA)) < 500 | as.numeric(colnames(data$spcA)) > 2450)

# data$DNspcA = data$spcA[, -edges]

# ...no cigar.


################################################################################
#                         PRINCIPAL COMPONENT ANALYSIS                         #
################################################################################
# reduce dimensionality with PCA
pcspec = prcomp(data$spcA) # raw data
pcspec2 = prcomp(data2$spcA) # SGf data
pcspec3 = prcomp(data2$spcAmovav) # SGf + SNV + movav
pcspec4 = prcomp(data2$spcARmovav) # SGf + SNV + movav + resample


summary(pcspec)[[6]][, 1:10] # good amount of variance explained by it...
summary(pcspec2)[[6]][, 1:10]
summary(pcspec3)[[6]][, 1:10]
summary(pcspec4)[[6]][, 1:10]

# scores
plot(pcspec$x[, 1],
     pcspec$x[, 2],
     xlab = "PC1",
     ylab = "PC2",
     pch = 16,
     col = rgb(0.5, 0.5, 0.5, alpha = 0.5))

# first loading
plot(colnames(data$spcA), pcspec$rotation[,1],
     type ="l",
     ylab ="Loading",
     xlab ="Wavelength (nm)",
     col = rgb(red = 1, green = 0, blue = 0, alpha = 1),
     ylim = c(-0.06, 0.06))

# second loading
lines(colnames(data$spcA), pcspec$rotation[,2],
      type ="l",
      ylab ="Loading",
      xlab ="Wavelength (nm)",
      col = rgb(red = 0, green = 0, blue = 1, alpha = 1))
# add a legend
legend("topright",
       legend = c("First loading", "Second loading"),
       lty = c(1, 1),
       col = c("red", "blue"))
abline(h=0)

# last principal component loading
plot(colnames(data$spcA), pcspec$rotation[,ncol(pcspec$rotation)],
     main = "Last PC loading (mostly noise)",
     type = "l",
     ylab = "Loading",
     xlab = "Wavelength (nm)",
     col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.5),
     ylim = c(-0.1, 0.1))

# OVERPLOTTED! 
# # biplot of the first two principal components scores and loadings
# biplot(pcspec,
#        col = c(rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 1),
#                rgb(red = 1, green = 0, blue = 0, alpha = 0.5)),
#        xlab = "PC 1",
#        ylab = "PC 2")
# abline(v=0, h=0)
# 

# spectral prediction domain
require(tripack)

# PCA scores triangulation
randTr = tri.mesh(pcspec$x[, 1], pcspec$x[, 2])
randTr2 = tri.mesh(pcspec2$x[, 1], pcspec2$x[, 2])
randTr3 = tri.mesh(pcspec3$x[, 1], pcspec3$x[, 2])
randTr4 = tri.mesh(pcspec4$x[, 1], pcspec4$x[, 2])

# nodes of the convex hull
randCH = convex.hull(randTr, plot.it = F)
randCH2 = convex.hull(randTr2, plot.it = F)
randCH3 = convex.hull(randTr3, plot.it = F)
randCH4 = convex.hull(randTr4, plot.it = F)

plot(pcspec$x[,1], pcspec$x[,2],
     xlab ="PC1",
     ylab ="PC2",
     # xlim = c(min(pcspec$x[,1:2]), max(pcspec$x[,1:2])),
     # ylim = c(min(pcspec3$x[,1:2]), max(pcspec3$x[,1:2])), # watch out for useless plot spaces
     ylim = c(-20, 20),
     pch = 16,
     col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.5))

lines(c(randCH$x, randCH$x[1]), c(randCH$y,randCH$y[1]),
      col="darkred",
      lwd=2)
#==============================================================================#
#           CHECKING IF NEW DATA ARE INSIDE THE PREDICTIVE SPACE               #
#==============================================================================#
new$spcA = as.matrix(new$spcA)
new_scores = sweep(new$spcA, 2, pcspec$center, "-") %*% pcspec$rotation

points(new_scores[,1], new_scores[,2],
     pch = 16,
     col = adjustcolor("blue", alpha.f = 0.4))

legend("topright",
       legend = c("Pristine polymer samples",
                  "Colored polymer samples",
                  "Convex hull delimitation"),
       col    = c("gray", "blue", "darkred"),
       pch    = c(16, 16, NA),
       lty    = c(NA, NA, 1),
       lwd    = c(NA, NA, 2))
  
require(splancs)
poly_coords = cbind(randCH$x, randCH$y)
poly_coords = rbind(poly_coords, poly_coords[1, ])

pts = as.matrix(new_scores[, 1:2])

inside = rep(FALSE, nrow(pts))
inside_idx = inpip(pts, poly_coords)  # returns indices

inside[inside_idx] = TRUE

# now 'inside' is logical
table(inside)

# plot outside points
points(pts[!inside, 1], pts[!inside, 2],
       pch = 4,
       col = adjustcolor("red", alpha.f = 0.8),
       lwd = 2)

require(dplyr)
new |> 
  filter(!inside) |> 
  select(SAMPLE_ID, SAMPLE_CODE, SIZE_INTERVALS_mm,
         MASS_mg, REPLICATE)

paste0(round((sum(inside == FALSE) * 100) / nrow(new), 2), "% ",
    "OF THE NEW SAMPLES ARE OUTSIDE THE PREDICTIVE DOMAIN")

## ggplot version =============================================================#
hull_x = c(randCH$x, randCH$x[1])
hull_y = c(randCH$y, randCH$y[1])

legend_levels = c("Pristine polymer samples",
                   "Colored polymer samples",
                   "Samples outside prediction domain",
                   "Convex hull delimitation")

ggplot() +
  geom_point(aes(x = pcspec$x[,1], y = pcspec$x[,2],
                 color = factor("Pristine polymer samples", levels = legend_levels)),
             size = 3, alpha = 0.5) +
  
  geom_point(aes(x = new_scores[,1], y = new_scores[,2],
                 color = factor("Colored polymer samples", levels = legend_levels)),
             size = 3, alpha = 0.6) +
  
  geom_point(aes(x = new_scores[!inside,1], y = new_scores[!inside,2],
                 color = factor("Samples outside prediction domain", levels = legend_levels)),
             shape = 4, size = 3, stroke = 1.5) +
  
  geom_path(aes(x = hull_x, y = hull_y,
                color = factor("Convex hull delimitation", levels = legend_levels)),
            linewidth = 1) +
  
  scale_color_manual(name = "",
                     values = c("Pristine polymer samples" = adjustcolor("steelblue"),
                                "Colored polymer samples" = adjustcolor("darkblue"),
                                "Samples outside prediction domain" = adjustcolor("red", 0.6),
                                "Convex hull delimitation" = "darkred")) +
  
  labs(x = "PC1", y = "PC2") +
  theme(legend.position = "top",
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1, 'cm'))
#==============================================================================#

require(ggplot2)
scores_df = as.data.frame(pcspec$x)

# loadings
loadings_df = as.data.frame(pcspec$rotation)

scores_df$Obs = rownames(scores_df)
loadings_df$Var = rownames(loadings_df)

# sum_abs_loadings_pc1 = abs(pcspec$rotation)
# important_bands = which(sum_abs_loadings_pc1 > quantile(sum_abs_loadings_pc1, 0.9))
# print(colnames(data$spcA)[important_bands])

# select variables (some of the wavelengths)
select_var = c(min(as.numeric(colnames(data$spcA))),
               median(as.numeric(colnames(data$spcA))),
               max(as.numeric(colnames(data$spcA))))

# filter loadings
# filt_loading = loadings_df[loadings_df$Var %in% select_var, ]

filt_loading = loadings_df[loadings_df$Var %in% c("350", "700", "1100", "2500"), ]

p = ggplot() +
  geom_point(data = scores_df,
             aes(x = PC1, y = PC2,
                 color = "Samples"),
             alpha = 0.5,
             size = 3) +
  labs(x = "PC1 (77.2%)", y = "PC2 (19.7%)") +
  scale_color_manual(name = NULL, values = c("Samples" = "steelblue", "Wavelengths" = "black")) +
  geom_segment(data = filt_loading, aes(x = 0, y = 0,
                                        xend = PC1*100,
                                        yend = PC2*100,
                                        color = "Wavelengths"),
               linewidth = 1,
               arrow = arrow(type = "open", length = unit(0.25, "cm"))) + 
  geom_text(data = filt_loading, colour = "black", size = 5, fontface = "bold",
            aes(x = PC1*100, y =PC2*100,
                label = Var),
            # fontface = "bold",
            hjust = .5,
            vjust = 2,
            nudge_x = c(0.4, 0.5, -0.5, 0),   # adjust for each wavelength
            nudge_y = c(0, 0.5, -0, 0)) +
  theme(axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        legend.position = "top",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1, 'cm'))

#==============================================================================#
require(factoextra)
# eigenvalues vs. number of dimensions
fviz_eig(pcspec, addlabels = T, geom = "bar",
         ncp = 3, main = "", ggtheme = theme_gray())


# spectra distribution (grouped by mass)
require(RColorBrewer)
p3 = fviz_pca_ind(pcspec, geom = "point",
             title = "",
             pointshape = 16,
             pointsize = 2,
             alpha.ind = 0.6,
             col.ind = as.factor(data$MASS_mg),
             # addEllipses = T,
             # ellipse.level = 0.5,
             ggtheme = theme_gray(),
             axes.linetype = "blank") +
  scale_fill_brewer(palette = "Spectral") +
  scale_color_brewer(palette = "Spectral") +
  labs(col = "Mass (mg)",
       x = "PC1",
       y = "PC2") + 
  guides(fill = "none") +
         # color = guide_legend(nrow = 1)) +
  theme(legend.position = "right",
        legend.direction = "vertical",
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1, 'cm'))
  # theme(legend.position = "top",
        # legend.direction = "horizontal")


# spectra distribution (grouped by size)
p2 = fviz_pca_ind(pcspec, geom = "point",
             title = "",
             pointshape = 16,
             pointsize = 2,
             alpha.ind = 0.8,
             addEllipses = F,
             col.ind = data$SIZE_INTERVALS_mm,
             ggtheme = theme_gray(),
             axes.linetype = "blank") +
  scale_fill_brewer(palette = "Spectral") +
  scale_color_brewer(palette = "Spectral") +
  labs(col = "Size intervals (mm)",
       x = "PC1",
       y = "PC2") + 
  guides(fill = "none") +
  theme(legend.position = "right",
        legend.direction = "vertical",
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1, 'cm'))


# spectra distribution (grouped by polymer)
p1 = fviz_pca_ind(pcspec, geom = "point",
             title = "",
             pointshape = 16,
             pointsize = 2,
             alpha.ind = 0.6,
             col.ind = data$POLYMER,
             addEllipses = F,
             # ellipse.level = 0.25,
             ggtheme = theme_gray(),
             axes.linetype = "blank") +
  scale_fill_brewer(palette = "Spectral") +
  scale_color_brewer(palette = "Spectral") +
  labs(col = "Polymer",
       x = "PC1",
       y = "PC2") + 
  guides(fill = "none") +
  theme(legend.position = "right",
        legend.direction = "vertical",
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1, 'cm'))

p4 = fviz_pca_ind(pcspec, geom = "point",
                  title = "",
                  pointshape = 16,
                  pointsize = 4,
                  alpha.ind = 0.6,
                  col.ind = as.factor(data$MASS_mg),
                  addEllipses = T,
                  # ellipse.level = 0.5,
                  ggtheme = theme_gray(),
                  axes.linetype = "blank") +
  scale_fill_brewer(palette = "Spectral") +
  scale_color_brewer(palette = "Spectral") +
  labs(col = "Mass (mg)",
       x = "PC1",
       y = "PC2") + 
  guides(fill = "none",
        color = guide_legend(nrow = 1, override.aes = list(size = 8))) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1, 'cm'))


fviz_pca_biplot(pcspec,
                geom.ind = "point",
                pointshape = 21,
                pointsize = 2.5,
                fill.ind = as.factor(data$MASS_mg),
                col.var = "darkred",
                select.var = list(name = select_var),
                arrowsize = 1) +
  labs(fill = "Plastic mass (mg)") +
  guides(color = "none")

p/(p1+p2+p3)/p4 +
  plot_annotation(tag_levels = 'a',
                  tag_prefix = '(', tag_suffix = ')') &
  theme(plot.tag = element_text(face = 'bold', size = 30))

################################################################################
#                           SPECTRA SIMILARITY ANALYSIS                        #
################################################################################
require(dplyr)
groups = data2 |> # deNoised_data.rds
  group_by(POLYMER, SIZE_CODE, SIZE_INTERVALS_mm, MASS_mg) |>
  summarise(mean_spectrum = list(colMeans(spcA)), .groups = "drop")

mean_spectra_matrix = do.call(rbind, groups$mean_spectrum)


require(resemble)
# compute Euclidian distance
EucD = f_diss(Xr = mean_spectra_matrix,
              Xu = mean_spectra_matrix,
              diss_method = "euclid",
              center = TRUE, scale = TRUE)


# compute Mahalanobis distance
PCA = pc_projection(mean_spectra_matrix,
                    pc_selection = list("cumvar", 0.99),
                    method = "pca",
                    center = TRUE, SCALE = FALSE)

mahD = f_diss(Xr = PCA$scores,
              Xu = PCA$scores,
              diss_method = "mahalanobis",
              center = FALSE, scale = FALSE)


# correlation similarity
scorr = cor(t(mean_spectra_matrix))
# compute similarity
cd1 = (1 - scorr)/2


# moving window for correlation similarity
mwcd = cor_diss(mean_spectra_matrix,
                mean_spectra_matrix,
                ws = 51, # must be an odd value
                center = FALSE,
                scale = FALSE)


# spectral angle mapper (dot-product cosine distance)
samD = f_diss(mean_spectra_matrix,
              mean_spectra_matrix,
              diss_method = "cosine",
              center = FALSE, scale = FALSE)

appr = c("EucD", "mahD", "cd1", "mwcd", "samD")

appr_dict = c(
  "EucD" = "Euclidean Distance",
  "mahD" = "Mahalanobis Distance",
  "cd1" = "Correlation Dissimilarity",
  "mwcd" = "Moving Window Correlation Dissimilarity",
  "samD" = "Spectral Angler Mapper")


# results = list()
# results2 = list()
results3 = list()

# for (mass in sort(unique(groups$MASS_mg), decreasing = T)) {
for (polymer in unique(groups$POLYMER)) {# NEED TO BE AWARE OF WHAT I WANT TO CONTROL
  # for (size in unique(groups$SIZE_CODE)) {
  for (mass in sort(unique(groups$MASS_mg), decreasing = T)) {# SUBSETTING SOME VARS FIX IT TO SEE VARIATON IN OTHER ONES
    
    # subset_idx = which(groups$MASS_mg == mass & groups$SIZE_CODE == size)
    # subset_idx = which(groups$POLYMER == polymer & groups$SIZE_CODE == size)
    subset_idx = which(groups$POLYMER == polymer & groups$MASS_mg == mass)
    if (length(subset_idx) < 2) next  # skip if less than two groups
    
    interval = unique(groups$SIZE_INTERVALS_mm[subset_idx])
    
    # cat("\n\n###### RESULTS FOR MASS:", mass, "mg | SIZE:", size,
    # cat("\n\n###### RESULTS FOR:", polymer,"| SIZE:", size,
    cat("\n\n###### RESULTS FOR:", polymer,"| MASS:", mass, "mg ######")
        # paste0("(", interval, " mm) ######"))
    
    # generate pairwise combinations
    idx_combos_subset = t(combn(subset_idx, 2))
    subset_results = list()
    
    for (name in appr) {
      cat("\n---------------", appr_dict[[name]], "---------------\n")
      approach_results = list()
      metric_matrix = get(name)
      
      for (k in seq_len(nrow(idx_combos_subset))) {
        i = idx_combos_subset[k, 1]
        j = idx_combos_subset[k, 2]
        
        label_i = paste0(groups$POLYMER[i], "_", groups$SIZE_CODE[i], groups$MASS_mg[i])
        label_j = paste0(groups$POLYMER[j], "_", groups$SIZE_CODE[j], groups$MASS_mg[j])
        comparison = paste(label_i, "VS", label_j)
        
        approach_results[[comparison]] = round(metric_matrix[i, j], 5)
      }
      
      subset_results[[appr_dict[name]]] = approach_results
      for (comp_name in names(approach_results)) {
        cat(comp_name, ":", approach_results[[comp_name]], "\n")
      }
    }
    # results[[paste0("Mass_", mass, "_Size_", size)]] = subset_results
    # results2[[paste0("Polymer_", polymer, "_Size_", size)]] = subset_results
    results3[[paste0("Polymer_", polymer, "_Mass_", mass)]] = subset_results
  }# WATCH OUT FOR RESULT STORAGE AS WELL
}

#==============================================================================#
# tidy up...

library(purrr)
library(dplyr)
library(tidyr)

# tidy_results = results |>
# tidy_results2 = results2 |>
tidy_results3 = results3 |>
  imap_dfr(function(context_data, context_name) {
    # For each metric in this context
    imap_dfr(context_data, function(matrix_data, metric_name) {
      # Convert matrix to data frame and add row names as a column
      df = as.data.frame(matrix_data)
      df$Row = rownames(df)  # Ensure row names are preserved
      
      df |> 
        pivot_longer(
          cols = -Row, 
          names_to = "Column", 
          values_to = "Value"
        ) |> 
        filter(Row < Column) |>   
        mutate(
          Context = context_name,
          Metric = metric_name
        ) |> 
        select(Context, Metric, Row, Column, Value)
    })
  })

## ranking dissimilarities
# tidy_results |>
# tidy_results2 |>
tidy_results3 |>
  group_by(Column, Metric) |>
  summarise(MetricMean = mean(Value), .groups = "drop") |>
  tidyr::pivot_wider(names_from = Metric, values_from = MetricMean) |>
  arrange(
    desc(`Mahalanobis Distance`),      # Primary: Best for finding outliers
    desc(`Euclidean Distance`),        # Secondary: Overall magnitude/shape difference
    desc(`Spectral Angler Mapper`),     # Tertiary: Pure shape difference (intensity-invariant)
    desc(`Correlation Dissimilarity`), # Quaternary: Global shape difference
    desc(`Moving Window Correlation Dissimilarity`) # Quinary: Localized shape differences
  ) |> 
  select(Column,
         `Mahalanobis Distance`,
         `Euclidean Distance`,
         `Spectral Angler Mapper`,
         `Correlation Dissimilarity`,
         `Moving Window Correlation Dissimilarity`)
#===============================================================================

# heatmap to help visualization
require(gplots)

labels = paste0(groups$POLYMER, "_", groups$SIZE_CODE, groups$MASS_mg)
rownames(mean_spectra_matrix) = labels

rownames(EucD) = colnames(EucD) = labels
heatmap.2(EucD,
          symm = TRUE,
          dendrogram = "none",
          Rowv = FALSE,
          Colv = FALSE,
          cexRow = 0.8,
          adjRow = c(0.2, 0.5),
          cexCol = 0.9,
          adjCol = c(0.8, 0.5),
          key = TRUE,
          keysize = 1,
          key.title = NA,
          key.par=list(mgp=c(1, 0.5, 0),
                       mar=c(5, 2, 1.8, 1)),
          key.xlab = "Euclidean Distance",
          trace = "none",
          density.info = "none",
          col = heat.colors(256, rev = T),
          margins = c(5, 5))

rownames(mahD) = colnames(mahD) = labels
heatmap.2(mahD,
          symm = TRUE,
          dendrogram = "none",
          Rowv = FALSE,
          Colv = FALSE,
          cexRow = 0.8,
          adjRow = c(0.2, 0.5),
          cexCol = 0.9,
          adjCol = c(0.8, 0.5),
          key = TRUE,
          keysize = 1,
          key.title = NA,
          key.par=list(mgp=c(1, 0.5, 0),
                       mar=c(5, 2, 1.8, 1)),
          key.xlab = "Mahalanobis Distance",
          trace = "none",
          density.info = "none",
          col = heat.colors(256, rev = T),
                            # , rev = TRUE), # it is called Mahalanobis DISTANCE...
          margins = c(5, 5))

rownames(cd1) = colnames(cd1) = labels
heatmap.2(cd1,
          symm = TRUE,
          dendrogram = "none",
          Rowv = FALSE,
          Colv = FALSE,
          cexRow = 0.8,
          adjRow = c(0.2, 0.5),
          cexCol = 0.9,
          adjCol = c(0.8, 0.5),
          key = TRUE,
          keysize = 1,
          key.title = NA,
          key.par=list(mgp=c(1, 0.5, 0),
                       mar=c(5, 2, 1.8, 1)),
          key.xlab = "Correlation Dissimilarity",
          trace = "none",
          density.info = "none",
          col = heat.colors(256, rev = T),
          margins = c(5, 5))

rownames(mwcd) = colnames(mwcd) = labels
heatmap.2(mwcd,
          symm = TRUE,
          dendrogram = "none",
          Rowv = FALSE,
          Colv = FALSE,
          cexRow = 0.8,
          adjRow = c(0.2, 0.5),
          cexCol = 0.9,
          adjCol = c(0.8, 0.5),
          key = TRUE,
          keysize = 1,
          key.title = NA,
          key.par=list(mgp=c(1, 0.5, 0),
                       mar=c(5, 2, 1.8, 1)),
          key.xlab = "Mov. Window Corr. Dissim.",
          trace = "none",
          density.info = "none",
          col = heat.colors(256, rev = T),
          margins = c(5, 5))

rownames(samD) = colnames(samD) = labels
heatmap.2(samD,
          symm = TRUE,
          dendrogram = "none",
          Rowv = FALSE,
          Colv = FALSE,
          cexRow = 0.8,
          adjRow = c(0.2, 0.5),
          cexCol = 0.9,
          adjCol = c(0.8, 0.5),
          key = TRUE,
          keysize = 1,
          key.title = NA,
          key.par=list(mgp=c(1, 0.5, 0),
                       mar=c(5, 2, 1.8, 1)),
          key.xlab = "Cosine Distance",
          trace = "none",
          density.info = "none",
          col = heat.colors(256, rev = T),
          margins = c(5, 5))


  #==============================================================================#
#                           SPECTRA SIMILARITY ANALYSIS                        #
#                              [pristine vs. colored]                          #
#==============================================================================#
groups_ref = data2 |> 
  group_by(POLYMER, SIZE_CODE, SIZE_INTERVALS_mm, MASS_mg) |> 
  summarise(mean_spectrum = list(colMeans(spcA)), .groups = "drop")

mean_spectra_ref = do.call(rbind, groups_ref$mean_spectrum)
rownames(mean_spectra_ref) = paste0(groups_ref$POLYMER, "_", 
                                    groups_ref$SIZE_CODE, 
                                    groups_ref$MASS_mg)


groups_new = new2 |> 
  group_by(POLYMER, SIZE_CODE, SIZE_INTERVALS_mm, MASS_mg) |> 
  summarise(mean_spectrum = list(colMeans(spcA)), .groups = "drop")

mean_spectra_new = do.call(rbind, groups_new$mean_spectrum)
rownames(mean_spectra_new) = paste0("ALL_",
                                    groups_new$SIZE_CODE,
                                    groups_new$MASS_mg)

# compute Euclidian distance (CROSS)
EucD_cross = f_diss(Xr = mean_spectra_ref,
                    Xu = mean_spectra_new,
                    diss_method = "euclid",
                    center = TRUE, scale = TRUE)


# compute Mahalanobis distance (CROSS)
PCA_ref = pc_projection(mean_spectra_ref,
                        pc_selection = list("cumvar", 0.99),
                        method = "pca",
                        center = TRUE, SCALE = FALSE)

X_new_centered = sweep(mean_spectra_new, 2, PCA_ref$center, FUN = "-")
X_new_scaled  = sweep(X_new_centered, 2, PCA_ref$scale, FUN = "/")

proj_new = X_new_scaled %*% t(PCA_ref$X_loadings)

mahD_cross = f_diss(Xr = PCA_ref$scores,
                    Xu = proj_new,
                    diss_method = "mahalanobis",
                    center = FALSE, scale = FALSE)


# correlation similarity (CROSS)
scorr_cross = cor(t(mean_spectra_ref), t(mean_spectra_new))
# compute similarity
cd1_cross = (1 - scorr_cross)/2


# moving window for correlation similarity (CROSS)
mwcd_cross = cor_diss(mean_spectra_ref,
                      mean_spectra_new,
                      ws = 51,
                      center = FALSE,
                      scale = FALSE)


# spectral angle mapper (CROSS)
samD_cross = f_diss(mean_spectra_ref,
                    mean_spectra_new,
                    center = FALSE, scale = FALSE)


appr_cross = c("EucD_cross", "mahD_cross", "cd1_cross",
               "mwcd_cross", "samD_cross")

appr_dict_cross = c(
  "EucD_cross" = "Euclidean Distance",
  "mahD_cross" = "Mahalanobis Distance",
  "cd1_cross" = "Correlation Dissimilarity",
  "mwcd_cross" = "Moving Window Correlation Dissimilarity",
  "samD_cross" = "Spectral Angler Mapper")


results_cross = list()

for (mass in sort(unique(groups_new$MASS_mg), decreasing = TRUE)) {
  for (size in unique(groups_new$SIZE_CODE)) {
    
    # find the EXTERNAL group name for this mass/size
    external_group_name = paste0("ALL_", size, mass)
    external_col_index = which(rownames(mean_spectra_new) == external_group_name)
    
    # find the INTERNAL group names for this mass/size
    internal_group_names = paste0(c("PE_", "PET_", "PP_", "PVC_"), size, mass)
    # find their row indices in the full mean_spectra_ref matrix
    internal_row_indices = which(rownames(mean_spectra_ref) %in% internal_group_names)
    
    # if either is not found, skip
    if (length(external_col_index) == 0 || length(internal_row_indices) == 0) next
    
    cat("\n\n###### CROSS RESULTS FOR MASS:", mass, "mg | SIZE:", size, "######\n")
    
    subset_results = list()
    
    for (name in appr_cross) {
      cat("\n---------------", appr_dict_cross[[name]], "---------------\n")
      metric_matrix = get(name)
      
      # for each internal polymer, get its dissimilarity to the external mixture
      for (i in seq_along(internal_row_indices)) {
        row_idx = internal_row_indices[i]
        internal_name = rownames(mean_spectra_ref)[row_idx]
        
        diss_value = round(metric_matrix[row_idx, external_col_index], 5)
        cat(internal_name, "VS", external_group_name, ":", diss_value, "\n")
        
        # store the result if needed
        subset_results[[appr_dict_cross[name]]][[paste0(internal_name, " VS ", external_group_name)]] = diss_value
      }
    }
    results_cross[[paste0("Mass_", mass, "_Size_", size)]] = subset_results
  }
}

#==============================================================================#
# tidy up, again...
require(purrr)
require(tidyr)

results_cross |> 
  imap_dfr(function(context_data, context_name) {
    
    imap_dfr(context_data, function(metric_vector, metric_name) {
      
      tibble(
        Column = names(metric_vector),
        Metric = metric_name,
        Value  = unlist(metric_vector)
      )
    })
  }) |> 
  pivot_wider(names_from = Metric, values_from = Value) |> 
  arrange(
    desc(`Mahalanobis Distance`),
    desc(`Euclidean Distance`),
    desc(`Spectral Angler Mapper`),
    desc(`Correlation Dissimilarity`),
    desc(`Moving Window Correlation Dissimilarity`)
  ) |> 
  select(Column,
         `Mahalanobis Distance`,
         `Euclidean Distance`,
         `Spectral Angler Mapper`,
         `Correlation Dissimilarity`,
         `Moving Window Correlation Dissimilarity`)
#===============================================================================

# heatmaps
labels_ref = paste0(groups_ref$POLYMER, "_", groups_ref$SIZE_CODE, "_", groups_ref$MASS_mg)
labels_new = paste0("ALL_", groups_new$SIZE_CODE, "_", groups_new$MASS_mg)

rownames(EucD_cross) = labels_ref
colnames(EucD_cross) = labels_new
heatmap.2(EucD_cross,
          symm = TRUE,
          dendrogram = "none",
          Rowv = FALSE,
          Colv = FALSE,
          cexRow = 0.8,
          adjRow = c(0.2, 0.5),
          cexCol = 0.9,
          adjCol = c(0.8, 0.5),
          key = TRUE,
          keysize = 1,
          key.title = NA,
          key.par=list(mgp=c(1, 0.5, 0),
                       mar=c(5, 2, 1.8, 1)),
          key.xlab = "Euclidean Distance",
          trace = "none",
          density.info = "none",
          col = heat.colors(256, rev = T),
          margins = c(8, 8))


rownames(mahD_cross) = labels_ref
colnames(mahD_cross) = labels_new
heatmap.2(mahD_cross,
          symm = TRUE,
          dendrogram = "none",
          Rowv = FALSE,
          Colv = FALSE,
          cexRow = 0.8,
          adjRow = c(0.2, 0.5),
          cexCol = 0.9,
          adjCol = c(0.8, 0.5),
          key = TRUE,
          keysize = 1,
          key.title = NA,
          key.par=list(mgp=c(1, 0.5, 0),
                       mar=c(5, 2, 1.8, 1)),
          key.xlab = "Mahalanobis Distance",
          trace = "none",
          density.info = "none",
          col = heat.colors(256, rev = T),
          margins = c(8, 8))

rownames(cd1_cross) = labels_ref
colnames(cd1_cross) = labels_new
heatmap.2(cd1_cross,
          symm = TRUE,
          dendrogram = "none",
          Rowv = FALSE,
          Colv = FALSE,
          cexRow = 0.8,
          adjRow = c(0.2, 0.5),
          cexCol = 0.9,
          adjCol = c(0.8, 0.5),
          key = TRUE,
          keysize = 1,
          key.title = NA,
          key.par=list(mgp=c(1, 0.5, 0),
                       mar=c(5, 2, 1.8, 1)),
          key.xlab = "Correlation Dissimilarity",
          trace = "none",
          density.info = "none",
          col = heat.colors(256, rev = T),
          margins = c(8, 8))

rownames(mwcd_cross) = labels_ref
colnames(mwcd_cross) = labels_new
heatmap.2(mwcd_cross,
          symm = TRUE,
          dendrogram = "none",
          Rowv = FALSE,
          Colv = FALSE,
          cexRow = 0.8,
          adjRow = c(0.2, 0.5),
          cexCol = 0.9,
          adjCol = c(0.8, 0.5),
          key = TRUE,
          keysize = 1,
          key.title = NA,
          key.par=list(mgp=c(1, 0.5, 0),
                       mar=c(5, 2, 1.8, 1)),
          key.xlab = "Mov. Window Corr. Dissim.",
          trace = "none",
          density.info = "none",
          col = heat.colors(256, rev = T),
          margins = c(8, 8))

rownames(samD_cross) = labels_ref
colnames(samD_cross) = labels_new
heatmap.2(samD_cross,
          symm = TRUE,
          dendrogram = "none",
          Rowv = FALSE,
          Colv = FALSE,
          cexRow = 0.8,
          adjRow = c(0.2, 0.5),
          cexCol = 0.9,
          adjCol = c(0.8, 0.5),
          key = TRUE,
          keysize = 1,
          key.title = NA,
          key.par=list(mgp=c(1, 0.5, 0),
                       mar=c(5, 2, 1.8, 1)),
          key.xlab = "Cosine Distance",
          trace = "none",
          density.info = "none",
          col = heat.colors(256, rev = T),
          margins = c(8, 8))

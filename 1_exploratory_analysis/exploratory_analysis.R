require(factoextra)
require(tripack)
require(splancs)
require(dplyr)
require(ggplot2)
require(RColorBrewer)
require(patchwork)
require(effectsize)
require(resemble)
require(purrr)
require(tidyr)
require(gplots)


# load data
## raw for PCA
pristine = readRDS("../raw_spectra/raw_pristine.rds")
colored = readRDS("../raw_spectra/raw_colored.rds")

## denoised for dissimilarity
pristine2 = readRDS("../preprocessed_data/pristine_denoised.rds")
colored2 = readRDS("../preprocessed_data/colored_denoised.rds")

# convert spectra to absorbance
pristine$spcA = log(1/pristine$spc)
colored$spcA = log(1/colored$spc)

pristine2$spcA = log(1/pristine2$spc)
colored2$spcA = log(1/colored2$spc)

################################################################################
#                         PRINCIPAL COMPONENT ANALYSIS                         #
################################################################################
pcspec = prcomp(pristine$spcA)
summary(pcspec)[[6]][, 1:10] # good amount of variance explained by it...

# eigenvalues vs. number of dimensions (screeplot)
fviz_eig(pcspec, addlabels = T, geom = "bar",
         ncp = 3, main = "", ggtheme = theme_gray())

# scores
plot(pcspec$x[, 1],
     pcspec$x[, 2],
     xlab = "PC1",
     ylab = "PC2",
     pch = 16,
     col = rgb(0.5, 0.5, 0.5, alpha = 0.5))

# first loading
plot(colnames(pristine$spcA), pcspec$rotation[,1],
     type ="l",
     ylab ="Loading",
     xlab ="Wavelength (nm)",
     col = rgb(red = 1, green = 0, blue = 0, alpha = 1),
     ylim = c(-0.06, 0.06))

# second loading
lines(colnames(pristine$spcA), pcspec$rotation[,2],
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
plot(colnames(pristine$spcA), pcspec$rotation[,ncol(pcspec$rotation)],
     main = "Last PC loading (mostly noise)",
     type = "l",
     ylab = "Loading",
     xlab = "Wavelength (nm)",
     col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.5),
     ylim = c(-0.1, 0.1))


#==============================================================================#
#                      CHECK SPECTRAL PREDICTION DOMAIN                        #
#==============================================================================#
# PCA scores triangulation
randTr = tri.mesh(pcspec$x[, 1], pcspec$x[, 2])

# nodes of the convex hull
randCH = convex.hull(randTr, plot.it = F)

plot(pcspec$x[,1], pcspec$x[,2],
     xlab ="PC1",
     ylab ="PC2",
     ylim = c(-20, 20),
     pch = 16,
     col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.5))

lines(c(randCH$x, randCH$x[1]), c(randCH$y,randCH$y[1]),
      col="darkred",
      lwd=2)


# how much of colored data are inside the predictive space?
colored$spcA = as.matrix(colored$spcA)
new_scores = sweep(colored$spcA, 2, pcspec$center, "-") %*% pcspec$rotation

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

colored |> 
  filter(!inside) |> 
  select(SAMPLE_ID, SAMPLE_CODE, SIZE_INTERVALS_mm,
         MASS_mg, REPLICATE)

paste0(round((sum(inside == FALSE) * 100) / nrow(colored), 2), "% ",
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


## check data structure =======================================================#

# extract scores and loadings
scores_df = as.data.frame(pcspec$x)
loadings_df = as.data.frame(pcspec$rotation)

scores_df$Obs = rownames(scores_df)
loadings_df$Var = rownames(loadings_df)

# select variables (some of the wavelengths)
select_var = c(min(as.numeric(colnames(pristine$spcA))),
               median(as.numeric(colnames(pristine$spcA))),
               max(as.numeric(colnames(pristine$spcA))))

# filter loadings
# filt_loading = loadings_df[loadings_df$Var %in% select_var, ]

# select representative wavelengths for Vis-NIR-SWIR ranges
filt_loading = loadings_df[loadings_df$Var %in% c("350", "700", "1100", "2500"), ]


# biplot
p = ggplot() +
  geom_point(data = scores_df,
             aes(x = PC1, y = PC2,
                 color = "Samples"),
             alpha = 0.5,
             size = 3) +
  labs(x = "PC1 (77.2%)", y = "PC2 (19.7%)") + # hardcoded here, can be pulled from `summary(pcspec)`
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

# spectra distribution (grouped by mass)
p3 = fviz_pca_ind(pcspec, geom = "point",
             title = "",
             pointshape = 16,
             pointsize = 2, # smaller is better for panel visualization...
             alpha.ind = 0.6,
             col.ind = as.factor(pristine$MASS_mg),
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
             col.ind = pristine$SIZE_INTERVALS_mm,
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
             col.ind = pristine$POLYMER,
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


# confidence ellipses (grouped by mass)
p4 = fviz_pca_ind(pcspec, geom = "point",
                  title = "",
                  pointshape = 16,
                  pointsize = 4,
                  alpha.ind = 0.6,
                  col.ind = as.factor(pristine$MASS_mg),
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


# compose panel
p/(p1+p2+p3)/p4 +
  plot_annotation(tag_levels = 'a',
                  tag_prefix = '(', tag_suffix = ')') &
  theme(plot.tag = element_text(face = 'bold', size = 20))


# BEYOND VISUAL INSPECTION...
## check MASS_mg correlation with PC1 and PC2
cor.test(pcspec$x[ ,1], pristine$MASS_mg)
cor.test(pcspec$x[ ,2], pristine$MASS_mg)

## AOV for POLYMER with PC1 and PC2
summary(aov(pcspec$x[ ,1] ~ pristine$POLYMER))
summary(aov(pcspec$x[ ,2] ~ pristine$POLYMER))

## size effect in PC1
eta_squared(aov(pcspec$x[ ,1] ~ pristine$MASS_mg))
eta_squared(aov(pcspec$x[ ,1] ~ pristine$POLYMER))

## size effect in PC2
eta_squared(aov(pcspec$x[ ,2] ~ pristine$MASS_mg))
eta_squared(aov(pcspec$x[ ,2] ~ pristine$POLYMER))


################################################################################
#                        SPECTRA DISSIMILARITY ANALYSIS                        #
################################################################################

# calculate mean absorbance spectra for pristine...
groups = pristine2 |>
  group_by(POLYMER, SIZE_CODE, SIZE_INTERVALS_mm, MASS_mg) |>
  summarise(mean_spectrum = list(colMeans(spcA)), .groups = "drop")

mean_spectra_matrix = do.call(rbind, groups$mean_spectrum)

# compute Euclidian distance
EucD = f_diss(Xr = mean_spectra_matrix,
              Xu = mean_spectra_matrix,
              diss_method = "euclid",
              center = TRUE, scale = TRUE)


# compute Mahalanobis distance
PCA = pc_projection(mean_spectra_matrix,
                    pc_selection = list("cumvar", 0.99),
                    method = "pca",
                    center = TRUE, scale = FALSE)

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
  "samD" = "Spectral Angle Mapper")


#==============================================================================#
#         CONFIG - change these to switch what is controlling vs. varying      #
#==============================================================================#
CONTROL_VARS = c("MASS_mg", "SIZE_CODE") # CONTEXT_i = c("MASS_mg", "SIZE_CODE")
                                         # CONTEXT_ii = c("POLYMER", "SIZE_CODE")
                                         # CONTEXT_iii = c("MASS_mg", "POLYMER")
VARY_LABEL = function(row) paste0(row$POLYMER, "_", row$SIZE_CODE, row$MASS_mg)


# run all approach comparisons for one subset of groups...
compare_subset = function(subset_idx, groups, appr, appr_dict) {
  idx_combos = t(combn(subset_idx, 2))
  
  map(set_names(appr), function(name) {
    metric_matrix = get(name)
    
    apply(idx_combos, 1, function(pair) {
      i = pair[1]; j = pair[2]
      label_i = VARY_LABEL(groups[i, ])
      label_j = VARY_LABEL(groups[j, ])
      round(metric_matrix[i, j], 5)
    }) |>
      set_names(apply(idx_combos, 1, function(pair) {
        paste(VARY_LABEL(groups[pair[1], ]), "vs", VARY_LABEL(groups[pair[2], ]))
      }))
  }) |>
    set_names(unlist(appr_dict[appr]))
}


# iterate over unique control-var combinations...
control_combos = groups |>
  select(all_of(CONTROL_VARS)) |>
  distinct() |>
  mutate(across(everything(), as.character))

results = pmap(control_combos, function(...) {
  vals = list(...)
  filters = map2(CONTROL_VARS, vals, ~ groups[[.x]] == .y)
  subset_idx = which(Reduce(`&`, filters))
  
  if (length(subset_idx) < 2) return(NULL)
  
  interval = unique(groups$SIZE_INTERVALS_mm[subset_idx])
  context_label = paste(map2_chr(CONTROL_VARS, vals, ~ paste0(.x, "_", .y)), collapse = "_")
  
  cat("\n\n######", paste(CONTROL_VARS, as.character(unlist(vals)), sep = ": ", collapse = " | "),
      paste0("(", interval, " mm) ######\n"))
  
  result = compare_subset(subset_idx, groups, appr, appr_dict)
  
  # print results
  iwalk(result, function(comparisons, metric_name) {
    cat("---------------", metric_name, "---------------\n")
    iwalk(comparisons, ~ cat(.y, ":", .x, "\n"))
  })
  
  result
}) |>
  set_names(apply(control_combos, 1, paste, collapse = "_")) |>
  compact()  # drop NULLs (skipped subsets)


# tidy up and rank...
tidy_results = imap_dfr(results, function(context_data, context_name) {
  imap_dfr(context_data, function(comparisons, metric_name) {
    tibble(
      Context    = context_name,
      Metric     = metric_name,
      Comparison = names(comparisons),
      Value      = unlist(comparisons)
    )
  })
})

METRIC_ORDER = c(
  "Mahalanobis Distance",
  "Euclidean Distance",
  "Spectral Angle Mapper",
  "Correlation Dissimilarity",
  "Moving Window Correlation Dissimilarity"
)

# ready to save...
tidy_results |>
  group_by(Comparison, Metric) |>
  summarise(MetricMean = mean(Value), .groups = "drop") |>
  pivot_wider(names_from = Metric, values_from = MetricMean) |>
  arrange(across(all_of(METRIC_ORDER), desc)) |>
  select(Comparison, all_of(METRIC_ORDER))
# |> write.csv(file = "CONTEXT_n.csv, row.names = F)

#==============================================================================#

# heatmap to help visualization
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
#                               pristine vs. colored                           #
#==============================================================================#
# calculate mean absorbance spectra for pristine...
groups_ref = pristine2 |> 
  group_by(POLYMER, SIZE_CODE, SIZE_INTERVALS_mm, MASS_mg) |> 
  summarise(mean_spectrum = list(colMeans(spcA)), .groups = "drop")

mean_spectra_ref = do.call(rbind, groups_ref$mean_spectrum)
rownames(mean_spectra_ref) = paste0(groups_ref$POLYMER, "_", 
                                    groups_ref$SIZE_CODE, 
                                    groups_ref$MASS_mg)

# calculate mean absorbance spectra for colored...
groups_new = colored2 |> 
  group_by(POLYMER, SIZE_CODE, SIZE_INTERVALS_mm, MASS_mg) |> 
  summarise(mean_spectrum = list(colMeans(spcA)), .groups = "drop")

mean_spectra_new = do.call(rbind, groups_new$mean_spectrum)
rownames(mean_spectra_new) = paste0("MIX_",
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
                        center = TRUE, scale = FALSE)

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
                    diss_method = "cosine",
                    center = FALSE, scale = FALSE)


appr_cross = c("EucD_cross", "mahD_cross", "cd1_cross",
               "mwcd_cross", "samD_cross")

appr_dict_cross = c(
  "EucD_cross" = "Euclidean Distance",
  "mahD_cross" = "Mahalanobis Distance",
  "cd1_cross" = "Correlation Dissimilarity",
  "mwcd_cross" = "Moving Window Correlation Dissimilarity",
  "samD_cross" = "Spectral Angle Mapper")


results_cross = list()

for (mass in sort(unique(groups_new$MASS_mg), decreasing = TRUE)) {
  for (size in unique(groups_new$SIZE_CODE)) {
    
    # find the EXTERNAL group name for this mass/size
    external_group_name = paste0("MIX_", size, mass)
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
        cat(internal_name, "vs", external_group_name, ":", diss_value, "\n")
        
        # store the result if needed
        subset_results[[appr_dict_cross[name]]][[paste0(internal_name, " vs ", external_group_name)]] = diss_value
      }
    }
    results_cross[[paste0("Mass_", mass, "_Size_", size)]] = subset_results
  }
}


# tidy up...
results_cross |> 
  imap_dfr(function(context_data, context_name) {
    
    imap_dfr(context_data, function(metric_vector, metric_name) {
      
      tibble(
        Comparison = names(metric_vector),
        Metric = metric_name,
        Value  = unlist(metric_vector)
      )
    })
  }) |> 
  pivot_wider(names_from = Metric, values_from = Value) |> 
  arrange(
    desc(`Mahalanobis Distance`),
    desc(`Euclidean Distance`),
    desc(`Spectral Angle Mapper`),
    desc(`Correlation Dissimilarity`),
    desc(`Moving Window Correlation Dissimilarity`)
  ) |> 
  select(Comparison,
         `Mahalanobis Distance`,
         `Euclidean Distance`,
         `Spectral Angle Mapper`,
         `Correlation Dissimilarity`,
         `Moving Window Correlation Dissimilarity`)
# |> write.csv(file = "CONTEXT_iv.csv, row.names = F)

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
          key.par = list(mgp = c(1, 0.5, 0),
                         mar = c(5, 2, 1.8, 1)),
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
          key.par = list(mgp = c(1, 0.5, 0),
                         mar = c(5, 2, 1.8, 1)),
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
          key.par =list (mgp = c(1, 0.5, 0),
                         mar = c(5, 2, 1.8, 1)),
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
          key.par = list(mgp = c(1, 0.5, 0),
                         mar = c(5, 2, 1.8, 1)),
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
          key.par = list(mgp = c(1, 0.5, 0),
                         mar = c(5, 2, 1.8, 1)),
          key.xlab = "Cosine Distance",
          trace = "none",
          density.info = "none",
          col = heat.colors(256, rev = T),
          margins = c(8, 8))
# repetitive because that's a lot of info and each may need different editing... 

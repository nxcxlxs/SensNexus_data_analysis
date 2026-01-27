require(ggplot2)
require(dplyr)
require(patchwork)

validation_results = data.frame(
  Model = rep(c("Full preprocessing",
                "Minimal preprocessing",
                "In-between preprocessing",
                "Raw data"),
              times = 4),  # 4 models × 4 validation types
  Validation_Type = rep(c("75/25 Hold-out",
                          "k-fold CV",
                          "External (naive)",
                          "External (tuned)"), each = 12),  # 3 metrics × 4 models = 12 per validation type
  Metric = rep(rep(c("ME", "RMSE", "R²"), times = 4), times = 4),  # 3 metrics per model
  Value = c(
    # 75/25 Hold-out validation
    -51.6718288, 368.3188416, 0.9259154,  # MODEL1
    34.4806711, 375.0805673, 0.9231703,   # MODEL2
    31.7480146, 368.6132237, 0.9257969,   # MODEL3
    29.2771195, 679.8320883, 0.7476034,   # MODEL4
    
    # k-fold CV (using your 4-fold CV values)
    -21.2562, 542.8593, 0.8391,           # MODEL1
    15.4975, 512.5971, 0.8565,            # MODEL2
    0.2660, 513.3000, 0.8561,             # MODEL3
    -24.8231, 573.7132, 0.8202,           # MODEL4
    
    # External "naive"
    -681.2096820, 1374.0183875, 0.0796633,    # MODEL1
    -791.5640335, 1523.9561530, -0.1321568,   # MODEL2 (negative R²)
    -622.7435488, 1287.4521536, 0.1919768,    # MODEL3
    -666.2196153, 1336.0768800, 0.1297891,    # MODEL4
    
    # External "tuned-on-external"
    -635.4124750, 1246.4591768, 0.2426132,    # MODEL1
    -731.5591977, 1392.2814677, 0.0550350,    # MODEL2
    -622.7435488, 1287.4521536, 0.1919768,    # MODEL3
    -666.2196153, 1336.0768800, 0.1297891     # MODEL4
  )
)

validation_results$Metric = factor(validation_results$Metric, 
                                   levels = c("ME", "RMSE", "R²"))

validation_results$Model = factor(validation_results$Model,
                                  levels = c("Raw data",
                                             "Minimal preprocessing",
                                             "In-between preprocessing",
                                             "Full preprocessing"))

validation_results$Validation_Type = factor(validation_results$Validation_Type,
                                            levels = c("75/25 Hold-out",
                                                       "k-fold CV",
                                                       "External (naive)",
                                                       "External (tuned)"))

# validation_results_filtered = validation_results |> 
#   filter(Validation_Type != "External (naive)")

ggplot(validation_results,
       aes(x = Model, y = Value, group = 1)) +
  geom_line(size = 1.2, color = "darkgreen") +
  geom_point(size = 1.5, color = "darkgreen") +
  geom_text(aes(label = round(Value, 2)),
            vjust = -1,
            size = 3, color = "darkgreen") +
  facet_grid(Metric ~ Validation_Type, scales = "free_y") +
  scale_y_continuous(position = "left") +
  labs(title = "Model Performance Trends",
       x = "Preprocessing Approach",
       y = "Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 14))

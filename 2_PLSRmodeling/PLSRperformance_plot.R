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
    66.6020032, 354.3433093, 0.9314309,     # MODEL1
    56.4854015, 430.3879971, 0.8988419,     # MODEL2
    79.5472914, 362.2250690, 0.9283465,     # MODEL3
    59.4212943, 427.6531738, 0.9001234,     # MODEL4
    
    # k-fold CV
    -5.8456, 429.0567, 0.8995,              # MODEL1
    0.2498, 451.2568, 0.8888,               # MODEL2
    3.2399, 418.9342, 0.9042,               # MODEL3
    1.7271, 454.4206, 0.8872,               # MODEL4
    
    # External "naive"
    -2650.4705204, 5719.8711412, -14.9490247,  # MODEL1
    -580.6311103, 918.6611902, 0.5885923,      # MODEL2
    -643.0621647, 1606.8055733, -0.2586017,    # MODEL3
    -667.5676297, 1627.9799779, -0.2919918,    # MODEL4
    
    # External "tuned-on-external"
    -229.0696, 695.7013, 0.7640569,         # MODEL1
    -187.2734, 628.5788, 0.807389,          # MODEL2
    -218.327, 711.8946, 0.7529454,          # MODEL3
    -570.8383, 762.579, 0.7165142           # MODEL4
  )
)

validation_results$Metric = factor(validation_results$Metric, 
                                    levels = c("ME", "RMSE", "R²"))

validation_results$Model = factor(validation_results$Model,
                                  levels = c("Raw data",
                                             "Minimal preprocessing",
                                             "In-between preprocessing",
                                             "Full preprocessing"))


# validation_results_filtered = validation_results |> 
#   filter(Validation_Type != "External (naive)")

ggplot(validation_results,
       aes(x = Model, y = Value, group = 1)) +
  geom_line(size = 1.2, color = "darkred") +
  geom_point(size = 1.5, color = "darkred") +
  geom_text(aes(label = round(Value, 2)),
            vjust = 1, hjust = 1.3,
            size = 3, color = "darkred") +
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


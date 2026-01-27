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
    -22.8830, 297.0873, 0.9518,     # MODEL1
    -13.0570, 307.6184, 0.9483,     # MODEL2
    -25.5752, 324.3925, 0.9425,     # MODEL3
    37.5988, 538.1751, 0.8418,      # MODEL4
    
    # k-fold CV
    -2.5077, 435.6986, 0.8963,      # MODEL1
    3.7972, 422.3500, 0.9026,       # MODEL2
    -1.8247, 434.2964, 0.8970,      # MODEL3
    4.9982, 426.7480, 0.9005,       # MODEL4
    
    # External "naive"
    -483.7766267, 1164.8982785, 0.3384881,  # MODEL1
    -528.8312800, 965.3670865, 0.5456959,   # MODEL2
    -520.3236000, 1181.3989517, 0.3196149,  # MODEL3
    -421.2462533, 802.4477389, 0.6860972,   # MODEL4
    
    # External "tuned-on-external"
    -471.0385067, 1093.2828086, 0.4173246,  # MODEL1
    -500.2920133, 912.4376097, 0.5941477,   # MODEL2
    -470.0807733, 1045.9738388, 0.4666611,  # MODEL3
    -377.7088267, 765.5373976, 0.7143104    # MODEL4
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
  geom_line(size = 1.2, color = "steelblue") +
  geom_point(size = 1.5, color = "steelblue") +
  geom_text(aes(label = round(Value, 2)),
            vjust = -1, hjust = 1,
            size = 3, color = "steelblue") +
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

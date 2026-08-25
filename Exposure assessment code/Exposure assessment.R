# Exposure assessment

library(mc2d)
library(ggplot2)
library(patchwork)

set.seed(123)

# Simulation settings
n_inner <- 1000
n_outer <- 500

ndvar(n_inner)
ndunc(n_outer)

# Input values
C_min <- 75
C_max <- 1100

M_shape <- 0.8760928
M_scale <- 43.276368

F_min <- 0.0018
F_max <- 0.0233

# Log-uniform distribution
rlogunif <- function(n, min, max) {
  exp(runif(n, log(min), log(max)))
}

# Input variability
C_retail_v <- mcstoc(
  rlogunif,
  type = "V",
  min = C_min,
  max = C_max
)

M_v <- mcstoc(
  rgamma,
  type = "V",
  shape = M_shape,
  scale = M_scale
)

F_pathogenic_v <- mcstoc(
  runif,
  type = "V",
  min = F_min,
  max = F_max
)

# Daily pathogenic dose
D_raw_v <- C_retail_v * M_v * F_pathogenic_v

# Model summary
exposure_model <- mc(
  C_retail_v,
  M_v,
  F_pathogenic_v,
  D_raw_v
)

print(exposure_model)

print(
  summary(
    exposure_model,
    probs = c(0.025, 0.50, 0.975)
  )
)

# Extract variability simulations
simulation_data <- data.frame(
  inner_id = seq_len(n_inner),
  C_retail = as.vector(unmc(C_retail_v, drop = TRUE)),
  M = as.vector(unmc(M_v, drop = TRUE)),
  F_pathogenic = as.vector(unmc(F_pathogenic_v, drop = TRUE)),
  D_raw = as.vector(unmc(D_raw_v, drop = TRUE))
)

# Save values for inspection
write.csv(
  simulation_data,
  "raw_shrimp_mc2d_simulation_data.csv",
  row.names = FALSE
)

# Plotting data
plot_data <- data.frame(
  Log10_dose = log10(simulation_data$D_raw)
)

plot_style <- theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 12)
  )

# Pathogenic-dose variability
dose_plot <- ggplot(plot_data, aes(Log10_dose)) +
  geom_histogram(
    binwidth = 0.10,
    boundary = 0,
    colour = "black",
    fill = "grey70"
  ) +
  labs(
    x = expression(
      Log[10] * " pathogenic dose (MPN per shrimp-consumption day)"
    ),
    y = "Frequency"
  ) +
  plot_style

print(dose_plot)

ggsave(
  "dose_variability_mc2d.png",
  dose_plot,
  width = 8,
  height = 6,
  dpi = 600,
  bg = "white"
)
#Data
shrimp <- read.csv("shrimp_consumption_days.csv")
x <- shrimp$shrimp_day_g
w <- shrimp$WTS_P * 630 / sum(shrimp$WTS_P)

# Lognormal
ln_fit <- optim(
  c(weighted.mean(log(x), w), log(sd(log(x)))), \(p) -sum(w * dlnorm(x, p[1], exp(p[2]), log = TRUE))
)

meanlog <- ln_fit$par[1]
sdlog <- exp(ln_fit$par[2])
ln_AIC <- 2 * ln_fit$value + 4

# Gamma
m <- weighted.mean(x, w)
v <- weighted.mean((x - m)^2, w)

gamma_fit <- optim(
  log(c(m^2 / v, m / v)), \(p) -sum(w * dgamma(x, exp(p[1]), rate = exp(p[2]), log = TRUE))
)

shape <- exp(gamma_fit$par[1])
scale <- 1 / exp(gamma_fit$par[2])
gamma_AIC <- 2 * gamma_fit$value + 4

results <- data.frame(
  distribution = c("Gamma", "Lognormal"),
  parameter_1 = c("shape", "meanlog"),
  estimate_1 = c(shape, meanlog),
  parameter_2 = c("scale", "sdlog"),
  estimate_2 = c(scale, sdlog),
  weighted_AIC = c(gamma_AIC, ln_AIC)
)

print(results, row.names = FALSE)

# Doxycycline Emax and EC50 estimation
library(dplyr)

# Data
dat <- read.csv("Antibiotic data with control.csv", check.names = FALSE)

# Starting values
rates <- dat %>%
  group_by(Concentration_mg_L) %>%
  summarise(rate = coef(lm(logCFU ~ Time_h))[["Time_h"]], .groups = "drop")

a_start <- mean(dat$logCFU[dat$Time_h == 0])
m0_start <- rates$rate[rates$Concentration_mg_L == 0]
Emax_start <- m0_start - min(rates$rate)

# Fit Emax model
emax_model <- nls(
  logCFU ~ a + (m0 - Emax * Concentration_mg_L /
                  (EC50 + Concentration_mg_L)) * Time_h,
  data = dat,
  start = list(a = a_start, m0 = m0_start,
               Emax = Emax_start, EC50 = 0.1)
)

summary(emax_model)

# Convert Emax to natural-log scale per day
pd_estimates <- coef(emax_model)[c("Emax", "EC50")]
pd_covariance <- vcov(emax_model)[c("Emax", "EC50"),
                                  c("Emax", "EC50")]

conversion <- c(24 * log(10), 1)
pd_estimates <- pd_estimates * conversion
pd_covariance <- diag(conversion) %*% pd_covariance %*% diag(conversion)
dimnames(pd_covariance) <- list(names(pd_estimates), names(pd_estimates))

# Parameter estimates, SDs and covariance
pd_estimates_with_sd <- cbind(
  Estimate = pd_estimates,
  SD = sqrt(diag(pd_covariance))
)

pd_estimates_with_sd
pd_covariance

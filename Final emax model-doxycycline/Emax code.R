# Doxycycline Emax and EC50 estimation

library(dplyr)

# Data
dat <- read.csv("Antibiotic data with control.csv", check.names = FALSE) 


# Starting values for the nls model
rates <- dat %>%
  group_by(Concentration_mg_L) %>%
  summarise(
    rate = coef(lm(logCFU ~ Time_h))[["Time_h"]],
    .groups = "drop"
  )

a_start <- mean(dat$logCFU[dat$Time_h == 0])

m0_start <- rates$rate[
  rates$Concentration_mg_L == 0
]

Emax_start <- m0_start - min(rates$rate)

EC50_start <- 0.1

# Fit the model
emax_model <- nls(
  logCFU ~ a + (m0 - Emax * Concentration_mg_L /
                  (EC50 + Concentration_mg_L)) * Time_h,
  data = dat,
  start = list(
    a = a_start,
    m0 = m0_start,
    Emax = Emax_start,
    EC50 = EC50_start
  )
)

summary(emax_model)

# Emax and EC50 estimates with SD
pd_sd <- sqrt(diag(pd_covariance))
pd_estimates_with_sd <- cbind(Estimate = pd_estimates, SD = pd_sd)

pd_estimates_with_sd
pd_covariance
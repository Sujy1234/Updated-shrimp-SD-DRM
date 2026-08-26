# SD-DRM model functions

# Numerical functions ----

# Stable calculation of log(1 - exp(-a))
log1mexpm <- function(a) {
  small <- a <= log(2); ans <- numeric(length(a))
  ans[small] <- log(-expm1(-a[small]))
  ans[!small] <- log1p(-exp(-a[!small]))
  ans
}

recycle_input <- function(x, n, name) {
  if (!(length(x) %in% c(1, n))) stop(paste0(name, " must have length 1 or ", n, "."))
  rep_len(x, n)
}


# Antibiotic-adjusted beta-Poisson parameters ----

# Convert r to the simple-death rate, add the PD effect, and refit r
fit_ab_effect <- function(mu, effect, tfs) {
  r_s <- -log1mexpm((mu + effect) * tfs)
  r_s <- pmin(pmax(r_s[is.finite(r_s)], .Machine$double.xmin), 1 - .Machine$double.eps)
  
  fit <- try(suppressWarnings(fitdistrplus::fitdist(r_s, "beta", method = "mle")), silent = TRUE)
  
  if (!inherits(fit, "try-error")) {
    pars <- unname(fit$estimate[c("shape1", "shape2")])
    if (all(is.finite(pars)) && all(pars > 0)) return(c(alpha_s = pars[1], beta_s = pars[2]))
  }
  
  m <- mean(r_s); v <- var(r_s); common <- m * (1 - m) / v - 1
  start <- if (!is.finite(common) || common <= 0) c(1, 1 / max(m, .Machine$double.xmin)) else pmax(c(m * common, (1 - m) * common), 1e-8)
  
  neg_loglik <- function(log_par) -sum(dbeta(r_s, exp(log_par[1]), exp(log_par[2]), log = TRUE))
  opt <- optim(log(start), neg_loglik, method = "BFGS", control = list(maxit = 10000))
  pars <- exp(opt$par)
  
  if (opt$convergence != 0 || any(!is.finite(pars)) || any(pars <= 0))
    stop(paste0("Beta MLE failed at PD effect = ", signif(effect, 6)))
  
  c(alpha_s = pars[1], beta_s = pars[2])
}


# Parameter lookup and interpolation ----

# Generate r once and refit alpha and beta across the PD-effect range
make_ab_lookup <- function(alpha, beta, Emax, EC50, C_max, tfs, n_grid = 250, seed = 0, nsim = 10000) {
  if (length(Emax) != length(EC50)) stop("Emax and EC50 must be paired vectors of equal length.")
  if (any(!is.finite(Emax)) || any(!is.finite(EC50)) || any(Emax < 0) || any(EC50 <= 0))
    stop("Invalid Emax or EC50 values.")
  
  set.seed(seed)
  r <- rbeta(nsim, alpha, beta)
  mu <- -log1mexpm(r) / tfs
  effect_max <- 1.001 * max(Emax * C_max / (EC50 + C_max))
  
  if (effect_max <= 0) return(data.frame(effect = 0, alpha_s = alpha, beta_s = beta))
  
  effect_grid <- c(0, exp(seq(log(effect_max * 1e-6), log(effect_max), length.out = n_grid - 1)))
  pars <- t(vapply(effect_grid[-1], function(x) fit_ab_effect(mu, x, tfs), numeric(2)))
  
  data.frame(effect = effect_grid,
             alpha_s = c(alpha, pars[, "alpha_s"]),
             beta_s = c(beta, pars[, "beta_s"]))
}

get_ab_from_lookup <- function(effect, ab_lookup) {
  upper_limit <- max(ab_lookup$effect)
  if (max(effect) > upper_limit * (1 + 1e-10))
    stop("The pharmacodynamic effect exceeds the lookup-table range.")
  
  effect <- pmin(pmax(effect, 0), upper_limit)
  
  data.frame(
    alpha_s = exp(approx(ab_lookup$effect, log(ab_lookup$alpha_s), xout = effect, rule = 2)$y),
    beta_s = exp(approx(ab_lookup$effect, log(ab_lookup$beta_s), xout = effect, rule = 2)$y)
  )
}


# SD-DRM illness risk ----

calc_sd_drm <- function(dose, fr, C, alpha, beta, Emax, EC50, tfs, ab_lookup) {
  n <- max(length(dose), length(fr), length(C), length(Emax), length(EC50))
  
  dose <- recycle_input(dose, n, "dose"); fr <- recycle_input(fr, n, "fr")
  C <- recycle_input(C, n, "C"); Emax <- recycle_input(Emax, n, "Emax")
  EC50 <- recycle_input(EC50, n, "EC50")
  
  effect <- ifelse(C == 0, 0, Emax * C / (EC50 + C))
  pars <- get_ab_from_lookup(effect, ab_lookup)
  
  Ns <- dose * (1 - fr); Nr <- dose * fr
  p_ext_s <- (1 + Ns / pars$beta_s)^(-pars$alpha_s)
  p_ext_r <- (1 + Nr / beta)^(-alpha)
  
  risk_total <- pmin(pmax(1 - p_ext_s * p_ext_r, 0), 1)
  risk_less_treatable <- pmin(pmax(1 - p_ext_r, 0), 1)
  risk_more_treatable <- pmin(pmax(p_ext_r * (1 - p_ext_s), 0), 1)
  delta_p_treat <- risk_more_treatable - risk_less_treatable
  status <- ifelse(delta_p_treat > 0, "More likely treatable", "Less likely treatable")
  
  data.frame(risk_total, delta_p_treat, status)
}
## =============================================================================
## Data generation
##
## Linear predictor (one measured covariate "climate", two confounders z1, z2):
##   eta_ij = beta0_j + beta_j * climate_i + delta1_j * z1_i + delta2_j * z2_i
## with
##   climate = 0.5 * X1 + N(0, 0.05^2)
##   z1 = -sX * X1 - eps1 ,  z2 = -sX * X1 + eps2
## as in Hui et al. (2024). X1 is the common (broad-scale) field, eps1 and eps2
## are the own parts of the confounders.
## =============================================================================

std <- function(v) (v - mean(v)) / sd(v)

#' Covariates for one dataset.
#' @param s A design row (list) with fields generator ("basis"/"matern"),
#'   x_scale ("coarse"/"fine"), sX, eps_sd.
#' @param S Output of setup_sites().
gen_covariates <- function(s, S) {
  N <- S$n_sites
  if (s$generator == "basis") {
    ## Hui et al.: coarse X uses phi_2..phi_5; fine X uses phi_11..phi_14
    if (s$x_scale == "coarse") {
      idx_X <- 2:5;   idx_e1 <- c(8, 12, 14); idx_e2 <- c(9, 13, 16)
    } else {
      idx_X <- 11:14; idx_e1 <- c(2, 4, 6);   idx_e2 <- c(3, 5, 7)
    }
    X1   <- drop(S$basis[, idx_X]  %*% rnorm(length(idx_X),  sd = 0.25))
    clim <- 0.5 * X1 + rnorm(N, sd = 0.05)
    eps1 <- drop(S$basis[, idx_e1] %*% rnorm(length(idx_e1), sd = s$eps_sd))
    eps2 <- drop(S$basis[, idx_e2] %*% rnorm(length(idx_e2), sd = s$eps_sd))
  } else {
    ## Matern (exponential) fields with the same marginal scales as the basis
    ## generator: sd(X1) ~ 0.5, sd(eps) ~ sqrt(3) * eps_sd. The own parts are
    ## residualised on X1 so that they are exactly orthogonal to it, as in the
    ## basis generator.
    L_X <- if (s$x_scale == "coarse") S$L_coarse else S$L_fine
    L_e <- if (s$x_scale == "coarse") S$L_fine   else S$L_coarse
    X1   <- 0.5 * std(drop(L_X %*% rnorm(N)))
    clim <- 0.5 * X1 + rnorm(N, sd = 0.05)
    own  <- function() sqrt(3) * s$eps_sd * std(resid(lm(drop(L_e %*% rnorm(N)) ~ X1)))
    eps1 <- own(); eps2 <- own()
  }
  z <- cbind(z1 = -s$sX * X1 - eps1, z2 = -s$sX * X1 + eps2)
  ## Regression of the confounders on the measured covariate (matrix A of Section 2)
  a <- coef(lm(z ~ clim))[2, ]
  list(climate = clim, z = z, X1 = X1, eps = cbind(eps1, eps2), a = a)
}

#' Species coefficients for one dataset (redrawn for every dataset).
#' beta_j (direct effect) and delta_j (responses to the confounders) are linked
#' through rho: rho is the correlation between beta_j and the species' combined
#' response to the confounders, (w1 + w2) / sqrt(2), where delta_jm = delta_mean + delta_sd * w_jm.
gen_coefficients <- function(s) {
  J  <- s$J
  w  <- matrix(rnorm(2 * J), J, 2)
  sc <- (w[, 1] + w[, 2]) / sqrt(2)
  es <- if (is.null(s$effect_scale)) 1 else s$effect_scale   # 0.5 for counts (see design.R)
  beta  <- es * (s$beta_mean + s$beta_sd * (s$rho * sc + sqrt(1 - s$rho^2) * rnorm(J)))
  beta0 <- rnorm(J, 0, sqrt(0.1))                          # intercepts as in Hui et al.
  delta <- es * (s$delta_mean + s$delta_sd * w)
  disp  <- if (s$family == "negbin") runif(J, 0, 5) else NULL   # as in Hui et al.
  list(beta0 = beta0, beta = beta, delta = delta, disp = disp)
}

#' Responses. Binary: Bernoulli-logit. Counts: NB2 with variance mu + phi * mu^2.
gen_responses <- function(s, cov, coefs) {
  N <- length(cov$climate); J <- length(coefs$beta)
  eta <- outer(rep(1, N), coefs$beta0) + outer(cov$climate, coefs$beta) + cov$z %*% t(coefs$delta)
  if (s$family == "binary") {
    y <- matrix(rbinom(N * J, 1, plogis(eta)), N, J)
  } else {
    mu <- exp(eta)
    y <- matrix(rnbinom(N * J, mu = mu, size = rep(1 / coefs$disp, each = N)), N, J)
  }
  y
}

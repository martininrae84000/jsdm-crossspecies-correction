## =============================================================================
## Cross-species correction: functions used by the application.
## Same computations as in the simulation code (../simulations/R/fit.R), written for
## any number of covariates P (column 1 of the coefficient matrices = intercept).
## The correction is a post-processing of a fitted SFA model: it uses only the
## estimated coefficients, the estimated loadings and their joint covariance.
## =============================================================================

#' Lower bounds on the diagonal of the loading matrix (identifiability), as in Hui et al.
loading_bounds <- function(obj, J, K) {
  blank <- matrix(0, J, K)
  blank[lower.tri(blank, diag = TRUE)] <- seq_len(sum(lower.tri(blank, diag = TRUE)))
  lower <- rep(-Inf, length(obj$par)); upper <- rep(Inf, length(obj$par))
  lower[grep("^loadings$", names(obj$par))][diag(blank)[seq_len(min(J, K))]] <- 1e-3
  list(lower = lower, upper = upper)
}

#' Joint covariance of the coefficients of covariate p (column p of the ADREPORTed
#' J x P matrix `bname`, stored column-major) and of the J x K loadings.
joint_cov_p <- function(sdr, bname, J, K, p) {
  rn <- names(sdr$value)
  ib <- which(rn == bname)[(p - 1) * J + 1:J]
  iL <- which(rn == "loadings_mat")[1:(J * K)]
  if (is.null(sdr$cov)) stop("sdreport() was run without the covariance of the reported quantities")
  sdr$cov[c(ib, iL), c(ib, iL)]
}

#' Cross-species correction for one covariate, with delta-method standard errors.
#' @param b  Coefficients of the covariate (length J), e.g. RSFA estimates.
#' @param L  Estimated loadings (J x K).
#' @param V  Joint covariance of c(b, vec(L)) (from joint_cov_p).
#' @return list(est, se, gamma, gamma_cov, cov_est)
cross_correct_R <- function(b, L, V) {
  J = length(b)
  K = ncol(L)
  ## h: from theta = c(coefficients, vec(loadings)) to c(corrected coefficients, gamma)
  h = function(theta) {
    bb = theta[1:J]                                   # coefficients (J)
    LL = matrix(theta[-(1:J)], J, K)                  # loadings (J x K)
    coef_reg = qr.coef(qr(cbind(1, LL)), bb)          # cross-species regression: intercept + K slopes
    gamma_hat = coef_reg[-1]                          # the K slopes (gamma_p, eq. 8)
    lambda_gamma = LL %*% gamma_hat                   # estimated bias of each species (J x 1)
    correct_val = bb - as.vector(lambda_gamma)        # corrected coefficients (eq. 9)
    return(c(correct_val, gamma_hat))                 # length J + K
  }
  theta = c(b, as.vector(L))
  val = h(theta)
  ## Jacobian of h by central finite differences, then delta method
  eps = 1e-6 * pmax(1, abs(theta))
  Jac = sapply(seq_along(theta), function(k) {
    e = numeric(length(theta))
    e[k] = eps[k]
    res = (h(theta + e) - h(theta - e)) / (2 * eps[k])
    return(res)
  })
  covh = Jac %*% V %*% t(Jac)
  list(est = val[1:J],
       se = sqrt(pmax(0, diag(covh)[1:J])),
       gamma = val[J + 1:K],
       gamma_cov = covh[J + 1:K, J + 1:K, drop = FALSE],
       cov_est = covh[1:J, 1:J])
}

#' Wald test of gamma = 0 (K degrees of freedom).
wald_gamma <- function(g, V) {
  g <- as.vector(g)
  stat <- tryCatch(drop(t(g) %*% solve(V) %*% g), error = function(e) NA_real_)
  c(stat = stat, df = length(g), p = if (is.na(stat)) NA_real_ else pchisq(stat, df = length(g), lower.tail = FALSE))
}

#' Diagnostics of condition C1 for a fitted loading matrix:
#'  - kappa: condition number of L' C L (large = loadings nearly collinear);
#'  - spread_to_se: for each factor, SD of the estimated loadings across species divided
#'    by their mean standard error (close to or below 1 = species hardly distinguishable,
#'    the correction cannot separate confounding from a common direct effect).
c1_diagnostics <- function(L, se_L) {
  Lc <- scale(L, center = TRUE, scale = FALSE)
  list(kappa = kappa(crossprod(Lc), exact = TRUE),
       spread_to_se = apply(L, 2, sd) / colMeans(se_L, na.rm = TRUE))
}

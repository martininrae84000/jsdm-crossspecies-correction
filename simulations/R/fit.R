## =============================================================================
## Model fitting with the TMB templates of Hui, Vu & Hooten (2024)
## (https://github.com/fhui28/spatialconfounding_JSDM, v0.6, GPL-3). Only change: the
## ADREPORT lines of the latent factors are commented out (not used here; they slow
## down sdreport()). Results are identical to the original files.
##   tmb/TMBINLA_SFA_binary.cpp  = simulation/binary_N400/TMBINLA_SFA.cpp
##   tmb/TMBINLA_SFA_negbin.cpp  = simulation/negbin_N400/TMBINLA_SFA.cpp
##   tmb/TMB_indFA_binary.cpp    = simulation/binary_N400/TMB_indFA.cpp
##   tmb/TMB_indFA_negbin.cpp    = simulation/negbin_N400/TMB_indFA.cpp
## Only the file names are changed (TMB names a compiled model after its file, and
## the binary and negative binomial versions must be loaded together).
## The cross-species correction is computed in R from the ADREPORTed quantities
## (betas, betas_rsr, loadings_mat) and their joint covariance from sdreport().
## The fitting steps (starting values, bounds, nlminb settings) follow Hui et al.
## The TMB libraries must have been compiled and loaded beforehand (load_tmb_models()).
## =============================================================================

tmb_models <- c("TMBINLA_SFA_binary", "TMBINLA_SFA_negbin", "TMB_indFA_binary", "TMB_indFA_negbin")

#' Compile (if needed) and load the four TMB models. Run from the project root.
load_tmb_models <- function(compile_if_needed = TRUE) {
  for (m in tmb_models) {
    lib <- TMB::dynlib(file.path("tmb", m))
    if (compile_if_needed && !file.exists(lib)) TMB::compile(file.path("tmb", paste0(m, ".cpp")))
    dyn.load(lib)
  }
  invisible(NULL)
}

## Lower bounds on the diagonal of the loading matrix (identifiability), as in Hui et al.
loading_bounds <- function(obj, J, K) {
  blank <- matrix(0, J, K)
  blank[lower.tri(blank, diag = TRUE)] <- seq_len(sum(lower.tri(blank, diag = TRUE)))
  lower <- rep(-Inf, length(obj$par)); upper <- rep(Inf, length(obj$par))
  lower[grep("^loadings$", names(obj$par))][diag(blank)[seq_len(min(J, K))]] <- 1e-3
  list(lower = lower, upper = upper)
}

## Estimates and standard errors of ADREPORTed quantities
extract_sdr <- function(sdr, names) {
  est <- as.list(sdr, "Estimate", report = TRUE)
  se  <- as.list(sdr, "Std. Error", report = TRUE)
  out <- lapply(names, function(n) list(est = est[[n]], se = se[[n]]))
  names(out) <- names
  out
}

#' Spatial factor analysis (SFA) and RSFA, model of Hui et al. (2024).
fit_sfa <- function(y, X, family, K, S, trace = 0) {
  J <- ncol(y); N <- nrow(y)
  XtXinv <- solve(crossprod(X))                         # (X'X)^{-1}, P x P
  ## Data passed to the C++ template (names must match its DATA_* declarations)
  data <- list(y = y, X = X, num_lv = K,
               OLSmatrix_transpose = X %*% XtXinv,        # X (X'X)^{-1}: gives the RSFA coefficients
               residualprojection = diag(N) - X %*% XtXinv %*% t(X),   # I - P_X: restricted factors
               n_i = S$mesh$n, meshidxloc = S$mesh$idx$loc - 1,       # mesh size; vertex of each site (0-based)
               G0 = S$fem$c0, G1 = S$fem$g1, G2 = S$fem$g2)           # SPDE finite-element matrices
  ## Starting values (names must match the PARAMETER_* declarations)
  pars <- list(betas = matrix(0, J, ncol(X)))              # species coefficients (fixed effects)
  if (family == "negbin") pars$log_dispparam <- numeric(J) # log NB dispersion, one per species
  pars$loadings <- numeric(sum(lower.tri(matrix(0, J, K), diag = TRUE)))   # lower-triangular loadings
  pars$lvs <- matrix(0, S$mesh$n, K)                       # latent fields at the mesh vertices
  pars$log_kappa <- rep(0, K)                              # log Matern range parameters
  ## The latent fields are random effects, integrated out by the Laplace approximation
  obj <- TMB::MakeADFun(data, pars, random = "lvs", DLL = paste0("TMBINLA_SFA_", family),
                        hessian = FALSE, silent = TRUE)
  b <- loading_bounds(obj, J, K)                           # diagonal loadings >= 1e-3
  t0 <- Sys.time()
  opt <- nlminb(obj$par, obj$fn, obj$gr, lower = b$lower, upper = b$upper,
                control = list(iter.max = 2000, eval.max = 2000, trace = trace))
  ## Standard errors and joint covariance of the ADREPORTed quantities
  ## (betas, betas_rsr, loadings_mat), at the optimum stored in obj$env$last.par.best
  sdr <- TMB::sdreport(obj)
  res <- extract_sdr(sdr, c("betas", "betas_rsr", "loadings_mat"))
  res$sdr <- sdr
  res$convergence <- opt$convergence
  res$pdHess <- sdr$pdHess
  res$time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  res
}

#' Independent (non-spatial) factor analysis, model of Hui et al. (2024).
fit_indfa <- function(y, X, family, K, trace = 0) {
  J <- ncol(y); N <- nrow(y)
  data <- list(y = y, X = X, num_lv = K)
  pars <- list(betas = matrix(0, J, ncol(X)))
  if (family == "negbin") pars$log_dispparam <- numeric(J)
  pars$loadings <- numeric(sum(lower.tri(matrix(0, J, K), diag = TRUE)))
  pars$lvs <- matrix(0, N, K)
  obj <- TMB::MakeADFun(data, pars, random = "lvs", DLL = paste0("TMB_indFA_", family),
                        hessian = FALSE, silent = TRUE)
  b <- loading_bounds(obj, J, K)
  t0 <- Sys.time()
  opt <- nlminb(obj$par, obj$fn, obj$gr, lower = b$lower, upper = b$upper,
                control = list(iter.max = 2000, eval.max = 2000, trace = trace))
  sdr <- TMB::sdreport(obj)
  res <- extract_sdr(sdr, c("betas", "loadings_mat"))
  res$sdr <- sdr
  res$convergence <- opt$convergence
  res$pdHess <- sdr$pdHess
  res$time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  res
}

#' Oracle: species-by-species GLM including the true confounders.
fit_oracle <- function(y, climate, z, family) {
  J <- ncol(y)
  out <- t(sapply(seq_len(J), function(j) {
    d <- data.frame(y = y[, j], climate = climate, z1 = z[, 1], z2 = z[, 2])
    m <- if (family == "binary") glm(y ~ climate + z1 + z2, family = binomial, data = d)
         else MASS::glm.nb(y ~ climate + z1 + z2, data = d)
    c(est = unname(coef(m)["climate"]), se = unname(sqrt(diag(vcov(m)))["climate"]))
  }))
  list(est = out[, "est"], se = out[, "se"])
}

#' Cross-species correction computed in R, with delta-method standard errors.
#' @param b Estimated slopes for one covariate (length J).
#' @param L Estimated loadings (J x K).
#' @param V Joint covariance matrix of c(b, vec(L)) (length J + J*K), from sdreport.
#' @return list(est, se, gamma, gamma_cov)
cross_correct_R <- function(b, L, V) {
  J <- length(b); K <- ncol(L)
  ## h: from theta = c(slopes, vec(loadings)) to c(corrected slopes, gamma)
  h <- function(theta) {
    bb <- theta[1:J]; LL <- matrix(theta[-(1:J)], J, K)
    g <- qr.coef(qr(cbind(1, LL)), bb)[-1]   # cross-species regression with intercept: K slopes (eq. 8)
    c(bb - drop(LL %*% g), g)                 # corrected slopes (eq. 9), then gamma
  }
  theta <- c(b, as.vector(L))
  val <- h(theta)
  ## Numerical Jacobian (central differences)
  eps <- 1e-6 * pmax(1, abs(theta))
  Jac <- sapply(seq_along(theta), function(k) {
    e <- numeric(length(theta)); e[k] <- eps[k]
    (h(theta + e) - h(theta - e)) / (2 * eps[k])
  })
  Vh <- Jac %*% V %*% t(Jac)                      # delta method
  list(est = val[1:J], se = sqrt(pmax(0, diag(Vh)[1:J])),
       gamma = val[J + 1:K], gamma_cov = Vh[J + 1:K, J + 1:K, drop = FALSE])
}

#' Joint covariance of the climate slopes (column 2 of the ADREPORTed matrix
#' named `bname`) and of the loadings, extracted from an sdreport object.
joint_cov <- function(sdr, bname, J, K) {
  rn <- names(sdr$value)
  ib <- which(rn == bname)[J + 1:J]                 # column 2 (column-major storage)
  iL <- which(rn == "loadings_mat")[1:(J * K)]
  sdr$cov[c(ib, iL), c(ib, iL)]
}

#' Wald test of gamma = 0 (K degrees of freedom) for the climate slope.
wald_gamma <- function(g, V) {
  g <- as.vector(g)
  stat <- tryCatch(drop(t(g) %*% solve(V) %*% g), error = function(e) NA_real_)
  c(stat = stat, p = if (is.na(stat)) NA_real_ else pchisq(stat, df = length(g), lower.tail = FALSE))
}

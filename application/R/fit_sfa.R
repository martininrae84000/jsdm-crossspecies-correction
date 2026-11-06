## =============================================================================
## SFA fit of Hui, Vu & Hooten (2024) with TMB, used by application_butterflies.R.
## Follows the "fit spatial factor analytic (SFA) model" section of
## application_butterflies/analysis.R of Hui et al.
## (https://github.com/fhui28/spatialconfounding_JSDM, v0.6, GPL-3), written as a
## function with the model ("random" or "fixed") and the number of factors K as arguments.
## =============================================================================

## -----------------------------------------------------------------------------
## fit_sfa: fits the SFA model of Hui et al. with TMB.
##   model: "random" (species coefficients as random effects) or "fixed"
##   K:     number of spatial latent factors
##   y:     response matrix (N sites x J species)
##   X:     covariate matrix (N x P), first column = intercept
##   mesh:  SPDE mesh with a vertex at each site (INLA or fmesher)
##   G:     list of the finite-element matrices G0, G1, G2 of the mesh
##   dll:   named vector of the compiled TMB libraries, c(random = ..., fixed = ...);
##          the libraries must be compiled and loaded beforehand
##          (see application_butterflies.R).
## Requires the package TMB and the function loading_bounds() of R/correction.R.
## Returns a list with the optimisation, the sdreport (estimates and covariances),
## the latent factors, the restricted factors and convergence diagnostics.
## -----------------------------------------------------------------------------
fit_sfa <- function(model, K, y, X, mesh, G, dll) {

  J <- ncol(y); P <- ncol(X)                       # numbers of species and of columns of X

  ## (X'X)^{-1}, P x P; crossprod(X) computes t(X) %*% X.
  XtXinv <- solve(crossprod(X))

  ## ---- Data passed to the C++ template ----
  ## Names must match the DATA_* declarations of the .cpp file.
  data <- list(y = y,                              # 0/1 responses: N sites x J species
               X = X,                              # covariates: N x P (column 1 = intercept)
               num_lv = K,                         # number of latent factors
               OLSmatrix_transpose = X %*% XtXinv, # X (X'X)^{-1}, N x P; in the .cpp:
                                                   # betas_rsr = betas + Lambda U' X (X'X)^{-1}
               residualprojection = matrix(0, 1, 1),   # not used: in Hui et al., the N x N matrix
                                                   # (I - P_X), removed from the template (memory)
               n_i = mesh$n,                       # number of mesh vertices
               meshidxloc = mesh$idx$loc - 1,      # vertex of each site; -1 because indices
                                                   # start at 0 in C++ and at 1 in R
               G0 = G$G0, G1 = G$G1, G2 = G$G2)    # SPDE finite-element matrices

  ## ---- Starting values of the parameters ----
  ## Names must match the PARAMETER_* declarations of the .cpp file. All start at zero.
  pars <- list(betas = matrix(0, J, P))            # species coefficients: J x P
  ## Model "random": beta_pj ~ N(alpha_p, sd_p^2) for each covariate p;
  ## alpha = community means, log_sds = log of the between-species standard deviations
  ## (the log keeps the standard deviations positive).
  if (model == "random") { pars$alpha <- numeric(P); pars$log_sds <- numeric(P) }
  ## Loadings: only the lower-triangular part (diagonal included) of the J x K matrix is
  ## estimated (identifiability constraint). It is stored as a vector; the .cpp rebuilds
  ## the matrix loadings_mat.
  pars$loadings <- numeric(sum(lower.tri(matrix(0, J, K), diag = TRUE)))
  pars$lvs <- matrix(0, mesh$n, K)                 # values of the K latent fields at the vertices
  pars$log_kappa <- rep(0, K)                      # log of the Matern range parameter, one per
                                                   # factor (range ~ sqrt(8) / kappa)

  ## ---- Random effects ----
  ## Integrated out by the Laplace approximation: the likelihood optimised is the marginal
  ## likelihood of the fixed parameters only. The latent fields are always random; in the
  ## "random" model, the species coefficients are random too.
  random <- if (model == "random") c("lvs", "betas") else "lvs"

  ## ---- Objective function ----
  ## MakeADFun builds obj$fn (minus the marginal log-likelihood) and obj$gr (its gradient,
  ## by automatic differentiation) from the compiled template dll[[model]].
  ## hessian = FALSE: no Hessian during the optimisation; silent = TRUE: no tracing output.
  obj <- MakeADFun(data, pars, random = random, DLL = dll[[model]], hessian = FALSE, silent = TRUE)

  ## ---- Bounds ----
  ## The diagonal loadings must be >= 1e-3 (positive, for identifiability);
  ## all other parameters are unbounded.
  b <- loading_bounds(obj, J, K)

  t0 <- Sys.time()                                 # start timing

  ## ---- Optimisation ----
  ## nlminb minimises obj$fn within the bounds; iter.max / eval.max = maximum numbers of
  ## iterations and evaluations; trace = 1 prints the progress at each iteration.
  ## During the optimisation, TMB stores the best point found (fixed parameters and modes
  ## of the random effects) in obj$env$last.par.best.
  opt <- nlminb(obj$par, obj$fn, obj$gr, lower = b$lower, upper = b$upper,
                control = list(iter.max = 10000, eval.max = 10000, trace = 1))

  ## ---- Uncertainty ----
  ## sdreport computes, at the optimum, the standard errors and the joint covariance of the
  ## parameters and of the ADREPORTed quantities of the .cpp (betas, betas_rsr, loadings_mat),
  ## from the Hessian (generalised delta method, which accounts for the random effects).
  sdr <- sdreport(obj)

  ## REPORTed quantities of the .cpp (values only, no standard errors), evaluated at the
  ## best point: here the latent factors at the sites, lvs_units (N x K).
  rep <- obj$report(obj$env$last.par.best)

  ## ---- Results ----
  list(model = model,                              # "random" or "fixed"
       K = K,
       opt = opt,                                  # raw output of nlminb
       sdr = sdr,                                  # raw output of sdreport (including sdr$cov)
       ## Estimates and standard errors as named lists:
       ## $betas (SFA), $betas_rsr (RSFA), $loadings_mat, ... (report = TRUE includes ADREPORTs).
       est = as.list(sdr, "Estimate", report = TRUE), se = as.list(sdr, "Std. Error", report = TRUE),
       lvs_units = rep$lvs_units,                  # latent factors U at the sites (N x K)
       ## Restricted factors (I - P_X) U: the part of U orthogonal to the covariates, as in
       ## Hui et al. Computed without forming the N x N matrix P_X: P_X U = X ((X'X)^{-1} (X'U)).
       lvs_rsr = rep$lvs_units - X %*% (XtXinv %*% crossprod(X, rep$lvs_units)),
       convergence = opt$convergence,              # 0 = the optimiser converged
       message = opt$message,                      # nlminb message, e.g. "relative convergence (4)"
       pdHess = sdr$pdHess,                        # TRUE if the Hessian is positive definite
       time_min = as.numeric(difftime(Sys.time(), t0, units = "mins")))   # duration in minutes
}

## =============================================================================
## Application: butterflies of Great Britain (presence-absence)
## Cross-species correction for spatial confounding, applied to the spatial factor
## analysis (SFA) of Hui, Vu & Hooten (2024).
##
## Data (not included): Wilkinson et al. (2019), https://zenodo.org/records/1452066,
## folder Datasets/Butterflies. Put these three files in ../data/ (see ../data/README.md):
##   Butterfly_PA.csv             presence-absence (sites x species)
##   Butterfly_LatLon_scaled.csv  site coordinates, rescaled
##   Butterfly_Cov.csv            site covariates
##
## The data preparation, the mesh and the model follow application_butterflies/analysis.R
## of Hui et al. (https://github.com/fhui28/spatialconfounding_JSDM, v0.6, GPL-3).
## Two changes to their TMB template, which do not affect the fit:
##   - the latent factors are REPORTed instead of ADREPORTed: their standard errors are
##     not needed, and with N = 2609 sites they make sdreport() very slow;
##   - the restricted factors lvs_rsr = (I - P_X) U are computed in R after the fit,
##     not inside the template: the dense N x N product needed several GB of memory.
## The templates are compiled with the TMBad framework (much less memory, see below).
##
## Run from the application/ folder:  Rscript application_butterflies.R
## (or source() it from RStudio after setwd() to this folder).
## Each fit is saved in results/ and is not refitted if the file already exists.
## =============================================================================

## Packages (startup messages hidden).
##   dplyr / tidyr: data handling; ggplot2: plots;
##   TMB: model fitting (Template Model Builder); Matrix: sparse matrices.
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(TMB); library(Matrix)
})
## Functions of the cross-species correction: loading_bounds, joint_cov_p,
## cross_correct_R, wald_gamma, c1_diagnostics.
source("R/correction.R")

## ---------------------------------------------------------------- settings ----
## Settings are read from environment variables, with a default value.
## Example, from R:  Sys.setenv(MODELS = "random,fixed", K_VALUES = "2,3")
##
## MODELS: which TMB model(s) to fit
##   "random": the model of Hui et al. for this case study, in which the species
##             coefficients are random effects, beta_pj ~ N(alpha_p, sd_p^2)
##             ("to stabilise fitting"); tmb/TMB_SFA_random.cpp
##   "fixed":  species coefficients as fixed effects, as in the simulations;
##             tmb/TMB_SFA_fixed.cpp
## K_VALUES: numbers of spatial factors (Hui et al.: 2). Adding 3 checks the sensitivity
## of the correction to K.
## strsplit(..., ",") splits a string such as "2,3" into the vector c("2", "3").
MODELS    <- strsplit(Sys.getenv("MODELS", "random"), ",")[[1]]
K_VALUES  <- as.integer(strsplit(Sys.getenv("K_VALUES", "2"), ",")[[1]])
## Folder that contains the three data files.
DATA_DIR  <- Sys.getenv("DATA_DIR", file.path("..", "data"))
## Hui et al. keep the species present at more than 150 sites (45 of 55).
MIN_PRESENCES <- 150
## DROP_DUPLICATES (default 1): in Butterfly_PA.csv, two pairs of columns are identical
## (small_blue = small_copper, common_blue = camberwell_beauty). A duplicated species would
## count twice in the loadings and in the cross-species regression; only the first column
## of each pair is kept. Set to 0 to keep both, as in Hui et al.
DROP_DUPLICATES <- Sys.getenv("DROP_DUPLICATES", "1") == "1"
## QUADRATIC (default 0): 1 adds the square of growing degree days (GDD) to the covariates.
## Species responses to climate are often unimodal; an omitted quadratic term would be
## captured by the latent factors, with loadings related to the species' climate responses,
## which would violate condition C3. Comparing the corrected effect of GDD with and without
## the quadratic term checks this.
QUADRATIC <- Sys.getenv("QUADRATIC", "0") == "1"
## Suffix added to the names of the result files, so that the variants are not mixed:
## "_dedup" when the duplicated species are removed, "_quad" with the quadratic term.
RUN_TAG <- paste0(if (DROP_DUPLICATES) "_dedup" else "", if (QUADRATIC) "_quad" else "")

## Create the results folder if needed (no warning if it already exists).
dir.create("results", showWarnings = FALSE)

## ------------------------------------------------------------------- data ----
## Check that the three files are there; otherwise stop with a message giving the
## working directory and where to download the data.
files_needed <- file.path(DATA_DIR, c("Butterfly_PA.csv", "Butterfly_LatLon_scaled.csv", "Butterfly_Cov.csv"))
if (!all(file.exists(files_needed))) stop(
  "Data files not found: ", paste(files_needed[!file.exists(files_needed)], collapse = ", "),
  "\nCurrent working directory: ", getwd(),
  "\nDownload them from https://zenodo.org/records/1452066 (folder Datasets/Butterflies) and put them in ",
  normalizePath(DATA_DIR, mustWork = FALSE), ", or set DATA_DIR to the folder that contains them.", call. = FALSE)
## Presence-absence: one row per site, one column per species (0/1).
resp_dat <- read.csv(file.path(DATA_DIR, "Butterfly_PA.csv"))
## colSums = number of sites where each species is present; keep those above the threshold.
resp_dat <- resp_dat[, colSums(resp_dat) > MIN_PRESENCES]
## duplicated() on the list of columns: TRUE for a column identical to an earlier one
## (the first column of each pair is FALSE, hence kept).
dup <- duplicated(as.list(resp_dat))
if (any(dup)) {
  ## For each duplicated column, find its "twin" (the first identical column) and print the pair.
  for (k in which(dup)) {
    twin <- names(resp_dat)[which(vapply(resp_dat, identical, logical(1), resp_dat[[k]]))[1]]
    cat("Identical columns:", twin, "=", names(resp_dat)[k], if (DROP_DUPLICATES) "-> dropped\n" else "-> kept\n")
  }
  ## Drop the duplicates if requested (the first column of each pair is kept).
  if (DROP_DUPLICATES) resp_dat <- resp_dat[, !dup]
}
## Covariates and site coordinates.
cov_dat  <- read.csv(file.path(DATA_DIR, "Butterfly_Cov.csv"))
longlat  <- read.csv(file.path(DATA_DIR, "Butterfly_LatLon_scaled.csv"))
## Bind coordinates (x, y) and covariates; the first column of the covariate file is a
## column of ones (intercept), removed here since model.matrix() adds it back.
cov_dat  <- cbind(longlat, cov_dat[, -1])
## Standardise the woodland percentages (mean 0, standard deviation 1).
## climate (GDD) is already standardised in the file.
cov_dat  <- cov_dat %>%
  mutate(blwood = as.vector(scale(blwood)),
         conwood = as.vector(scale(conwood)))

## Response matrix: N sites x J species.
y <- as.matrix(resp_dat)
## Covariate matrix X (N x P), with a column of ones for the intercept.
## climate = growing degree days (GDD), blwood = broadleaved woodland, conwood = coniferous
## woodland; with QUADRATIC, climate^2 is added (I() protects ^ inside the formula).
X <- if (QUADRATIC) model.matrix(~ climate + I(climate^2) + blwood + conwood, data = cov_dat) else
  model.matrix(~ climate + blwood + conwood, data = cov_dat)
## Covariate names, without the intercept.
covariates <- colnames(X)[-1]
## Readable labels for the plots.
cov_labels <- c(climate = "Growing degree days", `I(climate^2)` = "Growing degree days (squared)",
                blwood = "Broadleaved woodland", conwood = "Coniferous woodland")
## Dimensions: J species, N sites, P columns of X (intercept included).
J <- ncol(y); N <- nrow(y); P <- ncol(X)
cat("N =", N, "sites, J =", J, "species, covariates:", paste(covariates, collapse = ", "),
    "| run tag:", RUN_TAG, "\n")

## ------------------------------------------------------------------- mesh ----
## The spatial factors are Matern fields approximated by the SPDE approach
## (Lindgren et al. 2011), which needs a mesh (triangulation) of the domain.
## As in Hui et al., the mesh has a vertex at each site (refine = FALSE: no vertex added).
## INLA is used if installed, otherwise fmesher, which gives the same matrices.
## G0, G1, G2 are the finite-element matrices: in the TMB template, the precision matrix
## of each field is kappa^4 G0 + 2 kappa^2 G1 + G2, with kappa the range parameter.
Loc <- as.matrix(cov_dat[, c("x", "y")])
if (requireNamespace("INLA", quietly = TRUE)) {
  mesh <- INLA::inla.mesh.create(Loc, plot.delay = NULL, refine = FALSE)
  spde <- INLA::inla.spde2.matern(mesh, alpha = 2)   # alpha = 2: Matern smoothness nu = 1 in 2D
  G <- list(G0 = spde$param.inla$M0, G1 = spde$param.inla$M1, G2 = spde$param.inla$M2)
  cat("Mesh built with INLA:", mesh$n, "vertices\n")
} else {
  mesh <- fmesher::fm_rcdt_2d_inla(loc = Loc, refine = FALSE)
  fem <- fmesher::fm_fem(mesh, order = 2)             # c0, g1, g2 = M0, M1, M2 of INLA
  G <- list(G0 = fem$c0, G1 = fem$g1, G2 = fem$g2)
  cat("Mesh built with fmesher:", mesh$n, "vertices\n")
}

## ----------------------------------------------------------------- models ----
## Name of the compiled library (DLL) of each model.
dll <- c(random = "TMB_SFA_random", fixed = "TMB_SFA_fixed")
## The models are compiled with TMBad, an alternative automatic-differentiation engine of
## TMB. With the default engine (CppAD), building the sparse Hessian of the random effects
## for the 2609 sites needs more than 6 GB of memory; with TMBad, about 0.5 GB.
## A marker file (.tmbad) records that the library was compiled with TMBad, so that a
## library compiled earlier with CppAD is recompiled.
## Both libraries are compiled (if needed) and loaded, whatever MODELS contains.
for (m in dll) {
  marker <- file.path("tmb", paste0(m, ".tmbad"))
  ## Compile only if the library does not exist yet or was not compiled with TMBad.
  if (!file.exists(dynlib(file.path("tmb", m))) || !file.exists(marker)) {
    compile(file.path("tmb", paste0(m, ".cpp")), framework = "TMBad")
    file.create(marker)
  }
  ## Load the compiled library into the R session (unless already loaded).
  if (!m %in% names(getLoadedDLLs())) dyn.load(dynlib(file.path("tmb", m)))
}

## -----------------------------------------------------------------------------
## fit_sfa: fits the SFA model of Hui et al. with TMB.
##   model: "random" (species coefficients as random effects) or "fixed"
##   K:     number of spatial latent factors
## Uses objects created above: y, X, J, P, mesh, G, dll.
## Returns a list with the optimisation, the sdreport (estimates and covariances),
## the latent factors, the restricted factors and convergence diagnostics.
## -----------------------------------------------------------------------------
fit_sfa <- function(model, K) {

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

## Loop over models and values of K. Each fit is saved in an .rds file; if the file
## already exists, it is read instead of refitting (one fit takes 15 minutes or more).
fits <- list()
for (md in MODELS) for (K in K_VALUES) {
  f <- file.path("results", sprintf("fit_%s_K%d%s.rds", md, K, RUN_TAG))
  if (file.exists(f)) { fits[[basename(f)]] <- readRDS(f); next }
  cat("\n--- Fitting SFA, model", md, ", K =", K, "---\n")
  fit <- fit_sfa(md, K)
  saveRDS(fit, f)
  fits[[basename(f)]] <- fit
  cat("convergence:", fit$convergence, fit$message, "| pdHess:", fit$pdHess,
      "| time:", round(fit$time_min, 1), "min\n")
}

## ---------------------------------------------------------- post-process ----
## For each fit and each covariate: SFA, RSFA and corrected coefficients with their
## standard errors, community-mean effects, Wald test of gamma = 0, C1 diagnostics.
spp <- colnames(y)                                     # species names
coef_tab <- list(); cov_tab <- list(); c1_tab <- list() # tables filled in the loop
for (nm in names(fits)) {
  fit <- fits[[nm]]; K <- fit$K
  ## Estimated loadings (J x K) and their standard errors.
  L <- fit$est$loadings_mat; se_L <- fit$se$loadings_mat
  ## Diagnostics of condition C1 (heterogeneous responses of the species):
  ##   kappa = condition number of L'CL (large = loadings nearly collinear);
  ##   spread_to_se = spread of the loadings across species / their mean standard error,
  ##   for each factor (close to 1 or below = species hardly distinguishable).
  c1 <- c1_diagnostics(L, se_L)
  c1_tab[[nm]] <- data.frame(model = fit$model, K = K, convergence = fit$convergence,
                             pdHess = fit$pdHess, kappa_LCL = c1$kappa,
                             factor = paste0("LV", 1:K), spread_to_se = c1$spread_to_se)
  ## Loop over the covariates (p = 1 is the intercept, which is not corrected).
  for (p in 2:P) {
    v <- colnames(X)[p]
    ## SFA and RSFA coefficients of covariate p for the J species, with their standard errors.
    b_sfa <- fit$est$betas[, p];     se_sfa <- fit$se$betas[, p]
    b_rsr <- fit$est$betas_rsr[, p]; se_rsr <- fit$se$betas_rsr[, p]
    ## Joint covariance of the RSFA coefficients of covariate p and of the loadings,
    ## extracted from the sdreport: (J + J*K) x (J + J*K) matrix.
    V <- joint_cov_p(fit$sdr, "betas_rsr", J, K, p)
    ## Cross-species correction (equations 8-9 of the paper):
    ##   gamma = slopes of the regression of the RSFA coefficients on the loadings (with intercept),
    ##   est   = corrected coefficients = b_rsr - L gamma,
    ##   se, gamma_cov, cov_est = standard errors and covariances, obtained by propagating V
    ##   with the delta method (numerical Jacobian).
    cr <- cross_correct_R(b_rsr, L, V)
    ## Wald test of H0: gamma = 0 (K degrees of freedom):
    ## W = gamma' Var(gamma)^{-1} gamma, compared with a chi-square with K df.
    w  <- wald_gamma(cr$gamma, cr$gamma_cov)
    ## Invariance check (Proposition 1): the correction applied to the SFA coefficients
    ## must give the same corrected coefficients as when applied to the RSFA coefficients.
    cr_sfa <- cross_correct_R(b_sfa, L, joint_cov_p(fit$sdr, "betas", J, K, p))
    ## Species-level table: one row per species and method (SFA, RSFA, corrected).
    ## bias_hat = estimated confounding bias lambda_j' gamma.
    coef_tab[[paste(nm, v)]] <- bind_rows(
      data.frame(method = "SFA", est = b_sfa, se = se_sfa),
      data.frame(method = "RSFA", est = b_rsr, se = se_rsr),
      data.frame(method = "Corrected", est = cr$est, se = cr$se)) %>%
      mutate(model = fit$model, K = K, covariate = v, species = rep(spp, 3),
             bias_hat = rep(drop(L %*% cr$gamma), 3), .before = 1)
    ## Community-mean effects (mean over species). Variance of the mean of J correlated
    ## estimates: sum of all elements of their covariance matrix divided by J^2;
    ## hence se = sqrt(sum(V)) / J.
    Vb_rsr <- V[1:J, 1:J]                                       # covariance of the RSFA coefficients
    Vb_sfa <- joint_cov_p(fit$sdr, "betas", J, K, p)[1:J, 1:J]  # covariance of the SFA coefficients
    cov_tab[[paste(nm, v)]] <- data.frame(
      model = fit$model, K = K, covariate = v,
      mean_SFA = mean(b_sfa), se_mean_SFA = sqrt(sum(Vb_sfa)) / J,
      mean_RSFA = mean(b_rsr), se_mean_RSFA = sqrt(sum(Vb_rsr)) / J,
      mean_Corrected = mean(cr$est), se_mean_Corrected = sqrt(sum(cr$cov_est)) / J,
      ## Wald test: statistic, degrees of freedom, p-value.
      wald_stat = w[["stat"]], wald_df = w[["df"]], wald_p = w[["p"]],
      ## Estimated slope(s) gamma, one per factor, as a string "g1; g2".
      gamma = paste(sprintf("%.3f", cr$gamma), collapse = "; "),
      ## Proportion of species whose 95% interval excludes 0 (|z| > 1.96), for each method.
      prop_signif_SFA = mean(abs(b_sfa / se_sfa) > qnorm(0.975)),
      prop_signif_RSFA = mean(abs(b_rsr / se_rsr) > qnorm(0.975)),
      prop_signif_Corrected = mean(abs(cr$est / cr$se) > qnorm(0.975)),
      ## Maximum difference between the correction applied to SFA and to RSFA (~ 1e-14).
      invariance_maxdiff = max(abs(cr_sfa$est - cr$est)))
  }
}
## Stack the pieces into three tables and write them as CSV files.
coef_tab <- bind_rows(coef_tab); cov_tab <- bind_rows(cov_tab); c1_tab <- bind_rows(c1_tab)
write.csv(coef_tab, paste0("results/butterflies_coefficients", RUN_TAG, ".csv"), row.names = FALSE)
write.csv(cov_tab, paste0("results/butterflies_covariates", RUN_TAG, ".csv"), row.names = FALSE)
write.csv(c1_tab, paste0("results/butterflies_C1_diagnostics", RUN_TAG, ".csv"), row.names = FALSE)

## Print the main results in the console (3 significant digits).
options(width = 160)
cat("\n=== Community-mean effects, Wald test of gamma = 0, proportion of significant species ===\n")
print(cov_tab %>% mutate(across(where(is.double), ~ signif(.x, 3))), row.names = FALSE)
cat("\n=== Condition C1: spread of the loadings relative to their standard errors ===\n")
cat("(kappa_LCL large or spread_to_se close to 1: the correction is unreliable)\n")
print(c1_tab %>% mutate(across(where(is.double), ~ signif(.x, 3))), row.names = FALSE)

## ------------------------------------------------------------------ plots ----
## Colours of the three methods.
mcol <- c(SFA = "#9CA3AF", RSFA = "#C2410C", Corrected = "#1D4ED8")
for (nm in names(fits)) {
  ## tag = name of the fit without "fit_" and ".rds", used in the names of the PDF files.
  fit <- fits[[nm]]; tag <- sub("\\.rds$", "", sub("^fit_", "", nm))
  ## Rows of the species-level table for this fit, with readable labels.
  d <- coef_tab %>% filter(model == fit$model, K == fit$K) %>%
    mutate(covariate = factor(cov_labels[covariate], levels = cov_labels))

  ## 1. Cross-species plot: RSFA coefficient against the estimated bias lambda_j' gamma.
  ##    Under the model, the points lie around a line of slope 1; the corrected coefficient
  ##    is the vertical distance to that line (plus the intercept).
  d1 <- d %>% filter(method == "RSFA")
  p1 <- ggplot(d1, aes(bias_hat, est)) +
    geom_hline(yintercept = 0, colour = "grey85") +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#1D4ED8", linewidth = 0.6) +
    geom_point(colour = "#C2410C", size = 1.6) +
    facet_wrap(~ covariate, scales = "free") +           # one panel per covariate
    labs(x = expression("Estimated confounding bias  " * hat(lambda)[j]^T * hat(gamma)),
         y = "RSFA coefficient") +
    theme_bw(base_size = 10)
  ggsave(sprintf("results/butterflies_crossspecies_%s.pdf", tag), p1, width = 9, height = 3.4)

  ## 2. Caterpillar plot: SFA, RSFA and corrected coefficients with 95% intervals,
  ##    species ordered by their corrected estimate.
  d2 <- d %>%
    group_by(covariate) %>%
    ## Order of the species on the vertical axis: by corrected coefficient, within each panel.
    mutate(species = factor(species, levels = unique(species[method == "Corrected"][order(est[method == "Corrected"])]))) %>%
    ungroup()
  p2 <- ggplot(d2, aes(est, species, colour = method)) +
    geom_vline(xintercept = 0, colour = "grey70") +
    geom_pointrange(aes(xmin = est - 1.96 * se, xmax = est + 1.96 * se),
                    position = position_dodge(width = 0.6), size = 0.15) +   # three methods side by side
    scale_colour_manual(values = mcol) +
    facet_wrap(~ covariate, scales = "free_x") +
    labs(x = "Coefficient (95% interval)", y = NULL, colour = NULL) +
    theme_bw(base_size = 8) + theme(legend.position = "bottom")
  ggsave(sprintf("results/butterflies_caterpillar_%s.pdf", tag), p2, width = 10, height = 9)

  ## 3. Estimated loadings on the first two factors, one point (name) per species.
  if (fit$K >= 2) {
    Ld <- data.frame(species = spp, LV1 = fit$est$loadings_mat[, 1], LV2 = fit$est$loadings_mat[, 2])
    p3 <- ggplot(Ld, aes(LV1, LV2, label = species)) + geom_hline(yintercept = 0, colour = "grey85") +
      geom_vline(xintercept = 0, colour = "grey85") + geom_text(size = 2.4) + theme_bw(base_size = 10) +
      labs(title = "Estimated loadings (species)")
    ggsave(sprintf("results/butterflies_loadings_%s.pdf", tag), p3, width = 6, height = 5.5)
  }
}
cat("\nResults in results/: butterflies_coefficients", RUN_TAG, ".csv, butterflies_covariates", RUN_TAG,
    ".csv, butterflies_C1_diagnostics", RUN_TAG, ".csv and the PDF plots.\n", sep = "")
## The sensitivity analysis for condition C3 is in sensitivity_C3.R
## (it reads results/fit_random_K2_dedup.rds and results/butterflies_covariates_dedup.csv).

## =============================================================================
## Simulate one dataset for one scenario, fit all methods, return tidy results.
## =============================================================================

#' @param s Design row (list).
#' @param dataset Dataset index (also used to build the seed).
#' @param S Output of setup_sites().
run_one <- function(s, dataset, S) {
  ## Seed from seed_id (scenario number of an earlier version of the design, see design.R)
  set.seed(1e5 * (if (is.null(s$seed_id)) s$design_id else s$seed_id) + dataset)
  cov   <- gen_covariates(s, S)
  coefs <- gen_coefficients(s)
  y     <- gen_responses(s, cov, coefs)
  J <- s$J; K <- s$K_fit

  ## Truths for the climate slope
  truth_direct <- coefs$beta
  truth_total  <- coefs$beta + drop(coefs$delta %*% cov$a)

  X <- cbind(1, climate = cov$climate)
  out <- list(); diag <- list()

  ## --- SFA / RSFA / correction ------------------------------------------------
  sfa <- tryCatch(fit_sfa(y, X, s$family, K, S), error = function(e) e)
  if (!inherits(sfa, "error")) {
    L <- sfa$loadings_mat$est
    ## Correction computed in R from the RSFA slopes and loadings, delta-method SEs.
    ## crs applies it to the SFA slopes: identical result (Proposition 1), checked below.
    cr  <- cross_correct_R(sfa$betas_rsr$est[, 2], L, joint_cov(sfa$sdr, "betas_rsr", J, K))
    crs <- cross_correct_R(sfa$betas$est[, 2],     L, joint_cov(sfa$sdr, "betas", J, K))   # from SFA
    ok_sfa <- sfa$convergence == 0 && isTRUE(sfa$pdHess)
    out$SFA         <- data.frame(est = sfa$betas$est[, 2],     se = sfa$betas$se[, 2], ok = ok_sfa)
    out$RSFA        <- data.frame(est = sfa$betas_rsr$est[, 2], se = sfa$betas_rsr$se[, 2], ok = ok_sfa)
    out$Cross_RSFA  <- data.frame(est = cr$est, se = cr$se, ok = ok_sfa)
    w <- wald_gamma(cr$gamma, cr$gamma_cov)
    diag$sfa <- list(convergence = sfa$convergence, pdHess = sfa$pdHess, time = sfa$time,
                     invariance_maxdiff = max(abs(crs$est - cr$est)),
                     gamma = cr$gamma, wald_stat = w[["stat"]], wald_p = w[["p"]],
                     kappa_LCL = kappa(crossprod(scale(L, scale = FALSE))),
                     r2_delta_on_L = unname(sapply(summary(lm(coefs$delta ~ L)), `[[`, "r.squared")),
                     max_abs_loading = max(abs(L)),
                     loadings = L)
  } else diag$sfa <- list(error = conditionMessage(sfa))

  ## --- Independent factor analysis + correction ------------------------------
  if (isTRUE(s$fit_indfa)) {
    ind <- tryCatch(fit_indfa(y, X, s$family, K), error = function(e) e)
    if (!inherits(ind, "error")) {
      Li <- ind$loadings_mat$est
      cri <- cross_correct_R(ind$betas$est[, 2], Li, joint_cov(ind$sdr, "betas", J, K))
      ok_ind <- ind$convergence == 0 && isTRUE(ind$pdHess)
      out$Independent       <- data.frame(est = ind$betas$est[, 2], se = ind$betas$se[, 2], ok = ok_ind)
      out$Cross_Independent <- data.frame(est = cri$est, se = cri$se, ok = ok_ind)
      diag$indfa <- list(convergence = ind$convergence, pdHess = ind$pdHess, time = ind$time,
                         max_abs_loading = max(abs(Li)))
    } else diag$indfa <- list(error = conditionMessage(ind))
  }

  ## --- Spatial+ : SFA with the climate residualised on a spatial basis --------
  if (isTRUE(s$fit_spatialplus)) {
    k_sp <- s$spatialplus_k
    clim_res <- resid(lm(cov$climate ~ S$basis[, 2:k_sp]))
    Xsp <- cbind(1, climate = clim_res)
    sp <- tryCatch(fit_sfa(y, Xsp, s$family, K, S), error = function(e) e)
    if (!inherits(sp, "error")) {
      out$SpatialPlus <- data.frame(est = sp$betas$est[, 2], se = sp$betas$se[, 2],
                                    ok = sp$convergence == 0 && isTRUE(sp$pdHess))
      diag$spatialplus <- list(convergence = sp$convergence, pdHess = sp$pdHess, time = sp$time)
    } else diag$spatialplus <- list(error = conditionMessage(sp))
  }

  ## --- Oracle (true confounders observed) ---------------------------------
  orc <- tryCatch(fit_oracle(y, cov$climate, cov$z, s$family), error = function(e) e)
  if (!inherits(orc, "error")) out$Oracle <- data.frame(est = orc$est, se = orc$se, ok = TRUE)

  ## --- Tidy output --------------------------------------------------------
  res <- dplyr::bind_rows(lapply(names(out), function(m) {
    d <- out[[m]]
    data.frame(design_id = s$design_id, dataset = dataset, species = seq_len(J), method = m,
               est = d$est, se = d$se, converged = d$ok,
               truth_direct = truth_direct, truth_total = truth_total)
  }))
  list(results = res, diagnostics = diag,
       data_summary = list(a = cov$a, r2_z_on_clim = unname(sapply(summary(lm(cov$z ~ cov$climate)), `[[`, "r.squared"))))
}

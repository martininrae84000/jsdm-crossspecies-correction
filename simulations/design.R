## =============================================================================
## Simulation design: a reference scenario and one-at-a-time variations
## (the 8 scenarios of the paper, numbered 1 to 8).
## Each scenario answers one question and comes with an expected result (see below).
## In run_simulations.R, DESIGN_IDS can be "all" (or "core") or a list of numbers such as "1,2,5".
## =============================================================================

reference <- list(
  family          = "binary",   # "binary" or "negbin"
  sX              = 1,          # strength of confounding (Hui et al.: 0, 0.5, 1)
  delta_mean      = 1,          # mean response to the confounders (Hui et al.: 1)
  delta_sd        = 1,          # heterogeneity of responses to the confounders (Hui et al.: 0.1)
  eps_sd          = 0.5,        # own variation of the confounders (Hui et al.: 0.1)
  J               = 20,         # number of species (Hui et al.: 10)
  beta_mean       = 1,          # mean direct effect (Hui et al.: 1)
  beta_sd         = 0.5,        # heterogeneity of direct effects (Hui et al.: 0.1)
  rho             = 0,          # correlation between direct effects and responses to confounders
  effect_scale    = 1,          # multiplies all species effects (beta and delta: means and sds);
                                # 0.5 for counts, so that exp(eta) stays realistic (see below)
  x_scale         = "coarse",   # "coarse" or "fine" measured covariate (Hui et al.: both)
  generator       = "basis",    # "basis" (Hui et al.) or "matern"
  K_fit           = 2           # number of latent factors fitted (true number of confounders = 2)
)

#' @param label Short label (used in the CSV files).
#' @param set "core" (all scenarios of the paper).
#' @param paper_label Label used in the figures (English, for the paper).
#' @param seed_id Number used to build the random seed of each dataset
#'   (set.seed(1e5 * seed_id + dataset)). It keeps the scenario numbers of an earlier
#'   version of the design, so that the simulated data are exactly those analysed in
#'   the paper.
variant <- function(label, set = "core", paper_label = label, seed_id, ...) {
  s <- modifyList(reference, list(...))
  s$label <- label
  s$set <- set
  s$paper_label <- paper_label
  s$seed_id <- seed_id
  s
}

design_list <- list(
  ## ======================= Core scenarios (main text): 1 to 8 =======================
  ## 1 All conditions met: correction recovers the direct effect, RSFA the total effect.
  variant("reference", "core", "Reference", seed_id = 1),
  ## 2 No confounding: the correction should not harm (Wald test ~5% rejections).
  variant("sX = 0", "core", "No confounding", seed_id = 2, sX = 0),
  ## 3 C1 violated (homogeneous responses to the confounders): little or no gain expected.
  variant("delta_sd = 0.1", "core", "Homogeneous responses (C1)", seed_id = 4, delta_sd = 0.1),
  ## 4 C2 violated (confounders with little variation of their own): little or no gain expected.
  variant("eps_sd = 0.1", "core", "Weak residual variation (C2)", seed_id = 6, eps_sd = 0.1),
  ## 5 C3 violated: species-level bias about -beta_sd * rho * (s_j - mean(s)),
  ##   i.e. RMSE increased by about beta_sd * rho = 0.3.
  variant("rho = 0.6", "core", "Endogeneity (C3)", seed_id = 11, rho = 0.6),
  ## 6 Original design of Hui et al. (2024): C1 and C2 violated, no gain expected.
  variant("Hui et al., sX = 1", "core", "Design of Hui et al. (2024)", seed_id = 17,
          J = 10, delta_sd = 0.1, eps_sd = 0.1, beta_sd = 0.1),
  ## 7 Reference with counts. With effect_scale = 1, species with delta close to 3 get linear
  ##   predictors up to ~23 (counts of several billions) and the fits fail. All effects are
  ##   halved for counts (beta ~ N(0.5, 0.25^2), delta ~ N(0.5, 0.5^2)): conditions C1-C3 keep
  ##   the same relative strength.
  variant("NB reference", "core", "Counts (negative binomial)", seed_id = 19,
          family = "negbin", effect_scale = 0.5),
  ## 8 Reference with Matern (exponential) random fields instead of the basis functions of
  ##   Hui et al.: checks that the conclusions do not depend on the generator. With these
  ##   fields, SFA can be nearly as biased as RSFA for the direct effects.
  variant("Matern generator", "core", "Gaussian random fields", seed_id = 13, generator = "matern")
)
for (i in seq_along(design_list)) design_list[[i]]$design_id <- i

design_table <- dplyr::bind_rows(lapply(design_list, function(s) as.data.frame(s[c(
  "design_id", "seed_id", "set", "label", "paper_label", "family", "sX", "delta_sd", "eps_sd", "J",
  "beta_sd", "rho", "effect_scale", "x_scale", "generator", "K_fit")])))

#' Scenario numbers from a DESIGN_IDS string: "all", "core" or "1,2,5".
parse_design_ids <- function(x) {
  x <- trimws(x)
  if (x == "all") return(design_table$design_id)
  if (x == "core") return(design_table$design_id[design_table$set == x])
  as.integer(strsplit(x, ",")[[1]])
}

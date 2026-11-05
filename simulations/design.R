## =============================================================================
## Simulation design: a reference scenario and one-at-a-time variations.
##
## Scenarios are split in two sets:
##   - "core"       (8 scenarios): main text of the paper. Each one answers one
##                  question and comes with an expected result (see below).
##   - "supplement" (15 scenarios): further robustness checks, not reported in the paper.
## Core scenarios are numbered 1 to 8, supplementary ones 9 to 23. In run_simulations.R,
## DESIGN_IDS can be "core", "supplement", "all" or a list of numbers such as "1,2,16".
## =============================================================================

## Note on Spatial+. Spatial+ replaces the climate by its residual after regression on
## basis functions 2..spatialplus_k. With the basis-function generator, a coarse-scale climate
## is built on functions 2..5 and is therefore removed almost entirely (only the small iid
## noise remains), and a fine-scale climate is not removed at all: Spatial+ is degenerate in
## both cases and its fit may fail to converge. It is informative only with the Matern
## generator; it is fitted only in the supplementary scenario "Matern, fine-scale X".
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
  K_fit           = 2,          # number of latent factors fitted (true number of confounders = 2)
  fit_indfa       = TRUE,
  fit_spatialplus = FALSE,      # see note above; switched on in selected scenarios
  spatialplus_k   = 10          # Spatial+: climate residualised on basis functions 2..k
)

#' @param label Short label (used in the CSV files).
#' @param set "core" or "supplement".
#' @param paper_label Label used in the figures (English, for the paper).
#' @param seed_id Number used to build the random seed of each dataset
#'   (set.seed(1e5 * seed_id + dataset)). It keeps the scenario numbers of an earlier
#'   version of the design, so that the simulated data are exactly those analysed in
#'   the paper.
variant <- function(label, set = "supplement", paper_label = label, seed_id, ...) {
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
  ##   the same relative strength. Scenario 23 keeps the original scale of Hui et al.
  variant("NB reference", "core", "Counts (negative binomial)", seed_id = 19,
          family = "negbin", effect_scale = 0.5),
  ## 8 Reference with Matern (exponential) random fields instead of the basis functions of
  ##   Hui et al.: checks that the conclusions do not depend on the generator. With these
  ##   fields, SFA can be nearly as biased as RSFA for the direct effects.
  variant("Matern generator", "core", "Gaussian random fields", seed_id = 13, generator = "matern"),

  ## ===================== Supplementary scenarios: 9 to 23 =====================
  ## Intermediate values of C0, C1, C2, C3
  variant("sX = 0.5",           seed_id = 3,  sX = 0.5),
  variant("delta_sd = 0.5",     seed_id = 5,  delta_sd = 0.5),
  variant("eps_sd = 0.25",      seed_id = 7,  eps_sd = 0.25),
  variant("rho = 0.3",          seed_id = 10, rho = 0.3),
  ## Number of species and of factors
  variant("J = 10",             seed_id = 8,  J = 10),
  variant("J = 50",             seed_id = 9,  J = 50),
  variant("K = 1",              seed_id = 15, K_fit = 1),
  variant("K = 3",              seed_id = 16, K_fit = 3),
  ## Spatial scale and generator (Spatial+ is informative only with the Matern generator)
  variant("fine-scale X",       seed_id = 12, x_scale = "fine"),
  variant("Matern, fine-scale X", seed_id = 14, generator = "matern", x_scale = "fine", fit_spatialplus = TRUE),
  ## Design of Hui et al. without confounding
  variant("Hui et al., sX = 0", seed_id = 18, J = 10, delta_sd = 0.1, eps_sd = 0.1, beta_sd = 0.1, sX = 0),
  ## Further count scenarios
  variant("NB sX = 0",          seed_id = 20, family = "negbin", effect_scale = 0.5, sX = 0),
  variant("NB J = 10",          seed_id = 21, family = "negbin", effect_scale = 0.5, J = 10),
  variant("NB rho = 0.6",       seed_id = 22, family = "negbin", effect_scale = 0.5, rho = 0.6),
  variant("NB Hui et al.",      seed_id = 23, family = "negbin", J = 10, delta_sd = 0.1, eps_sd = 0.1, beta_sd = 0.1)
)
for (i in seq_along(design_list)) design_list[[i]]$design_id <- i

design_table <- dplyr::bind_rows(lapply(design_list, function(s) as.data.frame(s[c(
  "design_id", "seed_id", "set", "label", "paper_label", "family", "sX", "delta_sd", "eps_sd", "J",
  "beta_sd", "rho", "effect_scale", "x_scale", "generator", "K_fit", "fit_indfa", "fit_spatialplus")])))

#' Scenario numbers from a DESIGN_IDS string: "core", "supplement", "all" or "1,2,16".
parse_design_ids <- function(x) {
  x <- trimws(x)
  if (x == "all") return(design_table$design_id)
  if (x %in% c("core", "supplement")) return(design_table$design_id[design_table$set == x])
  as.integer(strsplit(x, ",")[[1]])
}

## Quick test: one dataset for the scenarios given in DESIGN_IDS (default 1), timings and checks.
suppressPackageStartupMessages({ library(TMB); library(dplyr); library(autoFRK); library(fmesher); library(MASS) })
source("R/setup.R"); source("R/generate.R"); source("R/fit.R"); source("R/one_dataset.R"); source("design.R")
load_tmb_models()
S <- setup_sites()
r0 <- function(x) if (is.null(x)) "-" else round(x)   # NULL when a model is not fitted
ids <- parse_design_ids(Sys.getenv("DESIGN_IDS", "1"))   # "1", "1,17", "core", ...
for (id in ids) {
  t0 <- Sys.time()
  r <- run_one(design_list[[id]], dataset = 1, S)
  cat("\n=== design", id, ":", design_list[[id]]$label, "| total time",
      round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s\n")
  cat("fit times (s): SFA", r0(r$diagnostics$sfa$time), " indFA", r0(r$diagnostics$indfa$time),
      " Spatial+", r0(r$diagnostics$spatialplus$time), "\n")
  cat("convergence SFA", r$diagnostics$sfa$convergence, " pdHess", r$diagnostics$sfa$pdHess,
      " | invariance max diff", signif(r$diagnostics$sfa$invariance_maxdiff, 3),
      " | Wald p", signif(r$diagnostics$sfa$wald_p, 3),
      " | R2(delta ~ L)", round(r$diagnostics$sfa$r2_delta_on_L, 2),
      " | max |loading| SFA", round(r$diagnostics$sfa$max_abs_loading, 1), " indFA", round(r$diagnostics$indfa$max_abs_loading, 1), "\n")
  print(r$results %>% group_by(method) %>%
          summarise(rmse_direct = sqrt(mean((est - truth_direct)^2)),
                    bias_direct = mean(est - truth_direct),
                    rmse_total = sqrt(mean((est - truth_total)^2)),
                    mean_se = mean(se), .groups = "drop") %>%
          mutate(across(where(is.numeric), ~ round(.x, 3))))
}

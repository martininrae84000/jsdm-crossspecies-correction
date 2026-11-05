## =============================================================================
## Summarise the simulation results.
##   Rscript summarise.R
##
## Methods are always compared on the SAME datasets: within a scenario, a dataset is
## used only if every method in the comparison converged on it (complete-case
## comparison, Morris et al. 2019, Stat. Med.). Otherwise a method that fails on the
## hardest datasets would look better than one that does not (e.g. the correction,
## which depends on SFA, against the oracle). The number of datasets used and
## excluded is reported for each scenario.
##
## Two comparisons are made:
##   - MAIN (main text): SFA, RSFA, Cross_RSFA (correction), Oracle. The first three
##     come from the same SFA fit, so the common datasets are those where SFA converged.
##   - ALL (supplement): every method fitted in the scenario (also Independent,
##     Cross_Independent, SpatialPlus). The common datasets are those where all of
##     them converged; with counts, the independent factor model often fails, so this
##     set can be small (see n_datasets).
##
## Outputs (folder results/):
##   main_figure.png / .pdf      RMSE of the direct effect, core scenarios, 4 methods,
##                               with +-2 Monte Carlo standard errors
##   main_rmse.csv               the same numbers as a table (scenario x method)
##   main_coverage.csv           correction: coverage of 95% intervals,
##                               rejection rate of the Wald test of gamma = 0
##   main_methods.csv            all criteria, MAIN comparison, all scenarios
##   summary_methods.csv         all criteria, ALL comparison, all scenarios
##   supplement_figure.png       RMSE of the direct effect, ALL comparison
##   summary_convergence.csv     proportion of converged fits per method and scenario
##                               (over all datasets) and number of datasets kept
##   summary_diagnostics.csv     SFA diagnostics (convergence, Hessian, invariance,
##                               Wald test, loadings)
##
## Options (environment variables):
##   STRICT=1    a fit counts as converged only if nlminb converged AND the Hessian is
##               positive definite (default: nlminb convergence only; see below)
##   KEEP_ALL=1  no exclusion at all (every dataset is used for every method)
##
## Why the default ignores the Hessian. With the lower-triangular loadings of Hui et al.
## (diagonal >= 1e-3), the first loading often sits at its bound when species 1 hardly
## loads on factor 1: nlminb converges but the Hessian is not positive definite. Point
## estimates are not affected (the correction depends on the column space of the
## loadings only); the proportion of such fits is reported (prop_pdHess).
## =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })
source("design.R")
options(width = 150)

strict   <- Sys.getenv("STRICT", "0") == "1"
keep_all <- Sys.getenv("KEEP_ALL", "0") == "1"

## -----------------------------------------------------------------------------
## Read the results
## -----------------------------------------------------------------------------
files <- list.files("results", pattern = "^dataset_.*\\.rds$", recursive = TRUE, full.names = TRUE)
runs <- lapply(files, readRDS)
ok_run <- vapply(runs, function(r) is.null(r$error), logical(1))
cat(sum(ok_run), "datasets read,", sum(!ok_run), "failed\n")
runs <- runs[ok_run]
if (!length(runs)) stop("No results found in results/")

## Convergence status of the fit each estimate comes from, rebuilt from the diagnostics
fit_status <- bind_rows(lapply(runs, function(r) {
  st <- function(d, methods) {
    if (is.null(d) || !is.null(d$error)) return(NULL)
    data.frame(method = methods, optimizer_ok = d$convergence == 0, pdHess = isTRUE(d$pdHess))
  }
  bind_rows(st(r$diagnostics$sfa, c("SFA", "RSFA", "Cross_RSFA")),
            st(r$diagnostics$indfa, c("Independent", "Cross_Independent")),
            st(r$diagnostics$spatialplus, "SpatialPlus"),
            data.frame(method = "Oracle", optimizer_ok = TRUE, pdHess = TRUE)) %>%
    mutate(design_id = r$results$design_id[1], dataset = r$results$dataset[1])
}))

res_all <- bind_rows(lapply(runs, `[[`, "results")) %>%
  dplyr::select(-any_of("converged")) %>%
  left_join(fit_status, by = c("design_id", "dataset", "method")) %>%
  mutate(optimizer_ok = coalesce(optimizer_ok, FALSE), pdHess = coalesce(pdHess, FALSE),
         ok = if (keep_all) TRUE else if (strict) optimizer_ok & pdHess else optimizer_ok) %>%
  left_join(design_table %>% dplyr::select(design_id, set, label, paper_label, family), by = "design_id")

## A method whose fit failed with an error on a dataset has no rows there: count it as not converged
n_runs <- res_all %>% distinct(design_id, dataset) %>% count(design_id, name = "n_runs")
expected <- res_all %>% distinct(design_id, method) %>%
  inner_join(res_all %>% distinct(design_id, dataset), by = "design_id", relationship = "many-to-many")
status <- expected %>%
  left_join(res_all %>% group_by(design_id, dataset, method) %>%
              summarise(ok = all(ok), optimizer_ok = all(optimizer_ok), pdHess = all(pdHess), .groups = "drop"),
            by = c("design_id", "dataset", "method")) %>%
  mutate(across(c(ok, optimizer_ok, pdHess), ~ coalesce(.x, FALSE)))

## Datasets on which every method of `methods` (among those fitted in the scenario) is ok
common_datasets <- function(methods) {
  status %>% filter(method %in% methods) %>%
    group_by(design_id, dataset) %>% summarise(all_ok = all(ok), .groups = "drop") %>%
    filter(all_ok) %>% dplyr::select(design_id, dataset)
}

## -----------------------------------------------------------------------------
## Performance criteria on a given set of datasets
## -----------------------------------------------------------------------------
z <- qnorm(0.975)
B_mc <- 200   # bootstrap replicates for the Monte Carlo standard error of the RMSE
mc_rmse <- function(err, dataset) {
  sse <- tapply(err^2, dataset, sum); n <- tapply(err, dataset, length)
  if (length(sse) < 2) return(NA_real_)
  sd(replicate(B_mc, { i <- sample(length(sse), replace = TRUE); sqrt(sum(sse[i]) / sum(n[i])) }))
}
criteria <- function(d) {
  set.seed(2024)
  d %>% group_by(design_id, set, label, paper_label, family, method) %>%
    summarise(
      n_datasets   = n_distinct(dataset),
      bias_direct  = mean(est - truth_direct),
      rmse_direct  = sqrt(mean((est - truth_direct)^2)),
      mcse_rmse_direct = mc_rmse(est - truth_direct, dataset),
      bias_total   = mean(est - truth_total),
      rmse_total   = sqrt(mean((est - truth_total)^2)),
      mean_se      = mean(se, na.rm = TRUE),
      cover_direct = mean(abs(est - truth_direct) <= z * se, na.rm = TRUE),
      cover_total  = mean(abs(est - truth_total)  <= z * se, na.rm = TRUE),
      .groups = "drop") %>%
    left_join(n_runs, by = "design_id") %>%
    mutate(n_excluded = n_runs - n_datasets) %>%
    arrange(design_id, method)
}

## -----------------------------------------------------------------------------
## MAIN comparison: SFA, RSFA, correction, oracle on common datasets
## -----------------------------------------------------------------------------
main_methods <- c(SFA = "SFA", RSFA = "RSFA", Cross_RSFA = "Cross-species correction",
                  Oracle = "Oracle (confounders observed)")
common_main <- common_datasets(names(main_methods))
res_main <- res_all %>% filter(method %in% names(main_methods)) %>% semi_join(common_main, by = c("design_id", "dataset"))
main_tab <- criteria(res_main)
write.csv(main_tab, "results/main_methods.csv", row.names = FALSE)

## -----------------------------------------------------------------------------
## ALL comparison: every method fitted in the scenario, on common datasets
## -----------------------------------------------------------------------------
common_all <- status %>% group_by(design_id, dataset) %>% summarise(all_ok = all(ok), .groups = "drop") %>%
  filter(all_ok) %>% dplyr::select(design_id, dataset)
res_cmp <- res_all %>% semi_join(common_all, by = c("design_id", "dataset"))
summary_methods <- criteria(res_cmp)
write.csv(summary_methods, "results/summary_methods.csv", row.names = FALSE)

## Convergence per method (over all datasets) and number of datasets kept
summary_conv <- status %>% group_by(design_id, method) %>%
  summarise(n_runs = n(), prop_converged = mean(optimizer_ok), prop_pdHess = mean(pdHess),
            prop_ok = mean(ok), .groups = "drop") %>%
  left_join(common_main %>% count(design_id, name = "n_kept_main"), by = "design_id") %>%
  left_join(common_all %>% count(design_id, name = "n_kept_all"), by = "design_id") %>%
  mutate(across(c(n_kept_main, n_kept_all), ~ coalesce(.x, 0L))) %>%
  left_join(design_table %>% dplyr::select(design_id, label), by = "design_id") %>%
  relocate(label, .after = design_id) %>% arrange(design_id, method)
write.csv(summary_conv, "results/summary_convergence.csv", row.names = FALSE)

## -----------------------------------------------------------------------------
## Diagnostics of the SFA fit. Convergence proportions over all datasets; the other
## diagnostics (Wald test, loadings) over the datasets of the MAIN comparison.
## -----------------------------------------------------------------------------
diag_sfa <- bind_rows(lapply(runs, function(r) {
  d <- r$diagnostics$sfa
  if (is.null(d) || !is.null(d$error)) return(NULL)
  data.frame(design_id = r$results$design_id[1], dataset = r$results$dataset[1],
             convergence = d$convergence, pdHess = d$pdHess, time_sfa = d$time,
             invariance_maxdiff = d$invariance_maxdiff, wald_p = d$wald_p,
             kappa_LCL = d$kappa_LCL, r2_delta_on_L_min = min(d$r2_delta_on_L))
}))
diag_conv <- diag_sfa %>% group_by(design_id) %>%
  summarise(n_runs = n(), prop_converged = mean(convergence == 0), prop_pdHess = mean(pdHess),
            mean_time_sfa_sec = mean(time_sfa), .groups = "drop")
diag_main <- diag_sfa %>% semi_join(common_main, by = c("design_id", "dataset")) %>% group_by(design_id) %>%
  summarise(n_kept = n(),
            max_invariance_diff = max(invariance_maxdiff),
            reject_gamma0_5pct  = mean(wald_p < 0.05, na.rm = TRUE),
            median_kappa_LCL    = median(kappa_LCL),
            median_r2_delta_L   = median(r2_delta_on_L_min), .groups = "drop")
summary_diag <- diag_conv %>% left_join(diag_main, by = "design_id") %>%
  left_join(design_table %>% dplyr::select(design_id, label), by = "design_id") %>% relocate(label, .after = design_id)
write.csv(summary_diag, "results/summary_diagnostics.csv", row.names = FALSE)

## -----------------------------------------------------------------------------
## Main-text outputs (core scenarios)
## -----------------------------------------------------------------------------
core <- main_tab %>% filter(set == "core")
if (nrow(core)) {
  core <- core %>%
    mutate(scenario = factor(paste0(design_id, ". ", paper_label),
                             levels = rev(paste0(design_table$design_id, ". ", design_table$paper_label))),
           Method = factor(main_methods[method], levels = main_methods))

  main_rmse <- core %>%
    mutate(cell = sprintf("%.2f (%.2f)", rmse_direct, mcse_rmse_direct)) %>%
    dplyr::select(design_id, paper_label, n_datasets, n_excluded, method, cell) %>%
    pivot_wider(names_from = method, values_from = cell) %>%
    dplyr::select(design_id, scenario = paper_label, n_datasets, n_excluded, any_of(names(main_methods))) %>%
    arrange(design_id)
  write.csv(main_rmse, "results/main_rmse.csv", row.names = FALSE)

  p <- ggplot(core, aes(x = rmse_direct, y = scenario, colour = Method, shape = Method)) +
    geom_errorbarh(aes(xmin = pmax(0, rmse_direct - 2 * mcse_rmse_direct), xmax = rmse_direct + 2 * mcse_rmse_direct),
                   height = 0, position = position_dodge(width = 0.6), alpha = 0.6) +
    geom_point(size = 2.8, position = position_dodge(width = 0.6)) +
    scale_colour_manual(values = c("#9CA3AF", "#C2410C", "#1D4ED8", "#111827")) +
    scale_shape_manual(values = c(16, 17, 15, 4)) +
    expand_limits(x = 0) +
    labs(x = "RMSE of the direct effect (±2 Monte Carlo standard errors)", y = NULL,
         colour = NULL, shape = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank()) +
    guides(colour = guide_legend(nrow = 2))
  ggsave("results/main_figure.png", p, width = 7.2, height = 4.6, dpi = 300, bg = "white")
  ggsave("results/main_figure.pdf", p, width = 7.2, height = 4.6)

  main_coverage <- core %>% filter(method == "Cross_RSFA") %>%
    dplyr::select(design_id, scenario = paper_label, n_datasets, coverage = cover_direct) %>%
    left_join(summary_diag %>% dplyr::select(design_id, reject_gamma0_5pct), by = "design_id") %>%
    arrange(design_id)
  write.csv(main_coverage, "results/main_coverage.csv", row.names = FALSE)

  cat("\n=== Main text: RMSE of the direct effect (Monte Carlo SE), core scenarios ===\n")
  cat("    (methods compared on the same datasets: those where the SFA fit is kept; n_excluded = datasets left out)\n")
  print(as.data.frame(main_rmse), row.names = FALSE)
  cat("\n=== Main text: correction, coverage of 95% intervals and Wald test of gamma = 0 ===\n")
  print(as.data.frame(main_coverage %>% mutate(across(where(is.double), ~ round(.x, 2)))), row.names = FALSE)
} else {
  cat("\nNo core scenario in the results: main-text outputs not produced.\n")
}

## -----------------------------------------------------------------------------
## Supplement figure (ALL comparison)
## -----------------------------------------------------------------------------
all_methods <- c(main_methods, Independent = "Independent FA",
                 Cross_Independent = "Cross-species correction (independent FA)",
                 SpatialPlus = "Spatial+")
supp <- summary_methods %>%
  mutate(scenario = factor(paste0(design_id, ". ", paper_label, " (n = ", n_datasets, ")"),
                           levels = rev(unique(paste0(design_id, ". ", paper_label, " (n = ", n_datasets, ")")[order(design_id)]))),
         Method = factor(all_methods[method], levels = all_methods))
if (nrow(supp)) {
  p2 <- ggplot(supp, aes(x = pmin(rmse_direct, 8), y = scenario, colour = Method, shape = Method)) +
    geom_point(size = 2.2, position = position_dodge(width = 0.6)) +
    scale_colour_manual(values = c("#9CA3AF", "#C2410C", "#1D4ED8", "#111827",
                                   "#A16207", "#0F766E", "#7C3AED"), drop = FALSE) +
    scale_shape_manual(values = c(16, 17, 15, 4, 1, 0, 5), drop = FALSE) +
    expand_limits(x = 0) +
    labs(x = "RMSE of the direct effect (truncated at 8); n = datasets where all methods converged",
         y = NULL, colour = NULL, shape = NULL) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
    guides(colour = guide_legend(nrow = 3))
  ggsave("results/supplement_figure.png", p2, width = 8.5,
         height = 2 + 0.45 * n_distinct(supp$design_id), dpi = 200, bg = "white", limitsize = FALSE)
}

## -----------------------------------------------------------------------------
## Warnings on the number of datasets kept
## -----------------------------------------------------------------------------
low <- summary_conv %>% distinct(design_id, label, n_runs, n_kept_main, n_kept_all) %>%
  filter(n_kept_main < 0.8 * n_runs | n_kept_all < 0.5 * n_runs)
if (nrow(low)) {
  cat("\nScenarios with many datasets excluded (see summary_convergence.csv):\n")
  print(as.data.frame(low), row.names = FALSE)
}
cat("\nOptions: STRICT =", strict, "| KEEP_ALL =", keep_all, "\n")
cat("Files written in results/: main_figure.png/.pdf, main_rmse.csv, main_coverage.csv, main_methods.csv,",
    "summary_methods.csv, supplement_figure.png, summary_convergence.csv, summary_diagnostics.csv\n")

## =============================================================================
## Figure of the simulation results in the paper (RMSE of the direct effect,
## core scenarios 1-8, four methods, +-2 Monte Carlo standard errors).
##
## Run from the simulations/ folder, after summarise.R:
##   Rscript make_figure.R
## Reads results/main_methods.csv (methods compared on the same datasets) and writes
## results/fig_sim_core.pdf and results/fig_sim_core.png.
## =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

tab <- read.csv(file.path("results", "main_methods.csv"))

## Scenario labels, in the order of the paper
lab <- c(`1` = "1. Reference", `2` = "2. No confounding", `3` = "3. Homogeneous responses (C1)",
         `4` = "4. Weak residual variation (C2)", `5` = "5. Endogeneity (C3)",
         `6` = "6. Design of Hui et al. (2024)", `7` = "7. Counts (negative binomial)",
         `8` = "8. Gaussian random fields")
## Methods shown, with their labels
meth <- c(SFA = "SFA", RSFA = "RSFA", Cross_RSFA = "Cross-species correction", Oracle = "Oracle")

d <- tab %>%
  filter(design_id %in% as.integer(names(lab)), method %in% names(meth)) %>%
  mutate(scen = factor(lab[as.character(design_id)], levels = rev(lab)),   # scenario 1 at the top
         Method = factor(meth[method], levels = meth))

xmax <- 10   # x-axis truncated at 10 (SFA and RSFA have very wide intervals in scenario 4)
p <- ggplot(d, aes(rmse_direct, scen, colour = Method, shape = Method)) +
  ## +-2 Monte Carlo standard errors of the RMSE (bootstrap over datasets, see summarise.R)
  geom_errorbarh(aes(xmin = pmax(0, rmse_direct - 2 * mcse_rmse_direct),
                     xmax = pmin(xmax, rmse_direct + 2 * mcse_rmse_direct)), height = 0,
                 position = position_dodge(width = 0.6), alpha = 0.6) +
  geom_point(size = 2.6, position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = c("#9CA3AF", "#C2410C", "#1D4ED8", "#111827")) +
  scale_shape_manual(values = c(16, 17, 15, 4)) +
  scale_x_continuous(limits = c(0, xmax), breaks = 0:5 * 2) +
  labs(x = expression("RMSE of the direct effect (" %+-% "2 Monte Carlo standard errors)"),
       y = NULL, colour = NULL, shape = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank())

ggsave(file.path("results", "fig_sim_core.pdf"), p, width = 7.2, height = 4.8)
ggsave(file.path("results", "fig_sim_core.png"), p, width = 7.2, height = 4.8, dpi = 200, bg = "white")
cat("Figure written to results/fig_sim_core.pdf and .png\n")

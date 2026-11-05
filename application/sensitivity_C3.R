## =============================================================================
## Sensitivity analysis for condition (C3), butterfly data.
## Main fit: species coefficients as random effects, K = 2, duplicated species removed.
##
## Run from the application/ folder, after application_butterflies.R:
##   Rscript sensitivity_C3.R [fit.rds] [covariates.csv]
## Defaults: results/fit_random_K2_dedup.rds and results/butterflies_covariates_dedup.csv.
## Writes results/butterflies_sensitivity_C3.csv and results/fig_butterflies_sensitivity_C3.pdf/.png.
##
## Principle (Section 3.5 of the paper). If (C3) fails, the direct effects contain a part
## linear in the loadings, beta_p = beta_cs_p + Lambda theta_p, and the corrected estimates
## miss it. The sensitivity parameter is r, the (multiple) correlation across species between
## the direct effects and the loadings, so that
##   ||C Lambda theta_p||^2 = ||C beta_cs_p||^2 * r^2 / (1 - r^2)   (C = centring matrix).
## For a given r, two directions for theta_p:
##   - "aligned":    theta_p proportional to gamma_hat_p (direct effects aligned with the
##                   estimated bias);
##   - "worst case": the direction that maximises |shift of the community mean|,
##                   |lbar' theta_p| = ||C Lambda theta_p|| * sqrt(lbar' (Lambda' C Lambda)^{-1} lbar).
## The community-mean direct effect is then mean(beta_cs_p) + lbar' theta_p.
## =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })

args <- commandArgs(TRUE)
fit_file <- if (length(args) >= 1) args[1] else file.path("results", "fit_random_K2_dedup.rds")
cov_file <- if (length(args) >= 2) args[2] else file.path("results", "butterflies_covariates_dedup.csv")

f  <- readRDS(fit_file)
L  <- f$est$loadings_mat                 # estimated loadings (J x K)
Bs <- f$est$betas_rsr                    # RSFA coefficients (J x P; column 1 = intercept)
J  <- nrow(L)

## Covariates of the main model (columns 2 to 4 of the coefficient matrix) and the standard
## errors of their corrected community means, read from the output of application_butterflies.R.
cov_codes <- c("climate", "blwood", "conwood")
cov_names <- c("Growing degree days", "Broadleaved woodland", "Coniferous woodland")
ct <- read.csv(cov_file)
se_mean <- sapply(cov_codes, function(v)
  ct$se_mean_Corrected[ct$model == f$model & ct$K == f$K & ct$covariate == v])

C    <- diag(J) - 1 / J                  # centring matrix over species
S    <- t(L) %*% C %*% L                 # Lambda' C Lambda (K x K)
lbar <- colMeans(L)                      # mean loadings over species
rs   <- seq(-0.6, 0.6, by = 0.01)        # values of the sensitivity parameter r

out <- list()
for (p in 1:3) {
  b    <- Bs[, p + 1]                    # RSFA coefficients of the covariate
  g    <- solve(S, t(L) %*% C %*% b)     # gamma_hat (equation 8)
  bcs  <- drop(b - L %*% g)              # corrected coefficients (equation 9)
  m0   <- mean(bcs)                      # corrected community mean (r = 0)
  v_cs <- sum((bcs - m0)^2)              # ||C beta_cs||^2
  ## Aligned direction: theta = c * gamma_hat, with c chosen so that the correlation equals r.
  ## Since C beta_cs is orthogonal to the columns of C Lambda,
  ## r = c ||C L g|| / sqrt(||C beta_cs||^2 + c^2 ||C L g||^2).
  q     <- drop(t(g) %*% S %*% g)        # ||C L g||^2
  c_al  <- rs * sqrt(v_cs) / (sqrt(q) * sqrt(1 - rs^2))
  sh_al <- c_al * sum(lbar * g)          # shift of the community mean: lbar' theta
  ## Least favourable direction (maximal shift of the community mean for a given r).
  sh_wc <- sign(rs) * sqrt(v_cs * rs^2 / (1 - rs^2) * drop(t(lbar) %*% solve(S, lbar)))
  out[[p]] <- data.frame(covariate = cov_names[p], r = rs, mean_corrected = m0,
                         aligned = m0 + sh_al, worst_case = m0 + sh_wc, se = se_mean[p])
}
d <- bind_rows(out)
write.csv(d, file.path("results", "butterflies_sensitivity_C3.csv"), row.names = FALSE)

## Values at a few r, printed in the console.
tab <- d %>% group_by(covariate) %>% summarise(
  mean_corrected = first(mean_corrected),
  aligned_r_0.3  = aligned[abs(r - 0.3) < 1e-9], aligned_r_m0.3 = aligned[abs(r + 0.3) < 1e-9],
  worst_r_0.3    = worst_case[abs(r - 0.3) < 1e-9], worst_r_m0.3 = worst_case[abs(r + 0.3) < 1e-9],
  .groups = "drop")
print(as.data.frame(tab), digits = 3)

## Figure: community-mean direct effect as a function of r, for both directions.
## Shaded band: +-1.96 standard errors of the corrected mean (r treated as known).
dl <- pivot_longer(d, c(aligned, worst_case), names_to = "direction")
dl$covariate <- factor(dl$covariate, cov_names)
dl$direction <- factor(dl$direction, c("aligned", "worst_case"),
                       c("Aligned with the estimated bias", "Least favourable direction"))
p <- ggplot(dl, aes(r, value, linetype = direction)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_ribbon(data = dl[dl$direction == "Aligned with the estimated bias", ],
              aes(ymin = value - 1.96 * se, ymax = value + 1.96 * se), fill = "#1D4ED8", alpha = 0.12) +
  geom_line(colour = "#1D4ED8") + facet_wrap(~ covariate, scales = "free_y") +
  labs(x = "Correlation across species between direct effects and loadings",
       y = "Community-mean direct effect", linetype = NULL) +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path("results", "fig_butterflies_sensitivity_C3.pdf"), p, width = 7.5, height = 3.2)
ggsave(file.path("results", "fig_butterflies_sensitivity_C3.png"), p, width = 7.5, height = 3.2, dpi = 200, bg = "white")
cat("Written: results/butterflies_sensitivity_C3.csv and results/fig_butterflies_sensitivity_C3.pdf/.png\n")

## =============================================================================
## Fixed objects shared by all scenarios: sites, spatial basis, SPDE mesh,
## Cholesky factors for the Matern (exponential) generator.
## =============================================================================

#' @param n_sites Number of sites sampled from the Great Britain butterfly coordinates
#'   (400 in Hui et al. 2024).
#' @param coords_file Path to Butterfly_LatLon_scaled.csv (see ../data/README.md).
#' @param range_coarse,range_fine Range parameters (exponential covariance) of the
#'   coarse- and fine-scale Gaussian fields used by the "matern" generator.
setup_sites <- function(n_sites = 400,
                        coords_file = file.path("..", "data", "Butterfly_LatLon_scaled.csv"),
                        range_coarse = 0.4, range_fine = 0.06) {
  coords_all <- read.csv(coords_file)
  set.seed(042024)                                   # same seed as Hui et al.
  coords <- dplyr::slice_sample(coords_all, n = n_sites)

  ## Multi-resolution thin-plate spline basis (as in Hui et al.)
  basis <- as.matrix(autoFRK::mrts(knot = coords[, c("x", "y")], k = 20))

  ## SPDE mesh for the spatial factors (fmesher instead of INLA; same matrices)
  mesh <- fmesher::fm_rcdt_2d_inla(loc = as.matrix(coords[, c("x", "y")]), refine = FALSE)
  fem  <- fmesher::fm_fem(mesh, order = 2)

  ## Cholesky factors for the Matern generator (exponential covariance = Matern 1/2)
  D <- as.matrix(dist(coords[, c("x", "y")]))
  chol_exp <- function(r) t(chol(exp(-D / r) + diag(1e-6, n_sites)))

  list(coords = coords, basis = basis,
       mesh = mesh, fem = fem,
       L_coarse = chol_exp(range_coarse), L_fine = chol_exp(range_fine),
       n_sites = n_sites)
}

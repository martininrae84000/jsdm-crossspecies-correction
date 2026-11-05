## =============================================================================
## Check the results folder against the current design.R.
##   source("check_results.R")      (from the project folder)
##
## 1. Lists the core scenarios of the current design.R (8 in the current version).
## 2. For every results/design_<id> folder: number of datasets done, and whether the
##    data stored in dataset_1.rds are those that the CURRENT design.R would simulate for
##    this scenario number (same true coefficients). "NO" means that the folder was
##    computed with an older numbering of the scenarios: do not use it as such.
## =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(autoFRK); library(fmesher) })
source("R/setup.R"); source("R/generate.R"); source("design.R")

core <- design_table %>% filter(set == "core") %>% dplyr::select(design_id, label)
cat("Core scenarios in the current design.R:", nrow(core), "\n")
print(as.data.frame(core), row.names = FALSE)

S <- setup_sites()
dirs <- list.dirs("results", recursive = FALSE)
if (!length(dirs)) stop("No results/design_<id> folder found")
chk <- bind_rows(lapply(dirs, function(d) {
  id <- as.integer(sub(".*design_", "", d))
  files <- list.files(d, pattern = "^dataset_.*\\.rds$", full.names = TRUE)
  ok <- NA
  f1 <- file.path(d, "dataset_1.rds")
  if (id <= length(design_list) && file.exists(f1)) {
    s <- design_list[[id]]
    set.seed(1e5 * (if (is.null(s$seed_id)) s$design_id else s$seed_id) + 1)
    cov <- gen_covariates(s, S); co <- gen_coefficients(s)
    r <- readRDS(f1)
    stored <- r$results$truth_direct[r$results$method == "Oracle"]
    ok <- length(stored) == length(co$beta) && isTRUE(all.equal(stored, co$beta))
  }
  data.frame(design_id = id,
             label = if (id <= length(design_list)) design_list[[id]]$label else "(not in design.R)",
             set = if (id <= length(design_list)) design_list[[id]]$set else NA,
             n_datasets = length(files),
             matches_current_design = ifelse(is.na(ok), "?", ifelse(ok, "yes", "NO")))
})) %>% arrange(design_id)
cat("\nResults folders:\n")
print(as.data.frame(chk), row.names = FALSE)
if (any(chk$matches_current_design == "NO"))
  cat("\nWARNING: folders marked NO were computed with another numbering of the scenarios.\n")
missing_core <- setdiff(core$design_id, chk$design_id)
if (length(missing_core)) cat("\nCore scenarios not computed yet:", paste(missing_core, collapse = ", "), "\n")

## =============================================================================
## Run the simulations.
##
## From the simulations/ folder:
##   Rscript run_simulations.R
## Settings (environment variables):
##   DESIGN_IDS    scenarios to run: "all" (the 8 scenarios of the paper; "core" is
##                 the same) or numbers such as "1,2,5"               (default "core")
##   NUM_DATASETS  datasets per scenario                               (default 100)
##   CORES         number of cores                                     (default: all - 1)
## Each dataset is saved in results/design_<id>/dataset_<k>.rds; existing files
## are skipped, so runs can be interrupted and resumed.
## =============================================================================

suppressPackageStartupMessages({
  library(TMB); library(dplyr); library(autoFRK); library(fmesher); library(MASS); library(parallel)
})

source("R/setup.R"); source("R/generate.R"); source("R/fit.R"); source("R/one_dataset.R")
source("design.R")

## Compile (once) and load the four TMB models of Hui et al. (2024)
load_tmb_models()

ids <- parse_design_ids(Sys.getenv("DESIGN_IDS", "core"))
n_datasets <- as.integer(Sys.getenv("NUM_DATASETS", "100"))
cores <- as.integer(Sys.getenv("CORES", max(1, parallel::detectCores() - 1)))

S <- setup_sites()

jobs <- expand.grid(dataset = seq_len(n_datasets), design_id = ids)
jobs$file <- file.path("results", paste0("design_", jobs$design_id), paste0("dataset_", jobs$dataset, ".rds"))
jobs <- jobs[!file.exists(jobs$file), ]
for (d in unique(dirname(jobs$file))) dir.create(d, recursive = TRUE, showWarnings = FALSE)
cat(nrow(jobs), "datasets to run on", cores, "cores\n")

run_job <- function(i) {
  s <- design_list[[jobs$design_id[i]]]
  r <- tryCatch(run_one(s, jobs$dataset[i], S), error = function(e) list(error = conditionMessage(e)))
  saveRDS(r, jobs$file[i])
  cat(format(Sys.time()), " design", jobs$design_id[i], " dataset", jobs$dataset[i],
      if (!is.null(r$error)) paste(" ERROR:", r$error) else "", "\n",
      file = "results/progress.log", append = TRUE)
  invisible(NULL)
}

## Parallel execution.
## - Linux / macOS: forked processes (mclapply).
## - Windows (or FORCE_PSOCK=1): a cluster of separate R sessions (PSOCK), since
##   mclapply cannot run in parallel on Windows. Each worker loads the packages,
##   the R functions and the compiled TMB models.
use_psock <- .Platform$OS.type == "windows" || Sys.getenv("FORCE_PSOCK", "0") == "1"
if (nrow(jobs) > 0) {
  if (cores == 1) {
    invisible(lapply(seq_len(nrow(jobs)), run_job))
  } else if (!use_psock) {
    invisible(parallel::mclapply(seq_len(nrow(jobs)), run_job, mc.cores = cores, mc.preschedule = FALSE))
  } else {
    cl <- parallel::makeCluster(cores)
    project_dir <- getwd()
    parallel::clusterExport(cl, c("jobs", "project_dir"))
    parallel::clusterEvalQ(cl, {
      setwd(project_dir)                 # workers start in the project folder
      suppressPackageStartupMessages({
        library(TMB); library(dplyr); library(autoFRK); library(fmesher); library(MASS)
      })
      source("R/setup.R"); source("R/generate.R"); source("R/fit.R"); source("R/one_dataset.R")
      source("design.R")
      load_tmb_models(compile_if_needed = FALSE)
      S <- setup_sites()
      NULL
    })
    parallel::clusterExport(cl, c("run_job"))
    invisible(parallel::parLapplyLB(cl, seq_len(nrow(jobs)), run_job))
    parallel::stopCluster(cl)
  }
}
cat("Done.\n")

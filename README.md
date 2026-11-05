# Cross-species correction for spatial confounding in joint species distribution models

Code for the simulations and the application of the paper
*Borrowing information across species to correct spatial confounding in joint species distribution models*
<!-- TODO: author, arXiv link and DOI once available -->.

The cross-species correction is a post-processing step applied to a fitted spatial factor analysis (SFA)
JSDM. For each covariate, the estimated species coefficients (unrestricted or restricted, RSFA) are
regressed on the estimated loadings across species, with an intercept, and the part explained by the
loadings is removed:

    gamma_hat_p = (Lambda' C Lambda)^{-1} Lambda' C beta*_p      (C: centring over species)
    beta_cs_p   = beta*_p - Lambda gamma_hat_p

Standard errors are obtained by the delta method from the joint covariance of the coefficients and the
loadings returned by TMB, and a Wald test of `gamma_p = 0` tests for confounding.

The models are fitted with the TMB templates of
[Hui, Vu & Hooten (2024)](https://github.com/fhui28/spatialconfounding_JSDM), *Spatial confounding in
joint species distribution models*, Methods in Ecology and Evolution 15, 1906–1921.

## Repository structure

```
data/            data files to download (not included), see data/README.md
simulations/     simulation study (Section 4 of the paper)
application/     butterflies of Great Britain (Section 5 of the paper)
```

### `simulations/`

| File | Content |
|---|---|
| `design.R` | Simulation design: 8 core scenarios (reported in the paper) and 15 supplementary ones |
| `run_simulations.R` | Runs the simulations (parallel; can be interrupted and resumed) |
| `summarise.R` | Performance criteria (RMSE, bias, coverage, Wald test), methods compared on the same datasets |
| `make_figure.R` | Figure of the paper (RMSE of the direct effect, scenarios 1–8) |
| `test_one.R` | Quick test on one dataset |
| `check_results.R` | Checks the `results/` folder against `design.R` |
| `R/setup.R` | Sites, spatial basis functions, SPDE mesh, Matérn generator |
| `R/generate.R` | Covariates, confounders, species coefficients, responses |
| `R/fit.R` | TMB fits (SFA, RSFA, independent factor analysis), oracle, cross-species correction, Wald test |
| `R/one_dataset.R` | One dataset: simulation, all methods, tidy results |
| `tmb/*.cpp` | TMB templates of Hui et al. (2024) for the simulations (see below) |

### `application/`

| File | Content |
|---|---|
| `application_butterflies.R` | Data preparation, SFA fits, cross-species correction, tables and plots |
| `sensitivity_C3.R` | Sensitivity analysis for condition (C3) |
| `R/correction.R` | Functions of the correction (joint covariance, correction, Wald test, diagnostics) |
| `tmb/TMB_SFA_random.cpp` | SFA with species coefficients as random effects (model of Hui et al. for this case study) |
| `tmb/TMB_SFA_fixed.cpp` | SFA with species coefficients as fixed effects |

## Requirements

R (≥ 4.3) with the packages `TMB`, `dplyr`, `tidyr`, `ggplot2`, `Matrix`, `MASS`, `autoFRK`, `fmesher`
(and optionally `INLA`, which gives the same mesh as `fmesher`). TMB needs a C++ compiler
(Rtools on Windows).

```r
install.packages(c("TMB", "dplyr", "tidyr", "ggplot2", "Matrix", "MASS", "autoFRK", "fmesher"))
```

## Data

The data are not included. Download the three files of the butterfly data of Wilkinson et al. (2019)
from [Zenodo (record 1452066)](https://zenodo.org/records/1452066), folder `Datasets/Butterflies`, and put
them in `data/` (see `data/README.md`). The simulations use only the site coordinates.

## Running the code

All scripts are run from their own folder (`simulations/` or `application/`).

**Simulations**

```sh
cd simulations
Rscript test_one.R                                 # quick check on one dataset
Rscript run_simulations.R                          # core scenarios, 100 datasets each
Rscript summarise.R                                # tables and figures in results/
Rscript make_figure.R                              # figure of the paper
```

Settings are passed as environment variables, e.g. `DESIGN_IDS=core NUM_DATASETS=100 CORES=8`
(`DESIGN_IDS` can be `core`, `supplement`, `all` or a list such as `1,2,5`). Each dataset is saved in
`results/design_<id>/dataset_<k>.rds`; existing files are skipped. With 100 datasets, the 8 core
scenarios take about one to two days on 8 cores.

**Application**

```sh
cd application
Rscript application_butterflies.R                              # main model: random coefficients, K = 2
MODELS=random,fixed K_VALUES=2,3 Rscript application_butterflies.R   # robustness (Table 4)
QUADRATIC=1 Rscript application_butterflies.R                  # quadratic effect of growing degree days
Rscript sensitivity_C3.R                                       # sensitivity analysis for (C3)
```

From R, use for instance `Sys.setenv(MODELS = "random,fixed", K_VALUES = "2,3")` before `source()`.
Each fit takes about 15 minutes or more and is saved in `results/`; it is not refitted if the file exists.
By default, the two pairs of species with identical records in the data file are reduced to one species
each (`DROP_DUPLICATES=0` keeps them, as in Hui et al.).

## Changes to the TMB templates of Hui et al. (2024)

The templates are those of https://github.com/fhui28/spatialconfounding_JSDM (v0.6). The likelihood,
the parameters and the reported coefficients and loadings are unchanged; only the reporting of the
latent factors was modified, to reduce computing time and memory:

| File here | Original file | Change |
|---|---|---|
| `simulations/tmb/TMBINLA_SFA_binary.cpp` | `simulation/binary_N400/TMBINLA_SFA.cpp` | `ADREPORT` of the latent factors disabled |
| `simulations/tmb/TMBINLA_SFA_negbin.cpp` | `simulation/negbin_N400/TMBINLA_SFA.cpp` | idem |
| `simulations/tmb/TMB_indFA_binary.cpp` | `simulation/binary_N400/TMB_indFA.cpp` | idem |
| `simulations/tmb/TMB_indFA_negbin.cpp` | `simulation/negbin_N400/TMB_indFA.cpp` | idem |
| `application/tmb/TMB_SFA_random.cpp` | `application_butterflies/TMB_SFA.cpp` | latent factors `REPORT`ed; restricted factors computed in R |
| `application/tmb/TMB_SFA_fixed.cpp` | `simulation/binary_N400/TMBINLA_SFA.cpp` | idem |

The file names differ because TMB names a compiled model after its file. The application templates are
compiled with the TMBad framework, which reduces the memory needed for the 2609 sites from more than
6 GB to about 0.5 GB.

## Licence

GNU General Public License v3 (see `LICENSE`), as the code of Hui et al. (2024) on which it builds.

## References

- Hui, F. K. C., Vu, Q. & Hooten, M. B. (2024). Spatial confounding in joint species distribution models.
  *Methods in Ecology and Evolution*, 15, 1906–1921. https://doi.org/10.1111/2041-210X.14420
- Wilkinson, D. P., Golding, N., Guillera-Arroita, G., Tingley, R. & McCarthy, M. A. (2019).
  A comparison of joint species distribution models for presence–absence data.
  *Methods in Ecology and Evolution*, 10, 198–211. https://doi.org/10.1111/2041-210X.13106
- Kristensen, K., Nielsen, A., Berg, C. W., Skaug, H. & Bell, B. M. (2016). TMB: automatic
  differentiation and Laplace approximation. *Journal of Statistical Software*, 70(5), 1–21.

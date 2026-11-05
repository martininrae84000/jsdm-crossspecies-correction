# Data

The data files are not included in this repository. Download them from the data archive of
Wilkinson et al. (2019), [Zenodo record 1452066](https://zenodo.org/records/1452066), folder
`Datasets/Butterflies`, and put the following three files in this folder:

| File | Content | Used by |
|---|---|---|
| `Butterfly_PA.csv` | presence–absence of 55 butterfly species at 2609 sites in Great Britain | application |
| `Butterfly_LatLon_scaled.csv` | site coordinates (rescaled), columns `x` and `y` | simulations and application |
| `Butterfly_Cov.csv` | site covariates, including growing degree days (`climate`), broadleaved (`blwood`) and coniferous (`conwood`) woodland | application |

The data originate from the Millennium Atlas of Butterflies in Britain and Ireland
(Asher et al. 2001) and were analysed by Ovaskainen et al. (2016), Wilkinson et al. (2019) and
Hui, Vu & Hooten (2024).

In `Butterfly_PA.csv`, two pairs of the 45 species with more than 150 presences have identical
records at every site (small blue and small copper; common blue and Camberwell beauty). By default,
`application/application_butterflies.R` keeps only the first column of each pair (43 species).

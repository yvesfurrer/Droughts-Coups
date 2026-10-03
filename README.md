# Droughts and Coups d'État in Africa

Replication materials for the master's seminar paper "Droughts and Coups d'État in Africa: A Quantitative Analysis of Coup Attempts from 1990 to 2023" (University of Lucerne, MA PPE, 2026).

## Research question

Droughts are climatological events that cannot be caused by coups, which rules out reverse causality and makes a causal reading plausible. Because the design cannot rule out time-varying confounders, the estimates are nevertheless reported as associations. 

## Main finding

Sustained, year-round drought (SPEI-12, December, lagged one year) is associated with a significantly higher probability of a coup attempt. The temporal scale of drought measurement is decisive: short seasonal droughts (SPEI-3, SPEI-6) show no association, while the annual scale carries a strong signal. The association is robust to a rare-events correction, country fixed effects, temporal-dependence controls, population-weighted drought exposure, and models without economic controls.

## Contents

- `droughts_coups_analysis.R`: full R script (data preparation, analysis, robustness checks, tables)
- `Output/`: generated tables (Tables 1-4 and appendix tables A1-A5 of the paper)

## Data sources (not included, freely available)

- Coup attempts: Powell & Thyne (2011), Coup d'État Dataset, country-year version (data and codebook), version V2024.10.30, http://www.uky.edu/~clthyn2/coup_data/home.html
- Drought (SPEI-3/6/12/24): SPEIbase v2.10, global NetCDF files, https://spei.csic.es/spei_database/
- Population weights (robustness check): Gridded Population of the World v4, population density 2020 (10 arc-minutes), downloaded automatically via the R package `geodata` on the first run
- World Bank, World Development Indicators, downloaded automatically via the R package `WDI` on the first run and cached locally, https://databank.worldbank.org/source/world-development-indicators
  - GDP per capita growth (annual %): NY.GDP.PCAP.KD.ZG
  - Military expenditure (% of general government expenditure): MS.MIL.XPND.ZS
  - Military expenditure (% of GDP, robustness check): MS.MIL.XPND.GD.ZS
  - Agriculture, forestry, and fishing, value added growth (annual %): NV.AGR.TOTL.KD.ZG
  - Agriculture, forestry, and fishing, value added (% of GDP): NV.AGR.TOTL.ZS

## Folder structure and expected file names

All paths in the script are relative to the repository folder, so the script must be run from there (or with the working directory set to it). The following files have to be placed manually before the first run:

```
Droughts-Coups/
├── droughts_coups_analysis.R
├── Data/
│   └── raw/
│       ├── powell_thyne_ccode_year.csv   # Powell & Thyne, country-year file, saved as CSV
│       ├── spei03.nc                     # SPEIbase, 3-month time scale
│       ├── spei06.nc                     # SPEIbase, 6-month time scale
│       ├── spei12.nc                     # SPEIbase, 12-month time scale
│       └── spei24.nc                     # SPEIbase, 24-month time scale
└── Output/                               # created by the script
```

The coup file must contain one row per country and year with the columns `ccode` (Correlates of War country code), `year` and `coup1` to `coup4` (1 = failed attempt, 2 = successful coup), as in the original country-year file.

All other files are created by the script:

- `Data/raw/wdi_raw.rds` and `Data/raw/wdi_mil_gdp.rds`: local cache of the World Bank downloads
- `Data/raw/population/pop/gpw_v4_population_density_rev11_2020_10m.tif`: population raster downloaded by `geodata`
- `Data/processed/panel.rds`: merged country-year panel with all derived variables
- `Data/processed/spei12_popweighted.rds`: population-weighted SPEI-12 series (robustness check)
- `Output/*.docx`: tables of the paper

## How to replicate

1. Download the coup data and the four SPEI files from the sources above and save them under the file names listed in the folder structure
2. Run the script with `RUN_DATA_PREP <- TRUE` (first run only). An internet connection is required for the World Bank API and the population raster
3. For later runs set `RUN_DATA_PREP <- FALSE`. The script then loads `Data/processed/panel.rds` and starts directly with the analysis
4. Tables are written to `Output/`. Additional checks that are not part of a table are printed to the console

## Software

R 4.4.2. Required packages are installed automatically in section 0 of the script. The `sessionInfo()` call at the end of the script prints the package versions used.

## Contact

Yves Furrer

Master's Programme Philosophy, Politics and Economics, University of Lucerne

yves.furrer@stud.unilu.ch

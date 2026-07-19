# Droughts and Coups in Africa
Replication materials for the master's seminar paper "Droughts and Political Instability in Africa" (University of Lucerne, MA PPE, 2026).

## Research question
What is the effect of droughts on coup attempts in Africa? The analysis covers 52 African countries, 1990 to 2023, using a country-year panel.

## Main finding
Sustained, year-round drought (SPEI-12, December, lagged one year) is associated with a significantly higher probability of a coup attempt. The temporal scale of drought measurement is decisive: short seasonal droughts (SPEI-3, SPEI-6) show no association, while the annual scale carries a strong signal. The association is robust to a rare-events correction, country fixed effects, temporal-dependence controls, population-weighted drought exposure, and models without economic controls.

## Contents
- `droughts_coups_analysis.R`: full R script (data preparation, analysis, robustness checks, tables)
- Tables 1-4 and appendix tables A1-A4 of the paper are written to `Output/` when the script is run

## Data sources (not included, freely available)
- Coup attempts: Powell & Thyne (2011), Coup d'État Dataset (data and codebook), http://www.uky.edu/~clthyn2/coup_data/home.html
- Drought (SPEI-3/6/12/24): SPEIbase v2.10, https://spei.csic.es/spei_database/
- Population weights (robustness check): Gridded Population of the World v4, downloaded automatically into `Data/raw/` via the R package `geodata` on the first run
- World Bank, World Development Indicators, https://databank.worldbank.org/source/world-development-indicators
  - GDP per capita growth (annual %): NY.GDP.PCAP.KD.ZG
  - Military expenditure (% of general government expenditure): MS.MIL.XPND.ZS
  - Agriculture, forestry, and fishing, value added growth (annual %): NV.AGR.TOTL.KD.ZG
  - Agriculture, forestry, and fishing, value added (% of GDP): NV.AGR.TOTL.ZS

## How to replicate
1. Download the coup and SPEI data from the sources above into `Data/raw/`
2. Adjust the file paths in section 1 of the script
3. Run the script with `RUN_DATA_PREP <- TRUE` (first run only; an internet connection is required for the World Bank API and the population raster)
4. Intermediate files (`Data/processed/panel.rds`, `Data/processed/spei12_popweighted.rds`) are created automatically; set `RUN_DATA_PREP <- FALSE` for later runs
5. Tables are written to `Output/`

## Software
R 4.4.2. Required packages are installed automatically in section 0 of the script. Package versions are documented by the `sessionInfo()` call at the end of the script.

## Contact
Yves Furrer,
yves.furrer@stud.unilu.ch,
Master's Program:
Philosophy, Politics and Economy,
University of Lucerne

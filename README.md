# Droughts and Coups in Africa

Replication materials for the master's seminar paper "Droughts and Political Instability in Africa" (University of Lucerne, MA PPE, 2026).

## Research question
What is the effect of droughts on coup attempts in Africa? The analysis covers 52 African countries, 1990 to 2023, using a country-year panel.

## Main finding
Sustained, year-round drought (SPEI-12, December, lagged one year) significantly raises the probability of a coup attempt. Short seasonal droughts (SPEI-3, SPEI-6) show no effect.

## Contents
- `droughts_coups_v5.R`: full R script (data preparation and analysis)

## Data sources (not included, freely available)
- Coup attempts: Powell & Thyne (2011), Coup d'État Dataset, http://www.uky.edu/~clthyn2/coup_data/home.html
- Drought (SPEI-3/6/12/24): SPEIbase v2.10, https://spei.csic.es/spei_database/
- World Bank, World Development Indicators, https://databank.worldbank.org/source/world-development-indicators
  - GDP per capita growth (annual %): NY.GDP.PCAP.KD.ZG
  - Military expenditure (% of general government expenditure): MS.MIL.XPND.ZS
  - Agriculture, forestry, and fishing, value added growth (annual %): NV.AGR.TOTL.KD.ZG
  - Agriculture, forestry, and fishing, value added (% of GDP): NV.AGR.TOTL.ZS
 
## How to replicate
1. Download the raw data from the sources above into `Data/raw/`
2. Adjust the file paths in section 1 of the script
3. Run the script with `RUN_DATA_PREP <- TRUE` (first run only)
4. Tables are written to `Output/`

## Software
R 4.4.2. Required packages are installed automatically in section 0 of the script.

## Contact
Yves Furrer
yves.furrer@stud.unilu.ch
Master's Programm:
Philosophy, Politics and Economy
University of Lucerne

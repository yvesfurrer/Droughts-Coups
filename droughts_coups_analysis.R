###############################################################################
# Masterseminarpaper: Droughts & Coups in Africa
# Author: Yves Furrer
#
# NOTE ON RE-RUNS: set RUN_DATA_PREP <- FALSE after the first successful
# run; the script then loads the saved panel and jumps to the analysis.
#
# This script was created with the assistance of Claude AI:
# Anthropic, model Claude Fable 5
# All code was reviewed and tested by the author.
###############################################################################

# ------------------------------------------------------------------ #
# 0. Packages
# ------------------------------------------------------------------ #
pkgs <- c("WDI", "terra", "sf", "rnaturalearth", "rnaturalearthdata",
          "ncdf4", "lubridate", "dplyr", "tidyr", "countrycode",
          "plm", "lmtest", "sandwich", "logistf", "modelsummary",
          "flextable", "pandoc", "survival",
          "geodata", "exactextractr")
new <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
if (length(new)) install.packages(new)

library(WDI); library(terra); library(sf); library(survival)
library(rnaturalearth); library(rnaturalearthdata)
library(ncdf4); library(lubridate); library(dplyr); library(tidyr)
library(countrycode); library(plm); library(lmtest); library(sandwich)
library(logistf); library(modelsummary); library(pandoc)
library(geodata); library(exactextractr); library(flextable)

# ------------------------------------------------------------------ #
# 1. Paths & parameters
#    All paths are relative to the project root; run the script from
#    the repository folder (or set the working directory there first).
# ------------------------------------------------------------------ #
path_coup   <- "Data/raw/powell_thyne_ccode_year.csv"
path_spei03 <- "Data/raw/spei03.nc"
path_spei06 <- "Data/raw/spei06.nc"
path_spei12 <- "Data/raw/spei12.nc"
path_spei24 <- "Data/raw/spei24.nc"
path_wdi    <- "Data/raw/wdi_raw.rds"       # local cache of WDI downloads
path_panel  <- "Data/processed/panel.rds"   # relative to getwd()

year_start <- 1990
year_end   <- 2023

options(timeout = 300)
RUN_DATA_PREP <- TRUE   # set to FALSE after the first successful run

if (RUN_DATA_PREP) {
  
  # ---------------------------------------------------------------- #
  # 2. World Bank indicators (cached locally after first download)
  # ---------------------------------------------------------------- #
  if (file.exists(path_wdi)) {
    wdi_raw <- readRDS(path_wdi)
    gdp <- wdi_raw$gdp; mil <- wdi_raw$mil
    agr <- wdi_raw$agr; agshare <- wdi_raw$agshare
  } else {
    gdp     <- WDI(country = "all", indicator = c(gdp_growth = "NY.GDP.PCAP.KD.ZG"),
                   start = year_start, end = year_end)
    mil     <- WDI(country = "all", indicator = c(mil_exp = "MS.MIL.XPND.ZS"),
                   start = year_start, end = year_end)
    agr     <- WDI(country = "all", indicator = c(agr_growth = "NV.AGR.TOTL.KD.ZG"),
                   start = year_start, end = year_end)
    agshare <- WDI(country = "all", indicator = c(agr_share = "NV.AGR.TOTL.ZS"),
                   start = year_start, end = year_end)
    dir.create(dirname(path_wdi), recursive = TRUE, showWarnings = FALSE)
    saveRDS(list(gdp = gdp, mil = mil, agr = agr, agshare = agshare), path_wdi)
  }
  
  clean_wdi <- function(d) d %>%
    filter(!is.na(iso3c), iso3c != "") %>%
    distinct(iso3c, year, .keep_all = TRUE)
  
  gdp <- clean_wdi(gdp); mil <- clean_wdi(mil)
  agr <- clean_wdi(agr); agshare <- clean_wdi(agshare)
  
  wdi <- gdp %>%
    left_join(mil     %>% select(iso3c, year, mil_exp),    by = c("iso3c", "year")) %>%
    left_join(agr     %>% select(iso3c, year, agr_growth), by = c("iso3c", "year")) %>%
    left_join(agshare %>% select(iso3c, year, agr_share),  by = c("iso3c", "year")) %>%
    # countrycode warns about WDI aggregate codes (AFE, ARB, SSF, WLD, ...);
    # these are not countries and are removed by the continent filter below
    mutate(continent = countrycode(iso3c, "iso3c", "continent")) %>%
    filter(continent == "Africa") %>%
    select(iso3c, country, year, gdp_growth, mil_exp, agr_growth, agr_share)
  
  # ---------------------------------------------------------------- #
  # 3. Coup data (Powell & Thyne, ccode/year format)
  # ---------------------------------------------------------------- #
  coup_raw <- read.csv(path_coup, stringsAsFactors = FALSE)
  
  coup <- coup_raw %>%
    mutate(across(starts_with("coup"), ~replace(., is.na(.), 0)),
           coup_attempt = as.integer(coup1 %in% c(1, 2) | coup2 %in% c(1, 2) |
                                       coup3 %in% c(1, 2) | coup4 %in% c(1, 2)),
           # countrycode warns about Gleditsch-Ward codes outside Africa and
           # about Somaliland (no ISO code); both are dropped by the filter below
           iso3c = countrycode(ccode, "gwn", "iso3c")) %>%
    filter(!is.na(iso3c), year >= year_start, year <= year_end) %>%
    group_by(iso3c, year) %>%
    summarise(coup_attempt = as.integer(any(coup_attempt == 1)), .groups = "drop")
  
  # ---------------------------------------------------------------- #
  # 4. SPEI (all four time scales): aggregate grid -> country-year
  #    spei12_dec : Dec value (Jan-Dec, full calendar year)  [MAIN IV]
  #    spei3_oct  : Oct value (Aug-Oct, acute shock)         [time scale]
  #    spei6_oct  : Oct value (May-Oct, growing season)      [time scale]
  #    spei6_mean : annual mean of monthly SPEI-6            [alt. measure]
  #    spei24_dec : Dec value (two full years)               [time scale]
  # ---------------------------------------------------------------- #
  africa <- ne_countries(scale = "medium", continent = "Africa",
                         returnclass = "sf")
  africa$iso3c <- countrycode(africa$admin, "country.name", "iso3c")
  africa <- africa[!is.na(africa$iso3c), ]
  africa_v <- vect(africa)
  
  extract_spei <- function(path) {
    r  <- rast(path)
    nc <- nc_open(path)
    t_val  <- ncvar_get(nc, "time")
    t_unit <- ncatt_get(nc, "time", "units")$value
    nc_close(nc)
    origin <- as.Date(sub(".*since\\s*", "", t_unit))
    dts <- if (grepl("month", t_unit)) origin %m+% months(as.integer(t_val))
    else origin + as.numeric(t_val)
    keep <- which(year(dts) >= year_start & year(dts) <= year_end)
    r <- r[[keep]]; dts <- dts[keep]
    r  <- crop(r, ext(africa_v))
    et <- terra::extract(r, africa_v, fun = mean, na.rm = TRUE)
    et$iso3c <- africa$iso3c
    et %>% select(-ID) %>%
      pivot_longer(-iso3c, names_to = "layer", values_to = "val") %>%
      group_by(iso3c) %>%
      mutate(date = rep(dts, times = 1)) %>%   # layer order == date order
      ungroup() %>%
      mutate(year = year(date), month = month(date)) %>%
      select(iso3c, year, month, val)
  }
  
  s03 <- extract_spei(path_spei03)
  s06 <- extract_spei(path_spei06)
  s12 <- extract_spei(path_spei12)
  s24 <- extract_spei(path_spei24)
  
  spei_long <- s12 %>% filter(month == 12) %>%
    select(iso3c, year, spei12_dec = val) %>%
    left_join(s06 %>% group_by(iso3c, year) %>%
                summarise(spei6_oct  = val[month == 10][1],
                          spei6_mean = mean(val, na.rm = TRUE),
                          .groups = "drop"),               by = c("iso3c", "year")) %>%
    left_join(s03 %>% filter(month == 10) %>%
                select(iso3c, year, spei3_oct = val),      by = c("iso3c", "year")) %>%
    left_join(s24 %>% filter(month == 12) %>%
                select(iso3c, year, spei24_dec = val),     by = c("iso3c", "year"))
  
  # ---------------------------------------------------------------- #
  # 5. Merge into country-year panel (base = WDI African countries)
  # ---------------------------------------------------------------- #
  panel <- wdi %>%
    left_join(coup,      by = c("iso3c", "year")) %>%
    left_join(spei_long, by = c("iso3c", "year")) %>%
    mutate(coup_attempt = replace(coup_attempt, is.na(coup_attempt), 0L))
  
  # ---------------------------------------------------------------- #
  # 6. Derived variables (ALL of them, so the saved panel is complete)
  # ---------------------------------------------------------------- #
  # 6a. drought indicators on the MAIN IV (SPEI-12, WMO thresholds)
  #     [WMO 2012, p. 4]: drought <= -1.0; moderate/severe/extreme
  panel <- panel %>%
    arrange(iso3c, year) %>%
    mutate(
      drought12 = as.integer(spei12_dec <= -1.0),
      drought12_cat = case_when(
        spei12_dec <= -2.0 ~ "extreme",
        spei12_dec <= -1.5 ~ "severe",
        spei12_dec <= -1.0 ~ "moderate",
        TRUE               ~ "none"),
      drought12_cat = factor(drought12_cat,
                             levels = c("none", "moderate", "severe", "extreme"))
    )
  
  # 6b. winsorised growth (1st/99th pct) to tame outliers
  q  <- quantile(panel$gdp_growth, c(.01, .99), na.rm = TRUE)
  qa <- quantile(panel$agr_growth, c(.01, .99), na.rm = TRUE)
  panel <- panel %>%
    mutate(gdp_growth_w = pmin(pmax(gdp_growth, q[1]),  q[2]),
           agr_growth_w = pmin(pmax(agr_growth, qa[1]), qa[2]))
  
  # 6c. onset + all lags (grouped by country, dplyr::lag!)
  panel <- panel %>%
    group_by(iso3c) %>%
    mutate(
      drought12_onset    = as.integer(drought12 == 1 & dplyr::lag(drought12, default = 0) == 0),
      drought12_onset_l1 = dplyr::lag(drought12_onset, 1),
      drought12_cat_l1   = dplyr::lag(drought12_cat, 1),
      spei12_dec_l1      = dplyr::lag(spei12_dec, 1),
      spei12_dec_l2      = dplyr::lag(spei12_dec, 2),
      spei12_dec_l3      = dplyr::lag(spei12_dec, 3),
      spei3_oct_l1       = dplyr::lag(spei3_oct, 1),
      spei6_oct_l1       = dplyr::lag(spei6_oct, 1),
      spei6_mean_l1      = dplyr::lag(spei6_mean, 1),
      spei24_dec_l1      = dplyr::lag(spei24_dec, 1),
      gdp_growth_w_l1    = dplyr::lag(gdp_growth_w, 1)
    ) %>%
    ungroup()
  
  # 6d. drought duration on SPEI-12 (consecutive years <= -1.0)
  panel <- panel %>%
    arrange(iso3c, year) %>%
    group_by(iso3c) %>%
    mutate(run = cumsum(drought12 == 0 | is.na(drought12)),
           drought12_dur = ifelse(drought12 == 1,
                                  ave(drought12, run, FUN = cumsum), 0),
           drought12_dur_l1 = dplyr::lag(drought12_dur, 1)) %>%
    ungroup() %>% select(-run)
  
  # 6e. agricultural dependence (country mean of agr_share, median split)
  panel <- panel %>%
    group_by(iso3c) %>%
    mutate(agr_share_mean = mean(agr_share, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(agr_dep = as.integer(agr_share_mean >
                                  median(agr_share_mean, na.rm = TRUE)))
  
  # 6f. decade dummies (avoid year-dummy separation in logit)
  panel$decade <- cut(panel$year, breaks = c(1989, 1999, 2009, 2019, 2023),
                      labels = c("1990s", "2000s", "2010s", "2020s"))
  
  # ---------------------------------------------------------------- #
  # 7. Quick checks & save
  # ---------------------------------------------------------------- #
  print(summary(panel))
  print(table(panel$drought12, useNA = "ifany"))
  print(table(panel$coup_attempt, useNA = "ifany"))
  print(panel %>% filter(is.na(spei12_dec)) %>% count(iso3c))  # expect MUS, SYC
  
  dir.create(dirname(path_panel), recursive = TRUE, showWarnings = FALSE)
  saveRDS(panel, path_panel)
  
} else {
  
  # fast re-entry: load the prepared panel from disk
  panel <- readRDS(path_panel)
  
}

###############################################################################
# Part B: Descriptives and Analysis (main IV: SPEI-12, December, t-1)
###############################################################################

# states enter the panel upon independence: a coup attempt is not defined
# for years in which the state did not exist (pre-1993 Eritrea is covered
# by Ethiopia, pre-2011 South Sudan by Sudan in the coup data); placed in
# Part B so the filter also applies when the cached panel is loaded
panel <- panel %>%
  filter(!(iso3c == "ERI" & year < 1993),
         !(iso3c == "SSD" & year < 2011))

dir.create("Output", showWarnings = FALSE)
pdat <- pdata.frame(panel, index = c("iso3c", "year"))

# ------------------------------------------------------------------ #
# 8. Temporal dependence: years since the last coup attempt
#     (Carter & Signorino 2010: peace years t, t^2, t^3 as controls)
#     dplyr::lag() ensures the counter refers to coups strictly BEFORE
#     year t; spells without a prior coup are left-censored at 1990
# ------------------------------------------------------------------ #
panel <- panel %>%
  arrange(iso3c, year) %>%
  group_by(iso3c) %>%
  mutate(
    .cy            = cummax(ifelse(coup_attempt == 1, year, -Inf)),
    .last_coup_pre = dplyr::lag(.cy),
    coup_years     = ifelse(is.finite(.last_coup_pre),
                            year - .last_coup_pre,
                            year - (year_start - 1)),
    cy2 = coup_years^2 / 100,   # scaled nuisance terms, not interpreted
    cy3 = coup_years^3 / 1000
  ) %>%
  ungroup() %>%
  select(-.cy, -.last_coup_pre)

# ------------------------------------------------------------------ #
# 9. Descriptive statistics (Appendix tables A3, A4)
# ------------------------------------------------------------------ #

# A4: country list with number of coup attempts (console check; MUS and SYC
#     are dropped so the count matches the 52-country analysis universe
#     of Chapter 3.2 and the Word table written in section 17)
panel %>% filter(!iso3c %in% c("MUS", "SYC")) %>%
  group_by(country) %>%
  summarise(coups = sum(coup_attempt)) %>%
  arrange(desc(coups)) %>% print(n = 60)

# A3: summary statistics of the main variables
#     The variables are renamed explicitly. as.numeric() strips the WDI
#     variable labels, which would otherwise override these names and
#     report the winsorised series under the label of the raw series.
a3 <- panel %>%
  transmute(
    `Coup attempt (0/1)`                         = coup_attempt,
    `SPEI-12 (December)`                         = as.numeric(spei12_dec),
    `GDP p.c. growth (annual %, wins.)`          = as.numeric(gdp_growth_w),
    `Military expenditure (% of gov. exp.)`      = as.numeric(mil_exp),
    `Agri. value added growth (annual %, wins.)` = as.numeric(agr_growth_w))

datasummary_skim(a3, output = "Output/tableA3_summary.docx")

# ------------------------------------------------------------------ #
# 10. H3: Drought and economic conditions (plausibility check)
#    Panel OLS, country + year FE, SE clustered by country
#    Note: lag() inside plm() formulas is plm's panel-aware lag.
# ------------------------------------------------------------------ #

# (i) GDP p.c. growth, contemporaneous
m_h3 <- plm(gdp_growth ~ spei12_dec,
            data = pdat, model = "within", effect = "twoways")
coeftest(m_h3, vcov = vcovHC(m_h3, cluster = "group", type = "HC1"))

# (ii) GDP p.c. growth, contemporaneous + lag
m_h3_both <- plm(gdp_growth ~ spei12_dec + lag(spei12_dec, 1),
                 data = pdat, model = "within", effect = "twoways")
coeftest(m_h3_both, vcov = vcovHC(m_h3_both, cluster = "group", type = "HC1"))

# (iii) robustness: winsorised GDP growth
m_h3_w <- plm(gdp_growth_w ~ spei12_dec + lag(spei12_dec, 1),
              data = pdat, model = "within", effect = "twoways")
coeftest(m_h3_w, vcov = vcovHC(m_h3_w, cluster = "group", type = "HC1"))

# (iv) robustness: excluding oil states
oil <- c("NGA","AGO","LBY","DZA","GAB","GNQ","COG","TCD","SSD")
pdat_ag <- pdata.frame(panel %>% filter(!iso3c %in% oil),
                       index = c("iso3c", "year"))
m_h3_ag <- plm(gdp_growth_w ~ spei12_dec + lag(spei12_dec, 1),
               data = pdat_ag, model = "within", effect = "twoways")
coeftest(m_h3_ag, vcov = vcovHC(m_h3_ag, cluster = "group", type = "HC1"))

# (v) agricultural value added growth (closest to the channel)
m_h3_agr <- plm(agr_growth_w ~ spei12_dec + lag(spei12_dec, 1),
                data = pdat, model = "within", effect = "twoways")
coeftest(m_h3_agr, vcov = vcovHC(m_h3_agr, cluster = "group", type = "HC1"))

# ------------------------------------------------------------------ #
# 11. H1: Drought and coup attempts (main effect, SPEI-12)
# ------------------------------------------------------------------ #

# (i) contemporaneous (timing within year ambiguous; t-1 is main spec)
m_h1 <- glm(coup_attempt ~ spei12_dec + gdp_growth_w + mil_exp + decade,
            data = panel, family = binomial)
coeftest(m_h1, vcov = vcovCL(m_h1, cluster = ~ iso3c, data = panel))

# (ii) lagged (t-1)  [MAIN MODEL]
m_h1_l1 <- glm(coup_attempt ~ spei12_dec_l1 + gdp_growth_w + mil_exp + decade,
               data = panel, family = binomial)
coeftest(m_h1_l1, vcov = vcovCL(m_h1_l1, cluster = ~ iso3c, data = panel))

# (iii) total effect without economic controls
#       gdp_growth_w is measured in the coup year and is itself affected by
#       coups (post-treatment concern); since the SPEI is exogenous, the
#       model without economic controls identifies the total effect
m_h1_nc <- glm(coup_attempt ~ spei12_dec_l1 + decade,
               data = panel, family = binomial)
coeftest(m_h1_nc, vcov = vcovCL(m_h1_nc, cluster = ~ iso3c, data = panel))

# ------------------------------------------------------------------ #
# 12. H2: Drought intensity (SPEI-12 categories, lagged)
# ------------------------------------------------------------------ #
m_h2 <- glm(coup_attempt ~ drought12_cat_l1 + gdp_growth_w + mil_exp + decade,
            data = panel, family = binomial)
coeftest(m_h2, vcov = vcovCL(m_h2, cluster = ~ iso3c, data = panel))

# Firth correction for the sparse intensity cells: the severe and extreme
# categories rest on few country-years, exactly the setting the penalised
# likelihood is designed for
mdat_h2 <- panel %>%
  select(coup_attempt, drought12_cat_l1, gdp_growth_w, mil_exp, decade, iso3c) %>%
  na.omit()
m_h2_firth <- logistf(coup_attempt ~ drought12_cat_l1 + gdp_growth_w + mil_exp + decade,
                      data = mdat_h2)
summary(m_h2_firth)

# ------------------------------------------------------------------ #
# 13. H4: Timing of the effect (SPEI-12)
# ------------------------------------------------------------------ #

# (i) drought onset in t-1
m_h4 <- glm(coup_attempt ~ drought12_onset_l1 + gdp_growth_w + mil_exp + decade,
            data = panel, family = binomial)
coeftest(m_h4, vcov = vcovCL(m_h4, cluster = ~ iso3c, data = panel))

# (ii) lag structure t to t-2 (main text)
m_h4_l2 <- glm(coup_attempt ~ spei12_dec + spei12_dec_l1 + spei12_dec_l2 +
                 gdp_growth_w + mil_exp + decade,
               data = panel, family = binomial)
coeftest(m_h4_l2, vcov = vcovCL(m_h4_l2, cluster = ~ iso3c, data = panel))

# (iii) lag structure t to t-3 (appendix)
m_h4_lags <- glm(coup_attempt ~ spei12_dec + spei12_dec_l1 + spei12_dec_l2 +
                   spei12_dec_l3 + gdp_growth_w + mil_exp + decade,
                 data = panel, family = binomial)
coeftest(m_h4_lags, vcov = vcovCL(m_h4_lags, cluster = ~ iso3c, data = panel))

# ------------------------------------------------------------------ #
# 14. Robustness checks (all on the main model, SPEI-12 t-1)
# ------------------------------------------------------------------ #

# (i) Firth logit (rare events bias correction)
mdat <- panel %>%
  select(coup_attempt, spei12_dec_l1, gdp_growth_w, mil_exp, decade, year, iso3c) %>%
  na.omit()
m_firth <- logistf(coup_attempt ~ spei12_dec_l1 + gdp_growth_w + mil_exp + decade,
                   data = mdat)
summary(m_firth)

# (i-b) no-controls baseline on the identical estimation sample (N = 1257)
m_h1_nc_cc <- glm(coup_attempt ~ spei12_dec_l1 + decade,
                  data = mdat, family = binomial)
coeftest(m_h1_nc_cc, vcov = vcovCL(m_h1_nc_cc, cluster = ~ iso3c, data = mdat))
cat("N =", nobs(m_h1_nc_cc),
    "| coups in sample:", sum(mdat$coup_attempt),
    "of", sum(panel$coup_attempt), "\n")
cat("coups in the no-controls sample:",
    sum(panel$coup_attempt[!is.na(panel$spei12_dec_l1)]), "\n")

# (ii) interaction: drought x agricultural dependence
m_int <- glm(coup_attempt ~ spei12_dec_l1 * agr_dep + gdp_growth_w + mil_exp + decade,
             data = panel, family = binomial)
coeftest(m_int, vcov = vcovCL(m_int, cluster = ~ iso3c, data = panel))

# (iii) conditional logit (country fixed effects)
m_clogit <- clogit(coup_attempt ~ spei12_dec_l1 + gdp_growth_w + mil_exp +
                     decade + strata(iso3c),
                   data = panel, method = "efron")
summary(m_clogit)

# (iii-b) effective (informative) sample of the conditional logit:
#         countries without variation in coup_attempt contribute a likelihood
#         factor of 1; survival::clogit keeps them in nobs(), so the
#         informative N must be computed and reported separately
vars_cl <- c("coup_attempt", "spei12_dec_l1", "gdp_growth_w", "mil_exp",
             "decade", "iso3c")
cc_cl  <- panel[complete.cases(panel[, vars_cl]), vars_cl]
inf_cl <- cc_cl %>%
  group_by(iso3c) %>%
  filter(sum(coup_attempt) > 0 & sum(coup_attempt) < dplyr::n()) %>%
  ungroup()
n_eff_clogit    <- nrow(inf_cl)
n_countries_inf <- dplyr::n_distinct(inf_cl$iso3c)
cat("Conditional logit: effective N =", n_eff_clogit,
    "| informative countries =", n_countries_inf, "\n")

# identity check: estimates must be identical on the informative sample
m_clogit_inf <- clogit(coup_attempt ~ spei12_dec_l1 + gdp_growth_w + mil_exp +
                         decade + strata(iso3c),
                       data = inf_cl, method = "efron")
summary(m_clogit_inf)

# (iv) Sub-Saharan Africa only
north_africa <- c("DZA", "EGY", "LBY", "MAR", "TUN")
panel_ssa <- panel %>% filter(!iso3c %in% north_africa)
m_ssa <- glm(coup_attempt ~ spei12_dec_l1 + gdp_growth_w + mil_exp + decade,
             data = panel_ssa, family = binomial)
coeftest(m_ssa, vcov = vcovCL(m_ssa, cluster = ~ iso3c, data = panel_ssa))

# (v) pre-existing economic downturn (lagged growth as control)
m_pre <- glm(coup_attempt ~ spei12_dec_l1 + gdp_growth_w + gdp_growth_w_l1 +
               mil_exp + decade, data = panel, family = binomial)
coeftest(m_pre, vcov = vcovCL(m_pre, cluster = ~ iso3c, data = panel))

# (vi) alternative measurement of the same construct: annual mean of SPEI-6
m_alt <- glm(coup_attempt ~ spei6_mean_l1 + gdp_growth_w + mil_exp + decade,
             data = panel, family = binomial)
coeftest(m_alt, vcov = vcovCL(m_alt, cluster = ~ iso3c, data = panel))

# (vii) drought duration (consecutive years, SPEI-12 based)
m_dur <- glm(coup_attempt ~ drought12_dur_l1 + gdp_growth_w + mil_exp + decade,
             data = panel, family = binomial)
coeftest(m_dur, vcov = vcovCL(m_dur, cluster = ~ iso3c, data = panel))

# (viii) temporal dependence / coup trap: peace-years polynomial
#        (Carter & Signorino 2010); the time terms are nuisance
#        parameters and are not interpreted
m_temp <- glm(coup_attempt ~ spei12_dec_l1 + gdp_growth_w + mil_exp + decade +
                coup_years + cy2 + cy3,
              data = panel, family = binomial)
coeftest(m_temp, vcov = vcovCL(m_temp, cluster = ~ iso3c, data = panel))

# (ix) population-weighted SPEI-12: grid cells weighted by 2020 population
#      (GPW v4 density via geodata, no login). Fixed weights avoid an
#      endogenous population response to drought and address the concern
#      that unweighted country means dilute exposure in large states.
path_pw <- "Data/processed/spei12_popweighted.rds"
if (file.exists(path_pw)) {
  spei12w <- readRDS(path_pw)
} else {
  # rebuild Africa polygons (Part A objects are not kept when RUN_DATA_PREP = FALSE)
  africa <- ne_countries(scale = "medium", continent = "Africa",
                         returnclass = "sf")
  africa$iso3c <- countrycode(africa$admin, "country.name", "iso3c")
  africa <- africa[!is.na(africa$iso3c), ]
  
  pop <- geodata::population(2020, res = 10, path = "Data/raw")
  r12 <- rast(path_spei12)
  nc  <- nc_open(path_spei12)
  tv  <- ncvar_get(nc, "time")
  tu  <- ncatt_get(nc, "time", "units")$value
  nc_close(nc)
  origin <- as.Date(sub(".*since\\s*", "", tu))
  dts <- if (grepl("month", tu)) origin %m+% months(as.integer(tv))
  else origin + as.numeric(tv)
  dec <- which(month(dts) == 12 & year(dts) >= year_start & year(dts) <= year_end)
  stopifnot(length(dec) == year_end - year_start + 1)
  r12  <- r12[[dec]]
  # 10' -> 30' density with na.rm (no NA bleed at coasts), then align to
  # the SPEI grid; default_weight = 0 treats remaining no-data cells as
  # weight zero instead of dropping the whole country
  pop05 <- terra::aggregate(pop, fact = 3, fun = "mean", na.rm = TRUE)
  popc  <- terra::resample(pop05, r12[[1]], method = "near") * cellSize(r12[[1]], unit = "km")
  
  wtab <- exact_extract(r12, africa, "weighted_mean", weights = popc,
                        default_weight = 0, progress = FALSE)
  wtab <- data.frame(iso3c = africa$iso3c, wtab)
  names(wtab)[-1] <- paste0("y", year_start:year_end)
  spei12w <- wtab %>%
    pivot_longer(-iso3c, names_to = "year", values_to = "spei12_dec_pw") %>%
    mutate(year = as.integer(sub("y", "", year)))
  saveRDS(spei12w, path_pw)
}

panel <- panel %>%
  select(-any_of(c("spei12_dec_pw", "spei12_pw_l1"))) %>%   # idempotent re-runs
  left_join(spei12w, by = c("iso3c", "year")) %>%
  group_by(iso3c) %>%
  mutate(spei12_pw_l1 = dplyr::lag(spei12_dec_pw, 1)) %>%
  ungroup()

m_pw <- glm(coup_attempt ~ spei12_pw_l1 + gdp_growth_w + mil_exp + decade,
            data = panel, family = binomial)
coeftest(m_pw, vcov = vcovCL(m_pw, cluster = ~ iso3c, data = panel))

# (x) two-way clustered SEs (country and year): droughts are spatially
#     correlated across neighbours and coups cluster in time, so
#     country-only clustering may understate the uncertainty
coeftest(m_h1_l1, vcov = vcovCL(m_h1_l1, cluster = ~ iso3c + year, data = panel))

# (xi) excluding the 2020s coup wave: years 1990-2019 only
panel_pre <- panel %>% filter(year <= 2019) %>% mutate(decade = droplevels(decade))
m_pre20 <- glm(coup_attempt ~ spei12_dec_l1 + gdp_growth_w + mil_exp + decade,
               data = panel_pre, family = binomial)
coeftest(m_pre20, vcov = vcovCL(m_pre20, cluster = ~ iso3c, data = panel_pre))
cat("coup attempts 1990-2019:", sum(panel_pre$coup_attempt), "of 77\n")

cc <- panel %>% filter(!is.na(spei12_dec_l1), !is.na(gdp_growth_w), !is.na(mil_exp))
cat("main estimation sample after state-membership filter: N =", nrow(cc), "\n")
cat("coup attempts in the estimation sample:",
    sum(cc$coup_attempt[cc$year <= 2019]), "in 1990-2019,",
    sum(cc$coup_attempt[cc$year >= 2020]), "in 2020-2023\n")

# (xii) linear year trend instead of decade dummies: the sample-wide drying
#       trend is not absorbed within decades (year centred at 2006)
m_trend <- glm(coup_attempt ~ spei12_dec_l1 + gdp_growth_w + mil_exp + I(year - 2006),
               data = panel, family = binomial)
coeftest(m_trend, vcov = vcovCL(m_trend, cluster = ~ iso3c, data = panel))

# (xiii) excluding the six countries of the post-2020 Sahel coup wave
wave <- c("MLI", "GIN", "BFA", "NER", "TCD", "SDN")
panel_nw <- panel %>% filter(!iso3c %in% wave)
m_nowave <- glm(coup_attempt ~ spei12_dec_l1 + gdp_growth_w + mil_exp + decade,
                data = panel_nw, family = binomial)
coeftest(m_nowave, vcov = vcovCL(m_nowave, cluster = ~ iso3c, data = panel_nw))
cat("coup attempts without the wave countries:", sum(panel_nw$coup_attempt), "of 77\n")

# (xiv) placebo: future drought (t+1) must not predict current coups;
#       estimated jointly with t-1, because drought persistence would let
#       an isolated lead pick up an echo of past dryness
panel <- panel %>%
  arrange(iso3c, year) %>%
  group_by(iso3c) %>%
  mutate(spei12_dec_f1 = dplyr::lead(spei12_dec, 1)) %>%
  ungroup()
m_placebo <- glm(coup_attempt ~ spei12_dec_l1 + spei12_dec_f1 +
                   gdp_growth_w + mil_exp + decade,
                 data = panel, family = binomial)
coeftest(m_placebo, vcov = vcovCL(m_placebo, cluster = ~ iso3c, data = panel))

# (xv) functional form: dry and wet components entered separately, so the
#      estimate is identified from the dry side rather than from
#      stabilising wet years (dry_l1 = drought magnitude, expected positive)
panel <- panel %>%
  mutate(dry_l1 = pmax(0, -spei12_dec_l1),
         wet_l1 = pmax(0,  spei12_dec_l1))
m_split <- glm(coup_attempt ~ dry_l1 + wet_l1 + gdp_growth_w + mil_exp + decade,
               data = panel, family = binomial)
coeftest(m_split, vcov = vcovCL(m_split, cluster = ~ iso3c, data = panel))

# (xvi) year fixed effects via Firth logit: with only 77 events, a full set
#       of year dummies gives roughly two events per parameter, and years
#       without a single coup attempt are perfectly predicted (quasi-
#       complete separation); the penalised likelihood handles both
print(table(factor(panel$year[panel$coup_attempt == 1],
                   levels = year_start:year_end)))  # zero-coup years visible
m_yearfe <- logistf(coup_attempt ~ spei12_dec_l1 + gdp_growth_w + mil_exp +
                      factor(year), data = mdat)
cat("Year-FE Firth: SPEI-12 t-1 =", round(coef(m_yearfe)["spei12_dec_l1"], 3),
    "| p =", signif(m_yearfe$prob["spei12_dec_l1"], 3), "\n")

# (xvii) military control and sample composition: listwise deletion costs
#        31 of the 77 coup years, almost all of them through military
#        expenditure as a share of government expenditure. The check
#        replaces that series with military expenditure as a share of GDP
#        (SIPRI via WDI, wider coverage), drops the military control
#        altogether, and re-estimates the reduced models on the two
#        restricted samples, so that the control and the sample can be
#        told apart. Console output only, no publication table is changed

# which series is the binding constraint among the coup years?
vars_main <- c("spei12_dec_l1", "gdp_growth_w", "mil_exp")
cat("missing values among the", sum(panel$coup_attempt), "coup years:\n")
print(colSums(is.na(panel[panel$coup_attempt == 1, vars_main])))

# alternative military control, cached like the other downloads
path_milgdp <- "Data/raw/wdi_mil_gdp.rds"
if (file.exists(path_milgdp)) {
  mil_gdp <- readRDS(path_milgdp)
} else {
  mil_gdp <- WDI(country = "all",
                 indicator = c(mil_gdp = "MS.MIL.XPND.GD.ZS"),
                 start = year_start, end = year_end) %>%
    filter(!is.na(iso3c), iso3c != "") %>%
    distinct(iso3c, year, .keep_all = TRUE) %>%
    select(iso3c, year, mil_gdp)
  dir.create(dirname(path_milgdp), recursive = TRUE, showWarnings = FALSE)
  saveRDS(mil_gdp, path_milgdp)
}

panel <- panel %>%
  select(-any_of("mil_gdp")) %>%              # idempotent re-runs
  left_join(mil_gdp, by = c("iso3c", "year"))

cat("coverage in the panel: mil_exp",
    round(100 * mean(!is.na(panel$mil_exp)), 1), "%, mil_gdp",
    round(100 * mean(!is.na(panel$mil_gdp)), 1), "%\n")

m_milgdp <- glm(coup_attempt ~ spei12_dec_l1 + gdp_growth_w + mil_gdp + decade,
                data = panel, family = binomial)
coeftest(m_milgdp, vcov = vcovCL(m_milgdp, cluster = ~ iso3c, data = panel))

m_nomil <- glm(coup_attempt ~ spei12_dec_l1 + gdp_growth_w + decade,
               data = panel, family = binomial)
coeftest(m_nomil, vcov = vcovCL(m_nomil, cluster = ~ iso3c, data = panel))

# same reduced model on the two restricted samples
smp_exp <- complete.cases(panel[, c("coup_attempt", "spei12_dec_l1",
                                    "gdp_growth_w", "mil_exp", "decade")])
smp_gdp <- complete.cases(panel[, c("coup_attempt", "spei12_dec_l1",
                                    "gdp_growth_w", "mil_gdp", "decade")])

m_nomil_sexp <- glm(coup_attempt ~ spei12_dec_l1 + gdp_growth_w + decade,
                    data = panel[smp_exp, ], family = binomial)
m_nomil_sgdp <- glm(coup_attempt ~ spei12_dec_l1 + gdp_growth_w + decade,
                    data = panel[smp_gdp, ], family = binomial)

# compact overview of the SPEI-12 t-1 coefficient across the five fits
spei_row <- function(m, lab, d) {
  ct <- coeftest(m, vcov = vcovCL(m, cluster = ~ iso3c, data = d))
  data.frame(Specification = lab,
             b  = round(ct["spei12_dec_l1", 1], 3),
             SE = round(ct["spei12_dec_l1", 2], 3),
             p  = round(ct["spei12_dec_l1", 4], 4),
             N  = nobs(m),
             `Coup years` = sum(model.frame(m)$coup_attempt),
             check.names = FALSE)
}

military_tab <- rbind(
  spei_row(m_h1_l1,      "Mil. exp. (% gov. exp.)",       panel),
  spei_row(m_nomil_sexp, "No mil. control, same sample",  panel[smp_exp, ]),
  spei_row(m_milgdp,     "Mil. exp. (% of GDP)",          panel),
  spei_row(m_nomil_sgdp, "No mil. control, same sample",  panel[smp_gdp, ]),
  spei_row(m_nomil,      "No mil. control, full sample",  panel))
print(military_tab, row.names = FALSE)

# ------------------------------------------------------------------ #
# 15. Time-scale profile: which temporal signature of drought matters?
#     SPEI-3 (acute) / SPEI-6 (season) / SPEI-12 (year) / SPEI-24 (2 yrs)
# ------------------------------------------------------------------ #
m_scale03 <- glm(coup_attempt ~ spei3_oct_l1 + gdp_growth_w + mil_exp + decade,
                 data = panel, family = binomial)
coeftest(m_scale03, vcov = vcovCL(m_scale03, cluster = ~ iso3c, data = panel))

m_scale06 <- glm(coup_attempt ~ spei6_oct_l1 + gdp_growth_w + mil_exp + decade,
                 data = panel, family = binomial)
coeftest(m_scale06, vcov = vcovCL(m_scale06, cluster = ~ iso3c, data = panel))

# SPEI-12 (Dec), t-1 = m_h1_l1 (main model, estimated above)

m_scale24 <- glm(coup_attempt ~ spei24_dec_l1 + gdp_growth_w + mil_exp + decade,
                 data = panel, family = binomial)
coeftest(m_scale24, vcov = vcovCL(m_scale24, cluster = ~ iso3c, data = panel))

# inference on the time-scale contrast: comparing significance levels is not
# a test of the coefficient difference (Gelman & Stern 2006); clustered
# pairs bootstrap over countries for the SPEI-12 vs. SPEI-3 difference
set.seed(1990)
B <- 1000
ctry <- unique(panel$iso3c)
boot_diff <- rep(NA_real_, B)
for (b in seq_len(B)) {
  draw <- sample(ctry, length(ctry), replace = TRUE)
  bd   <- dplyr::bind_rows(lapply(draw, function(cc) panel[panel$iso3c == cc, ]))
  b12  <- tryCatch(coef(glm(coup_attempt ~ spei12_dec_l1 + gdp_growth_w +
                              mil_exp + decade, data = bd, family = binomial))["spei12_dec_l1"],
                   error = function(e) NA_real_)
  b03  <- tryCatch(coef(glm(coup_attempt ~ spei3_oct_l1 + gdp_growth_w +
                              mil_exp + decade, data = bd, family = binomial))["spei3_oct_l1"],
                   error = function(e) NA_real_)
  boot_diff[b] <- b12 - b03
}
boot_diff <- boot_diff[is.finite(boot_diff)]
cat("bootstrap draws used:", length(boot_diff), "of", B, "\n")
cat("difference SPEI-12 minus SPEI-3 (point estimate):",
    round(unname(coef(m_h1_l1)["spei12_dec_l1"] - coef(m_scale03)["spei3_oct_l1"]), 3), "\n")
cat("bootstrap 95% CI:", round(unname(quantile(boot_diff, c(0.025, 0.975))), 3), "\n")
cat("share of draws with difference >= 0:", round(mean(boot_diff >= 0), 3), "\n")

# ------------------------------------------------------------------ #
# 16. Publication tables (modelsummary -> Word)
# ------------------------------------------------------------------ #

# GOF rows: "N" and "Log-likelihood" (logit tables), "N" and "R2" (plm)
gof_logit <- list(
  list(raw = "nobs",   clean = "N",              fmt = 0),
  list(raw = "logLik", clean = "Log-likelihood", fmt = 3))
gof_plm <- list(
  list(raw = "nobs",      clean = "N",  fmt = 0),
  list(raw = "r.squared", clean = "R2", fmt = 3))

# Table 1: H1, H2, H4 (drought & coups, SPEI-12)
cm_coup <- c("spei12_dec"                = "SPEI-12 (Dec), t",
             "spei12_dec_l1"             = "SPEI-12 (Dec), t-1",
             "spei12_dec_l2"             = "SPEI-12 (Dec), t-2",
             "drought12_cat_l1moderate"  = "Moderate drought, t-1",
             "drought12_cat_l1severe"    = "Severe drought, t-1",
             "drought12_cat_l1extreme"   = "Extreme drought, t-1",
             "drought12_onset_l1"        = "Drought onset, t-1",
             "gdp_growth_w"              = "GDP p.c. growth (wins.)",
             "mil_exp"                   = "Military expenditure")

modelsummary(
  list("H1: SPEI t" = m_h1, "H1: SPEI t-1" = m_h1_l1,
       "H1: No controls" = m_h1_nc, "H2: Intensity" = m_h2,
       "H4: Onset" = m_h4, "H4: Lags" = m_h4_l2),
  vcov = ~ iso3c, coef_map = cm_coup, stars = c("+" = .1, "*" = .05, "**" = .01, "***" = .001),
  gof_map = gof_logit,
  notes = "Logit. Decade dummies included. SE clustered by country.",
  output = "Output/table1_coups.docx")

# Table 2: robustness (incl. pre-existing downturn and pop-weighted SPEI)
cm_rob <- c("spei12_dec_l1"          = "SPEI-12 (Dec), t-1",
            "spei6_mean_l1"          = "SPEI-6 (ann. mean), t-1",
            "spei12_pw_l1"           = "SPEI-12 pop.-weighted, t-1",
            "drought12_dur_l1"       = "Drought duration (yrs), t-1",
            "agr_dep"                = "Agri. dependence",
            "spei12_dec_l1:agr_dep"  = "SPEI-12 t-1 x Agri. dep.",
            "gdp_growth_w"           = "GDP p.c. growth (wins.)",
            "gdp_growth_w_l1"        = "GDP p.c. growth, t-1 (wins.)",
            "mil_exp"                = "Military expenditure")

modelsummary(
  list("Firth logit" = m_firth, "Interaction" = m_int,
       "Cond. logit (FE)" = m_clogit, "Sub-Saharan" = m_ssa,
       "Pre-existing downturn" = m_pre, "Alt. measure" = m_alt,
       "Pop. weighted" = m_pw, "Duration" = m_dur),
  vcov = list(NULL, ~ iso3c, NULL, ~ iso3c, ~ iso3c, ~ iso3c, ~ iso3c, ~ iso3c),
  coef_map = cm_rob, stars = c("+" = .1, "*" = .05, "**" = .01, "***" = .001),
  gof_map = gof_logit,
  notes = paste0("Col. 1: penalised ML (Firth) with profile penalised-likelihood ",
                 "inference, SE not clustered. Col. 3: conditional logit (Efron ",
                 "approximation), country strata, model-based SE. Countries without ",
                 "any coup attempt do not contribute to the conditional likelihood ",
                 "(informative sample: ", n_eff_clogit, " country-years in ",
                 n_countries_inf, " countries). All other columns: SE clustered ",
                 "by country. Decade dummies included. Col. 7 uses the ",
                 "population-weighted SPEI-12 (grid cells weighted by 2020 population)."),
  output = "Output/table2_robustness.docx")

# Table 3: time-scale profile (SPEI-3 / -6 / -12 / -24, each lagged t-1)
cm_scale <- c("spei3_oct_l1"  = "SPEI-3 (Oct), t-1",
              "spei6_oct_l1"  = "SPEI-6 (Oct), t-1",
              "spei12_dec_l1" = "SPEI-12 (Dec), t-1",
              "spei24_dec_l1" = "SPEI-24 (Dec), t-1",
              "gdp_growth_w"  = "GDP p.c. growth (wins.)",
              "mil_exp"       = "Military expenditure")

modelsummary(
  list("SPEI-3" = m_scale03, "SPEI-6" = m_scale06,
       "SPEI-12" = m_h1_l1, "SPEI-24" = m_scale24),
  vcov = ~ iso3c, coef_map = cm_scale, stars = c("+" = .1, "*" = .05, "**" = .01, "***" = .001),
  gof_map = gof_logit,
  notes = "Logit. Decade dummies included. SE clustered by country.",
  output = "Output/table3_timescales.docx")

# Table 4: plausibility check, economic channel (H3)
cm_h3 <- c("spei12_dec"         = "SPEI-12, t",
           "lag(spei12_dec, 1)" = "SPEI-12, t-1")

modelsummary(
  list("GDP growth" = m_h3, "GDP growth (with lag)" = m_h3_both,
       "GDP growth (wins.)" = m_h3_w,
       "GDP growth excl. oil (wins.)" = m_h3_ag,
       "Agri. growth (wins.)" = m_h3_agr),
  vcov = function(x) vcovHC(x, cluster = "group", type = "HC1"),
  coef_map = cm_h3, stars = c("+" = .1, "*" = .05, "**" = .01, "***" = .001),
  gof_map = gof_plm,
  notes = paste0("Country and year fixed effects. SE clustered by country. ",
                 "Col. 4 excludes the nine major oil producers (Algeria, ",
                 "Angola, Chad, the Republic of Congo, Equatorial Guinea, ",
                 "Gabon, Libya, Nigeria and South Sudan)."),
  output = "Output/table4_h3.docx")

# ------------------------------------------------------------------ #
# 17. Appendix tables
# ------------------------------------------------------------------ #

# A1: extended lag structure (t to t-3)
modelsummary(
  list("Lags to t-3" = m_h4_lags),
  vcov = ~ iso3c,
  coef_map = c("spei12_dec"    = "SPEI-12 (Dec), t",
               "spei12_dec_l1" = "SPEI-12 (Dec), t-1",
               "spei12_dec_l2" = "SPEI-12 (Dec), t-2",
               "spei12_dec_l3" = "SPEI-12 (Dec), t-3",
               "gdp_growth_w"  = "GDP p.c. growth (wins.)",
               "mil_exp"       = "Military expenditure"),
  stars = c("+" = .1, "*" = .05, "**" = .01, "***" = .001), gof_map = gof_logit,
  notes = "Logit. Decade dummies included. SE clustered by country.",
  output = "Output/tableA1_lags.docx")

# A2: temporal dependence (coup trap) robustness, full results
cm_temp <- c("spei12_dec_l1" = "SPEI-12 (Dec), t-1",
             "coup_years"    = "Years since last coup",
             "cy2"           = "Years since last coup^2 / 100",
             "cy3"           = "Years since last coup^3 / 1000",
             "gdp_growth_w"  = "GDP p.c. growth (wins.)",
             "mil_exp"       = "Military expenditure")

modelsummary(
  list("Temporal dependence" = m_temp),
  vcov = ~ iso3c, coef_map = cm_temp, stars = c("+" = .1, "*" = .05, "**" = .01, "***" = .001),
  gof_map = gof_logit,
  notes = "Logit. Decade dummies included. SE clustered by country. Peace-years polynomial follows Carter & Signorino (2010).",
  output = "Output/tableA2_temporal.docx")

# A5: identification and functional-form checks
cm_id <- c("spei12_dec_l1"  = "SPEI-12 (Dec), t-1",
           "spei12_dec_f1"  = "SPEI-12 (Dec), t+1 (placebo)",
           "dry_l1"         = "Dry component, t-1",
           "wet_l1"         = "Wet component, t-1",
           "I(year - 2006)" = "Linear year trend",
           "gdp_growth_w"   = "GDP p.c. growth (wins.)",
           "mil_exp"        = "Military expenditure")

modelsummary(
  list("Two-way SE" = m_h1_l1, "1990-2019" = m_pre20, "Linear trend" = m_trend,
       "Excl. wave" = m_nowave, "Placebo lead" = m_placebo, "Dry/wet split" = m_split),
  vcov = list(~ iso3c + year, ~ iso3c, ~ iso3c, ~ iso3c, ~ iso3c, ~ iso3c),
  coef_map = cm_id, stars = c("+" = .1, "*" = .05, "**" = .01, "***" = .001),
  gof_map = gof_logit,
  notes = paste0("Logit. Decade dummies included except Col. 3 (linear year ",
                 "trend). SE clustered by country. Col. 1 clustered by country ",
                 "and year. Dry component = |min(SPEI, 0)|, wet component = ",
                 "max(SPEI, 0), so that a positive dry and a negative wet ",
                 "coefficient both indicate higher coup risk under drier ",
                 "conditions."),
  output = "Output/tableA5_identification.docx")

# A4: country list as Word table (excl. MUS/SYC, not covered by SPEIbase,
#     so the list matches the 52-country analysis universe of Chapter 3.2)
country_tab <- panel %>%
  filter(!iso3c %in% c("MUS", "SYC")) %>%
  group_by(Country = country) %>%
  summarise(`Coup years` = sum(coup_attempt)) %>%
  arrange(desc(`Coup years`))
save_as_docx(flextable(country_tab), path = "Output/tableA4_countries.docx")


# ------------------------------------------------------------------ #
# 18. Session info (reproducibility)
# ------------------------------------------------------------------ #
sessionInfo()

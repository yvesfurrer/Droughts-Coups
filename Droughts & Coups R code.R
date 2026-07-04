###############################################################################
# Masterseminarpaper: Droughts & Coups in Africa
# Author: Yves Furrer
# Main IV: SPEI-12 (December) = full calendar-year water balance
# Part A: Data preparation (coups, SPEI 3/6/12/24, WDI)
# Part B: Descriptives, Analysis (H1-H4), Robustness, Time scales, Tables
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
          "flextable", "pandoc", "survival")
new <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
if (length(new)) install.packages(new)

library(WDI); library(terra); library(sf); library(survival)
library(rnaturalearth); library(rnaturalearthdata)
library(ncdf4); library(lubridate); library(dplyr); library(tidyr)
library(countrycode); library(plm); library(lmtest); library(sandwich)
library(logistf); library(modelsummary); library(pandoc)

# CAUTION: plm masks dplyr::lag. Inside mutate() we therefore ALWAYS
# write dplyr::lag() explicitly. plm's own lag() is still used inside
# plm() formulas, where it is the correct (panel-aware) one.

# ------------------------------------------------------------------ #
# 1. Paths & parameters
# ------------------------------------------------------------------ #
path_coup   <- "C:/Users/yfurr/Documents/1 - UNILU/Masterseminararbeit 1/R/Data/raw/powell_thyne_ccode_year.csv"
path_spei03 <- "C:/Users/yfurr/Documents/1 - UNILU/Masterseminararbeit 1/R/Data/raw/spei03.nc"
path_spei06 <- "C:/Users/yfurr/Documents/1 - UNILU/Masterseminararbeit 1/R/Data/raw/spei06.nc"
path_spei12 <- "C:/Users/yfurr/Documents/1 - UNILU/Masterseminararbeit 1/R/Data/raw/spei12.nc"
path_spei24 <- "C:/Users/yfurr/Documents/1 - UNILU/Masterseminararbeit 1/R/Data/raw/spei24.nc"
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

dir.create("Output", showWarnings = FALSE)
pdat <- pdata.frame(panel, index = c("iso3c", "year"))

# ------------------------------------------------------------------ #
# 8. Descriptive statistics (Appendix tables A3, A4)
# ------------------------------------------------------------------ #

# A3: country list with number of coup attempts
panel %>% group_by(country) %>%
  summarise(coups = sum(coup_attempt)) %>%
  arrange(desc(coups)) %>% print(n = 60)

# A4: summary statistics of the main variables
datasummary_skim(panel %>% select(coup_attempt, spei12_dec, gdp_growth_w,
                                  mil_exp, agr_growth_w),
                 output = "Output/tableA2_summary.docx")

# ------------------------------------------------------------------ #
# 9. H3: Drought and economic conditions (plausibility check)
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
# 10. H1: Drought and coup attempts (main effect, SPEI-12)
# ------------------------------------------------------------------ #

# (i) contemporaneous (timing within year ambiguous; t-1 is main spec)
m_h1 <- glm(coup_attempt ~ spei12_dec + gdp_growth_w + mil_exp + decade,
            data = panel, family = binomial)
coeftest(m_h1, vcov = vcovCL(m_h1, cluster = ~ iso3c, data = panel))

# (ii) lagged (t-1)  [MAIN MODEL]
m_h1_l1 <- glm(coup_attempt ~ spei12_dec_l1 + gdp_growth_w + mil_exp + decade,
               data = panel, family = binomial)
coeftest(m_h1_l1, vcov = vcovCL(m_h1_l1, cluster = ~ iso3c, data = panel))

# ------------------------------------------------------------------ #
# 11. H2: Drought intensity (SPEI-12 categories, lagged)
# ------------------------------------------------------------------ #
m_h2 <- glm(coup_attempt ~ drought12_cat_l1 + gdp_growth_w + mil_exp + decade,
            data = panel, family = binomial)
coeftest(m_h2, vcov = vcovCL(m_h2, cluster = ~ iso3c, data = panel))

# ------------------------------------------------------------------ #
# 12. H4: Timing of the effect (SPEI-12)
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
# 13. Robustness checks (all on the main model, SPEI-12 t-1)
# ------------------------------------------------------------------ #

# (i) Firth logit (rare events bias correction)
mdat <- panel %>%
  select(coup_attempt, spei12_dec_l1, gdp_growth_w, mil_exp, decade, iso3c) %>%
  na.omit()
m_firth <- logistf(coup_attempt ~ spei12_dec_l1 + gdp_growth_w + mil_exp + decade,
                   data = mdat)
summary(m_firth)

# (ii) interaction: drought x agricultural dependence
m_int <- glm(coup_attempt ~ spei12_dec_l1 * agr_dep + gdp_growth_w + mil_exp + decade,
             data = panel, family = binomial)
coeftest(m_int, vcov = vcovCL(m_int, cluster = ~ iso3c, data = panel))

# (iii) conditional logit (country fixed effects)
m_clogit <- clogit(coup_attempt ~ spei12_dec_l1 + gdp_growth_w + mil_exp +
                     decade + strata(iso3c),
                   data = panel, method = "efron")
summary(m_clogit)

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

# ------------------------------------------------------------------ #
# 14. Time-scale profile: which temporal signature of drought matters?
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

# ------------------------------------------------------------------ #
# 15. Publication tables (modelsummary -> Word)
# ------------------------------------------------------------------ #

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
  list("H1: t" = m_h1, "H1: t-1" = m_h1_l1, "H2: intensity" = m_h2,
       "H4: onset" = m_h4, "H4: lags" = m_h4_l2),
  vcov = ~ iso3c, coef_map = cm_coup, stars = TRUE,
  gof_omit = "IC|RMSE|Std.Errors",
  notes = "Logit. Decade dummies included. SE clustered by country.",
  output = "Output/table1_coups.docx")

# Table 2: robustness (incl. pre-existing downturn)
cm_rob <- c("spei12_dec_l1"          = "SPEI-12 (Dec), t-1",
            "spei6_mean_l1"          = "SPEI-6 (ann. mean), t-1",
            "drought12_dur_l1"       = "Drought duration (yrs), t-1",
            "agr_dep"                = "Agri. dependence",
            "spei12_dec_l1:agr_dep"  = "SPEI-12 t-1 x Agri. dep.",
            "gdp_growth_w"           = "GDP p.c. growth (wins.)",
            "gdp_growth_w_l1"        = "GDP p.c. growth, t-1 (wins.)",
            "mil_exp"                = "Military expenditure")

modelsummary(
  list("Firth logit" = m_firth, "Interaction" = m_int,
       "Cond. logit (FE)" = m_clogit, "SSA only" = m_ssa,
       "Pre-exist. downturn" = m_pre, "Alt. measure" = m_alt,
       "Duration" = m_dur),
  vcov = list(NULL, ~ iso3c, NULL, ~ iso3c, ~ iso3c, ~ iso3c, ~ iso3c),
  coef_map = cm_rob, stars = TRUE,
  gof_omit = "IC|RMSE|Std.Errors",
  notes = "Col. 1: penalized ML (Firth). Col. 3: conditional logit, country strata. Decade dummies included.",
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
  vcov = ~ iso3c, coef_map = cm_scale, stars = TRUE,
  gof_omit = "IC|RMSE|Std.Errors",
  notes = "Logit. Decade dummies included. SE clustered by country.",
  output = "Output/table3_timescales.docx")

# Table 4: plausibility check, economic channel (H3)
cm_h3 <- c("spei12_dec"         = "SPEI-12 (Dec), t",
           "lag(spei12_dec, 1)" = "SPEI-12 (Dec), t-1")

modelsummary(
  list("(1) GDP growth" = m_h3, "(2) GDP growth" = m_h3_both,
       "(3) GDP growth (wins.)" = m_h3_w, "(4) Agri. growth (wins.)" = m_h3_agr),
  vcov = function(x) vcovHC(x, cluster = "group", type = "HC1"),
  coef_map = cm_h3, stars = TRUE,
  gof_omit = "IC|RMSE|Adj|Within|FE|Std.Errors",
  notes = "Country and year fixed effects. SE clustered by country.",
  output = "Output/table4_h3.docx")

# ------------------------------------------------------------------ #
# 15. Appendix tables
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
  stars = TRUE, gof_omit = "IC|RMSE|Std.Errors",
  notes = "Logit. Decade dummies included. SE clustered by country.",
  output = "Output/tableA1_lags.docx")

# A3: country list as Word table
library(flextable)
country_tab <- panel %>% group_by(country) %>%
  summarise(`Coup attempts` = sum(coup_attempt)) %>%
  arrange(desc(`Coup attempts`))
save_as_docx(flextable(country_tab), path = "Output/tableA3_countries.docx")

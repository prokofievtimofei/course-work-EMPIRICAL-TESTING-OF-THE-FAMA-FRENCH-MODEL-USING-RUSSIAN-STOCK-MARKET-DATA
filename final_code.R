# ================================================================
#  Fama-French 5-Factor Model — Russian market
#  R_i - R_f = α + β(R_m-R_f) + s*SMB + h*HML + r*RMW + c*CMA
# ================================================================

# ================================================================
#  Input data format:

#    monthly_prices.csv    — wide: rows = dates, columns = tickers

#    fundament.csv       — long: date | symbol | totalAssets |
#                          totalStockholdersEquity | revenue | costOfRevenue |
#                          operatingIncome | sellingGeneralAndAdministrativeExpenses |
#                           interestExpense | marketCapitalization

#    mkt_factor.csv      — long: TRADEDATE,MKT

#    rf_monthly.csv      — long | RF_month (key rate of the Central Bank, % per annum)
# ================================================================

library(tidyverse)
library(lubridate)
library(broom)

# ================================================================
# PATHS AND PARAMETERS
# ================================================================


# Paths to data files
PATH_PRICES    <- "/Users/prokofiev/Yandex.Disk.localized/Finance/Fama-French model/final/monthly_prices.csv"
PATH_FUNDAMENT <- "/Users/prokofiev/Yandex.Disk.localized/Finance/Fama-French model/final/fundament.csv"
PATH_MKT       <- "/Users/prokofiev/Yandex.Disk.localized/Finance/Fama-French model/final/mkt_factor.csv"
PATH_RF        <- "/Users/prokofiev/Yandex.Disk.localized/Finance/Fama-French model/final/rf_monthly.csv"


# Thresholds for 2x3 sorting (Fama-French 2015)
LOW_PCTILE  <- 0.30
HIGH_PCTILE <- 0.70


# ================================================================
# 1. DATA LOADING
# ================================================================


# 1.1 Loading prices, converting from wide to long
# VTBR: 5000:1 reverse split in July 2024
prices <- read_csv(PATH_PRICES) %>%
  rename(date = TRADEDATE) %>%
  mutate(date = as.Date(date)) %>%
  pivot_longer(-date, names_to = "ticker", values_to = "price") %>%
  filter(!is.na(price), price > 0) %>%
  mutate(
    price = if_else(ticker == "VTBR" & date < as.Date("2024-07-01"),
                    price * 5000, price),
    ym = floor_date(date, "month"))


## 1.2 Monthly returns
monthly_ret <- prices %>%
  group_by(ticker, ym) %>%
  slice_max(date, n = 1) %>%
  ungroup() %>%
  arrange(ticker, ym) %>%
  group_by(ticker) %>%
  mutate(ret = price / dplyr::lag(price) - 1) %>%
  ungroup() %>%
  filter(!is.na(ret)) %>%
  select(date = ym, ticker, ret, price)

# Winsorization of monthly returns at the level of 1% and 99% 
monthly_ret <- monthly_ret %>%
  group_by(date) %>%
  mutate(ret = pmin(pmax(ret, quantile(ret, 0.01, na.rm=TRUE)),
                    quantile(ret, 0.99, na.rm=TRUE))) %>%
  ungroup()


## 1.3 Fundamental data
fundament <- read_csv(PATH_FUNDAMENT) %>%
  mutate(date = as.Date(date), year = year(date)) %>%
  rename(ticker = symbol) %>%
  filter(totalStockholdersEquity > 0) %>%
  mutate(
    # OP = (Revenue - COGS - SG&A - Interest) / Book Equity  (following Fama-French 2015)
    op_ff = (revenue - costOfRevenue -
               sellingGeneralAndAdministrativeExpenses -
               interestExpense) / totalStockholdersEquity
  ) %>%
  arrange(ticker, year) %>%
  group_by(ticker) %>%
  mutate(ta_lag = dplyr::lag(totalAssets),
         inv    = (totalAssets - ta_lag) / ta_lag) %>%
  ungroup() %>%
  filter(is.finite(op_ff), is.finite(inv))


## 1.4 Market Factor (IMOEX)
mkt <- read_csv(PATH_MKT) %>%
  rename(date = TRADEDATE) %>%
  mutate(date = floor_date(as.Date(date), "month")) %>%
  group_by(date) %>% slice_tail(n = 1) %>% ungroup()


## 1.5 Risk-free rate (key rate of the Central Bank, monthly)
rf <- read_csv(PATH_RF) %>%
  rename(date = Date, rf_raw = RF_month) %>%
  mutate(date = floor_date(as.Date(date), "month"),
         rf = rf_raw) %>%
  select(date, rf) %>%
  group_by(date) %>% slice_tail(n = 1) %>% ungroup()


# ================================================================
# 2. Market capitalization at the end of June (for sorting)
# Number of shares = MC December(t-1) / Price December(t-1)
# MC June(t) = Number of shares × Price June(t)
# ================================================================


price_dec <- monthly_ret %>%
  filter(month(date) == 12) %>%
  mutate(year = year(date)) %>%
  select(ticker, year, price_dec = price)

price_june <- monthly_ret %>%
  filter(month(date) == 6) %>%
  mutate(sort_year = year(date)) %>%
  select(ticker, sort_year, price_june = price)

mcap_june <- fundament %>%
  select(ticker, year, mktcap_dec = marketCapitalization) %>%
  inner_join(price_dec, by = c("ticker", "year")) %>%
  mutate(shares    = mktcap_dec / price_dec,
         sort_year = year + 1) %>%
  inner_join(price_june, by = c("ticker", "sort_year")) %>%
  mutate(mktcap_june = shares * price_june) %>%
  filter(mktcap_june > 0, !is.na(mktcap_june)) %>%
  select(ticker, sort_year, mktcap_june)


# ================================================================
# 3. CHARACTERISTICS FOR SORTING (B/M, OP, INV)
# ================================================================


char_data <- mcap_june %>%
  left_join(
    fundament %>%
      mutate(sort_year = year + 1) %>%
      select(ticker, sort_year,
             book_equity = totalStockholdersEquity, op_ff, inv),
    by = c("ticker", "sort_year")
  ) %>%
  mutate(bm = book_equity / mktcap_june,
         op = op_ff) %>%
  select(-op_ff) %>%
  filter(!is.na(bm), !is.na(op), !is.na(inv),
         bm > 0, is.finite(op), is.finite(inv))


# ================================================================
# 4. ASSIGNING 2x3 PORTFOLIOS (Size x Characteristic)
# ================================================================


assign_2x3 <- function(df, char_col) {
  df %>%
    group_by(sort_year) %>%
    mutate(
      size_grp = if_else(mktcap_june <= median(mktcap_june), "S", "B"),
      q_lo     = quantile(.data[[char_col]], LOW_PCTILE,  na.rm = TRUE),
      q_hi     = quantile(.data[[char_col]], HIGH_PCTILE, na.rm = TRUE),
      char_grp = case_when(
        .data[[char_col]] <= q_lo ~ "L",
        .data[[char_col]] <= q_hi ~ "N",
        TRUE                      ~ "H"
      ),
      portfolio = paste0(size_grp, char_grp)
    ) %>%
    ungroup() %>%
    select(ticker, sort_year, mktcap_june, portfolio)}

ports_bm  <- assign_2x3(char_data, "bm")
ports_op  <- assign_2x3(char_data, "op")
ports_inv <- assign_2x3(char_data, "inv")


# ================================================================
# 5. VALUE-WEIGHTED PORTFOLIO RETURN
# ================================================================


vw_ret <- function(assignments) {
  map_dfr(unique(assignments$sort_year), function(yr) {
    monthly_ret %>%
      filter(date >= make_date(yr, 7, 1),
             date <= make_date(yr + 1, 6, 30)) %>%
      inner_join(assignments %>% filter(sort_year == yr), by = "ticker") %>%
      group_by(date, portfolio) %>%
      summarise(r = weighted.mean(ret, mktcap_june, na.rm = TRUE),
                .groups = "drop")})}

ret_bm  <- vw_ret(ports_bm)
ret_op  <- vw_ret(ports_op)
ret_inv <- vw_ret(ports_inv)


# ================================================================
# 6. CONSTRUCTING FF5 FACTORS
# ================================================================


w <- function(df, suffix) {
  df %>%
    pivot_wider(names_from = portfolio, values_from = r) %>%
    rename_with(~paste0(., suffix), -date)}

factors <- w(ret_bm, "_bm") %>%
  left_join(w(ret_op,  "_op"),  by = "date") %>%
  left_join(w(ret_inv, "_inv"), by = "date") %>%
  mutate(
    SMB = ((SL_bm + SN_bm + SH_bm) / 3 - (BL_bm + BN_bm + BH_bm) / 3 +
             (SL_op + SN_op + SH_op) / 3 - (BL_op + BN_op + BH_op) / 3 +
             (SL_inv + SN_inv + SH_inv)/ 3 - (BL_inv + BN_inv + BH_inv)/ 3) / 3,
    HML = (SH_bm + BH_bm) / 2 - (SL_bm + BL_bm) / 2,
    RMW = (SH_op + BH_op) / 2 - (SL_op + BL_op) / 2,
    CMA = (SL_inv + BL_inv)/ 2 - (SH_inv + BH_inv)/ 2
  ) %>%
  select(date, SMB, HML, RMW, CMA) %>%
  left_join(mkt, by = "date") %>%
  left_join(rf,  by = "date") %>%
  filter(!is.na(MKT), !is.na(rf))


# ================================================================
# 7. DESCRIPTIVE STATISTICS
# ================================================================


fac_stats <- factors %>%
  select(MKT, SMB, HML, RMW, CMA) %>%
  summarise(across(everything(), list(
    mean   = ~mean(.)   * 100,
    sd     = ~sd(.)     * 100,
    t_stat = ~mean(.) / (sd(.) / sqrt(n())),
    sharpe = ~mean(.) / sd(.) * sqrt(12)
  ), .names = "{.col}__{.fn}")) %>%
  pivot_longer(everything(), names_to = c("Factor", "Stat"), names_sep = "__") %>%
  pivot_wider(names_from = Stat, values_from = value) %>%
  mutate(across(-Factor, ~round(., 3)))


cat("\nFactors (% per month)\n"); print(fac_stats)
cat("\nCorrelations\n"); print(round(cor(factors %>% select(MKT, SMB, HML, RMW, CMA)), 3))


# ================================================================
# 8. TEST PORTFOLIOS: 25 (5×5 Size × B/M)
# ================================================================


ports_25 <- char_data %>%
  group_by(sort_year) %>%
  mutate(port25 = paste0("S", ntile(mktcap_june, 5),
                         "B", ntile(bm, 5))) %>%
  ungroup()


test_ret <- map_dfr(unique(ports_25$sort_year), function(yr) {
  monthly_ret %>%
    filter(date >= make_date(yr, 7, 1),
           date <= make_date(yr + 1, 6, 30)) %>%
    inner_join(ports_25 %>% filter(sort_year == yr) %>%
                 select(ticker, port25, mktcap_june), by = "ticker") %>%
    group_by(date, port25) %>%
    summarise(ret = weighted.mean(ret, mktcap_june, na.rm = TRUE),
              .groups = "drop")})


reg_data <- test_ret %>%
  left_join(factors, by = "date") %>%
  mutate(excess_ret = ret - rf) %>%
  filter(!is.na(excess_ret))


# ================================================================
# 9. FF5 REGRESSIONS FOR 25 TEST PORTFOLIOS
# ================================================================


results <- reg_data %>%
  group_by(port25) %>%
  group_modify(~{
    fit <- lm(excess_ret ~ MKT + SMB + HML + RMW + CMA, data = .x)
    bind_cols(tidy(fit),
              tibble(adj_r2 = glance(fit)$adj.r.squared, n = nrow(.x)))
  }) %>%
  ungroup()

alphas <- results %>%
  filter(term == "(Intercept)") %>%
  transmute(
    Port   = port25,
    Alpha  = round(estimate * 100, 3),
    t      = round(statistic, 2),
    p      = round(p.value, 4),
    adj_R2 = round(adj_r2, 3),
    `*`    = case_when(p.value < .01 ~ "***", 
                       p.value < .05 ~ "**",
                       p.value < .10 ~ "*",   
                       TRUE          ~ ""))

print(arrange(alphas, Port), n = 25)
cat(sprintf("\nMean|α| = %.3f%% significant (5%%) = %d/25\n", mean(abs(alphas$Alpha)), sum(alphas$p < .05)))


# ================================================================
# 9b. ROBUSTNESS: EQUAL-WEIGHTED PORTFOLIOS
# ================================================================


ew_ret <- function(assignments) {
  map_dfr(unique(assignments$sort_year), function(yr) {
    monthly_ret %>%
      filter(date >= make_date(yr, 7, 1),
             date <= make_date(yr + 1, 6, 30)) %>%
      inner_join(assignments %>% filter(sort_year == yr), by = "ticker") %>%
      group_by(date, portfolio) %>%
      summarise(r = mean(ret, na.rm = TRUE),   # equal weight
                .groups = "drop")})}


ew_ret_bm  <- ew_ret(ports_bm)
ew_ret_op  <- ew_ret(ports_op)
ew_ret_inv <- ew_ret(ports_inv)


factors_ew <- w(ew_ret_bm, "_bm") %>%
  left_join(w(ew_ret_op,  "_op"),  by = "date") %>%
  left_join(w(ew_ret_inv, "_inv"), by = "date") %>%
  mutate(
    SMB = ((SL_bm + SN_bm + SH_bm) / 3 - (BL_bm + BN_bm + BH_bm) / 3 +
             (SL_op + SN_op + SH_op) / 3 - (BL_op + BN_op + BH_op) / 3 +
             (SL_inv + SN_inv + SH_inv)/ 3 - (BL_inv + BN_inv + BH_inv)/ 3) / 3,
    HML = (SH_bm + BH_bm) / 2 - (SL_bm + BL_bm) / 2,
    RMW = (SH_op + BH_op) / 2 - (SL_op + BL_op) / 2,
    CMA = (SL_inv + BL_inv)/ 2 - (SH_inv + BH_inv)/ 2
  ) %>%
  select(date, SMB, HML, RMW, CMA) %>%
  left_join(mkt, by = "date") %>%
  left_join(rf,  by = "date") %>%
  filter(!is.na(MKT), !is.na(rf))


## EW Test Portfolios
test_ret_ew <- map_dfr(unique(ports_25$sort_year), function(yr) {
  monthly_ret %>%
    filter(date >= make_date(yr, 7, 1),
           date <= make_date(yr + 1, 6, 30)) %>%
    inner_join(ports_25 %>% filter(sort_year == yr) %>%
                 select(ticker, port25, mktcap_june), by = "ticker") %>%
    group_by(date, port25) %>%
    summarise(ret = mean(ret, na.rm = TRUE),
              .groups = "drop")})

reg_data_ew <- test_ret_ew %>%
  left_join(factors_ew, by = "date") %>%
  mutate(excess_ret = ret - rf) %>%
  filter(!is.na(excess_ret))

results_ew <- reg_data_ew %>%
  group_by(port25) %>%
  group_modify(~{
    fit <- lm(excess_ret ~ MKT + SMB + HML + RMW + CMA, data = .x)
    bind_cols(tidy(fit),
              tibble(adj_r2 = glance(fit)$adj.r.squared, n = nrow(.x)))
  }) %>%
  ungroup()

alphas_ew <- results_ew %>%
  filter(term == "(Intercept)") %>%
  transmute(
    Port   = port25,
    Alpha  = round(estimate * 100, 3),
    t      = round(statistic, 2),
    p      = round(p.value, 4),
    adj_R2 = round(adj_r2, 3),
    `*`    = case_when(p.value < .01 ~ "***", 
                       p.value < .05 ~ "**",
                       p.value < .10 ~ "*",   
                       TRUE          ~ ""))

cat("\nRobustness: Equal-weighted alphas (% month)\n"); print(arrange(alphas_ew, Port), n = 25)
cat(sprintf("\nMean|α| = %.3f%% significant (5%%) = %d/25\n", mean(abs(alphas_ew$Alpha)), sum(alphas_ew$p < .05)))


## Comparison of VW vs. EW
comparison <- alphas %>%
  select(Port, Alpha_VW = Alpha, t_VW = t) %>%
  left_join(alphas_ew %>% select(Port, Alpha_EW = Alpha, t_EW = t), by = "Port")

cat("\nComparison VW vs EW\n"); print(comparison, n = 25)


# ================================================================
# 10. GRS-TEST
# ================================================================


grs <- function(E, F_) {
  T <- nrow(E); N <- ncol(E); K <- ncol(F_)
  a <- numeric(N); U <- matrix(NA, T, N)
  for (i in 1:N) {
    m <- lm(E[, i] ~ F_)
    a[i]    <- coef(m)[1]
    U[, i]  <- resid(m)
  }
  sh2  <- as.numeric(t(colMeans(F_)) %*% solve(cov(F_)) %*% colMeans(F_))
  stat <- (T/N) * ((T-N-K)/(T-K-1)) *
    as.numeric(t(a) %*% solve(cov(U)) %*% a) / (1 + sh2)
  list(stat      = round(stat, 4),
       pval      = round(pf(stat, N, T-N-K, lower.tail = FALSE), 4),
       mean_abs_a = round(mean(abs(a)) * 100, 4))}

ok_dates <- reg_data %>%
  group_by(date) %>%
  filter(n_distinct(port25) == 25) %>%
  pull(date) %>% unique() %>% sort()

E_mat <- reg_data %>%
  filter(date %in% ok_dates) %>%
  select(date, port25, excess_ret) %>%
  pivot_wider(names_from = port25, values_from = excess_ret) %>%
  arrange(date) %>% select(-date) %>% as.matrix()

F_mat <- reg_data %>%
  filter(date %in% ok_dates, port25 == unique(reg_data$port25)[1]) %>%
  arrange(date) %>%
  select(MKT, SMB, HML, RMW, CMA) %>% as.matrix()

g <- grs(E_mat, F_mat)

cat(sprintf("\nGRS-test\nF = %.4f   p = %.4f   mean|α| = %.4f%%\n", g$stat, g$pval, g$mean_abs_a))

# ================================================================
# 11 SUBPERIOD ANALYSIS
# ================================================================


run_period <- function(label, date_from, date_to) {
  f <- factors %>% filter(date >= date_from, date < date_to)
  cat(sprintf("\nPeriod: %s (%d months)\n", label, nrow(f)))
  
  stats <- f %>%
    select(MKT, SMB, HML, RMW, CMA) %>%
    summarise(across(everything(), list(
      mean   = ~round(mean(.) * 100, 3),
      t_stat = ~round(mean(.) / (sd(.) / sqrt(n())), 2)
    ), .names = "{.col}__{.fn}")) %>%
    pivot_longer(everything(), names_to = c("Factor","Stat"), names_sep="__") %>%
    pivot_wider(names_from=Stat, values_from=value)
  print(stats)}


run_period("Until 2022",    as.Date("2016-07-01"), as.Date("2022-02-01"))
run_period("After 2022", as.Date("2022-02-01"), as.Date("2026-03-01"))


# ================================================================
# 11b. GRS TEST — FULL SAMPLE
# ================================================================


run_grs_full <- function(label, date_from, date_to, n_ports = 25) {
  
  cat(sprintf("\nGRS Test: %s\n", label))
  
  rd_sub <- reg_data %>%
    filter(date >= date_from, date < date_to)
  
  ok_dates_sub <- rd_sub %>%
    group_by(date) %>%
    filter(n_distinct(port25) == n_ports) %>%
    pull(date) %>% unique() %>% sort()
  
  cat(sprintf("Actual period: %s – %s\n", min(ok_dates_sub), max(ok_dates_sub)))
  
  T <- length(ok_dates_sub)
  N <- n_ports
  K <- 5
  
  cat(sprintf("T = %d | N = %d | K = %d | df = %d\n", T, N, K, T - N - K))
  
  if (T - N - K <= 0) {
    cat("ERROR: insufficient degrees of freedom (T - N - K <= 0)\n")
    return(invisible(NULL))
  }
  
  E_sub <- rd_sub %>%
    filter(date %in% ok_dates_sub) %>%
    select(date, port25, excess_ret) %>%
    pivot_wider(names_from = port25, values_from = excess_ret) %>%
    arrange(date) %>% select(-date) %>% as.matrix()
  
  F_sub <- rd_sub %>%
    filter(date %in% ok_dates_sub, port25 == unique(rd_sub$port25)[1]) %>%
    arrange(date) %>%
    select(MKT, SMB, HML, RMW, CMA) %>% as.matrix()
  
  g_sub <- grs(E_sub, F_sub)
  
  cat(sprintf("GRS F-statistic : %.4f\n", g_sub$stat))
  cat(sprintf("p-value         : %.4f\n", g_sub$pval))
  cat(sprintf("Mean|alpha|     : %.4f%% per month\n", g_sub$mean_abs_a))
  
  return(invisible(g_sub))}


run_grs_full("Full sample (Jul 2016 – Feb 2026)", as.Date("2016-07-01"), as.Date("2026-03-01"))


# ================================================================
# 11c. GRS TEST — SUBPERIODS (N=22, excluding chronically empty portfolios)
# ================================================================
# S2B1, S4B4, S5B5 excluded: unpopulated before Jul 2019.
# Retaining all 25 portfolios yields df = T - N - K = 31 - 25 - 5 = 1
# for the pre-sanctions subperiod, rendering the GRS statistic unreliable.
# Reducing to N=22 raises df to 40, sufficient for stable inference.


exclude_ports <- c("S2B1", "S4B4", "S5B5")

run_grs_22 <- function(label, date_from, date_to) {
  
  cat(sprintf("\nGRS Test (N=22): %s\n", label))
  
  rd_sub <- reg_data %>%
    filter(date >= date_from, date < date_to,
           !port25 %in% exclude_ports)
  
  ok_dates_sub <- rd_sub %>%
    group_by(date) %>%
    filter(n_distinct(port25) == 22) %>%
    pull(date) %>% unique() %>% sort()
  
  cat(sprintf("Actual period: %s – %s\n", min(ok_dates_sub), max(ok_dates_sub)))
  
  T <- length(ok_dates_sub)
  N <- 22
  K <- 5
  
  cat(sprintf("T = %d | N = %d | K = %d | df = %d\n", T, N, K, T - N - K))
  
  if (T - N - K <= 0) {
    cat("ERROR: insufficient degrees of freedom\n")
    return(invisible(NULL))
  }
  
  E_sub <- rd_sub %>%
    filter(date %in% ok_dates_sub) %>%
    select(date, port25, excess_ret) %>%
    pivot_wider(names_from = port25, values_from = excess_ret) %>%
    arrange(date) %>% select(-date) %>% as.matrix()
  
  F_sub <- rd_sub %>%
    filter(date %in% ok_dates_sub,
           port25 == unique(rd_sub$port25)[1]) %>%
    arrange(date) %>%
    select(MKT, SMB, HML, RMW, CMA) %>% as.matrix()
  
  g_sub <- grs(E_sub, F_sub)
  
  cat(sprintf("GRS F-statistic : %.4f\n", g_sub$stat))
  cat(sprintf("p-value         : %.4f\n", g_sub$pval))
  cat(sprintf("Mean|alpha|     : %.4f%% per month\n", g_sub$mean_abs_a))
  
  return(invisible(g_sub))}


run_grs_22("Pre-sanctions  (Jul 2016 – Jan 2022)", as.Date("2016-07-01"), as.Date("2022-02-01"))
run_grs_22("Post-sanctions (Feb 2022 – Feb 2026)", as.Date("2022-02-01"), as.Date("2026-03-01"))
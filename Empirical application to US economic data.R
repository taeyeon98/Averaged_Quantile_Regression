################################################################################
## 0. Preparations -------------------------------------------------------------
################################################################################
library(lubridate)
library(dplyr)
library(nowcasting)     # PRTDB, Bpanel, nowcast, etc.
library(quantreg)
library(Metrics)        # rmse
library(moments)        # skewness, kurtosis

## Data ------------------------------------------------------------------------
kocsv <- read.csv("C:\\Users\\US_economic_data.csv", header = TRUE)
kocsv$X <- as.Date(kocsv$X)

## (A) End‑date sequences ------------------------------------------------------
## Rolling‑window end‑date sequences (12 different windows)
enddate_list <- list(
  # ① Post‑dot‑com low‑volatility period (4 yr)
  "2001Q1–2004Q4" = seq(ymd("2001-01-01"), ymd("2005-01-01"), by = "quarter"),

  # ② Low‑rate housing boom (3 yr)
  "2003Q1–2005Q4" = seq(ymd("2003-01-01"), ymd("2006-01-01"), by = "quarter"),

  # ③ Sub‑prime warning to pre‑crisis (2.5 yr)
  "2005Q3–2007Q4" = seq(ymd("2005-07-01"), ymd("2008-01-01"), by = "quarter"),

  # ④ Global Financial Crisis core (high kurtosis, 2 yr)
  "2007Q3–2009Q2" = seq(ymd("2007-07-01"), ymd("2009-07-01"), by = "quarter"),

  # ⑤ Immediate post‑crisis + QE1 (3 yr)
  "2008Q3–2011Q2" = seq(ymd("2008-07-01"), ymd("2011-07-01"), by = "quarter"),

  # ⑥ Trump tax‑cut & trade dispute (3 yr)
  "2016Q1–2018Q4" = seq(ymd("2016-01-01"), ymd("2019-01-01"), by = "quarter"),

  # ⑦ Two years just before COVID‑19
  "2018Q1–2019Q4" = seq(ymd("2018-01-01"), ymd("2020-01-01"), by = "quarter"),

  # ⑧ Full 2000s decade (10 yr window)
  "2001Q1–2010Q4" = seq(ymd("2001-01-01"), ymd("2011-01-01"), by = "quarter"),

  # ⑨ Most recent 11 yr (taper‑tantrum → pandemic → inflation surge)
  "2013Q1–2023Q4" = seq(ymd("2013-01-01"), ymd("2024-01-01"), by = "quarter"),

  # ⑩ Pre‑crisis through QE1 completion (5 yr)
  "2006Q1–2010Q4" = seq(ymd("2006-01-01"), ymd("2011-01-01"), by = "quarter"),

  # ⑪ Crisis core + early aftermath (high kurtosis, 3.5 yr)
  "2007Q3–2010Q4" = seq(ymd("2007-07-01"), ymd("2011-01-01"), by = "quarter"),

  # ⑫ US‑China trade‑war window (high kurtosis, 3.5 yr)
  "2017Q1–2020Q2" = seq(ymd("2017-01-01"), ymd("2020-07-01"), by = "quarter")
)

train_year <- 4   # training‑window length (years)

################################################################################
## 1. Core computation function (original logic kept) --------------------------
################################################################################
calc_metrics <- function(enddate_vec) {

  leng <- length(enddate_vec)
  MSE_lm  <- numeric(leng)
  MSE_QR  <- numeric(leng)
  MSE_AQR <- numeric(leng)
  skew <- numeric(leng)
  kurt <- numeric(leng)

  for (i in seq_len(leng)) {

    ## tryCatch: skip this step if an error/NaN occurs
    tryCatch({
      q           <- enddate_vec[i]
      start_date  <- q %m-% years(train_year)
      start_year  <- year(start_date)
      start_mon   <- month(start_date)
      start_qua   <- (start_mon + 2) / 3

      ## ---------- Split X and CPI -------------------------------------------
      kocsv_x <- kocsv |>
        select(-CPI) |>
        filter(X >= start_date & X < q)

      cpi_vec <- kocsv |>
        filter(X >= start_date & X < q) |>
        pull(CPI)
      if (length(cpi_vec) < 12) stop("insufficient sample")

      cpi_q <- cpi_vec[seq(1, length(cpi_vec), by = 3)]
      cpi_q[length(cpi_q)] <- NA        # target to forecast

      ## ---------- Factor nowcasting -----------------------------------------
      base   <- ts(kocsv_x[,-1], start = c(start_year, start_mon), frequency = 12)
      trans  <- rep(3, ncol(base))
      delay  <- rep(0, ncol(base))

      GDP_qtr <- ts(cpi_q, start = c(start_year, start_qua), frequency = 4)

      BRGDP   <- list(base = base, trans = trans, delay = delay)
      vintage <- PRTDB(mts = BRGDP$base, delay = BRGDP$delay,
                       vintage = "2019-10-01")
      base_v  <- window(vintage, start = c(start_year, start_mon), frequency = 12)
      x_panel <- Bpanel(base = base_v, trans = BRGDP$trans)

      y_month <- diff(diff(GDP_qtr, 4))
      y_month <- qtr2month(y_month)

      data_all  <- cbind(y_month, x_panel)
      frequency <- c(4, rep(12, ncol(x_panel)))

      now_res <- nowcast(
        formula   = y_month ~ .,
        data      = data_all,
        r = 2, q = 2, p = 1, method = "2s",
        frequency = frequency
      )

      factors_ts  <- now_res$factors$dynamic_factors
      factors_qtr <- month2qtr(factors_ts)

      ## ---------- Regression data (1‑quarter‑ahead) -------------------------
      y_target <- diff(diff(GDP_qtr, 4))
      y_vec    <- as.vector(y_target)

      f1 <- factors_qtr[,1]; f1_lag <- f1[6:(length(f1) - 5)]
      f2 <- factors_qtr[,2]; f2_lag <- f2[6:(length(f2) - 5)]
      y_lag <- y_vec[1:(length(y_vec) - 1)]

      regress_df <- data.frame(y_lag, f1_lag, f2_lag)

      ## ---------- (A) OLS ----------------------------------------------------
      lm_fit <- lm(y_lag ~ f1_lag + f2_lag, data = regress_df)
      lm_pred_full <- lm_fit$coefficients[1] +
        lm_fit$coefficients[2]*f1 +
        lm_fit$coefficients[3]*f2
      pred_OLS <- tail(lm_pred_full, 5)[1]

      ## ---------- (B) QR / AQR ----------------------------------------------
      tau_vec <- c(0.25, 0.50, 0.75)
      qr_pred_each <- sapply(tau_vec, function(tau) {
        rq_fit <- rq(y_lag ~ f1_lag + f2_lag, data = regress_df, tau = tau)
        rq_pred_full <- rq_fit$coefficients[1] +
          rq_fit$coefficients[2]*f1 +
          rq_fit$coefficients[3]*f2
        tail(rq_pred_full, 5)[1]
      })
      pred_QR  <- qr_pred_each[2]             # τ = 0.50
      pred_AQR <- mean(qr_pred_each)          # equal‑weight average

      ## ---------- Store metrics ---------------------------------------------
      y_real_full <- diff(diff(ts(cpi_vec, frequency = 4), 4))
      y_real <- tail(as.vector(y_real_full), 5)[1]

      MSE_lm[i]  <- rmse(y_real, pred_OLS)
      MSE_QR[i]  <- rmse(y_real, pred_QR)
      MSE_AQR[i] <- rmse(y_real, pred_AQR)

      skew[i] <- skewness(y_lag, na.rm = TRUE)
      kurt[i] <- kurtosis(y_lag, na.rm = TRUE)

    }, error = function(e) {
      ## If an error occurs, store NA (ignored in mean calculations)
      MSE_lm[i]  <<- NA;  MSE_QR[i]  <<- NA;  MSE_AQR[i] <<- NA
      skew[i]    <<- NA;  kurt[i]    <<- NA
    })
  } # end for i

  list(
    mean_abs_skew = mean(abs(skew), na.rm = TRUE),
    mean_kurt     = mean(kurt,      na.rm = TRUE),
    RMSE_OLS      = mean(MSE_lm,    na.rm = TRUE),
    RMSE_QR       = mean(MSE_QR,    na.rm = TRUE),
    RMSE_AQR      = mean(MSE_AQR,   na.rm = TRUE)
  )
}

################################################################################
## 2. Compute & summarise results per window -----------------------------------
################################################################################
results_all <- lapply(enddate_list, calc_metrics)

summary_tbl <- do.call(
  rbind,
  Map(function(name, res) {
    data.frame(
      Window        = name,
      mean_abs_skew = round(res$mean_abs_skew, 3),
      mean_kurtosis = round(res$mean_kurt, 3),
      RMSE_OLS      = round(res$RMSE_OLS, 3),
      RMSE_QR       = round(res$RMSE_QR, 3),
      RMSE_AQR      = round(res$RMSE_AQR, 3),
      check.names   = FALSE
    )
  },
  names(enddate_list), results_all)
)

print(summary_tbl, row.names = FALSE)

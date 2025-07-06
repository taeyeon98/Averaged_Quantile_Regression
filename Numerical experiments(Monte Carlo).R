################################################################################
## (1) Monte Carlo 
################################################################################
library(doParallel)
library(foreach)
library(lubridate)
library(dplyr)
library(astsa)
library(quantreg)
library(Metrics)
library(mvtnorm)
library(EnvStats)
library(moments)
library(VGAM)
library(forecast)
library(TSPred)
library(ggplot2)
library(nowcasting)  
library(stabledist)  
library(knitr)     

nCores <- parallel::detectCores()
cl <- makeCluster(nCores)
registerDoParallel(cl)

alphaVec <- seq(0.4, 2.0, by=0.2)
betaVec  <- c(-1, -0.5, 0, 0.5, 1)

stable_scenarios <- list()
cnt <- 1
for(a in alphaVec){
  for(b in betaVec){
    stable_scenarios[[cnt]] <- list(
      name   = paste0("stable_a",a,"_b",b),
      distType="stable",
      alpha=a,
      beta=b,
      gamma=1,
      delta=0
    )
    cnt <- cnt+1
  }
}
all_scenarios <- stable_scenarios

simOne <- function(w, scenario) {
 
  distType <- scenario$distType
  alpha    <- scenario$alpha
  beta     <- scenario$beta
  gamma    <- scenario$gamma
  delta    <- scenario$delta
 
  N <- 10
  T <- 300
  rho <- 0.9
  d   <- 0.5
  tau <- 0.5
  u   <- 0.1
 
  FM <- matrix(0, nrow=T, ncol=2)
  for (i in 2:T) {
    u1 <- rnorm(1, 0, sqrt(1 - rho^2))
    u2 <- rnorm(1, 0, sqrt(1 - rho^2))
    FM[i,1] <- rho*FM[i-1,1] + u1
    FM[i,2] <- rho*FM[i-1,2] + u2  
  }
 
  LM <- matrix(rnorm(N*2, 0, 1), nrow=N, ncol=2)
 
  # idiosyncratic 
  alpha_M <- numeric(N)
  beta_M  <- numeric(N)
  for (i in 1:N) {
    beta_M[i]  <- runif(1, u, 1-u)
    b          <- beta_M[i]
    alpha_M[i] <- b*sum(LM[i,]^2)/(1 - b)
  }
 
  TM <- matrix(0, nrow=N, ncol=N)
  for (i in 1:N) {
    for (j in 1:N) {
      TM[i,j] <- sqrt(alpha_M[i]*alpha_M[j])*(0.5^abs(i-j))*(1 - d^2)
    }
  }
 
  # idiosyncratic error
  eps_M <- matrix(0, nrow=N, ncol=T)
  for (i in 2:T) {
    V_M <- mvtnorm::rmvnorm(1, mean=rep(0,N), sigma=TM)
    eps_M[, i] <- as.vector(V_M) + d * eps_M[, i-1]
  }
 
  CPI <- numeric(T)
  for(i in 1:T){
    if(distType == "stable"){
      ranu <- rstable(1, alpha=alpha, beta=beta, gamma=gamma, delta=delta)
      CPI[i] <- FM[i,1] + FM[i,2] + ranu
    } else {
      stop("Only stable scenario is used.")
    }
  }
 
  X_M <- matrix(0, nrow=N, ncol=T)
  for (i in 1:T) {
    X_M[,i] <- LM %*% FM[i,] + eps_M[,i]
  }
 
  X <- lubridate::ymd("2001-01-01") + months(0:(T-1))
  mon_sam_2 <- cbind(X, CPI, t(X_M))
  mon_sam   <- as.data.frame(mon_sam_2)
 
  kocsv      <- mon_sam
  kocsv$X    <- as.Date(kocsv$X, origin="1970-01-01")
 
  enddate_local <- seq(lubridate::ymd('2007-01-01'),
                       lubridate::ymd('2020-01-01'),
                       by='years')
  q <- enddate_local[10]
  train_year <- 4
  start_date <- q %m-% years(train_year)
  start_year <- as.numeric(substring(start_date, 1,4))
 
  kocsv2 <- kocsv %>%
    dplyr::select(-CPI) %>%
    dplyr::filter(X >= start_date & X < q)
 
  cpi_vec <- kocsv %>%
    dplyr::filter(X >= start_date & X < q) %>%
    dplyr::pull(CPI)
 
  cpi_std <- (cpi_vec - mean(cpi_vec)) / sd(cpi_vec)
  cpi2    <- cpi_std[seq(1,length(cpi_std), by=3)]
 
  cpi3 <- cpi2
 
  GDP_qtr <- ts(cpi3, start=c(start_year,1), frequency=4)
  ppp     <- as.vector(GDP_qtr)
 
  base <- ts(kocsv2[,-1], start=c(start_year,1), frequency=12)
  trans <- rep(0,N)
  delay <- rep(0,N)
 
  GDP_qtr_ts <- ts(ppp, start=c(start_year,1), frequency=4)
  BRGDP   <- list(base=base, trans=trans, delay=delay)
 
  vintage <- PRTDB(mts=BRGDP$base, delay=BRGDP$delay, vintage="2019-10-01")
  base    <- window(vintage, start=c(start_year,1), frequency=12)
  x       <- Bpanel(base=base, trans=BRGDP$trans)
 
  y        <- qtr2month(GDP_qtr_ts)
  data     <- cbind(y,x)
  frequency<- c(4, rep(12, ncol(x)))
 
  now11 <- nowcast(
    formula=y ~ ., data=data, r=2, q=2, p=1, method="2s",
    frequency=frequency
  )
 
  fatoresTS  <- now11$factors$dynamic_factors
  fatoresTRI <- month2qtr(fatoresTS)
  ppp_len    <- length(ppp)
 
  www  <- fatoresTRI[,1]
  www2 <- www[1:(length(www)-5)]
  wwww <- fatoresTRI[,2]
  www3 <- wwww[1:(length(wwww)-5)]
 
  www1 <- ppp[1:(ppp_len - 1)]
  dop  <- data.frame(www1, www2, www3)
 
  fit <- lm(www1 ~ www2 + www3, data=dop)
  final <- fit$coefficients[1] +
    fit$coefficients[2]*www +
    fit$coefficients[3]*wwww
  predict_final_lm <- final[length(final)]
 
  tau_vec <- c(0.25, 0.50, 0.75)
  pred_vals <- numeric(length(tau_vec))
 
  for(i in seq_along(tau_vec)){
    rqfit <- rq(www1 ~ www2 + www3, data=dop, tau=tau_vec[i])
    final_i <- rqfit$coefficients[1] +
      rqfit$coefficients[2]*www +
      rqfit$coefficients[3]*wwww
    pred_vals[i] <- final_i[length(final_i)]
  }
 
  # Q1 ~ Q4
  w1 <- c(1/3,1/3,1/3)
  predict_final_q1 <- sum(pred_vals * w1)
 
  w2 <- c(0.6,0.3,0.1)
  predict_final_q2 <- sum(pred_vals * w2)
 
  w3 <- c(0.1,0.3,0.6)
  predict_final_q3 <- sum(pred_vals * w3)
 
  w4_raw <- c(1,2,1)
  w4 <- w4_raw / sum(w4_raw)
  predict_final_q4 <- sum(pred_vals * w4)
 
  # Q5, Q6, Q7
  predict_final_q5 <- pred_vals[2]
  predict_final_q6 <- pred_vals[1]
  predict_final_q7 <- pred_vals[3]
 
  # ARIMA, Naive
  arima_data <- ppp[1:(ppp_len - 1)]
  arima_for  <- TSPred::arimapred(arima_data, n.ahead=1)
  predict_final_ar <- as.vector(arima_for)[1]
 
  predict_final_ba <- ppp[ppp_len - 1]
  real <- ppp[ppp_len]
 
  # MSE
  mse_ols <- mse(predict_final_lm, real)
  mse_q1  <- mse(predict_final_q1, real)
  mse_q2  <- mse(predict_final_q2, real)
  mse_q3  <- mse(predict_final_q3, real)
  mse_q4  <- mse(predict_final_q4, real)
  mse_q5  <- mse(predict_final_q5, real)
  mse_q6  <- mse(predict_final_q6, real)
  mse_q7  <- mse(predict_final_q7, real)
  mse_ar  <- mse(predict_final_ar, real)
  mse_nv  <- mse(predict_final_ba, real)
 
  # Skewness, Kurtosis
  out_var  <- var(www1)
  out_skew <- skewness(www1) 
  out_kurt <- kurtosis(www1) 
 
  c(mse_ols, mse_q1, mse_q2, mse_q3, mse_q4, mse_q5, mse_q6, mse_q7,
    mse_ar, mse_nv, out_var, out_skew, out_kurt)
}

## 2) Monte Carlo
results_table <- data.frame()
registerDoParallel(cl)

for(sce in all_scenarios){
 
  cat("\n=== Running scenario:", sce$name, " ===\n")
  mon <- 3000
 
  res <- foreach(w = 1:mon, .combine=rbind,
                 .errorhandling='remove',
                 .packages=c("lubridate","dplyr","astsa","quantreg","Metrics",
                             "mvtnorm","EnvStats","moments","VGAM","forecast",
                             "TSPred","ggplot2","nowcasting","stabledist")) %dopar% {
    simOne(w, sce)
  }
 
  cm <- colMeans(res, na.rm=TRUE)
  sdev <- apply(res, 2, sd)
  nRep <- nrow(res)
  se   <- sdev / sqrt(nRep)
 
  # simOne() 
  # 1:MSE_OLS, 2:MSE_Q1, 3:MSE_Q2, 4:MSE_Q3, 5:MSE_Q4,
  # 6:MSE_Q5, 7:MSE_Q6, 8:MSE_Q7, 9:MSE_AR, 10:MSE_NV,
  # 11:Var, 12:Skew, 13:Kurt
 
  mse_ols    <- cm[1];  mse_ols_se <- se[1]
  mse_q1     <- cm[2];  mse_q1_se  <- se[2]
  mse_q2     <- cm[3];  mse_q2_se  <- se[3]
  mse_q3     <- cm[4];  mse_q3_se  <- se[4]
  mse_q4     <- cm[5];  mse_q4_se  <- se[5]
  mse_q5     <- cm[6];  mse_q5_se  <- se[6]
  mse_q6     <- cm[7];  mse_q6_se  <- se[7]
  mse_q7     <- cm[8];  mse_q7_se  <- se[8]
  mse_ar     <- cm[9];  mse_ar_se  <- se[9]
  mse_nv     <- cm[10]; mse_nv_se  <- se[10]
 
  mvar       <- cm[11]; var_se     <- se[11]
  mskew      <- cm[12]; skew_se    <- se[12]  
  mkurt      <- cm[13]; kurt_se    <- se[13]  
 
  newrow <- data.frame(
    Scenario    = sce$name,
    alpha       = sce$alpha,
    beta        = sce$beta,
   
    MSE_OLS     = mse_ols,
    MSE_OLS_SE  = mse_ols_se,
    MSE_Q1      = mse_q1,
    MSE_Q1_SE   = mse_q1_se,
    MSE_Q2      = mse_q2,
    MSE_Q2_SE   = mse_q2_se,
    MSE_Q3      = mse_q3,
    MSE_Q3_SE   = mse_q3_se,
    MSE_Q4      = mse_q4,
    MSE_Q4_SE   = mse_q4_se,
    MSE_Q5      = mse_q5,
    MSE_Q5_SE   = mse_q5_se,
    MSE_Q6      = mse_q6,
    MSE_Q6_SE   = mse_q6_se,
    MSE_Q7      = mse_q7,
    MSE_Q7_SE   = mse_q7_se,
   
    MSE_AR      = mse_ar,
    MSE_AR_SE   = mse_ar_se,
    MSE_NV      = mse_nv,
    MSE_NV_SE   = mse_nv_se,
   
    Var         = mvar,
    Var_SE      = var_se,
    Skew        = mskew,
    Skew_SE     = skew_se, 
    Kurt        = mkurt,
    Kurt_SE     = kurt_se  
  )
 
  results_table <- rbind(results_table, newrow)
}

stopCluster(cl)

cat("\n===== 최종 결과 =====\n")
print(results_table)

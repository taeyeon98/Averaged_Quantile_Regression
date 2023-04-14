#clean environment 
remove(list=ls())

#required packages
library(dplyr)
library(nowcasting)
#Package 'nowcasting' was removed from the CRAN repository. Available versions can be obtained from the archive.
library(astsa)
library(lubridate)
library(quantreg)
library(Metrics)
library(mvtnorm)
library(EnvStats)
library(progress)
library(moments)
library(VGAM)
library(forecast)
library(TSPred)
library(ggplot2)


#Monte_Carlo

#Number of repetitions
mon<-2000


err1<-rep(0,mon)
err2<-rep(0,mon)
err3<-rep(0,mon)
err4<-rep(0,mon)
err5<-rep(0,mon)
err6<-rep(0,mon)
err7<-rep(0,mon)

skew<-0
kurt<-0
vari<-0

pb <- progress_bar$new(total = mon)
enddate<-seq(ymd('2007-01-01'),ymd('2020-01-01'), by = 'years')


for (w in 1:mon) {
  
  #Parameter setting related to data generation.
  N<-10
  T<-300
  rho<-0.9
  d<-0.5
  tau<-0.5
  u<-0.1
  
  #Data generation
  LM<-matrix(, nrow = N, ncol = 2)
  FM<-matrix(,nrow = T, ncol = 2)
  eps_M<-matrix(,nrow = N, ncol = T)
  X_M<-matrix(,nrow = N, ncol = T)
  TM<-matrix(,nrow = N, ncol = N)
  CPI<-rep(0, T)
  
  FM[1,]<-c(0,0)
  
  for (i in 2:T) {
    u1<-rnorm(1,0,sqrt(1-rho^2))
    u2<-rnorm(1,0,sqrt(1-rho^2))
    FM[i,1]<-rho*FM[i-1,1]+u1
    FM[i,2]<-rho*FM[i-1,2]+u2  
  }
  
  for (i in 1:N) {
    LM[i,1]<-rnorm(1,0,1)
    LM[i,2]<-rnorm(1,0,1)  
  }
  
  alpha_M<-rep(0,N)
  beta_M<-rep(0,N)
  
  for (i in 1:N) {
    beta_M[i]<-runif(1,u,1-u)
    b<-beta_M[i]
    alpha_M[i]<-b*sum((LM[i,])^2)/(1-b)
  }
  
  
  for (i in 1:N) {
    for (j in 1:N) {
      TM[i,j]<-sqrt(alpha_M[i]*alpha_M[j])*(tau^abs(i-j))*(1-d^2)
    }
  }
  
  
  eps_M[,1]<-rep(0,N)
  for (i in 2:T) {
    V_M2<-rmvnorm(1,mean=rep(0,N), sigma=TM)
    V_M<-as.vector(V_M2)
    eps_t<-as.vector(eps_M[,i-1])
    eps_p<-V_M+d*eps_t
    eps_M[,i]<-eps_p
  }
  
  #error term setting
  #In addition to the following distribution, the code 'rnorm' can also be applied in the normal distribution and 'rt' in the t-distribution.
  for (i in 1:T) {
    ranu<-rnormMix(1, mean1 = 0, sd1 = 2, mean2 = 2, sd2 = 1, p.mix = 0.8)
    CPI[i]<-FM[i,1]+FM[i,2]+ranu
  }
  
  
  for (i in 1:T) {
    X_M[,i]<-LM%*%FM[i,]+eps_M[,i]
  }
  
  X<-ymd("2001-01-01") + months(0:(T-1))
  mon_sam_2<-cbind(X, CPI,t(X_M))
  mon_sam<-as.data.frame(mon_sam_2)
  
  kocsv<-mon_sam
  kocsv$X<-as.Date(kocsv$X,origin = "1970-01-01")
  
  q<-enddate[10]
  train_year<-4
  styear<-as.numeric(substring(q, 3,4))
  start_date<-q %m-% years(train_year)
  start_year<-as.numeric(substring(start_date, 1,4))
  
  
  kocsv2<-kocsv%>%
    select(-CPI)%>%
    filter(X>=start_date &X<q)
  
  cpi<-kocsv%>%
    filter(X>=start_date &X<q)%>%
    select(CPI)
  
  cpi<-as.vector(cpi$CPI)
  
  cpi<-(cpi-mean(cpi))/sd(cpi)
  cpi2<-cpi[seq(1, length(cpi), by=3)]
  cpi3<-replace(cpi2, length(cpi2),NA)
  
  base <- ts(kocsv2[,-1], start=c(start_year, 1), frequency=12)
  trans<-rep(0,N)
  delay<-rep(0,N)
  
  GDP_qtr<-ts(cpi3,start = c(start_year,1), frequency=4)
  BRGDP<-list(base=base,trans=trans,delay=delay)
  
  vintage <- PRTDB(mts = BRGDP$base, delay = BRGDP$delay, vintage = "2019-10-01")
  base <- window(vintage, start = c(start_year,1), frequency = 12)
  x <- Bpanel(base = base, trans = BRGDP$trans)
  
  y <- qtr2month(GDP_qtr)
  
  data <- cbind(y,x)
  frequency <- c(4,rep(12,ncol(x)))
  
  #Dynamic Factor Extraction
  now11 <- nowcast(formula = y~., data = data, r = 2, q = 2 , p = 1, method = "2s",
                   frequency = frequency)
  
  fatoresTS<-now11$factors$dynamic_factors
  fatoresTRI <- month2qtr(fatoresTS)
  ppp<-as.vector(GDP_qtr)
  
  www<-fatoresTRI[,1]
  www2<-www[1:(length(www)-5)]
  wwww<-fatoresTRI[,2]
  www3<-wwww[1:(length(www)-5)]
  www1<-ppp[1:(length(ppp)-1)]
  dop<-as.data.frame(cbind(www1, www2, www3))
  fit<-lm(www1~www2+www3, data = dop)
  
  #OLS approach
  final<-fit$coefficients[1]+(fit$coefficients[2])*(www)+(fit$coefficients[3])*(wwww)
  predict_final_lm<-(as.vector(final))[(length(final)-4):length(final)]
  
  
  #Quantile approaches
  #The following three models are divided according to the weights set in quantile combination.
  
  qvector<-rep(0,5)
  for (j in 1:2) {
    rqfit<-rq(www1~www2+www3, data = dop, tau=j/3)
    final_2<-rqfit$coefficients[1]+(rqfit$coefficients[2])*(www)+(rqfit$coefficients[3])*(wwww)
    predict_final_2<-(as.vector(final_2))[(length(final_2)-4):length(final_2)]
    qvector<-qvector+predict_final_2
  }
  predict_final_qu<-qvector/2
  
  qvector_2<-rep(0,5)
  for (j in 1:3) {
    rqfit<-rq(www1~www2+www3, data = dop, tau=j/4)
    final_2<-rqfit$coefficients[1]+(rqfit$coefficients[2])*(www)+(rqfit$coefficients[3])*(wwww)
    predict_final_2<-(as.vector(final_2))[(length(final_2)-4):length(final_2)]
    qvector_2<-qvector_2+predict_final_2
  }
  predict_final_qu_2<-qvector_2/3
  
  
  qvector_3<-rep(0,5)
  for (j in 1:4) {
    rqfit<-rq(www1~www2+www3, data = dop, tau=j/5)
    final_2<-rqfit$coefficients[1]+(rqfit$coefficients[2])*(www)+(rqfit$coefficients[3])*(wwww)
    predict_final_2<-(as.vector(final_2))[(length(final_2)-4):length(final_2)]
    qvector_3<-qvector_3+predict_final_2
  }
  predict_final_qu_3<-qvector_3/4
  
  qvector_4<-rep(0,5)
  for (j in 1:4) {
    rqfit<-rq(www1~www2+www3, data = dop, tau=j/5)
    final_2<-rqfit$coefficients[1]+(rqfit$coefficients[2])*(www)+(rqfit$coefficients[3])*(wwww)
    predict_final_2<-(as.vector(final_2))[(length(final_2)-4):length(final_2)]
    qvector_4<-qvector_4+predict_final_2*(2.5-abs(2.5-j))
  }
  predict_final_qu_4<-qvector_4/6
  
  
  cpi4<-kocsv%>%
    filter(X>=start_date &X<enddate[11])%>%
    select(CPI)
  
  cpi5<-as.vector(cpi4$CPI)
  cpi5<-(cpi5-mean(cpi5))/sd(cpi5)
  cpi6<-cpi5[seq(1, length(cpi5), by=3)]
  GDP_qtr_2<-ts(cpi6,start = c(start_year,1), frequency=4)
  ppp2<-as.vector(GDP_qtr_2)
  #Baseline model 
  real<-ppp2[(length(ppp2)-4):length(ppp2)]
  
  #Arima model 
  arima_data<-ppp2[1:(length(ppp2)-5)]
  arima_for<-as.vector(arimapred(arima_data, n.ahead = 5))
  
  
  predict_final_lm<-predict_final_lm[1]
  predict_final_qu<-predict_final_qu[1]
  real<-real[1]
  predict_final_ar<-arima_for[1]
  predict_final_ba<-ppp2[length(ppp2)-5]
  
  err1[w]<-mse(predict_final_lm,real)
  err2[w]<-mse(predict_final_qu,real)
  err3[w]<-mse(predict_final_qu_2,real)
  err4[w]<-mse(predict_final_qu_3,real)
  err5[w]<-mse(predict_final_qu_4,real)
  err6[w]<-mse(predict_final_ar, real)
  err7[w]<-mse(predict_final_ba, real)
  
  vari<-vari+var(www1)
  skew<-skew+abs(skewness(www1))
  kurt<-kurt+kurtosis(www1)
  pb$tick()
  
}

#Skewness, variance, and kurtosis of the data used in model fitting
vari<-vari/mon
sk<-skew/mon
ku<-kurt/mon

#Ratio of MSE to OLS
sum(err2)/sum(err1)
sum(err3)/sum(err1)
sum(err4)/sum(err1)
sum(err5)/sum(err1)
sum(err6)/sum(err1)
sum(err7)/sum(err1)

sk
ku
vari

output1<-c(N, 50, 2, rho, d, tau, u, sk, vari, ku,
           sum(err2)/sum(err1), sum(err3)/sum(err1), sum(err4)/sum(err1), 
           sum(err5)/sum(err1), sum(err6)/sum(err1), sum(err7)/sum(err1))

output1

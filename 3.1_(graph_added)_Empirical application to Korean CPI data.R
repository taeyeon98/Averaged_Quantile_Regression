
#The CSV file can be obtained by downloading the 'inflation data' file at the main page. 
kocsv<-read.csv("C:\\Users\\inflation_data.csv", header = TRUE)
kocsv$X<-as.Date(kocsv$X)


#Initial value setting 
enddate<-seq(ymd('2007-01-01'),ymd('2020-01-01'), by = 'quarter')
enddate[15]
leng<-15
train_year<-4

MSE_lm<-rep(0,leng)
MSE_qu_1<-rep(0,leng)
MSE_qu_2<-rep(0,leng)
MSE_qu_3<-rep(0,leng)
MSE_ba<-rep(0,leng)
MSE_sa<-rep(0,leng)
skew<-rep(0,leng)


#One quarter ahead
for (i in 1:leng) {
  q<-enddate[i]
  styear<-as.numeric(substring(q, 3,4))
  start_date<-q %m-% years(train_year)
  start_year<-as.numeric(substring(start_date, 1,4))
  start_mon<-as.numeric(substring(start_date, 6,7))
  start_qua<-(start_mon+2)/3
  
  start_date_2<-start_date %m-% months(3)
  styear_2<-as.numeric(substring(start_date_2, 3,4))
  start_year_2<-as.numeric(substring(start_date_2, 1,4))
  start_mon_2<-as.numeric(substring(start_date_2, 6,7))
  start_qua_2<-(start_mon_2+2)/3
  
  kocsv2<-kocsv%>%
    select(-CPI)%>%
    filter(X>=start_date &X<q)
  
  cpi<-kocsv%>%
    filter(X>=start_date &X<q)%>%
    select(CPI)
  
  cpi<-as.vector(cpi$CPI)
  cpi2<-cpi[seq(1, length(cpi), by=3)]
  cpi3<-replace(cpi2, length(cpi2),NA)
  
  base <- ts(kocsv2[,-1], start=c(start_year,start_mon ), frequency=12)
  trans<-rep(3,14)
  delay<-rep(0,14)
  
  GDP_qtr<-ts(cpi3,start = c(start_year,start_qua), frequency=4)
  BRGDP<-list(base=base,trans=trans,delay=delay)
  
  vintage <- PRTDB(mts = BRGDP$base, delay = BRGDP$delay, vintage = "2019-10-01")
  base <- window(vintage, start = c(start_year,start_mon), frequency = 12)
  x <- Bpanel(base = base, trans = BRGDP$trans)
  
  y <- diff(diff(GDP_qtr,4))
  y <- qtr2month(y)
  
  data <- cbind(y,x)
  frequency <- c(4,rep(12,ncol(x)))
  
  now11 <- nowcast(formula = y~., data = data, r = 2, q = 2 , p = 1, method = "2s",
                   frequency = frequency)
  
  fatoresTS<-now11$factors$dynamic_factors
  fatoresTRI <- month2qtr(fatoresTS)
  ppp<-diff(diff(GDP_qtr,4))
  ppp<-as.vector(ppp)
  
  www<-fatoresTRI[,1]
  www2<-www[6:(length(www)-5)]
  wwww<-fatoresTRI[,2]
  www3<-wwww[6:(length(www)-5)]
  www1<-ppp[1:(length(ppp)-1)]
  dop<-as.data.frame(cbind(www1, www2, www3))
  fit<-lm(www1~www2+www3, data = dop)
  
  final<-fit$coefficients[1]+(fit$coefficients[2])*(www)+(fit$coefficients[3])*(wwww)
  predict_final_lm<-(as.vector(final))[(length(final)-4):length(final)]
  
  #quantile reg
  qvector_1<-rep(0,5)
  for (j in 1:3) {
    rqfit<-rq(www1~www2+www3, data = dop, tau=j/4)
    final_2<-rqfit$coefficients[1]+(rqfit$coefficients[2])*(www)+(rqfit$coefficients[3])*(wwww)
    predict_final_2<-(as.vector(final_2))[(length(final_2)-4):length(final_2)]
    qvector_1<-qvector_1+predict_final_2
  }
  predict_final_qu_1<-qvector_1/3
  
  qvector_2<-rep(0,5)
  for (j in 1:4) {
    rqfit_2<-rq(www1~www2+www3, data = dop, tau=j/5)
    final_3<-rqfit_2$coefficients[1]+(rqfit_2$coefficients[2])*(www)+(rqfit_2$coefficients[3])*(wwww)
    predict_final_3<-(as.vector(final_3))[(length(final_3)-4):length(final_3)]
    qvector_2<-qvector_2+predict_final_3
  }
  predict_final_qu_2<-qvector_2/4
  
  qvector_3<-rep(0,5)
  for (j in 1:4) {
    rqfit_3<-rq(www1~www2+www3, data = dop, tau=j/5)
    final_4<-rqfit_3$coefficients[1]+(rqfit_3$coefficients[2])*(www)+(rqfit_3$coefficients[3])*(wwww)
    predict_final_4<-(as.vector(final_3))[(length(final_3)-4):length(final_3)]
    qvector_3<-qvector_3+predict_final_4*(2.5-abs(2.5-j))
  }
  predict_final_qu_3<-qvector_3/6
  
  cpi4<-kocsv%>%
    filter(X>=start_date)%>%
    select(X, CPI)
  
  cpi5<-as.vector(cpi4$CPI)
  cpi6<-cpi5[seq(1, length(cpi5), by=3)]
  
  GDP_qtr_2<-ts(cpi6,start = c(start_year,start_qua), frequency=4)
  
  y2 <- diff(diff(GDP_qtr_2,4))
  
  lowyear<-4+start_year_2
  upyear<-5+start_year_2
  real<-as.vector(window(y2, start=c(lowyear, start_qua_2), end=c(upyear, start_qua_2)))
  
  GDP_qtr_3<-ts(lag(cpi6),start = c(start_year,start_qua), frequency=4)
  
  y3 <- diff(diff(GDP_qtr_3,4))
  
  baseline<-as.vector(window(y3,start=c(lowyear, start_qua_2), end=c(upyear, start_qua_2)))
  
  cpi7<-as.vector(diff(diff(GDP_qtr,4)))
  sarima_data<-cpi7[1:(length(cpi7)-1)]
  sarima_for<-as.vector(arimapred(sarima_data, n.ahead = 5))
  
  MSE_lm[i]<-mse(real[1], predict_final_lm[1])
  MSE_qu_1[i]<-mse(real[1], predict_final_qu_1[1])
  MSE_qu_2[i]<-mse(real[1], predict_final_qu_2[1])
  MSE_qu_3[i]<-mse(real[1], predict_final_qu_3[1])
  MSE_ba[i]<-mse(real[1], baseline[1])
  MSE_sa[i]<-mse(real[1],sarima_for[1])
  skew[i]<-skewness(www1)
}

mean(abs(skew))
sum(MSE_qu_1)/sum(MSE_lm) 
sum(MSE_qu_2)/sum(MSE_lm) 
sum(MSE_sa)/sum(MSE_lm)
sum(MSE_ba)/sum(MSE_lm) 

MSE_lm[i]<-mse(real[1], predict_final_lm[1])
MSE_qu_1[i]<-mse(real[1], predict_final_qu_1[1])
MSE_qu_2[i]<-mse(real[1], predict_final_qu_2[1])
MSE_qu_3[i]<-mse(real[1], predict_final_qu_3[1])
MSE_ba[i]<-mse(real[1], baseline[1])
MSE_sa[i]<-mse(real[1],sarima_for[1])




#Three quarters ahead
for (i in 1:leng) {
  q<-enddate[i]
  styear<-as.numeric(substring(q, 3,4))
  start_date<-q %m-% years(train_year)
  start_year<-as.numeric(substring(start_date, 1,4))
  start_mon<-as.numeric(substring(start_date, 6,7))
  start_qua<-(start_mon+2)/3
  
  start_date_2<-start_date %m-% months(3)
  styear_2<-as.numeric(substring(start_date_2, 3,4))
  start_year_2<-as.numeric(substring(start_date_2, 1,4))
  start_mon_2<-as.numeric(substring(start_date_2, 6,7))
  start_qua_2<-(start_mon_2+2)/3
  
  kocsv2<-kocsv%>%
    select(-CPI)%>%
    filter(X>=start_date &X<q)
  
  cpi<-kocsv%>%
    filter(X>=start_date &X<q)%>%
    select(CPI)
  
  cpi<-as.vector(cpi$CPI)
  cpi2<-cpi[seq(1, length(cpi), by=3)]
  cpi3<-replace(cpi2, length(cpi2),NA)
  
  base <- ts(kocsv2[,-1], start=c(start_year,start_mon ), frequency=12)
  trans<-rep(3,14)
  delay<-rep(0,14)
  
  GDP_qtr<-ts(cpi3,start = c(start_year,start_qua), frequency=4)
  BRGDP<-list(base=base,trans=trans,delay=delay)
  
  vintage <- PRTDB(mts = BRGDP$base, delay = BRGDP$delay, vintage = "2019-10-01")
  base <- window(vintage, start = c(start_year,start_mon), frequency = 12)
  x <- Bpanel(base = base, trans = BRGDP$trans)

  y <- diff(diff(GDP_qtr,4))
  y <- qtr2month(y)
  
  data <- cbind(y,x)
  frequency <- c(4,rep(12,ncol(x)))
  
  now11 <- nowcast(formula = y~., data = data, r = 2, q = 2 , p = 1, method = "2s",
                   frequency = frequency)
  
  fatoresTS<-now11$factors$dynamic_factors
  fatoresTRI <- month2qtr(fatoresTS)
  ppp<-diff(diff(GDP_qtr,4))
  ppp<-as.vector(ppp)
  
  www<-fatoresTRI[,1]
  www2<-www[6:(length(www)-5)]
  wwww<-fatoresTRI[,2]
  www3<-wwww[6:(length(www)-5)]
  www1<-ppp[1:(length(ppp)-1)]
  dop<-as.data.frame(cbind(www1, www2, www3))
  fit<-lm(www1~www2+www3, data = dop)
  
  final<-fit$coefficients[1]+(fit$coefficients[2])*(www)+(fit$coefficients[3])*(wwww)
  predict_final_lm<-(as.vector(final))[(length(final)-4):length(final)]
  
  #quantile reg
  qvector_1<-rep(0,5)
  for (j in 1:3) {
    rqfit<-rq(www1~www2+www3, data = dop, tau=j/4)
    final_2<-rqfit$coefficients[1]+(rqfit$coefficients[2])*(www)+(rqfit$coefficients[3])*(wwww)
    predict_final_2<-(as.vector(final_2))[(length(final_2)-4):length(final_2)]
    qvector_1<-qvector_1+predict_final_2
  }
  predict_final_qu_1<-qvector_1/3
  
  qvector_2<-rep(0,5)
  for (j in 1:4) {
    rqfit_2<-rq(www1~www2+www3, data = dop, tau=j/5)
    final_3<-rqfit_2$coefficients[1]+(rqfit_2$coefficients[2])*(www)+(rqfit_2$coefficients[3])*(wwww)
    predict_final_3<-(as.vector(final_3))[(length(final_3)-4):length(final_3)]
    qvector_2<-qvector_2+predict_final_3
  }
  predict_final_qu_2<-qvector_2/4
  
  qvector_3<-rep(0,5)
  for (j in 1:4) {
    rqfit_3<-rq(www1~www2+www3, data = dop, tau=j/5)
    final_4<-rqfit_3$coefficients[1]+(rqfit_3$coefficients[2])*(www)+(rqfit_3$coefficients[3])*(wwww)
    predict_final_4<-(as.vector(final_3))[(length(final_3)-4):length(final_3)]
    qvector_3<-qvector_3+predict_final_4*(2.5-abs(2.5-j))
  }
  predict_final_qu_3<-qvector_3/6
  
  cpi4<-kocsv%>%
    filter(X>=start_date)%>%
    select(X, CPI)
  
  cpi5<-as.vector(cpi4$CPI)
  cpi6<-cpi5[seq(1, length(cpi5), by=3)]
  
  GDP_qtr_2<-ts(cpi6,start = c(start_year,start_qua), frequency=4)
  
  y2 <- diff(diff(GDP_qtr_2,4))
  
  lowyear<-4+start_year_2
  upyear<-5+start_year_2
  real<-as.vector(window(y2, start=c(lowyear, start_qua_2), end=c(upyear, start_qua_2)))
  
  GDP_qtr_3<-ts(lag(cpi6),start = c(start_year,start_qua), frequency=4)
  
  y3 <- diff(diff(GDP_qtr_3,4))
  
  baseline<-as.vector(window(y3,start=c(lowyear, start_qua_2), end=c(upyear, start_qua_2)))
  
  cpi7<-as.vector(diff(diff(GDP_qtr,4)))
  sarima_data<-cpi7[1:(length(cpi7)-1)]
  sarima_for<-as.vector(arimapred(sarima_data, n.ahead = 5))
  
  MSE_lm[i]<-mse(real[1:3], predict_final_lm[1:3])
  MSE_qu_1[i]<-mse(real[1:3], predict_final_qu_1[1:3])
  MSE_qu_2[i]<-mse(real[1:3], predict_final_qu_2[1:3])
  MSE_qu_3[i]<-mse(real[1:3], predict_final_qu_3[1:3])
  MSE_ba[i]<-mse(real[1:3], baseline[1:3])
  MSE_sa[i]<-mse(real[1],sarima_for[1:3])
  skew[i]<-skewness(www1)  
  
  cpi77<-cpi5[seq(1, length(cpi5), by=3)]
  cpi78<-cpi77[1:18]
  GDP_qtr_23<-ts(cpi78,start = c(start_year,start_qua), frequency=4)
  y23 <- diff(diff(GDP_qtr_23,4))
  c(rep(NA, 10),predict_final_qu_2[1:3])
  y_24<-ts(c(rep(NA, 10),predict_final_qu_2[1:3]),start = c(start_year,start_qua+5), frequency=4)
  plot(y23, xlab="Time(quarter)", ylab="YOY percent change in CPI",ylim=c(-2,2), type='b')
  par(new=TRUE)
  plot(y_24, xlab="Time(quarter)", ylab="YOY percent change in CPI", ylim=c(-2,2), type='o', col='red',pch = 16)
  legend(x = "topright",          
         legend = c("TRUE", "Prediction"),  
         lty = c(1, 1),          
         col = c("Black", "Red"),      
         cex = 0.9,
         lwd = 2)  
  

}

mean(abs(skew))
sum(MSE_qu_1)/sum(MSE_lm) 
sum(MSE_qu_2)/sum(MSE_lm) 
sum(MSE_sa)/sum(MSE_lm)
sum(MSE_ba)/sum(MSE_lm) 




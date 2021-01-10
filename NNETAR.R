

### Library
library(forecast)
library(timetk)
library(sweep)
library(forecast)
library(tidyquant)
library(tidyverse)
library(LINselect)
library(tsoutliers)
library(expsmooth)
library(fma)
library(vars)

################################################################  NNETAR Parameters  ##############################################
L <-10 #short term forecast horizon
RPT <-1000 #Repeats
decay_N<-.01  #Low Decay
S_N<-10
PATHS<-100

######################################### Find Outliers in the Data Set  ######################################################

get_outliers<-function(myts_TGT) {
  

  
  
  ## Create Outliers Regressors for ARIMAX
  ## Two type of outliers Level Shift (LS) and Temprory Change (TC)
  
  n<-length(myts_TGT)
  ### AO
  resAO<- tso(y = myts_TGT, types = c("AO"), discard.method = "bottom-up", tsmethod = "auto.arima",
              args.tsmethod = list(allowdrift = FALSE, ic = "bic"))
  
  
  A<-resAO$outliers$ind
  A[is.na(A)]<-0
  
  ### TC
  
  resTC<- tso(y = myts_TGT, types = c("TC"),discard.method = "bottom-up", tsmethod = "auto.arima",
              args.tsmethod = list(allowdrift = FALSE, ic = "bic"))
  
  
  T<-resTC$outliers$ind
  T[is.na(T)]<-0
  
  ### LS
  
  resLS<- tso(y = myts_TGT, types = c("LS"),discard.method = "bottom-up", tsmethod = "auto.arima",
              args.tsmethod = list(allowdrift = FALSE, ic = "bic"))
  
  
  L<-resLS$outliers$ind
  L[is.na(L)]<-0
  
  ### Models
  
  mo.ls <- outliers("LS", L)
  ls <- outliers.effects(mo.ls, n)
  
  mo.tc <- outliers("TC", T)
  tc <- outliers.effects(mo.tc, n)
  
  mo.ao <- outliers("AO", A)
  ao <- outliers.effects(mo.ao, n)
  
  df<-data.frame(ao,tc,ls)
  
  return(df)
  
}

### Use Function
outlier<-get_outliers(lynx)

### Bind Together As Exogenous Matrix
outlier <-cbind(log(lynx),outlier)  #LOg Target Value

### Make Dataframe
df <-data.frame(outlier)

my_reg <-  lm(df$log.lynx. ~ ., data = df)
tmp <-exp(my_reg$fitted.values)  #convert back from log


### XREG Using Auto.Arima Function
fit <- auto.arima(tmp)
XREG1 <-fit$fitted  #fitted Values
XREG2 <-data.frame(forecast(fit,L))
XREG2 <-XREG2$Point.Forecast

###################################################################  AUto-Arima With Multiple XREG##################################################
### XREG Using Auto.Arima Function
fit <- auto.arima(lynx)
XREG1 <-fit$fitted  #fitted Values
XREG2 <-data.frame(forecast(fit,L))
XREG2 <-XREG2$Point.Forecast

### NNETAR Forecast
#TGT_N <-as.numeric(lynx)  #Convert Dependant Variable into Numeric

TGT_N <-lynx #Time Series as Target Variable

### Fit NNETAR
fit_NNETAR <- nnetar(TGT_N, size=S_N, repeats = RPT, xreg = XREG1,lambda = "auto", decay = decay_N)
sweep::sw_glance(fit_NNETAR)

### Take Standard Deviations from Residuals
res_sd <- sd(fit_NNETAR$residuals, na.rm=TRUE)

sims <- rnorm(L*PATHS, mean=0, sd=res_sd)  #Build Intervals

fcast_NNETAR <- forecast(fit_NNETAR, h=L, PI=TRUE,  xreg = XREG2,  npaths=PATHS, innov=sims)  #Forecast
autoplot(fcast_NNETAR$mean)

fcast_NNETAR <-data.frame(fcast_NNETAR)
fcast_NNETAR <-as.numeric(fcast_NNETAR$Point.Forecast)

acc<-accuracy(fit_NNETAR)
MAPE_Arima_MR <-acc[5]  #MAPE



###################################################################  AUto-Arima ##################################################

### XREG Using Auto.Arima Function
fit <- auto.arima(lynx)
XREG1 <-fit$fitted  #fitted Values
XREG2 <-data.frame(forecast(fit,L))
XREG2 <-XREG2$Point.Forecast

### NNETAR Forecast
#TGT_N <-as.numeric(lynx)  #Convert Dependant Variable into Numeric

TGT_N <-lynx #Time Series as Target Variable

### Fit NNETAR
fit_NNETAR <- nnetar(TGT_N, size=S_N, repeats = RPT, xreg = XREG1,lambda = "auto", decay = decay_N)
sweep::sw_glance(fit_NNETAR)

### Take Standard Deviations from Residuals
res_sd <- sd(fit_NNETAR$residuals, na.rm=TRUE)

sims <- rnorm(L*PATHS, mean=0, sd=res_sd)  #Build Intervals

fcast_NNETAR <- forecast(fit_NNETAR, h=L, PI=TRUE,  xreg = XREG2,  npaths=PATHS, innov=sims)  #Forecast
autoplot(fcast_NNETAR$mean)

fcast_NNETAR <-data.frame(fcast_NNETAR)
fcast_NNETAR <-as.numeric(fcast_NNETAR$Point.Forecast)

acc<-accuracy(fit_NNETAR)
MAPE_Arima <-acc[5]  #MAPE

###################################################################  TBATS XREG  ##################################################


### XREG Using TBATS
fit <- tbats(lynx)
XREG1 <-fit$fitted  #fitted Values
XREG2 <-data.frame(forecast(fit,L))
XREG2 <-XREG2$Point.Forecast


### NNETAR Forecast
#TGT_N <-as.numeric(lynx)  #Convert Dependant Variable into Numeric

TGT_N <-lynx #Time Series as Target Variable

### Fit NNETAR

fit_NNETAR <- nnetar(TGT_N, size=S_N, repeats = RPT, xreg = XREG1,lambda = "auto", decay = decay_N)
sweep::sw_glance(fit_NNETAR)

### Take Standard Deviations from Residuals
res_sd <- sd(fit_NNETAR$residuals, na.rm=TRUE)

sims <- rnorm(L*PATHS, mean=0, sd=res_sd)  #Build Intervals

fcast_NNETAR <- forecast(fit_NNETAR, h=L, PI=TRUE,  xreg = XREG2,  npaths=PATHS, innov=sims)  #Forecast
autoplot(fcast_NNETAR$mean)

fcast_NNETAR <-data.frame(fcast_NNETAR)
fcast_NNETAR <-as.numeric(fcast_NNETAR$Point.Forecast)

acc<-accuracy(fit_NNETAR)
MAPE_TBATS <-acc[5]  #MAPE

###############################################################    XREG Using ETS #######################################################
fit <- ets(lynx)
XREG1 <-fit$fitted  #fitted Values
XREG2 <-data.frame(forecast(fit,L))
XREG2 <-XREG2$Point.Forecast


### NNETAR Forecast
#TGT_N <-as.numeric(lynx)  #Convert Dependant Variable into Numeric

TGT_N <-lynx #Time Series as Target Variable

### Fit NNETAR

fit_NNETAR <- nnetar(TGT_N, size=S_N, repeats = RPT, xreg = XREG1,lambda = "auto", decay = decay_N)
sweep::sw_glance(fit_NNETAR)

### Take Standard Deviations from Residuals
res_sd <- sd(fit_NNETAR$residuals, na.rm=TRUE)

sims <- rnorm(L*PATHS, mean=0, sd=res_sd)  #Build Intervals

fcast_NNETAR <- forecast(fit_NNETAR, h=L, PI=TRUE,  xreg = XREG2,  npaths=PATHS, innov=sims)  #Forecast
autoplot(fcast_NNETAR$mean)

fcast_NNETAR <-data.frame(fcast_NNETAR)
fcast_NNETAR <-as.numeric(fcast_NNETAR$Point.Forecast)

acc<-accuracy(fit_NNETAR)
MAPE_ETS <-acc[5]  #MAPE


###################################################################  No XREG (Pure) ##################################################

### NNETAR Forecast
#TGT_N <-as.numeric(lynx)  #Convert Dependant Variable into Numeric

TGT_N <-lynx #Time Series as Target Variable

### Fit NNETAR

fit_NNETAR <- nnetar(TGT_N, size=S_N, repeats = RPT, lambda = "auto", decay = decay_N)
sweep::sw_glance(fit_NNETAR)

### Take Standard Deviations from Residuals
res_sd <- sd(fit_NNETAR$residuals, na.rm=TRUE)

sims <- rnorm(L*PATHS, mean=0, sd=res_sd)  #Build Intervals

fcast_NNETAR <- forecast(fit_NNETAR, h=L, PI=TRUE,  npaths=PATHS, innov=sims)  #Forecast
autoplot(fcast_NNETAR$mean)

fcast_NNETAR <-data.frame(fcast_NNETAR)
fcast_NNETAR <-as.numeric(fcast_NNETAR$Point.Forecast)

acc<-accuracy(fit_NNETAR)
MAPE_NO_XREG <-acc[5]  #MAPE

########################################################  No XREG With Log Norm Interval ##################################################

### NNETAR Forecast
#TGT_N <-as.numeric(lynx)  #Convert Dependant Variable into Numeric

TGT_N <-lynx #Time Series as Target Variable

### Fit NNETAR

fit_NNETAR <- nnetar(TGT_N, size=S_N, repeats = RPT, lambda = "auto", decay = decay_N)
sweep::sw_glance(fit_NNETAR)

### Take Standard Deviations from Residuals
res_mean <- mean(fit_NNETAR$residuals, na.rm=TRUE)

sims <- rpois(L*PATHS, res_mean)  #Build Intervals

fcast_NNETAR <- forecast(fit_NNETAR, h=L, PI=TRUE,  npaths=PATHS, innov=sims)  #Forecast
autoplot(fcast_NNETAR$mean)

fcast_NNETAR <-data.frame(fcast_NNETAR)
fcast_NNETAR <-as.numeric(fcast_NNETAR$Point.Forecast)

acc<-accuracy(fit_NNETAR)
MAPE_RPOIS <-acc[5]  #MAPE




### Results
accuracy_results <-data.frame(MAPE_NO_XREG, MAPE_TBATS, MAPE_ETS, MAPE_RPOIS,  MAPE_Arima, MAPE_Arima_MR )
colnames(accuracy_results)<-c("No_XREG", "TBATS", "ETS", "RPOIS", "Arima", "Matrix_XREG")
accuracy_results

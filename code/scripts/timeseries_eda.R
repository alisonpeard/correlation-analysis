rm(list=ls())
library(dplyr)
library(glarma)
library(gamlss)
library(lubridate)
library(ggplot2)

decompose.column <- function(df, column){
  ts <- ts(df[[column]], frequency = 12)
  ts.decomp <- stl(ts, s.window='periodic')$time.series
  colnames(ts.decomp) <- c(paste0(column, '.seasonal'), paste0(column, '.trend'), paste0(column, '.remainder'))
  return(data.frame(ts.decomp))
}
lag.column <- function(df, column, k = 1){
  df[[paste0(column, '.lag', k)]] = dplyr::lag(df[[column]], k)
  return(df)
}
cm.binary <- function(x, y){
  tp <- sum(x * y)
  fp <- sum((1 - x) * y)
  fn <- sum(x * (1 - y))
  tn <- sum((1 - x) * (1 - y))
  cm <- matrix(c(tp, fp, fn, tn), nrow=2, ncol=2)
  dimnames(cm) <- list(c("T", "F"), c("T", "F"))
  return(cm)
}

# load data
SCENARIO <- 'nf'
df.all <- read.csv(paste0("/Users/alison/Documents/RAPID/correlation-analysis/davids_results/data_results/", SCENARIO, "/full_timeseries/ts_with_levels.csv"))
df.all <- na.omit(df.all)
RZ_ID = 117 # London
df.all <- df.all[df.all$RZ_ID == RZ_ID,]

# training subset by ensemble
ENSEMBLE <- paste0(toupper(SCENARIO), '1')
df <- df.all[df.all$ensemble == ENSEMBLE,]
df$n <- lubridate::days_in_month(df$Date)
n <- nrow(df)
colnames(df)

# testing subset by ensemble
ENSEMBLE <- paste0(toupper(SCENARIO), '4')
df.test <- df.all[df.all$ensemble == ENSEMBLE,]
df.test$n <- lubridate::days_in_month(df.test$Date)

# analysis bits 
if(TRUE){
  par(mfrow=c(1,2))
  hist(df$LoS, breaks = 100, freq=FALSE)
  par(mfrow=c(1,1))
  hist(df[df$LoS > 0,]$LoS, breaks = 100, freq=FALSE)
  perc.zeros <- round(sum(as.numeric(df$LoS == 0)) / nrow(df), 6)
  print(paste0(perc.zeros * 100, "% zeros in LoS."))
  
  lambda <- 10
  mu <- 10
  size <- 1
  x <- as.numeric(1:max(df$LoS))
  y.pois <- dpois(x, lambda)
  y.nb <- dnbinom(x, mu=mu, size=size)
  y.zip <- dZIP(x, 1000, 0.98)
  lines(x, y.pois, col='red', lwd=2)
  lines(x, y.nb, col='blue', lwd=2)
  lines(x, y.zip, col='green', lwd=2)
  
} # basic EDA
if(FALSE){
  par(mfrow=c(3,1))
  plot(df$LoS, type='l')
  plot(df$LoS.binary, type='l')
  plot(df$si24, type='l')
} # plot LoS
if(FALSE){
  par(mfrow=c(3,1))
  method='min'
  plot(df$LoS, type='l')
  plot(rank(df$LoS, ties.method=method), type='l')
  plot(rank(df$si12, ties.method=method), type='l')
} # plot ranked Los
if(FALSE){
  par(mfrow=c(2,2))
  acf(df$LoS) # not stationary
  acf(diff(df$LoS)) # stationary, cuts off at lag 2
  pacf(df$LoS) # exponential decay
  pacf(diff(df$LoS)) # exponential decay, best guess MA(3)
} # ACF
if(FALSE){
  par(mfrow=c(2,2))
  method='min'
  acf(rank(df$LoS, ties.method=method))
  pacf(rank(df$LoS, ties.method=method))
  acf(diff(rank(df$LoS, ties.method=method)))
  pacf(diff(rank(df$LoS, ties.method=method)))
}# Spearman ACF
if(TRUE){
  ymin <- -.15
  ymax <- .1
  method = 'min'
  par(mfrow=c(3, 3))
  ccf(rank(df$ep_total, ties.method=method), rank(df$LoS, ties.method=method), ylim=c(ymin, ymax), main='Total EP')
  ccf(rank(df$anomaly_mean, ties.method=method), rank(df$LoS, ties.method=method), ylim=c(ymin, ymax), main='Mean anomaly')
  ccf(rank(df$anomaly_q50, ties.method=method), rank(df$LoS, ties.method=method), ylim=c(ymin, ymax), main='Q50 anomaly')
  
  ccf(rank(df$deficit_mean, ties.method=method), rank(df$LoS, ties.method=method), ylim=c(ymin, ymax), main='Mean deficit')
  ccf(rank(df$deficit_q50, ties.method=method), rank(df$LoS, ties.method=method), ylim=c(ymin, ymax), main='Q50 deficit')
  ccf(rank(df$deficit_q90, ties.method=method), rank(df$LoS, ties.method=method), ylim=c(ymin, ymax), main='Q90 deficit')
  
  ccf(rank(df$si6, ties.method=method), rank(df$LoS, ties.method=method), ylim=c(ymin, ymax), main='SI6')
  ccf(rank(df$si12, ties.method=method), rank(df$LoS, ties.method=method), ylim=c(ymin, ymax), main='SSI12')
  ccf(rank(df$si24, ties.method=method), rank(df$LoS, ties.method=method), ylim=c(ymin, ymax), main='SI24')
}# Spearman CCF
if(FALSE){
  par(mfrow=c(2, 2))
  method <- 'min'
  plot(df$ep_total, df$LoS)
  plot(rank(df$ep_total, ties.method=method), rank(df$LoS, ties.method=method))
  plot(df$si24, df$LoS)
  plot(rank(df$si24, ties.method=method), rank(df$LoS, ties.method=method))
} # Scatter plots (ranked)
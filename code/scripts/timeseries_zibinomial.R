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

if(TRUE){
  INDICATORS <- c('si6', 'si12', 'si24', 'ep_total', 'anomaly_q50')
  df.model <- na.omit(df[, c('LoS_l0', INDICATORS, 'n')])
  for(INDICATOR in INDICATORS){
    df.model <- cbind(df.model, decompose.column(df.model, INDICATOR))
    INDICATORS <- c(INDICATORS, paste0(INDICATOR, '.trend'))
    INDICATORS <- c(INDICATORS, paste0(INDICATOR, '.seasonal'))
    INDICATORS <- c(INDICATORS, paste0(INDICATOR, '.remainder'))
  }
  for(INDICATOR in INDICATORS){
    for(k in 1:24){
      df.model <- lag.column(df.model, INDICATOR, k)
      INDICATORS <- c(INDICATORS, paste0(INDICATOR, '.lag', k))
    }
  }
  # binomial MLE
  df.model$MLE <- df.model$LoS_l0 / df.model$n
  mu.hat <- mean(df.model[df.model$LoS_l0 > 0,]$MLE) # MLE for binomial
  sigma.hat <- mean(df.model$LoS_l0 > 0) # starting estimate for bernoulli
  
  # need a matrix response for this
  plot(df.model$LoS_l0, type="l")
  hist(df.model$LoS_l0)
  y <- cbind(df.model$LoS_l0, df.model$n - df.model$LoS_l0)
  df.model$y <- y
  
  df.sub <- na.omit(df.model[, c('y', INDICATORS)])
  
  sigma.INDICATORS <- c('anomaly_q50.trend')
  sigma.regressors <- paste(c(sigma.INDICATORS, paste(sigma.INDICATORS, '.lag', c(1, 3, 12, 24), sep="")), collapse=" + ")
  sigma.formula <- as.formula(paste0('~ ', sigma.regressors)); print(sigma.formula)
  
  bin.mod <- gamlss(data = df.sub,
                    formula = y ~ .,
                    sigma.formula = sigma.formula,
                    mu.start = p.hat,
                    sigma.start = sigma.hat,
                    family = "ZIBI",
                    method=RS(20))
  summary(bin.mod)
  # view fit
  mu <- bin.mod$mu.fv
  sigma <- bin.mod$sigma.fv
  bd <- df.model$n
  
  q50 <- qZIBI(0.5, bd, mu, sigma)
  lower <- qZIBI(0.05, bd, mu, sigma)
  upper <- qZIBI(0.95, bd, mu, sigma)
  
  # get dates for plotting
  index <- rownames(df.sub)
  date <- as.Date(df.all[index,]$Date)
  
  par(mfrow=c(2, 1))
  plot(date, lower, type='l', lty=3, col='red', ylim=c(0, 31))
  lines(date, upper, type='l', lty=3, col='red')
  lines(date, q50, type='l', col='blue')
  points(date, df.sub$y[,1], col='black', pch=16)
  plot(date, df.sub$y[,1], type='l')
  summary(bin.mod)
  
  #ggplot
  df.sub$q50 <- q50
  df.sub$lower <- lower
  df.sub$upper <- upper
  df.sub$date <- date
  
  # CI polygons
  Qlower <- subset(df.sub, select=c('date', 'lower'))
  Qupper <- subset(df.sub, select=c('date', 'upper'))
  names(Qlower) <- c('x', 'y'); Qlower <- Qlower[order(Qlower$x),]
  names(Qupper) <- c('x', 'y'); Qupper <- Qupper[order(Qupper$x, decreasing=TRUE),]
  conf.int <- rbind(Qlower, Qupper)
  
  ggplot(df.sub) + theme_bw() + 
    geom_polygon(data=conf.int, aes(x=x, y=y), fill='lightblue', alpha=0.5) +
    geom_line(aes(y=q50, x=date), col='blue') + 
    geom_point(aes(y=y[,1], x=date), col='black', pch=16) + 
    xlab("Year") + ylab("LoS Days") + 
    ggtitle('Zero-inflated binomial model')
    
  ## assess residuals for persistence
  #res <- bin.mod$residuals
  #par(mfrow=c(1, 3))
  #hist(res)
  #acf(res)
  #pacf(res)
} # GAMLSS zero-inflated binomial (ZIBI) model


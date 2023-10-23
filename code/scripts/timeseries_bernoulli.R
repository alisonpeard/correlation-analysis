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

# GAMLSS models
if(FALSE){
  INDICATORS <- c('si6', 'si12', 'si24', 'ep_total','anomaly_mean', 'anomaly_q50')
  df.model <- na.omit(df[, c('LoS_l0', 'RZ_ID', INDICATORS)])
  df.model.test <- na.omit(df.test[, c('LoS_l0', 'RZ_ID', INDICATORS)])
  df.model$LoS.binary <- as.numeric(df.model$LoS_l0 > 0)
  df.model.test$LoS.binary <- as.numeric(df.model.test$LoS_l0 > 0)
  p.hat <- mean(df.model$LoS.binary)
  
  if(TRUE){
    for(INDICATOR in INDICATORS){
      df.model <- cbind(df.model, decompose.column(df.model, INDICATOR))
      df.model.test <- cbind(df.model.test, decompose.column(df.model.test, INDICATOR))
      INDICATORS <- c(INDICATORS, paste0(INDICATOR, '.trend'))
      INDICATORS <- c(INDICATORS, paste0(INDICATOR, '.seasonal'))
      INDICATORS <- c(INDICATORS, paste0(INDICATOR, '.remainder'))
    }
  }
  
  par(mfrow=c(1,1))
  plot(df.model$si6, type='l')
  for(INDICATOR in INDICATORS){
    for(k in 1:24){
      df.model <- lag.column(df.model, INDICATOR, k)
      df.model.test <- lag.column(df.model.test, INDICATOR, k)
      INDICATORS <- c(INDICATORS, paste0(INDICATOR, '.lag', k))
      lines(df.model[[paste0('si6.lag', k)]], col='blue')
    }
  }
  
  y <- cbind(df.model$LoS.binary, 1 - df.model$LoS.binary)
  y.test <- cbind(df.model.test$LoS.binary, 1 - df.model.test$LoS.binary)
  df.model$y <- y
  df.model.test$y <- y.test
  df.sub <- na.omit(df.model[, c('y', INDICATORS)])
  df.sub.test <- na.omit(df.model.test[, c('y', INDICATORS)])
  
  if(FALSE){
    par(mfrow=c(2,2))
    hist(df.model$LoS_l0, main='London NF LoS days for single ensemble', xlab='LoS days')
    hist(df.model$LoS.binary, breaks=2, xaxt='n', main='London NF LoS occurence for single ensemble', xlab='Yes or no')
    axis(side = 1, at = c(0, 0.25, 0.75, 1), labels = c("", "No", "Yes", ""))
    plot(df.model$LoS_l0, xlab='Months', ylab='Number of days', type='l', main='London NF LoS days for single ensemble')
    plot(df.model$LoS.binary,xlab='Months', ylab='Yes or no', type='l', main='London NF LoS occurence for single ensemble')
    
    par(mfrow=c(2,1))
    plot(df.sub$y[,1], type='l', main='train')
    plot(df.sub.test$y[,1], type='l', main='test')
  }
  
  index <- rownames(df.sub)
  date <- as.Date(df.all[index,]$Date)
  
  if(FALSE){
    look.at <- 'ep_total'
    par(mfrow=c(3, 1))
    plot(date, df.sub[[paste0(look.at, '.trend')]], type='l', xlab='date', ylab='trend', main='Q50 anomaly decomposed')
    plot(date, df.sub[[paste0(look.at, '.seasonal')]], type='l', xlab='date', ylab='seasonal')
    plot(date, df.sub[[paste0(look.at, '.remainder')]], type='l', xlab='date', ylab='remainder')
  }
  if(FALSE){
    look.at <- 'anomaly_q50.trend'
    par(mfrow=c(4,1))
    ccf(rank(df.sub[[look.at]]), rank(df.sub$y[,1]), main=paste0('CCF x=rank(', look.at, '), y=rank(binary LoS)'))
    plot(date, df.sub$y[,1], type='l', ylab='LoS binary')
    plot(date, df.sub[[look.at]], type='l', ylab=look.at)
    plot(df.sub[[look.at]], df.sub$y[,1], pch=16, xlab=look.at, ylab='LoS binary', main='scatter plot')
  }
  
  INDICATORS <- c('si6.trend', 'anomaly_q50.trend')
  regressors <- paste(c(INDICATORS, paste(INDICATORS, '.lag', c(1, 3, 12, 24), sep="")), collapse=" + ")
  formula <- as.formula(paste0('y ~ ', regressors)); print(formula)
  bin.mod <- gamlss(formula,
                    data = df.sub,
                    mu.start = p.hat,
                    family = "BI",
                    method=RS(20))
  summary(bin.mod)
  
  mu <- bin.mod$mu.fv
  bd <- 1
  q50 <- qBI(0.5, bd, mu)
  lower <- qBI(0.025, bd, mu)
  upper <- qBI(0.975, bd, mu)
  
  df.sub$q50 <- q50
  df.sub$lower <- lower
  df.sub$upper <- upper
  df.sub$date <- date
  
  Qlower <- subset(df.sub, select=c('date', 'lower'))
  Qupper <- subset(df.sub, select=c('date', 'upper'))
  names(Qlower) <- c('x', 'y'); Qlower <- Qlower[order(Qlower$x),]
  names(Qupper) <- c('x', 'y'); Qupper <- Qupper[order(Qupper$x, decreasing=TRUE),]
  conf.int <- rbind(Qlower, Qupper)
  
  ggplot(df.sub) + theme_bw() + 
    geom_polygon(data=conf.int, aes(x=x, y=y), fill='lightblue', alpha=0.5) +
    geom_line(aes(y=q50, x=date), col='blue') + 
    geom_point(aes(y=y[,1], x=date), col='black', pch=16) + 
    xlab("Year") + ylab("LoS occurs (yes/no)") + 
    ggtitle('Binary binomial model for restrictions')
  
  # test on test data
  mu.pred <- predict(bin.mod, newdata=df.sub.test, what="mu", type='response')
  index <- rownames(df.sub.test)
  date <- as.Date(df.all[index,]$Date)
  
  df.sub.test$q50 <- qBI(0.5, bd, mu.pred)
  df.sub.test$lower <- qBI(0.025, bd, mu.pred)
  df.sub.test$upper <- qBI(0.975, bd, mu.pred)
  df.sub.test$date <- date
  
  Qlower <- subset(df.sub.test, select=c('date', 'lower'))
  Qupper <- subset(df.sub.test, select=c('date', 'upper'))
  names(Qlower) <- c('x', 'y'); Qlower <- Qlower[order(Qlower$x),]
  names(Qupper) <- c('x', 'y'); Qupper <- Qupper[order(Qupper$x, decreasing=TRUE),]
  conf.int <- rbind(Qlower, Qupper)
  
  ggplot(df.sub.test) + theme_bw() + 
    geom_polygon(data=conf.int, aes(x=x, y=y), fill='lightblue', alpha=0.5) +
    geom_line(aes(y=q50, x=date), col='blue') + 
    geom_point(aes(y=y[,1], x=date), col='black', pch=16) + 
    xlab("Year") + ylab("LoS occurs (yes/no)") + 
    ggtitle('Binary binomial predictions on test set')
  
  # Evaluation criteria
  y.true <- df.sub.test$y[,1]
  y.pred <- df.sub.test$q50
  p <- cbind(mu.pred, 1 - mu.pred)
  n <- length(mu.pred)
  score.bce <-  - sum(df.sub.test$y * log(p)) / n # want at least less than 0.2
  score.cm <- cm.binary(y.true, y.pred)
  score.accuracy <- (score.cm[1,1] + score.cm[2,2]) / sum(score.cm)
  score.precision <- score.cm[1,1] / (score.cm[1,1] + score.cm[2,1])
  score.recall <- score.cm[1,1] / (score.cm[1,1] + score.cm[1,2])
  score.f1_score <- 2 * score.precision * score.recall / (score.precision + score.recall)
  
  print(paste0('Binary cross entropy: ', round(score.bce, 4)))
  print(paste0('Precision: ', round(score.precision, 6)))
  print(paste0('Recall: ', round(score.recall, 6)))
  print(paste0('F1-score: ', round(score.f1_score, 6)))
  
} # GAMLSS to predict binary - use anomaly_q50.trend lags c(1,3,12,24)!
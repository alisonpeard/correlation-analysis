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
ENSEMBLE <- paste0(toupper(SCENARIO), '2')
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

# GAMLSS models
if(TRUE){
  INDICATORS <- c('si6', 'si12', 'si24', 'anomaly_mean', 'anomaly_q50')
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
      lines(df.model[[paste0('si6.lag', k)]], col='blue', alpha=.5)
    }
  }
  
  y <- cbind(df.model$LoS.binary, 1 - df.model$LoS.binary)
  y.test <- cbind(df.model.test$LoS.binary, 1 - df.model.test$LoS.binary)
  df.model$y <- y
  df.model.test$y <- y.test
  df.sub <- na.omit(df.model[, c('y', INDICATORS)])
  df.sub.test <- na.omit(df.model.test[, c('y', INDICATORS)])
  
  par(mfrow=c(2,1))
  plot(df.sub$y[,1], type='l', main='train')
  plot(df.sub.test$y[,1], type='l', main='test')
  
  par(mfrow=c(4,1))
  look.at <- 'anomaly_q50.trend'
  ccf(rank(df.sub[[look.at]]), rank(df.sub$y[,1]))
  plot(df.sub$y[,1], type='l')
  plot(df.sub[[look.at]])
  plot(df.sub[[look.at]], df.sub$y[,1], pch=16)
  
  INDICATORS <- c('anomaly_q50.trend')
  regressors <- paste(c(INDICATORS, paste(INDICATORS, '.lag', c(1, 3, 12, 24), sep="")), collapse=" + ")
  formula <- as.formula(paste0('y ~ ', regressors))
  bin.mod <- gamlss(formula,
                    data = df.sub,
                    mu.start = p.hat,
                    family = "BI",
                    method=RS(20))
  summary(bin.mod)
  
  mu <- bin.mod$mu.fv
  bd <- 1
  q50 <- qBI(0.5, bd, mu)
  lower <- qBI(0.05, bd, mu)
  upper <- qBI(0.95, bd, mu)
  
  index <- rownames(df.sub)
  date <- as.Date(df.all[index,]$Date)
  
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
  df.sub.test$lower <- qBI(0.05, bd, mu.pred)
  df.sub.test$upper <- qBI(0.95, bd, mu.pred)
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
} # GAMLSS to predict binary - getting there!





# OTHER MODELS
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
    for(k in 1:4){
      df.model <- lag.column(df.model, INDICATOR, k)
      INDICATORS <- c(INDICATORS, paste0(INDICATOR, '.lag', k))
    }
  }
  # binomial MLE
  df.model$MLE <- df.model$LoS_l0 / df.model$n
  df.model <- df.model[df.model$LoS_l0 > 0,] # take only positives to avoid ZIBI
  p.hat <- mean(df.model$MLE) # MLE for binomial
  
  # need a matrix response for this
  par(mfrow=c(2, 1))
  plot(df.model$LoS_l0, type="l", ylab='LoS days', xlab='Month')
  hist(df.model$LoS_l0, xlab="Days", main='No LoS days per month')
  
  y <- cbind(df.model$LoS_l0, df.model$n - df.model$LoS_l0)
  df.model$y <- y
  df.sub <- df.model[, c('y', INDICATORS)]
  
  bin.mod <- gamlss(data = df.sub,
                    formula = y ~ .,
                    mu.start = p.hat,
                    family = "BI",
                    method=RS(20))
  
  # view fit
  mu <- bin.mod$mu.fv
  bd <- df.model$n
  q50 <- qBI(0.5, bd, mu)
  lower <- qBI(0.05, bd, mu)
  upper <- qBI(0.95, bd, mu)
  
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
    ggtitle('Binomial model for positive counts')
  
  #par(mfrow=c(2, 1))
  #plot(q10, type='l', lty=3, col='red')
  #lines(q90, type='l', lty=3, col='red')
  #lines(q50, type='l', col='blue')
  #points(df.model$LoS_l0, col='black', pch=16)
  #plot(df.model$LoS_l0, col='black', type='l')
  
  #lines(df$LoS_l0, col='black')
  #plot(df.model$y[,1], type='l')
  #lines(fitted.mu * df.model$n, type='l', col='blue')
  #summary(bin.mod)

  ## assess residuals for persistence
  #res <- bin.mod$residuals
  #par(mfrow=c(1, 3))
  #hist(res)
  #acf(res)
  #pacf(res)
  
  ## worm plot
  #wp(bin.mod)
} # GAMLSS binomial model
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
    for(k in 1:4){
      df.model <- lag.column(df.model, INDICATOR, k)
      INDICATORS <- c(INDICATORS, paste0(INDICATOR, '.lag', k))
    }
  }
  # binomial MLE
  df.model$MLE <- df.model$LoS_l0 / df.model$n
  p.hat <- mean(df.model[df.model$LoS_l0 > 0,]$MLE) # MLE for binomial
  
  perc.zeros <- mean(df.model$LoS == 0)
  
  # need a matrix response for this
  plot(df.model$LoS_l0, type="l")
  hist(df.model$LoS_l0)
  y <- cbind(df.model$LoS_l0, df.model$n - df.model$LoS_l0)
  df.model$y <- y
  
  df.sub <- df.model[, c('y', INDICATORS)]
  bin.mod <- gamlss(data = df.sub,
                    formula = y ~ .,
                    sigma.formula = ~.,
                    mu.start = p.hat,
                    sigma.start = perc.zeros,
                    family = "ZIBI",
                    method=RS(20))
  
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
# OLD STUFF
if(FALSE){
  if(TRUE){
    # fit a Poisson to the > 0 values
    # fit a ZIP to both with sigma.start=%zeros
    # look at ACF, PACF, and CCF to determine GLARMA lags
    
    # fit the model
    INDICATORS <- c('si6', 'si12', 'si24', 'ep_total', 'anomaly_q50')
    df.model <- na.omit(df[, c('LoS', INDICATORS)])
    for(INDICATOR in INDICATORS){
      df.model <- cbind(df.model, decompose.column(df.model, INDICATOR))
    }
    for(INDICATOR in INDICATORS){
      for(k in 1:4){
        df.model <- lag.column(df.model, INDICATOR, k)
      }
    }
    gam.mod <- gamlss(data = df.model,
                      formula = LoS ~ .,
                      sigma.start = perc.zeros,
                      sigma.fix = TRUE,
                      family = "ZIP",
                      method=RS(5))
    
    # view fit
    fitted.mu <- gam.mod$mu.fv
    fitted.sigma <- gam.mod$sigma.fv
    head(fitted.mu)
    head(fitted.sigma)
    
    par(mfrow=c(2, 1))
    plot(fitted.mu, type='l', col='blue')
    lines(df$LoS, col='black')
    plot(df$LoS, type='l')
    summary(gam.mod)
    
    par(mfrow=c(1, 1))
    plot(df[,INDICATORS])
    
    # assess residuals for persistence
    res <- gam.mod$residuals
    par(mfrow=c(1, 3))
    hist(res)
    acf(res)
    pacf(res)
  } # GAMLSS Poi/NB model (no lags)
}

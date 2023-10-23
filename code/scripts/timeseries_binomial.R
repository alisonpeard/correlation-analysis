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
print(paste0("Number of LoS days in train: ", sum(df$LoS_l0 > 0)))

# testing subset by ensemble]
if(TRUE){
  ENSEMBLE <- paste0(toupper(SCENARIO), '83')
  df.test <- df.all[df.all$ensemble == ENSEMBLE,]
  df.test$n <- lubridate::days_in_month(df.test$Date)
  print(paste0("Number of LoS days in test: ", sum(df.test$LoS_l0 > 0)))
}


if(TRUE){
  if(TRUE){
    INDICATORS <- c('si6', 'si12', 'si24', 'ep_total', 'anomaly_q50')
    df.model <- na.omit(df[, c('LoS_l0', INDICATORS, 'n')])
    df.model.test <- na.omit(df.test[, c('LoS_l0', INDICATORS, 'n')])
    
    # calculate lags and decomposition BEFORE removing zeros
    for(INDICATOR in INDICATORS){
      df.model <- cbind(df.model, decompose.column(df.model, INDICATOR))
      df.model.test <- cbind(df.model.test, decompose.column(df.model.test, INDICATOR))
      INDICATORS <- c(INDICATORS, paste0(INDICATOR, '.trend'))
      INDICATORS <- c(INDICATORS, paste0(INDICATOR, '.seasonal'))
      INDICATORS <- c(INDICATORS, paste0(INDICATOR, '.remainder'))
    }
    for(INDICATOR in INDICATORS){
      for(k in 1:24){
        df.model <- lag.column(df.model, INDICATOR, k)
        df.model.test <- lag.column(df.model.test, INDICATOR, k)
        INDICATORS <- c(INDICATORS, paste0(INDICATOR, '.lag', k))
      }
    }
    
    # binomial MLE
    df.model$MLE <- df.model$LoS_l0 / df.model$n
    df.model <- df.model[df.model$LoS_l0 > 0,] # take only positives to avoid ZIBI
    df.model.test <- df.model.test[df.model.test$LoS_l0 > 0,]
    p.hat <- mean(df.model$MLE) # MLE for binomial
    
    # need a matrix response for this
    par(mfrow=c(2, 1))
    plot(df.model$LoS_l0, type="l", ylab='LoS days', xlab='Month')
    hist(df.model$LoS_l0, xlab="Days", main='No LoS days per month')
    
    y <- cbind(df.model$LoS_l0, df.model$n - df.model$LoS_l0)
    y.test <- cbind(df.model.test$LoS_l0, df.model.test$n - df.model.test$LoS_l0)
    df.model$y <- y
    df.model.test$y <- y.test
    df.sub <- na.omit(df.model[, c('y', INDICATORS, 'n')])
    df.sub.test <- na.omit(df.model.test[, c('y', INDICATORS, 'n')])
    
    INDICATORS <- c('anomaly_q50.trend')
    cols <- c(INDICATORS, paste(INDICATORS, '.lag', c(1, 3, 12, 24), sep=""))
    regressors <- paste(cols, collapse=" + ")
    formula <- as.formula(paste0('y ~ ', regressors)); print(formula)
    bin.mod <- gamlss(data = df.sub,
                      formula = formula,
                      mu.start = p.hat,
                      family = "BI",
                      method=RS(20))
    summary(bin.mod)
  }
  if(TRUE){
    # view fit
    mu <- bin.mod$mu.fv
    bd <- df.model$n
    q50 <- qBI(0.5, bd, mu)
    lower <- qBI(0.05, bd, mu)
    upper <- qBI(0.95, bd, mu)
    
    # get dates for plotting
    index <- rownames(df.sub)
    date <- as.Date(df.all[index,]$Date)
    all.dates <- data.frame(date=seq(from=min(date), to=max(date), by="month"))
    
    #ggplot
    df.sub$q50 <- q50
    df.sub$lower <- lower
    df.sub$upper <- upper
    df.sub$date <- date
    
    # NAs for missing dates
    df.plot <- data.frame(date = date, q50=q50, lower=lower, upper=upper, y=y[,1])
    df.plot <- merge(all.dates, df.plot, by = "date", all.x = TRUE)
    
    # CI polygons
    par(mfrow=c(3,1))
    Qlower <- subset(df.plot, select=c('date', 'lower'))
    Qupper <- subset(df.plot, select=c('date', 'upper'))
    names(Qlower) <- c('x', 'y'); Qlower <- Qlower[order(Qlower$x),]
    names(Qupper) <- c('x', 'y'); Qupper <- Qupper[order(Qupper$x, decreasing=TRUE),]
    rle.lower <- rle(!is.na(Qlower$y))
    rle.upper <- rle(!is.na(Qupper$y))
    Qlower$group <- as.factor(rep(seq_along(rle.lower$values), rle.lower$lengths))
    Qupper$group <- as.factor(rep(rev(seq_along(rev(rle.upper$values))), rle.upper$lengths)) # this line is tricky
    conf.int <- na.omit(rbind(Qlower, Qupper))
    conf.int <- conf.int[order(conf.int$group),]
    
    ggplot(df.plot) + theme_bw() + 
      geom_polygon(data=conf.int, aes(x=x, y=y, group=group), fill='lightblue', alpha=0.5) +
      geom_path(data=conf.int, aes(x=x, y=y, group=group), col='lightblue', alpha=0.5) +
      geom_line(aes(y=q50, x=date), col='blue') +
      geom_point(aes(y=y, x=date), col='black', pch=16) +
      xlab("Year") + ylab("LoS Days") + 
      ggtitle('Binomial model for positive counts')
  }
  if(TRUE){
    # predictions
    mu.pred <- predict(bin.mod, newdata=df.sub.test, what="mu", type='response')
    bd <- df.sub.test$n
    q50 <- qBI(0.5, bd, mu.pred)
    lower <- qBI(0.05, bd, mu.pred)
    upper <- qBI(0.95, bd, mu.pred)
    
    index <- rownames(df.sub.test)
    date <- as.Date(df.all[index,]$Date)
    all.dates <- data.frame(date=seq(from=min(date), to=max(date), by="month"))
    
    df.sub.test$q50 <- q50
    df.sub.test$lower <- lower
    df.sub.test$upper <- upper
    df.sub.test$date <- date
    
    # NAs for missing dates
    df.plot.test <- data.frame(date = date, q50=q50, lower=lower, upper=upper, y=df.sub.test$y[,1])
    df.plot.test <- merge(all.dates, df.plot.test, by = "date", all.x = TRUE)
    
    # CI polygons
    par(mfrow=c(3,1))
    Qlower <- subset(df.plot.test, select=c('date', 'lower'))
    Qupper <- subset(df.plot.test, select=c('date', 'upper'))
    names(Qlower) <- c('x', 'y'); Qlower <- Qlower[order(Qlower$x),]
    names(Qupper) <- c('x', 'y'); Qupper <- Qupper[order(Qupper$x, decreasing=TRUE),]
    rle.lower <- rle(!is.na(Qlower$y))
    rle.upper <- rle(!is.na(Qupper$y))
    Qlower$group <- as.factor(rep(seq_along(rle.lower$values), rle.lower$lengths))
    Qupper$group <- as.factor(rep(rev(seq_along(rev(rle.upper$values))), rle.upper$lengths)) # this line is tricky
    conf.int <- na.omit(rbind(Qlower, Qupper))
    conf.int <- conf.int[order(conf.int$group),]
    
    ggplot(df.plot.test) + theme_bw() + 
      geom_polygon(data=conf.int, aes(x=x, y=y, group=group), fill='lightblue', alpha=0.5) +
      geom_path(data=conf.int, aes(x=x, y=y, group=group), col='lightblue', alpha=0.5) +
      geom_line(aes(y=q50, x=date), col='blue') +
      geom_point(aes(y=q50, x=date), col='blue', pch=20, cex=.2) +
      geom_point(aes(y=y, x=date), col='black', pch=20) +
      xlab("Year") + ylab("LoS Days") + 
      ggtitle('Binomial predictions on test set')
    
    score.rmse <- sqrt(mean((df.sub.test$y[,1] - df.sub.test$q50)^2))
    print(paste0("RMSE: ", round(score.rmse, 6)))
  }

  

} # GAMLSS binomial model, 6)

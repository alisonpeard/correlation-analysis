library(zoo)
library(forecast) # tsdisplay, stl
library(arrow) # parquet reading
library(tseries) # adf.test
library(lubridate)
library(glarma)

indicator <- "ep_total"
buffer <- 25

# load data
df1 <- read_parquet("/Users/alison/Documents/RAPID/correlation-analysis/data_results/nf/full_timeseries/wrz_1/NF4.parquet")
df2 <- read_parquet("/Users/alison/Documents/RAPID/correlation-analysis/data_results/nf/full_timeseries/wrz_1/NF2.parquet")
df3 <- read_parquet("/Users/alison/Documents/RAPID/correlation-analysis/data_results/nf/full_timeseries/wrz_1/NF3.parquet")


# process dataframe: concatenate the ensembles to help fit
format.date <- function(date_string) {
  date <- as.Date(date_string)
  year <- as.integer(format(date, "%Y"))
  month <- as.integer(format(date, "%m"))
  day <- as.integer(format(date, "%d"))
  return(c(year, month, day))
}
make_ts <- function(df,indicator, buffer, delta=c(0, 0, 0)){
  # process df
  df <- df[(df$buffer==buffer),]
  df$date <- make_date(year=df$Year, month=df$Month)
  all.dates <- seq.Date(min(df$date), max(df$date), "month")
  df <- merge(x=data.frame(date=all.dates), y=df, all.x=TRUE)
  
  # get key dates
  start.date <- format.date(df[1, "date"]) + delta
  end.date <- format.date(df[dim(df)[1], "date"]) + delta
  
  df$date <- df$date + years(delta[1]) + months(delta[2]) + days(delta[3])
  # fill in any NaNs
  zoo.values <- zoo(df[, c(indicator, "LoS")], df$date)
  approx.values <- na.approx(zoo.values)
  
  return(list(ts=approx.values, start=start.date, end=end.date))
}
concat_ts <- function(dfs, indicator, buffer){
  delta = c(0, 0, 0)
  res <- make_ts(dfs[[1]],indicator, buffer)
  start <- res$start
  end <- res$end
  ts <- res$ts
  delta <- end - start + c(0, 0, 1)
  
  for(df in dfs[2:length(dfs)]){
    res <- make_ts(df, indicator, buffer, delta)
    ts <- rbind(ts, res$ts)
    end <- res$end
    delta <- end - start + c(0, 0, 1)
  }
  
  ts <- ts(ts, start=start, end=end, frequency=12)
  return(list(ts=ts, start=start, end=end))
}
res <- concat_ts(dfs=list(df1, df2, df3), indicator=indicator, buffer=buffer)
ts <- res$ts

tsdisplay(ts[,indicator]) # lags don't die out but oscillate => trend/seasonality
tsdisplay(ts[,'LoS']) # lags die out so ok stationary

adf.test(ts[,indicator]) # p-value < 0.01, alternative hypothesis: stationary
adf.test(ts[,"LoS"]) # p-value < 0.01, alternative hypothesis: stationary, only for residuals

# cumulative periodograms => trend and periodicity
par(mfrow=c(2, 2))
cpgram(ts[,"LoS"])
cpgram(diff(ts[,"LoS"]))
cpgram(ts[,indicator])
cpgram(diff(diff(ts[,indicator], 12)))

diffs <- function(x, diffs=c(1), start=res$start, end=res$end){
  for(d in diffs){
    x <- diff(x, d)
  }
  return(x)
}
d1 <- diffs(ts[,indicator], c(1, 12))
ts.d <- ts[1:length(d1),]
ts.d[, indicator] <- d1
ts.d[, "LoS"] <- ts[1:length(d1), "LoS"]
ts.d <- ts(ts.d, start=res$start, end=res$end, frequency=12)

# Autopersistence functions Startz (2008)
autopersistence <- function(x, maxlags=1){
  x <- as.numeric(x > 0) # make sure its binary
  apf.0 <- numeric(maxlags)
  apf.1 <- numeric(maxlags)
  n <- length(x)
  for(i in 1:maxlags){
    xk <- x[(1+i):n]
    x0 <- x[1:(n-i)]
    
    # conditional probabilities
    one.zero <- sum((xk - x0) == 1)  # (1, 0)-pairs
    zero.zero <- sum((2 - (xk + x0) == 2)) # (0, 0)-pairs
    one.one <- sum((2 - (xk + x0) == 0)) # (1, 1)-pairs
    zero.one <- sum((xk - x0) == -1) # (0, 1)-pairs
    
    apf.0[i] <- one.zero / (one.zero + zero.zero)
    apf.1[i] <- one.one / (one.one + zero.one)
  }
  return(list(apf0=apf.0, apf1=apf.1))
}
apfs <- autopersistence(x, maxlags=36)
par(mfrow=c(1, 1))
plot(1:36, apfs$apf0, type="l", col='blue', ylim=c(.3, .7), xlab="lags", ylab="probability")
lines(1:36, apfs$apf1, col='green')
title("Empirical Autopersistence Graphs")
legend(16, 0.7, legend=c("APG0 P(1|0)", "APG1 P(1|1)"),
       col=c("blue", "green"), lty=c(1, 1), cex=0.8)


# ACDF and PACFs
par(mfrow=c(2,2))
acf(ts.d[,indicator])
acf(ts.d[,"LoS"])
pacf(ts.d[,indicator])
pacf(ts.d[,"LoS"])
# What do these suggest? That both look like moving average processes.

par(mfrow=c(1, 1))
plot(decompose(ts.d[,indicator]))
decompose(ts[,indicator])$seasonal

par(mfrow=c(2, 1))
plot(ts[,"LoS"])
plot(ts.d[,indicator])
y <- as.integer(ts[,"LoS"])
X <- as.numeric(ts.d[,indicator])
X <- cbind(intersect=rep(1, length(X)), variable=X)
y.bin <- as.numeric(y>0)
y.bin <- matrix(nrow=length(y), ncol=2)
y.bin[,1] <- as.numeric(y>0)
y.bin[,2] <- as.numeric(y==0)

barma.mod <- glarmaBinomialIdentity(y.bin, X, phiLags=c(1), thetaLags=c(1), delta=c(1,1,1))

poisson.mod <- glarma(y, X, type="Poi", thetaLags=c(12), method="NR")
plot(poisson.mod)
summary(poisson.mod)

plot(approx.values[,indicator1])
#tsdisplay(approx.values$LoS, lag=12)
#tsdisplay(approx.values[, indicator1], lag=12)
#tsdisplay(approx.values[, indicator2], lag=12)

# check stationarity



# individual timeseries

plot(ts)

ts <- cbind(ts, event=ts(as.numeric(ts[,'LoS'] > 0), start=start.date, end=end.date, frequency=12))
colnames(ts) <- c("LoS", indicator1, indicator2, 'event')

plot(ts)

# create differenced vectors
diff.ts <- function(x, start=start.date, end=end.date, frequency=12, d=1){
  x <- c(numeric(d), x[1: length(x)-d])
  x <- ts(x, start=start, end=end, frequency=frequency)
  print(head(x))
  print(length(x))
  return(x)
}

diffs.ts <- function(x, start=start.date, end=end.date, frequency=12, ds=c(1, 1, 12)){
  for(d in ds){
    print(length(x))
    x <- diff.ts(x, start=start, end=end, frequency=frequency, d=d)
    print(length(x))
  }
  return(x)
}

y <- diffs.ts(ts[,indicator1], c(1, 1))

ts <- cbind(ts, diff.ts(diff.ts(df[,indicator1])), diff.ts(diff.ts(df[,indicator1])))
colnames(ts) <- c("LoS", indicator1, indicator2, 'event', paste(indicator1, ".d"), paste(indicator2, ".d"))
plot(ts)

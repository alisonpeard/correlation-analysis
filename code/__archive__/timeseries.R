library(zoo)
library(arrow)
library(stats)
library(forecast)
library(lubridate)
library(TSstudio)


# load data
indicator <- "q50_deficit_total"
df <- read_parquet("/Users/alison/Documents/RAPID/RAPID Input/processed_data/hist/full_timeseries/wah_36.parquet")
df <- df[(df$buffer==25),]
df$date <- make_date(year=df$Year, month=df$Month)

all.dates <- seq.Date(min(df$date), max(df$date), "month")
df <- merge(x=data.frame(date=all.dates), y=df, all.x=TRUE)

zoo.values <- zoo(df[, c(indicator, "LoS")], df$date)
approx.values <- na.approx(zoo.values)
plot(approx.values[,indicator])

tsdisplay(approx.values$LoS, lag=12)
tsdisplay(approx.values[, indicator], lag=12)

indicator.ts <- ts(approx.values[, indicator], start=c(1962,2), end=c(2015,12), frequency=12)
los.ts <- ts(approx.values$LoS, start=c(1962,2), end=c(2015,12), frequency=12)
# some plots
plot(stl(indicator.ts, s.window="periodic"))
plot(stl(los.ts, s.window="periodic"))

# split into train and validation sets
indicator.tr = window(indicator.ts, start=c(1962,2), end=c(2012,12))
indicator.valid = window(indicator.ts, start=c(2012,12),end=c(2015,12))
los.tr = window(los.ts, start=c(1962,2), end=c(2012,12))
los.valid = window(los.ts, start=c(2012,12),end=c(2015,12))

# combine
train <- cbind(indicator=indicator.tr, los=los.tr)
valid <- cbind(indicator=indicator.valid, los=los.valid)

# Vector Autoregression
library(vars)

VARselect(train, lag.max=12, type="const")[["selection"]]
# https://www.r-bloggers.com/2020/02/testing-for-a-causal-effect-with-2-time-series/
# https://otexts.com/fpp2/VAR.html
var = VAR(train, lag.max=12, type="const")
serial.test(var, lags.pt=12, type="PT.asymptotic") # h0: residuals white noise
coefficients(var)
causality(var, cause="indicator")

library(ggplot2)
forecast(var) %>%
  autoplot() + xlab("Year")

# compare
predictions <- predict(var, n.ahead=12)
los.preds <- ts(predictions$fcst$los[,1], start=c(2012, 12), end=c(2015,12), frequency=12)
plot(los.preds)
lines(los.valid, col='red')

# trend check
Phi = matrix(c(coefficients(var)$indicator[1:2,1],coefficients(var)$los[1:2,1]),2,2)
eigen(Phi)
plot(train)

# try with differencing
var.diff = VAR(diff(train), lag.max=12, type = "const")
coefficients(var.diff)
causality(var.diff, cause="indicator")

# OLD STUFF
#####
indicator.ts <- ts(as.vector(df[,indicator]), start=c(1962,2), end=c(1987,12), frequency=12)


los.ts <- ts(as.vector(df$LoS), start=c(1962, 2), end=c(1991, 12), frequency=12)
anomaly.ts <- ts(as.vector(df$q50_anomaly), start=c(1962, 2), end=c(1991, 12), frequency=12)

# have a look
plot(rain.ts)
plot(stl(rain.ts, s.window="periodic"))
tsdisplay(rain.ts)
tsdisplay(diff(rain.ts))  # ACF and PACF to choose p & q
tsdisplay(diff(rain.ts, lag=12))

plot(los.ts)
plot(stl(los.ts, s.window="periodic"))
tsdisplay(los.ts)

plot(anomaly.ts)
plot(stl(anomaly.ts, s.window="periodic"))
tsdisplay(anomaly.ts)

par(mfrow=c(3, 1))
plot(rain.ts)
plot(los.ts)
plot(anomaly.ts)


# split into train and validation sets
rain.tr = window(rain.ts, start=c(1962,2), end=c(1990,12))
rain.valid = window(rain.ts, start=c(1990,12),end=c(1991,12))

# fit a model and do some basic predicting
fit <- arima(rain.tr, order=c(0,1,8), seasonal=list(order=c(11,0,11), period=12))


tsdiag(fit)
# residuals look normal(ish) but a little seasonal
# ACF of residuals shows no autocorrelations remaining
# no significant p-values in LB plot mean residuals are random
cpgram(residuals(fit)) # looks good


rain.pred = predict(fit, n.ahead=12)

plot(rain.ts)
lines(rain.pred$pred, col="blue", lwd=2)
lines(rain.pred$pred + 2*rain.pred$se, col="blue", lwd=1, lty="dashed")
lines(rain.pred$pred - 2*rain.pred$se, col="blue", lwd=1, lty="dashed")

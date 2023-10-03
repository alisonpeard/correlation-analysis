library(zoo)
library(arrow)
library(stats)
library(forecast)
library(lubridate)
library(TSstudio)
library(dynlm)

# load data
indicator <- "ep_total"
buffer <- 25

df <- read_parquet("/Users/alison/Documents/RAPID/correlation-analysis/data_results/nf/full_timeseries/wrz_1/NF1.parquet")
df <- df[(df$buffer==buffer),]
df$date <- make_date(year=df$Year, month=df$Month)

all.dates <- seq.Date(min(df$date), max(df$date), "month")
df <- merge(x=data.frame(date=all.dates), y=df, all.x=TRUE)

zoo.values <- zoo(df[, c(indicator1, indicator2, "LoS")], df$date)
approx.values <- na.approx(zoo.values)
plot(approx.values[,indicator])

tsdisplay(approx.values$LoS, lag=12)
tsdisplay(approx.values[, indicator1], lag=12)
tsdisplay(approx.values[, indicator2], lag=12)

#ts <- ts(approx.values[, c('LoS', indicator1)], start=c(1962,2), end=c(2015,12), frequency=12)
indicator1.ts <- ts(approx.values[, indicator1], start=c(1962,2), end=c(2015,12), frequency=12)
indicator2.ts <- ts(approx.values[, indicator2], start=c(1962,2), end=c(2015,12), frequency=12)
los.ts <- ts(approx.values$LoS, start=c(1962,2), end=c(2015,12), frequency=12)
event.ts <- ts(as.numeric(los.ts > 0), start=c(1962,2), end=c(2015,12), frequency=12)

# look for trends/seasonality
plot(stl(indicator1.ts, s.window='periodic'))
plot(stl(diff(diff(indicator1.ts)), s.window="periodic")) # stationary but seasonal
plot(stl(diff(diff(diff(indicator1.ts)), 12), s.window="periodic"))

indicator1.ts.detrended <- diff(diff(diff(indicator1.ts)), 12)
indicator2.ts.detrended <- diff(diff(diff(indicator2.ts)), 12)
los.ts.detrended <- ts(los.ts[1:length(indicator1.ts.detrended)], start=c(1963,4), end=c(2015,12), frequency=12)
event.ts.detrended <- ts(event.ts[1:length(indicator1.ts.detrended)], start=c(1963,4), end=c(2015,12), frequency=12)

# binomial model to predict event
indicator1.ts.l6 <- stats::lag(indicator1.ts, k=-2)
event.ts.l1 <- ts(c(0, event.ts[1:length(event.ts) - 1]), start=c(1962,2), end=c(2015,12), frequency=12) # custom lag(?)
binomial.model <- glm(event.ts ~ event.ts.l1, family="binomial")

#binomial.model <- glm(event.ts.detrended ~ indicator1.ts.detrended + indicator2.ts.detrended, family="binomial")
#binomial.model <- glm(event.ts ~ stats::lag(event.ts, k=-1), family="binomial")
binomial.beta <- binomial.model$coefficients
summary(binomial.model)
paste0("Model rank: ", binomial.model$rank)
pred <- predict(binomial.model, indicator1.ts.detrended + indicator2.ts.detrended, type="response")
pred <- ts(pred, start=c(1963,3), end=c(2015,12), frequency=12)

# plot results
plot(event.ts)
lines(pred, col="red")

# poisson model to predict #LoS days in a month, λ = exp(θx) for Poisson
poisson.model <- glm(los.ts.detrended ~ indicator1.ts.detrended + indicator2.ts.detrended, family="poisson")  # NA for lagged values
poisson.beta <- poisson.model$coefficients
summary(poisson.model)
paste0("Model rank: ", poisson.model$rank)
poisson.pred <- predict(poisson.model, indicator1.ts.detrended + indicator2.ts.detrended, type="response")
poisson.pred <- ts(poisson.pred, start=c(1962,2), end=c(2015,12), frequency=12)

# plot results
plot(los.ts.detrended)
lines(poisson.pred, col="red")

# some plots
plot(stl(indicator.ts, s.window="periodic"))
plot(stl(los.ts, s.window="periodic"))


# OLDER STUF
#######
# split into train and validation sets
indicator.tr = window(indicator.ts, start=c(1962,2), end=c(2012,12))
indicator.valid = window(indicator.ts, start=c(2012,12),end=c(2015,12))
los.tr = window(los.ts, start=c(1962,2), end=c(2012,12))
los.valid = window(los.ts, start=c(2012,12),end=c(2015,12))

# combine
train <- cbind(indicator=indicator.tr, los=los.tr)
valid <- cbind(indicator=indicator.valid, los=los.valid)

####### Vector Autoregression
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

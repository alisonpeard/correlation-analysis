# GLARMA Example

### Example from Davis, Dunsmuir Wang (1999)
## MA(1,2,5), Pearson Residuals, Fisher Scoring

data(Polio)
y <- Polio[, 2]
X <- as.matrix(Polio[, 3:8])
colnames(Polio)

par(mfrow=c(3, 3))
plot(Polio[, 2], type='l')
plot(Polio[, 3], type='l')
plot(Polio[, 4], type='l')
plot(Polio[, 5], type='l')
plot(Polio[, 6], type='l')
plot(Polio[, 7], type='l')
plot(Polio[, 8], type='l')

glarmamod <- glarma(y, X, thetaLags = c(1,2,5), type = "Poi", method = "FS",
                    residuals = "Pearson", maxit = 100, grad = 1e-6)

glarmamod
summary(glarmamod)


## Score Type (GAS)  Residuals, Newton Raphson
## Note: Newton Raphson fails to converge from GLM initial estimates.
## Setting up the initial estimates by ourselves
init.delta <- glarmamod$delta
beta <- init.delta[1:6]
thetaInit <- init.delta[7:9]

glarmamod <- glarma(y, X, beta = beta, thetaLags = c(1, 2, 5),
                    thetaInit = thetaInit, type ="Poi", method = "NR",
                    residuals = "Score", maxit = 100, grad = 1e-6)
glarmamod
summary(glarmamod)

## AR(1,5), Pearson Residuals, Fisher Scoring
glarmamod <- glarma(y, X, phiLags = c(1, 5), type = "Poi", method = "FS",
                    residuals = "Pearson", maxit = 100, grad = 1e-6)
glarmamod
summary(glarmamod)

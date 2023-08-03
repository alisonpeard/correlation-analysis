library(gamlss)
library(ggplot2)
data <- read.csv("/Users/alison/Documents/RAPID/correlation-analysis/data_results/for_r/sample_ts.csv")
colnames(data)
data <- data[, c('date', 'LoS', 'ep_total', 'q50_anomaly_total', 'q50_deficit_total')]
data$cumulative.month <- 1:nrow(data)

model <- gamlss(data=data, formula=LoS ~ ep_total + q50_deficit_total + cumulative.month, trace=FALSE, family='PO')
summary(model)

model_mu <- fitted(model, what='mu')
data$model_Q50 <- qPO(0.50, model_mu)

ggplot(data) + theme_bw()+
  geom_point(aes(y= LoS, x = cumulative.month), col = "black") +
  geom_line(aes(y = model_q50, x = cumulative.month), col = "red")

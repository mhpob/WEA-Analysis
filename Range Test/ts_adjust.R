data <- readRDS('data and imports/rangetest_logit_data.rds')
data <- d_probs
names(data) <- tolower(gsub(' ', '_', names(data)))

# Remove outliers before modeling
library(tsoutliers)

d50_wea <- ts(na.omit(data[data$array == 'MD WEA', 'd50']))
tso_wea <- tso(d50_wea, types = 'AO',
               tsmethod = 'auto.arima',
               # cval = 5,
               maxit.iloop = 10, maxit.oloop = 10)
plot(tso_wea)

wea_adj <- tso_wea$yadj
tso_wea2 <- tso(wea_adj, types = 'AO',
               tsmethod = 'auto.arima',
               # cval = 5,
               maxit.iloop = 10, maxit.oloop = 10)


f <- forecast::tsclean(wea_adj)


forecast::findfrequency(wea_adj)

wea_adj <- ts(wea_adj, frequency = 9)
plot(stl(wea_adj, s.window = 'p'))

forecast::ndiffs(wea_adj)
d_wea_adj <- diff(wea_adj)
forecast::ndiffs(d_wea_adj)



d50_wea <- ts(na.omit(data[data$array == 'MD WEA', 'd50']), frequency = 9)
stl_wea <- stl(d50_wea, s.window = 'p')



d50_inner <- ts(na.omit(data[data$array == 'Inner', 'd50']))
tso_inner <- tso(d50_inner, types = 'IO', maxit.iloop = 10, maxit.oloop = 20)
plot(tso_inner)


c('IO', 'AO', 'LS', 'TC', 'SLS')



library(mgcv)

k <- gam(d50 ~ s(average_noise) +
           s(average_temperature) +
           array, data = data, method = 'REML', select = T)
summary(k)


library(forecast)
k <- seasadj(stl(ts(na.omit(data[data$array == 'Inner', 'd50']), frequency = 8), s.window = 'per'))

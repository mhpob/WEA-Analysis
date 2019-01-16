library(dplyr)
#' ---
#' title: "Removing outliers from D50 time series"
#' author: Mike O'Brien
#' ---
#' ```{r setup, echo = FALSE}
#' knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file())
#' ```
data <- readRDS('data and imports/rangetest_logit_data.rds')
names(data) <- tolower(gsub(' ', '_', names(data)))

# Remove outliers ----
# The WEA array is much more noisy than
# Two methods: tsoutliers::tso ('intervention analysis') and
# forecast::tsclean (Friedman's "super smoother").
d50_wea <- ts(na.omit(data[data$array == 'MD WEA', 'd50']))

tso_wea <- tsoutliers::tso(d50_wea, types = 'AO',
               tsmethod = 'auto.arima',
               maxit.iloop = 10, maxit.oloop = 10)
plot(tso_wea)

wea_adj <- tso_wea$yadj
f <- forecast::tsclean(d50_wea)

plot(wea_adj, type = 'n')
lines(d50_wea, lwd = 5)
lines(wea_adj, lwd = 3, col = 'gray')
lines(f, col = 'red')

forecast::findfrequency(wea_adj)

# Need to remove the two ridiculously incorrect points before trying to use
# unsupervised outlier removal.
data <- data %>%
  group_by(array) %>%
  mutate_at(vars(grep('d\\d\\d', names(.))),
            funs(ifelse((. > (4 * lead(.)) & . > (4 * lag(.))) |
                          . <= 0,
                        ((lag(.) + lead(.)) / 2), .)))




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

library(ggplot2); library(dplyr); library(changepoint)
# https://tdhock.github.io/change-tutorial/RK-CptWorkshop.html
# Needs changepoint version >2.2.4 or variance estimation will break with PELT
data <- readRDS('data and imports/rangetest_no_outliers.rds')
wea_all <- data[data$array == 'MD WEA',]
inner_all <- data[data$array == 'Inner',]

# D50 ----
# MD WEA
cp_d50wea <- cpt.meanvar(data[data$array == 'MD WEA',][['d50_adj']],
                      method = 'PELT', penalty = 'CROPS',
                      pen.value = c(15, 500))

plot(cp_d50wea, diagnostic = T)
abline(v = 2, col = 'blue')

plot(cp_d50wea, ncpts = 2)
abline(v = which(wea_all$date == '2018-04-11')) #tending date
abline(v = which(wea_all$date == '2018-08-09')) #tending date
cpwea_d50_dates <- wea_all[na.omit(
  cpts.full(cp_d50wea)[which(
    apply(cpts.full(cp_d50wea), 1,
          function(x) length(na.omit(x)))
    == 2),]), 'date'] #this number is the number of changepoints

list(mean = param.est(param(cp_d50wea, ncpts = 2))[[1]],
     stdev = sqrt(param.est(param(cp_d50wea, ncpts = 2))[[2]]))

# Inner
cp_d50inner <- cpt.meanvar(data[data$array == 'Inner',][['d50_adj']],
                        method = 'PELT', penalty = 'CROPS',
                        pen.value = c(15, 500))
plot(cp_d50inner, diagnostic = T)
abline(v = 3, col = 'blue')

plot(cp_d50inner, ncpts = 3)
abline(v = which(wea_all$date == '2018-04-11')) #tending date
abline(v = which(wea_all$date == '2018-08-09')) #tending date

cpinn_d50_dates <- inner_all[na.omit(
  cpts.full(cp_d50inner)[which(
    apply(cpts.full(cp_d50inner), 1,
          function(x) length(na.omit(x)))
    == 3),]), 'date']

list(mean = param.est(param(cp_d50inner, ncpts = 3))[[1]],
     stdev = sqrt(param.est(param(cp_d50inner, ncpts = 3))[[2]]))

# dT ----
# MD WEA
cp_dtwea <- cpt.meanvar(na.omit(data[data$array == 'MD WEA', 'dt']),
                        method = 'PELT', penalty = 'CROPS',
                        pen.value = c(20, 500))
plot(cp_dtwea, diagnostic = T)
abline(v = 3, col = 'blue')

plot(cp_dtwea, ncpts = 3)
cpwea_dt_dates <- wea_all[na.omit(
  cpts.full(cp_dtwea)[which(
    apply(cpts.full(cp_dtwea), 1,
          function(x) length(na.omit(x)))
    == 3),]), 'date']

list(mean = param.est(param(cp_dtwea, ncpts = 3))[[1]],
     stdev = sqrt(param.est(param(cp_dtwea, ncpts = 3))[[2]]))


# Inner
cp_dtinn <- cpt.meanvar(na.omit(data[data$array == 'Inner', 'dt']),
                        method = 'PELT', penalty = 'CROPS',
                        pen.value = c(25, 500))
plot(cp_dtinn, diagnostic = T)
abline(v = 1, col = 'blue')

plot(cp_dtinn, ncpts = 1)
cpinn_dt_dates <- inner_all[na.omit(
  cpts.full(cp_dtinn)[which(
    apply(cpts.full(cp_dtinn), 1,
          function(x) length(na.omit(x)))
    == 1),]), 'date']

list(mean = param.est(param(cp_dtinn, ncpts = 1))[[1]],
     stdev = sqrt(param.est(param(cp_dtinn, ncpts = 1))[[2]]))

# Noise ----
# MD WEA
cp_noisewea <- cpt.meanvar(na.omit(data[data$array == 'MD WEA', 'Average noise']),
                           method = 'PELT', penalty = 'CROPS',
                           # test.stat = 'Gamma',
                           pen.value = c(0, 500))
plot(cp_noisewea, diagnostic = T)
abline(v = 2, col = 'blue')

plot(cp_noisewea, ncpts = 2)
cpwea_noise_dates <- wea_all[na.omit(
  cpts.full(cp_noisewea)[which(
    apply(cpts.full(cp_noisewea), 1,
          function(x) length(na.omit(x)))
    == 2),]), 'date']

list(mean = param.est(param(cp_noisewea, ncpts = 2))[[1]],
     stdev = sqrt(param.est(param(cp_noisewea, ncpts = 2))[[2]]))


# Inner
cp_noiseinn <- cpt.meanvar(na.omit(data[data$array == 'Inner', 'Average noise']),
                        method = 'PELT', penalty = 'CROPS',
                        # test.stat = 'Gamma',
                        pen.value = c(0.1, 500))
plot(cp_noiseinn, diagnostic = T)
abline(v = 1, col = 'blue')

plot(cp_noiseinn, ncpts = 1)
cpinn_noise_dates <- inner_all[na.omit(
  cpts.full(cp_noiseinn)[which(
    apply(cpts.full(cp_noiseinn), 1,
          function(x) length(na.omit(x)))
    == 1),]), 'date']

list(mean = param.est(param(cp_noiseinn, ncpts = 1))[[1]],
     stdev = sqrt(param.est(param(cp_noiseinn, ncpts = 1))[[2]]))

# Temperature ----
# MD WEA
cp_tempwea <- cpt.meanvar(na.omit(data[data$array == 'MD WEA', 'Average temperature']),
                           method = 'PELT', penalty = 'CROPS',
                           pen.value = c(0, 500))
plot(cp_tempwea, diagnostic = T)
abline(v = 3, col = 'blue')

plot(cp_tempwea, ncpts = 3)
cpwea_temp_dates <- wea_all[na.omit(
  cpts.full(cp_tempwea)[which(
    apply(cpts.full(cp_tempwea), 1,
          function(x) length(na.omit(x)))
    == 3),]), 'date']


# Inner
cp_tempinn <- cpt.meanvar(na.omit(data[data$array == 'Inner', 'Average temperature']),
                           method = 'PELT', penalty = 'CROPS',
                           # test.stat = 'Gamma',
                           pen.value = c(0.1, 500))
plot(cp_tempinn, diagnostic = T)
abline(v = 4, col = 'blue')

plot(cp_tempinn, ncpts = 4)
cpinn_temp_dates <- inner_all[na.omit(
  cpts.full(cp_tempinn)[which(
    apply(cpts.full(cp_tempinn), 1,
          function(x) length(na.omit(x)))
    == 4),]), 'date']

# plotting ----
cp_df <- data.frame(
  start = c(min(wea_all$date), pull(cpwea_d50_dates),
            min(inner_all$date), pull(cpinn_d50_dates)),
  end = c(pull(cpwea_d50_dates), max(wea_all$date),
          pull(cpinn_d50_dates), max(inner_all$date)),
  mean = c(param.est(param(cp_d50wea, ncpts = 2))[[1]],
           param.est(param(cp_d50inner, ncpts = 3))[[1]]),
  stdev = c(sqrt(param.est(param(cp_d50wea, ncpts = 2))[[2]]),
            sqrt(param.est(param(cp_d50inner, ncpts = 3))[[2]])),
  array = c(rep('MD WEA', 3), rep('Inner', 4)))


ggplot() +
  geom_ribbon(data = data, aes(x = date, ymin = d95_adj, ymax = d5_adj),
              fill = 'lightgray')+
  geom_line(data = data, aes(x = date, y = d50_adj), lwd = 1) +
  geom_segment(data = cp_df, aes(x = start, xend = end, y = mean, yend = mean),
               col = 'red', lwd = 1) +
  labs(x = NULL, y = 'D50') +
  scale_x_date(date_breaks = 'month', date_labels = '%b') +
  scale_y_continuous(breaks = seq(0, 3000, 500)) +
  facet_wrap(~ array, ncol = 1, scales = 'free_y') +
  theme_bw()


cp_df <- data.frame(
  variable = c(rep('▲T', 6),
               rep('D50', 4),
               rep('Average noise', 5),
               rep('Average temperature', 9)),
  array = c(rep('MD WEA', 4), rep('Inner', 2),
            rep('MD WEA', 2), rep('Inner', 2),
            rep('MD WEA', 3), rep('Inner', 2),
            rep('MD WEA', 4), rep('Inner', 5)),
  start = c(min(data$date), cpwea_dt_dates,
            min(data$date), cpinn_dt_dates,
            min(data$date), cpwea_d50_dates,
            min(data$date), cpinn_d50_dates,
            min(data$date), cpwea_noise_dates,
            min(data$date), cpinn_noise_dates,
            min(data$date), cpwea_temp_dates,
            min(data$date), cpinn_temp_dates),
  end = c(cpwea_dt_dates, max(data$date),
          cpinn_dt_dates, max(data$date),
          cpwea_d50_dates, max(data$date),
          cpinn_d50_dates, max(data$date),
          cpwea_noise_dates, max(data$date),
          cpinn_noise_dates, max(data$date),
          cpwea_temp_dates, max(data$date),
          cpinn_temp_dates, max(data$date)),
  mean = c(param.est(param(cp_dtwea, ncpts = 3))[[1]],
           param.est(param(cp_dtinn, ncpts = 1))[[1]],
           param.est(param(cp_d50wea, ncpts = 1))[[1]],
           param.est(param(cp_d50inner, ncpts = 1))[[1]],
           param.est(param(cp_noisewea, ncpts = 2))[[1]],
           param.est(param(cp_noiseinn, ncpts = 1))[[1]],
           param.est(param(cp_tempwea, ncpts = 3))[[1]],
           param.est(param(cp_tempinn, ncpts = 4))[[1]]),
  stdev = c(sqrt(param.est(param(cp_dtwea, ncpts = 3))[[2]]),
            sqrt(param.est(param(cp_dtinn, ncpts = 1))[[2]]),
            sqrt(param.est(param(cp_d50wea, ncpts = 1))[[2]]),
            sqrt(param.est(param(cp_d50inner, ncpts = 1))[[2]]),
            sqrt(param.est(param(cp_noisewea, ncpts = 2))[[2]]),
            sqrt(param.est(param(cp_noiseinn, ncpts = 1))[[2]]),
            sqrt(param.est(param(cp_tempwea, ncpts = 3))[[2]]),
            sqrt(param.est(param(cp_tempinn, ncpts = 4))[[2]])))


# ggplot() +
#   geom_line(data = data, aes(x = date, y = dt), lwd = 1) +
#   geom_segment(data = cp_df, aes(x = start, xend = end, y = mean, yend = mean),
#                col = 'red', lwd = 1) +
#   labs(x = NULL, y = expression(Delta*T~'('*degree*'C)')) +
#   scale_x_date(date_breaks = 'month', date_labels = '%b') +
#   facet_wrap(~ array, ncol = 1, scales = 'free_y') +
#   theme_bw()


j <- select(data, array, date, `Average noise`,  `Average temperature`, dt, D50) %>%
  tidyr::gather(key = 'variable', value = 'value',
                `Average noise`,  `Average temperature`, dt, D50) %>%
  mutate(variable = ifelse(variable == 'dt', '▲T', variable),
         variable = factor(variable, levels = c('D50', 'Average temperature',
                                                'Average noise', '▲T'),
                           ordered= T))

ggplot() +
  geom_line(data = j, aes(x = date, y = value)) +
  geom_segment(data = cp_df, aes(x = start, xend = end, y = mean, yend = mean),
               col = 'red', lwd = 1) +
  facet_grid(variable ~ array, scale = 'free_y') +
  theme_bw()

ggplot() +
  geom_line(data = filter(j, grepl('noise|T', variable)),
            aes(x = date, y = value)) +
  geom_segment(data = filter(cp_df, grepl('noise|T', variable)),
               aes(x = start, xend = end, y = mean, yend = mean),
               col = 'red', lwd = 1) +
  facet_grid(variable ~ array, scale = 'free_y') +
  labs(x = NULL, y = 'Value') +
  theme_bw()

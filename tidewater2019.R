library(dplyr)
data <- readRDS('data and imports/rangetest_logit_binary.RDS') %>%
  mutate(
    array = ifelse(array == 'MD WEA', 'Outer', 'Inner'),
         array = as.factor(array),
         date_f = as.factor(date)) %>%
  rename(noise = `Average noise`,
         tilt = `Tilt angle`,
         bwt = `Average temperature`)

library(lme4)
m1 <- glmer(cbind(success, fail) ~ distance + (0+distance | date_f:array),
            data = data, family = 'binomial')

coefs <- coef(m1)[[1]]

plot.ts(-coefs[grepl('Inner', row.names(coefs)), 1] /
          coefs[grepl('Inner', row.names(coefs)), 2])

plot.ts(-coefs[grepl('Outer', row.names(coefs)), 1] /
          coefs[grepl('Outer', row.names(coefs)), 2])


d50 <- data.frame(array = 'Inner',
                  d50 = -coefs[grepl('Inner', row.names(coefs)), 1] /
                    coefs[grepl('Inner', row.names(coefs)), 2],
                  date = unique(data$date))
d50 <- rbind(d50,
             data.frame(array = 'Outer',
                        d50 = -coefs[grepl('Outer', row.names(coefs)), 1] /
                          coefs[grepl('Outer', row.names(coefs)), 2],
             date = unique(data$date)))

library(ggplot2)

ggplot() + geom_line(data = d50, aes(x = date, y = d50)) +
  facet_wrap(~array, nrow = 2)




# extract random effects and conditional vcov matrix
j <- ranef(m1, condVar = T)
# pull out variances
f <- attr(j[[1]], 'postVar')


new_vc <- vcov(m1)
new_vc[2,2] <- f[1]

param_reps <- MASS::mvrnorm(100, c(fixef(m1)[1], j[[1]][1,1]), new_vc)



ci.pred <- function(n, obj){
  r_efs <- ranef(obj, condVar = T)
  r_vars <- attr(r_efs[[1]], 'postVar')

  coefs <- coef(obj)[[1]]

  vc <- vcov(obj)
  sd_estimate <- rep(NULL, nrow(coef(obj)[[1]]))

  for(i in seq(1, nrow(coefs), 1)){
    vc[2, 2] <- vcov(obj)[2,2] + r_vars[i]
    param_reps <- as.data.frame(
      MASS::mvrnorm(n, as.matrix(coefs[i,]), vc)
    )
    param_reps$estimate <- -param_reps[1]/param_reps[2]
    sd_estimate[i] <- sqrt(var(param_reps$estimate))
  }

  mod_d50 <- -coefs[,1]/coefs[,2]


  out <- data.frame(d50 = mod_d50,
           uci = mod_d50 + 1.96 * sd_estimate,
           lci = mod_d50 - 1.96 * sd_estimate)
  out
}

k <- ci.pred(1000, m1)


library(ggplot2)
k <- as.data.frame(k)
k$array <- rep(c('Inner', 'Outer'))
k$date <- rep(seq.Date(min(data$date), max(data$date), 'day'), each = 2)


ggplot() + geom_ribbon(data = k, aes(x = date, ymin = lci, ymax = uci),
                       fill = 'gray') +
  geom_line(data = k, aes(x = date, y = d50), lwd = 1) +
  facet_wrap(~array, ncol = 1) +
  labs(x = NULL, y = 'D50 (m)') +
  theme_bw() +
  theme(text = element_text(size=20))

## Changepoints ----

library(changepoint)
cp_d50out <- cpt.meanvar(k[k$array == 'Outer',][['d50']],
                         method = 'PELT', penalty = 'CROPS',
                         pen.value = c(15, 500))
cpout_d50_dates <- k[k$array == 'Outer',][na.omit(
  cpts.full(cp_d50out)[which(
    apply(cpts.full(cp_d50out), 1,
          function(x) length(na.omit(x)))
    == 2),]), 'date'] #this number is the number of changepoints

list(mean = param.est(param(cp_d50out, ncpts = 2))[[1]],
     stdev = sqrt(param.est(param(cp_d50out, ncpts = 2))[[2]]))


cp_d50inn <- cpt.meanvar(k[k$array == 'Inner',][['d50']],
                         method = 'PELT', penalty = 'CROPS',
                         pen.value = c(15, 500))
cpinn_d50_dates <- k[k$array == 'Inner',][na.omit(
  cpts.full(cp_d50inn)[which(
    apply(cpts.full(cp_d50inn), 1,
          function(x) length(na.omit(x)))
    == 2),]), 'date'] #this number is the number of changepoints

list(mean = param.est(param(cp_d50inn, ncpts = 2))[[1]],
     stdev = sqrt(param.est(param(cp_d50inn, ncpts = 2))[[2]]))

cp_df <- data.frame(
  start = c(min(k[k$array == 'Outer', 'date']), cpout_d50_dates,
            min(k[k$array == 'Inner', 'date']), cpinn_d50_dates),
  end = c(cpout_d50_dates, max(k[k$array == 'Outer', 'date']),
          cpinn_d50_dates, max(k[k$array == 'Inner', 'date'])),
  mean = c(param.est(param(cp_d50out, ncpts = 2))[[1]],
           param.est(param(cp_d50inner, ncpts = 2))[[1]]),
  stdev = c(sqrt(param.est(param(cp_d50out, ncpts = 2))[[2]]),
            sqrt(param.est(param(cp_d50inner, ncpts = 2))[[2]])),
  array = c(rep('Outer', 3), rep('Inner', 3)))

ggplot() + geom_ribbon(data = k, aes(x = date, ymin = lci, ymax = uci),
                       fill = 'gray') +
  geom_line(data = k, aes(x = date, y = d50), lwd = 1) +
  geom_segment(data = cp_df, aes(x = start, xend = end, y = mean, yend = mean),
               col = 'red', lwd = 1) +
  facet_wrap(~array, ncol = 1) +
  labs(x = NULL, y = 'D50 (m)') +
  theme_bw() +
  theme(text = element_text(size = 20))

## DelT ----
library(ggplot2)
ggplot() + geom_line(data = data, aes(x = date, y = dt, color = array),
                     lwd = 1.5) +
  scale_y_reverse() +
  labs(color = 'Array', x = NULL, y = expression(Delta*T~'(°C)')) +
  theme_bw() +
  theme(text = element_text(size = 20),
        legend.position = c(0.9, 0.2))


## Modeling ----
library(dplyr)
data <- left_join(data, d50) %>%
  distinct(date, array, .keep_all = T)

mod_data <- data[, names(data) %in% c('date', 'array', 'd50', 'noise',
                                      'tilt', 'dt', 'wdir', 'wspd', 'wvht',
                                      'dpd', 'apd', 'mwd', 'pres')]
mod_data <- mod_data[complete.cases(mod_data),]

library(mgcv)
m1 <- gam(d50 ~ s(dt, array, bs = 're') + s(dt, bs = 'ts') +
             s(noise, array, bs = 're') + s(noise, bs = 'ts') +
             s(tilt, array, bs = 're') + s(tilt, bs = 'ts') +
             s(array, bs = 're') +
             s(wdir, bs = 'ts')+
             s(wspd, bs = 'ts') +
             s(wvht, bs = 'ts') +
             s(dpd, bs = 'ts') +
             s(mwd, bs = 'ts') +
             s(pres, bs = 'ts'),
           data = mod_data,
           family = Gamma(),
           method = 'REML',
           verbosePQL = F)

summary(m1)

m1b <- gam(d50 ~ s(dt, bs = 'ts') +
            s(noise, bs = 'ts') +
            s(tilt, bs = 'ts') +
            s(array, bs = 're') +
            s(wdir, bs = 'ts')+
            s(wspd, bs = 'ts') +
            s(wvht, bs = 'ts') +
            s(dpd, bs = 'ts') +
            s(mwd, bs = 'ts') +
            s(pres, bs = 'ts'),
          data = mod_data,
          family = Gamma(),
          method = 'REML',
          verbosePQL = F)

summary(m1b)
AIC(m1b)

m2 <- gam(d50 ~ s(dt, array, bs = 're') + s(dt, bs = 'ts') +
             s(noise, array, bs = 're') + s(noise, bs = 'ts') +
             s(tilt, array, bs = 're') + s(tilt, bs = 'ts') +
             s(array, bs = 're') +
             s(wspd, bs = 'ts') +
             s(wvht, bs = 'ts') +
             s(dpd, bs = 'ts') +
             s(mwd, bs = 'ts') +
             s(pres, bs = 'ts'),
           data = mod_data,
           family = Gamma(),
           method = 'REML',
           verbosePQL = F)

summary(m2)

m3 <- gam(d50 ~ s(dt, bs = 'ts') +
            s(noise, array, bs = 're') + s(noise, bs = 'ts') +
            s(tilt, array, bs = 're') + s(tilt, bs = 'ts') +
            s(array, bs = 're') +
            s(wspd, bs = 'ts') +
            s(wvht, bs = 'ts') +
            s(dpd, bs = 'ts') +
            s(mwd, bs = 'ts') +
            s(pres, bs = 'ts'),
          data = mod_data,
          family = Gamma(),
          method = 'REML',
          verbosePQL = F)

summary(m3)

m4 <- gam(d50 ~ s(dt, bs = 'ts') +
            s(noise, array, bs = 're') + s(noise, bs = 'ts') +
            s(tilt, array, bs = 're') + s(tilt, bs = 'ts') +
            s(array, bs = 're') +
            s(wspd, bs = 'ts') +
            s(wvht, bs = 'ts') +
            s(dpd, bs = 'ts') +
            s(mwd, bs = 'ts'),
          data = mod_data,
          family = Gamma(),
          method = 'REML',
          verbosePQL = F)

summary(m4)

m5 <- gam(d50 ~ s(dt, bs = 'ts') +
            s(noise, array, bs = 're') + s(noise, bs = 'ts') +
            s(array, bs = 're') +
            s(wspd, bs = 'ts') +
            s(wvht, bs = 'ts') +
            s(dpd, bs = 'ts') +
            s(mwd, bs = 'ts'),
          data = mod_data,
          family = Gamma(),
          method = 'REML',
          verbosePQL = F)

summary(m5)

m6 <- gam(d50 ~ s(dt, bs = 'ts') +
            s(noise, array, bs = 're') + s(noise, bs = 'ts') +
            s(array, bs = 're') +
            s(wspd, bs = 'ts') +
            s(wvht, bs = 'ts') +
            s(mwd, bs = 'ts'),
          data = mod_data,
          family = Gamma(),
          method = 'REML',
          verbosePQL = F)

summary(m6)

m7 <- gam(d50 ~ s(dt, bs = 'ts') +
            s(noise, array, bs = 're') + s(noise, bs = 'ts') +
            s(array, bs = 're') +
            wspd +
            wvht +
            s(mwd, bs = 'ts'),
          data = mod_data,
          family = Gamma(),
          method = 'REML',
          verbosePQL = F)

summary(m7)

m8 <- gam(d50 ~ s(dt, bs = 'ts') +
            s(noise, bs = 'ts') +
            s(array, bs = 're') +
            wspd +
            wvht +
            s(mwd, bs = 'ts'),
          data = mod_data,
          family = Gamma(),
          method = 'REML',
          verbosePQL = F)

summary(m8)
AIC(m8)


## Start with only a global random intercept
m1b <- gam(d50 ~ s(dt, bs = 'ts') +
             s(noise, bs = 'ts') +
             s(tilt, bs = 'ts') +
             s(array, bs = 're') +
             s(wdir, bs = 'ts')+
             s(wspd, bs = 'ts') +
             s(wvht, bs = 'ts') +
             s(dpd, bs = 'ts') +
             s(mwd, bs = 'ts') +
             s(pres, bs = 'ts'),
           data = mod_data,
           family = Gamma(),
           method = 'REML',
           verbosePQL = F)

summary(m1b)
AIC(m1b)

m2b <- gam(d50 ~ s(dt, bs = 'ts') +
             s(noise, bs = 'ts') +
             s(tilt, bs = 'ts') +
             s(array, bs = 're') +
             s(wspd, bs = 'ts') +
             s(wvht, bs = 'ts') +
             s(dpd, bs = 'ts') +
             s(mwd, bs = 'ts') +
             s(pres, bs = 'ts'),
           data = mod_data,
           family = Gamma(),
           method = 'REML',
           verbosePQL = F)

summary(m2b)
AIC(m2b)

m3b <- gam(d50 ~ s(dt, bs = 'ts') +
             s(noise, bs = 'ts') +
             s(array, bs = 're') +
             s(wspd, bs = 'ts') +
             s(wvht, bs = 'ts') +
             s(dpd, bs = 'ts') +
             s(mwd, bs = 'ts') +
             s(pres, bs = 'ts'),
           data = mod_data,
           family = Gamma(),
           method = 'REML',
           verbosePQL = F)

summary(m3b)
AIC(m3b)


m4b <- gam(d50 ~ s(dt, bs = 'ts') +
             s(noise, bs = 'ts') +
             s(array, bs = 're') +
             s(wspd, bs = 'ts') +
             s(wvht, bs = 'ts') +
             s(dpd, bs = 'ts') +
             s(mwd, bs = 'ts'),
           data = mod_data,
           family = Gamma(),
           method = 'REML',
           verbosePQL = F)

summary(m4b)
AIC(m4b)

m5b <- gam(d50 ~ s(dt, bs = 'ts') +
             s(noise, bs = 'ts') +
             s(array, bs = 're') +
             s(tilt, bs = 'ts') +
             s(wspd, bs = 'ts') +
             s(wvht, bs = 'ts') +
             s(dpd, bs = 'ts') +
             s(mwd, bs = 'ts'),
           data = mod_data,
           family = Gamma(),
           method = 'REML',
           verbosePQL = F)

summary(m5b)
AIC(m5b)





m1 <- gam(d50 ~ s(dt, array, bs = 're') + s(dt, bs = 'ts') +
            s(noise, array, bs = 're') + s(noise, bs = 'ts') +
            s(tilt, array, bs = 're') + s(tilt, bs = 'ts') +
            s(array, bs = 're') +
            s(wdir, bs = 'ts')+
            s(wspd, bs = 'ts') +
            s(wvht, bs = 'ts') +
            s(dpd, bs = 'ts') +
            s(mwd, bs = 'ts') +
            s(pres, bs = 'ts'),
          data = mod_data,
          family = Gamma(),
          method = 'REML',
          verbosePQL = F, na.action = 'na.fail')

# cl <- parallel::makeCluster(parallel::detectCores() - 2)
# parallel::clusterExport(cl, 'mod_data')
# parallel::clusterEvalQ(cl, library(mgcv))
# parallel::clusterEvalQ(cl, options(na.action = "na.fail"))
# k <- pdredge(m1, evaluate = T, fixed = ~ s(array, bs = 're'), cluster = cl)

model.sel(k)[1:15]

m_sel1c <- gam(d50 ~
                s(dt, array, bs = 're') + s(dt, bs = 'ts') +
                s(noise, array, bs = 're') + s(noise, bs = 'ts') +
                s(tilt, array, bs = 're') + s(tilt, bs = 'ts') +
                s(wspd, bs = 'ts') +
                s(wvht, bs = 'ts') +
                s(mwd, bs = 'ts'),
              data = mod_data,
              family = Gamma(),
              method = 'REML')

m_sel1b <- gam(d50 ~ s(array, bs = 're') +
                s(dt, array, bs = 're') + s(dt, bs = 'ts') +
                s(noise, array, bs = 're') + s(noise, bs = 'ts') +
                s(tilt, array, bs = 're') + s(tilt, bs = 'ts') +
                s(wspd, bs = 'ts') +
                wvht +
                s(mwd, bs = 'ts'),
              data = mod_data,
              family = Gamma(),
              method = 'REML')

m_sel2 <- gam(d50 ~ s(array, bs = 're') +
                s(dt, bs = 'ts') +
                s(noise, array, bs = 're') + s(noise, bs = 'ts') +
                s(tilt, array, bs = 're') + s(tilt, bs = 'ts') +
                s(wspd, bs = 'ts') +
                s(wvht, bs = 'ts') +
                s(mwd, bs = 'ts'),
              data = mod_data,
              family = Gamma(),
              method = 'REML')

m_sel2b <- gam(d50 ~
                 s(dt, bs = 'ts') +
                 s(noise, array, bs = 're') + s(noise, bs = 'ts') +
                 s(tilt, array, bs = 're') + s(tilt, bs = 'ts') +
                 s(wspd, bs = 'ts') +
                 wvht +
                 s(mwd, bs = 'ts'),
               data = mod_data,
               family = Gamma(),
               method = 'REML')
m_sel2c <- gam(d50 ~ s(array, bs = 're') +
                 s(dt, bs = 'ts') +
                 s(noise, array, bs = 're') + s(noise, bs = 'ts') +
                 s(tilt, array, bs = 're') + s(tilt, bs = 'ts') +
                 wspd +
                 wvht +
                 s(mwd, bs = 'ts'),
               data = mod_data,
               family = Gamma(),
               method = 'REML')

m_sel2d <- gam(d50 ~ s(array, bs = 're') +
                 s(dt, bs = 'ts') +
                 s(noise, array, bs = 're') + s(noise, bs = 'ts') +
                 s(tilt, array, bs = 're') + s(tilt, bs = 'ts') +
                 wspd +
                 wvht +
                 te(mwd, wdir, bs = 'ts'),
               data = mod_data,
               family = Gamma(),
               method = 'REML')

## Choose m_sel2b
dt_p <- data.frame(dt = seq(min(mod_data$dt), max(mod_data$dt), 0.1),
                   array = 'Inner',
                   noise = median(mod_data$noise),
                   tilt = median(mod_data$tilt),
                   wspd = mean(mod_data$wspd),
                   wvht = mean(mod_data$wvht),
                   mwd = mean(mod_data$mwd))
pred <- predict(m_sel2b, dt_p, type = 'response')
pred_se <- predict(m_sel2b, dt_p, type = 'link', se = T)
pred <- data.frame(dt = dt_p$dt,
                   pred,
                   UCI = 1/(pred_se$fit - 1.96*pred_se$se.fit),
                   LCI = 1/(pred_se$fit + 1.96*pred_se$se.fit))

library(ggplot2)
ggplot(data = pred) +
  geom_ribbon(aes(x = dt, ymin = LCI, ymax = UCI),
              fill = 'gray') +
  geom_line(aes(x = dt, y = pred)) +
  geom_rug() +
  labs(x = expression(Delta*T~'(°C)'), y = 'D50 (m)') +
  theme_bw() +
  theme(text = element_text(size = 20))



dt_p <- data.frame(dt = 0,
                   array = 'Inner',
                   noise = seq(min(mod_data$noise), max(mod_data$noise), 1),
                   tilt = median(mod_data$tilt),
                   wspd = mean(mod_data$wspd),
                   wvht = mean(mod_data$wvht),
                   mwd = mean(mod_data$mwd))
pred <- predict(m_sel2b, dt_p, type = 'response')
pred_se <- predict(m_sel2b, dt_p, type = 'link', se = T)
pred <- data.frame(var = dt_p$noise,
                   pred,
                   UCI = 1/(pred_se$fit - 1.96*pred_se$se.fit),
                   LCI = 1/(pred_se$fit + 1.96*pred_se$se.fit))

library(ggplot2)
ggplot(data = pred) +
  geom_ribbon(aes(x = var, ymin = LCI, ymax = UCI),
              fill = 'gray') +
  geom_line(aes(x = var, y = pred)) +
  labs(x = 'Noise (mV)', y = 'D50 (m)') +
  ylim(0, 1500)+
  theme_bw() +
  theme(text = element_text(size = 20))



p_mod <- predict(m_sel1c, type = 'response',
                 # next line removes random effects (anything with 'array')
                 exclude = grep('array', row.names(summary(m_sel1c)$s.table),
                                value = T))

p_mod_se <- predict(m_sel1c, type = 'link', se = T,
                    exclude = grep('array', row.names(summary(m_sel1c)$s.table),
                                   value = T))

library(dplyr)
p_mod <- cbind(data.frame(mod_data), data.frame(response_fit = p_mod), p_mod_se)
p_mod <- p_mod %>%
  # This calculates the lower and upper CIs and transforms them to the response scale
  mutate(LCI = 1/(fit + 2 * se.fit),
         UCI = 1/(fit - 2 * se.fit))

p <- p_mod %>%
  select(date, array, Observed = d50, Estimated = response_fit, LCI, UCI) %>%
  # use the time period where all predictors are available
  filter(date <= '2018-08-19') %>%
  tidyr::gather(key = 'key', value = 'val', -date, -array, -LCI, -UCI)

library(ggplot2)
ggplot() +
  geom_ribbon(data = p[p$key == 'Estimated',],
              aes(x = date, ymin = LCI, ymax =  UCI),
              fill = 'gray', alpha = 0.5) +
  geom_line(data = p[p$key %in% c('Observed', 'Estimated'),],
            aes(x = date, y = val, linetype = key)) +
  labs(x=NULL, y = 'Estimated D50', linetype = NULL) +
  facet_wrap(~array, nrow = 2) +
  theme_bw()





## Ran out of time. have to go with previously-fitted model
mod7 <- gam(d50 ~ s(dt, bs = 'ts') +
              s(noise, array, bs = 're') + s(noise, bs = 'ts') +
              s(tilt, bs = 'ts') +
              wspd +
              s(wvht, bs = 'ts') +
              dpd,
            data = mod_data,
            family = Gamma(),
            method = 'REML')

summary(mod7)

p_mod <- predict(mod7, type = 'response',
                 # next line removes random effects (anything with 'array')
                 exclude = grep('array', row.names(summary(mod7)$s.table),
                                value = T))

p_mod_se <- predict(mod7, type = 'link', se = T,
                    exclude = grep('array', row.names(summary(mod7)$s.table),
                                   value = T))
library(dplyr)
p_mod <- cbind(data.frame(mod_data), data.frame(response_fit = p_mod), p_mod_se)
p_mod <- p_mod %>%
  # This calculates the lower and upper CIs and transforms them to the response scale
  mutate(LCI = 1/(fit + 2 * se.fit),
         UCI = 1/(fit - 2 * se.fit))

p <- p_mod %>%
  select(date, array, Estimated = d50, Predicted = response_fit, LCI, UCI) %>%
  # use the time period where all predictors are available
  filter(date <= '2018-08-19') %>%
  tidyr::gather(key = 'key', value = 'val', -date, -array, -LCI, -UCI)

library(ggplot2)
ggplot() +
  geom_ribbon(data = p[p$key == 'Predicted',],
              aes(x = date, ymin = LCI, ymax =  UCI),
              fill = 'gray', alpha = 0.5) +
  geom_line(data = p[p$key %in% c('Estimated', 'Predicted'),],
            aes(x = date, y = val, linetype = key)) +
  labs(x=NULL, y = 'D50 (m)', linetype = NULL) +
  facet_wrap(~array, nrow = 2) +
  theme_bw()+
  theme(text = element_text(size = 20),
        legend.position = c(0.2, 0.9))





dt_p <- data.frame(dt = median(mod_data$dt),
                   array = 'Inner',
                   noise = seq(min(mod_data$noise),max(mod_data$noise),length =100),
                   tilt = median(mod_data$tilt),
                   wspd = median(mod_data$wspd),
                   wvht = median(mod_data$wvht),
                   dpd = median(mod_data$dpd))
pred <- predict(mod7, dt_p, type = 'response',
                exclude = grep('array', row.names(summary(m_sel2b)$s.table),
                               value = T))
pred_se <- predict(mod7, dt_p, type = 'link', se = T,
                   exclude = grep('array', row.names(summary(mod7)$s.table),
                                  value = T))
pred <- data.frame(var = dt_p$noise,
                   pred,
                   UCI = 1/(pred_se$fit - 1.96*pred_se$se.fit),
                   LCI = 1/(pred_se$fit + 1.96*pred_se$se.fit))

library(ggplot2)
ggplot(data = pred) +
  geom_ribbon(aes(x = var, ymin = LCI, ymax = UCI),
              fill = 'gray') +
  geom_line(aes(x = var, y = pred)) +
  geom_rug(data = mod_data, aes(x= noise)) +
  labs(x ='Noise (mV)', y = 'D50 (m)') +
  theme_bw() +
  theme(text = element_text(size = 20))


data <- readRDS('data and imports/rangetest_logit_binary.RDS')

pl_dat <- filter(data, date >= '2018-09-01')
m1 <- glm(cbind(success, fail) ~ distance + array,
            data = pl_dat, family = 'binomial')

new.dat <- data.frame(distance = seq(0, 800, 10),
                      array = rep(c('Inner', 'MD WEA'),
                                  each = length(seq(0,800,10) *2)))
pred <- predict(m1, new.dat, type = 'response',
             se.fit = T)

pred <- cbind(new.dat,
           data.frame(fit = exp(pred$fit)/(1+exp(pred$fit)),
                      lci = exp(pred$fit - 1.96 * pred$se.fit)/(1+exp(pred$fit- 1.96*pred$se.fit)),
                      uci = exp(pred$fit+1.96*pred$se.fit)/(1+exp(pred$fit+1.96*pred$se.fit))))


pred <- cbind(new.dat,
              data.frame(fit = pred$fit,
                         lci = pred$fit- 1.96*pred$se.fit,
                         uci = pred$fit+1.96*pred$se.fit))

# coefs <- coef(m1)[[1]]
#
# m1_fort <- fortify(m1) %>%
#   mutate(fit = exp(.fitted)/(1 + exp(.fitted)))
# k <- predict(m1, data.frame(date = unique()))
#
# pl_dat <- filter(data, date >= '2018-09-01')
# k <- predict(m1, data.frame(date = rep(unique(pl_dat$date),
#                                        times = 2),
#                             distance = seq(0, 800, 10),
#                             array = rep(c('Inner', 'MD WEA'),
#                                         each = nrow(pl_dat)/2)))

library(ggplot2)
ggplot() +
  geom_hline(aes(yintercept = 0.5)) +
  geom_boxplot(data = pl_dat,
               aes(x = distance, y = success/(success + fail),
               group = interaction(distance, array),
                   fill = array))+
  scale_fill_manual(values = c('#f8766d', '#00ba38')) +
  geom_smooth(data = pl_dat, aes(x = distance, y = success/(success + fail),
                                 lty = array),
              method = 'glm', method.args = list(family = 'binomial'))+


  geom_ribbon(data = pred, aes(x = distance, ymin = lci, ymax = uci,
                               group = array),
              fill = 'grey') +
  geom_line(data = pred, aes(x = distance, y = fit, color = array)) +
  scale_color_manual(values = c('red', 'green')) +
  # facet_wrap(~ array, nrow = 2) +
  labs(x = 'Distance', y = 'Frequency of detection',
       color = 'Array', fill = 'Array') +
  theme_bw() +
  theme(legend.justification = c(0, 0), legend.position = c(0.01, 0))






library(vegan)
v_dat <- data[, names(data) %in% c('noise', 'bwt',
                                   'tilt', 'sst', 'dt', 'wdir','wspd',
                                   'gst','wvht','dpd','apd','mwd','pres','atmp',
                                   'wtmp', 'd50')]
v_dat <- v_dat[complete.cases(v_dat),]

pca <- rda(v_dat, scale = T)
biplot(pca, type = c('text', 'points'), scaling = -1)


library(mgcv)

m1 <- gam(d50 ~ s(bwt, bs = 'ts') + s(noise, bs = 'ts') + s(tilt, bs = 'ts') +
            s(array, bs = 're'),
          data = data,
          family = Gamma,
          method = 'REML')

m2 <- gam(d50 ~ s(bwt, bs = 'ts') + s(noise, bs = 'ts') +
            s(array, bs = 're'),
          data = data,
          family = Gamma,
          method = 'REML')


m3 <- gam(d50 ~ s(bwt, bs = 'ts') +  s(bwt, array, bs = 're') +
            s(noise, bs = 'ts') + s(noise, array, bs = 're') +
            s(tilt, bs = 'ts') + s(tilt, array, bs = 're') +
            s(array, bs = 're'),
          data = data,
          family = Gamma,
          method = 'REML')
# AIC(m3)
# [1] 9184.208

m4 <- gam(d50 ~ s(bwt, bs = 'ts') +  s(bwt, array, bs = 're') +
            s(noise, bs = 'ts') + s(noise, array, bs = 're') +
            s(array, bs = 're'),
          data = data,
          family = Gamma,
          method = 'REML')
plot(m4)


#choose m2
p_mod <- predict(m2, type = 'link', se = T,
                    exclude = grep('array', row.names(summary(m2)$s.table),
                                   value = T))
library(dplyr)
p_mod <- cbind(data, data.frame(p_mod))
p_mod <- p_mod %>%
  # This calculates the lower and upper CIs and transforms them to the response scale
  mutate(Predicted = 1/fit,
         LCI = 1/(fit + 2 * se.fit),
         UCI = 1/(fit - 2 * se.fit))

p <- p_mod %>%
  select(date, array, Estimated = d50, Predicted, LCI, UCI) %>%
  # use the time period where all predictors are available
  # filter(date <= '2018-08-19') %>%
  tidyr::gather(key = 'key', value = 'val', -date, -array, -LCI, -UCI)

library(ggplot2)
ggplot() +
  geom_ribbon(data = p[p$key == 'Predicted',],
              aes(x = date, ymin = LCI, ymax =  UCI),
              fill = 'gray', alpha = 0.5) +
  geom_line(data = p[p$key %in% c('Estimated', 'Predicted'),],
            aes(x = date, y = val, linetype = key)) +
  scale_linetype_manual(values=c("dashed", "solid")) +
  labs(x=NULL, y = 'D50 (m)', linetype = NULL) +
  facet_wrap(~array, nrow = 2) +
  theme_bw()+
  theme(text = element_text(size = 20),
        legend.position = c(0.2, 0.9))



bwt_p <- data.frame(bwt = seq(min(data$bwt),max(data$bwt),length =100),
                   noise = median(data$noise),
                   array = 'Inner')
pred <- predict(m2, bwt_p, type = 'link',
                exclude = grep('array', row.names(summary(m2)$s.table),
                               value = T),
                se = T)
pred <- data.frame(var = bwt_p$bwt,
                   pred = 1/pred$fit,
                   UCI = 1/(pred$fit - 2*pred$se.fit),
                   LCI = 1/(pred$fit + 2*pred$se.fit))
library(ggplot2)
ggplot() +
  geom_ribbon(data = pred, aes(x = var, ymin = LCI, ymax = UCI),
              fill = 'pink') +
  geom_point(data = data, aes(x = bwt, y = d50)) +
  geom_line(data = pred, aes(x = var, y = pred), lwd = 1, color = 'red') +
  geom_rug(data = data, aes(x= bwt)) +
  labs(x ='Bottom Water Temperature (°C)', y = 'D50 (m)') +
  theme_bw() +
  theme(text = element_text(size = 20))



noise_p <- data.frame(bwt = mean(data$bwt),
                    noise = seq(min(data$noise),max(data$noise),length =100),
                    array = 'Inner')
pred <- predict(m2, noise_p, type = 'link',
                exclude = grep('array', row.names(summary(m2)$s.table),
                               value = T),
                se = T)
pred <- data.frame(var = noise_p$noise,
                   pred = 1/pred$fit,
                   UCI = 1/(pred$fit - 2*pred$se.fit),
                   LCI = 1/(pred$fit + 2*pred$se.fit))
library(ggplot2)
ggplot(data = pred) +
  geom_point(data =data, aes(x = noise, y = d50)) +
  geom_ribbon(aes(x = var, ymin = LCI, ymax = UCI),
              fill = 'pink') +
  geom_point(data =data, aes(x = noise, y = d50)) +
  geom_line(aes(x = var, y = pred), color = 'red', lwd = 1) +
  geom_rug(data = data, aes(x = noise)) +

  labs(x ='Ambient Noise (mV)', y = 'D50 (m)') +
  theme_bw() +
  theme(text = element_text(size = 20))

# ggplot(data = pred) +
#   geom_point(data = data, aes(x = noise, y = d50)) +
#   geom_line(aes(x = var, y = LCI), linetype = 'dashed', color = 'red', lwd = 1) +
#   geom_line(aes(x = var, y = UCI), linetype = 'dashed', color = 'red', lwd = 1) +
#   geom_line(aes(x = var, y = pred), color = 'red', lwd = 1) +
#   geom_rug(data = data, aes(x= noise)) +
#
#   labs(x ='Ambient Noise (mV)', y = 'D50 (m)') +
#   theme_bw() +
#   theme(text = element_text(size = 20))




# Next up ----
m1 <- gam(d50 ~
            s(sst, bs = 'ts') +
            s(mwd, bs = 'ts') +
            s(wspd, bs = 'ts') +
            s(apd, bs = 'ts') +
            s(array, bs = 're'),
          data = data,
          family = Gamma,
          method = 'REML')

m2 <- gam(d50 ~
            s(sst, bs = 'ts') +
            s(wspd, bs = 'ts') +
            s(apd, bs = 'ts') +
            s(array, bs = 're'),
          data = data,
          family = Gamma,
          method = 'REML')



sst_p <- data.frame(sst = seq(min(data$sst),max(data$sst),length =100),
                      wspd = mean(data$wspd, na.rm = T),
                      apd = mean(data$apd, na.rm = T),
                      array = 'Inner')
pred <- predict(m2, sst_p, type = 'link',
                exclude = grep('array', row.names(summary(m2)$s.table),
                               value = T),
                se = T)
pred <- data.frame(var = sst_p$sst,
                   pred = 1/pred$fit,
                   UCI = 1/(pred$fit - 2*pred$se.fit),
                   LCI = 1/(pred$fit + 2*pred$se.fit))
library(ggplot2)
ggplot(data = pred) +
  geom_ribbon(aes(x = var, ymin = LCI, ymax = UCI),
              fill = 'pink') +
  # geom_point(data =data, aes(x = sst, y = d50)) +
  geom_line(aes(x = var, y = pred), color = 'red', lwd = 1) +
  geom_rug(data = data, aes(x = sst)) +
  labs(x ='Sea Surface Temperature (°C)', y = 'D50 (m)') +
  theme_bw() +
  theme(text = element_text(size = 20))

wspd_p <- data.frame(sst = mean(data$sst, na.rm = T),
                    wspd = seq(min(data$wspd, na.rm = T),max(data$wspd, na.rm = T),length =100),
                    apd = mean(data$apd, na.rm = T),
                    array = 'Inner')
pred <- predict(m2, wspd_p, type = 'link',
                exclude = grep('array', row.names(summary(m2)$s.table),
                               value = T),
                se = T)
pred <- data.frame(var = wspd_p$wspd,
                   pred = 1/pred$fit,
                   UCI = 1/(pred$fit - 2*pred$se.fit),
                   LCI = 1/(pred$fit + 2*pred$se.fit))
library(ggplot2)
ggplot() +
  geom_ribbon(data = pred, aes(x = var, ymin = LCI, ymax = UCI),
              fill = 'pink') +
  # geom_point(data = data, aes(x = wspd, y = d50)) +
  geom_line(data = pred, aes(x = var, y = pred), color = 'red', lwd = 1) +
  geom_rug(data = data, aes(x = wspd)) +
  labs(x ='Wind Speed (m/s)', y = 'D50 (m)') +
  theme_bw() +
  theme(text = element_text(size = 20))

apd_p <- data.frame(sst = mean(data$sst, na.rm = T),
                    wspd = mean(data$wspd, na.rm = T),
                    apd = seq(min(data$apd, na.rm = T),max(data$apd, na.rm = T),length =100),
                    array = 'Inner')
pred <- predict(m2, apd_p, type = 'link',
                exclude = grep('array', row.names(summary(m2)$s.table),
                               value = T),
                se = T)
pred <- data.frame(var = apd_p$apd,
                   pred = 1/pred$fit,
                   UCI = 1/(pred$fit - 2*pred$se.fit),
                   LCI = 1/(pred$fit + 2*pred$se.fit))
library(ggplot2)
ggplot() +
  geom_ribbon(data = pred, aes(x = var, ymin = LCI, ymax = UCI),
              fill = 'pink') +
  # geom_point(data = data, aes(x = apd, y = d50)) +
  geom_line(data = pred, aes(x = var, y = pred), color = 'red', lwd = 1) +
  geom_rug(data = data, aes(x = apd)) +
  labs(x ='Average Wave Period (s)', y = 'D50 (m)') +
  theme_bw() +
  theme(text = element_text(size = 20))


p_mod <- predict(m2, type = 'link', se = T,
                 exclude = grep('array', row.names(summary(m2)$s.table),
                                value = T))
library(dplyr)
p_mod <- cbind(data[complete.cases(
  data[, names(data) %in% c('sst','wspd', 'apd')]
  ),], data.frame(p_mod))
p_mod <- p_mod %>%
  filter(date <= '2018-09-30') %>%
  # This calculates the lower and upper CIs and transforms them to the response scale
  mutate(Predicted = 1/fit,
         LCI = 1/(fit + 2 * se.fit),
         UCI = 1/(fit - 2 * se.fit))

p <- p_mod %>%
  select(date, array, Estimated = d50, Predicted, LCI, UCI) %>%
  # use the time period where all predictors are available
  # filter(date <= '2018-08-19') %>%
  tidyr::gather(key = 'key', value = 'val', -date, -array, -LCI, -UCI)

library(ggplot2)
ggplot() +
  geom_ribbon(data = p[p$key == 'Predicted',],
              aes(x = date, ymin = LCI, ymax =  UCI),
              fill = 'gray', alpha = 0.5) +
  geom_line(data = p[p$key %in% c('Estimated', 'Predicted'),],
            aes(x = date, y = val, linetype = key)) +
  scale_linetype_manual(values=c("dashed", "solid")) +
  labs(x=NULL, y = 'D50 (m)', linetype = NULL) +
  facet_wrap(~array, nrow = 2) +
  theme_bw()+
  theme(text = element_text(size = 20),
        legend.position = c(0.2, 0.9))

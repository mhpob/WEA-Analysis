library(dplyr)
library(mgcv)

data <- readRDS('data and imports/rangetest_logit_binary.RDS') %>%
  mutate(array = as.factor(array),
         date_f = as.factor(date)) %>%
  rename(noise = `Average noise`,
         tilt = `Tilt angle`,
         bwt = `Average temperature`)

m <- gam(cbind(success, fail) ~ distance +
            s(array, bs = 're'),
          data = data,
          family = 'binomial',
          method = 'ML', gamma = 1.4)

m_n <- gam(cbind(success, fail) ~ distance +
           s(array, bs = 're') +
           s(noise, bs = 'tp'),
         data = data,
         family = 'binomial',
         method = 'ML', gamma = 1.4)

m_t <- gam(cbind(success, fail) ~ distance +
            s(array, bs = 're') +
            s(dt, bs = 'tp'),
          data = data,
          family = 'binomial',
          method = 'ML', gamma = 1.4)

m_nt <- gam(cbind(success, fail) ~ distance +
            s(array, bs = 're') +
            s(dt, bs = 'tp') +
            s(noise, bs = 'tp'),
          data = data,
          family = 'binomial', gamma = 2,
          method = 'REML')
summary(m_nt)
plot(m_nt, select = 3, trans = m_nt$family$linkinv)

m_nr <- gam(cbind(success, fail) ~ distance +
              s(array, bs = 're') +
              s(dt, bs = 'tp') +
              s(noise, bs = 're'),
            data = data,
            family = 'binomial',
            method = 'ML')

plot_smooth(m_nr, 'dt', transform = m_nt$family$linkinv,
            cond = list(distance = 900),
            rug = T, ylab = 'Detection Efficiency',
            ylim = c(0, 1))
plot_smooth(m_nr, 'distance', transform = m_nt$family$linkinv, rm.ranef = T,
            cond = list(dt = 3),
            rug = T, ylab = 'Detection Efficiency',
            ylim = c(0, 1))

library(itsadug)

plot_smooth(m_nt, 'distance', transform = m_nt$family$linkinv, rm.ranef = T,
            cond = list(dt = 0, noise = 230),
            rug = T, xlab = 'Distance (m)', ylab = 'Detection Efficiency',
            ylim = c(0, 1), xlim = c(0, 1500))


plot_smooth(m_nt, 'distance', transform = m_nt$family$linkinv, rm.ranef = T,
            cond = list(dt = 3, noise = 230),
            rug = T, xlab = 'Distance (m)', ylab = 'Detection Efficiency',
            ylim = c(0, 1), xlim = c(0, 1500))

plot_smooth(m_nt, 'distance', transform = m_nt$family$linkinv, rm.ranef = T,
            cond = list(dt = 5, noise = 230),
            rug = T, xlab = 'Distance (m)', ylab = 'Detection Efficiency',
            ylim = c(0, 1), xlim = c(0, 1500))



plot_smooth(m_nt, 'noise', transform = m_nt$family$linkinv, rm.ranef = T,
            cond = list(distance = 800, dt = 0),
            rug = T, xlab = 'Ambient noise (mV)', ylab = 'Detection Efficiency',
            ylim = c(0, 1))
plot_smooth(m_nt, 'noise', transform = m_nt$family$linkinv, rm.ranef = T,
            cond = list(distance = 800, dt = 3),
            rug = T, xlab = 'Ambient noise (mV)', ylab = 'Detection Efficiency',
            ylim = c(0, 1))
plot_smooth(m_nt, 'noise', transform = m_nt$family$linkinv, rm.ranef = T,
            cond = list(distance = 800, dt = 5),
            rug = T, xlab = 'Ambient noise (mV)', ylab = 'Detection Efficiency',
            ylim = c(0, 1))

plot_smooth(m_nt, 'dt', transform = m_nt$family$linkinv, rm.ranef = T,
            cond = list(distance = 600, noise = 230),
            rug = T, xlab = 'dT (C)', ylab = 'Detection Efficiency',
            ylim = c(0, 1), v0 = -1.65, h0 = 0.5)


plot_smooth(m_nt, 'dt', transform = m_nt$family$linkinv, rm.ranef = T,
            cond = list(distance = 1000, noise = 230),
            rug = T, xlab = 'dT (C)', ylab = 'Detection Efficiency',
            ylim = c(0, 1), v0 = 2.2)


pred_data <- data.frame(distance = seq(0, 1500, 10),
                        noise = 230,
                        dt = 3,
                        array = 'Inner')
pred_se <- predict(m_nt, pred_data,
                   exclude = grep('array', row.names(summary(m_nt)$s.table),
                                  value = T),
                   type = 'link', se = T)
pred <- data.frame(distance = pred_data$distance,
                   pred = m_nt$family$linkinv(pred_se$fit),
                   uci = m_nt$family$linkinv(pred_se$fit + 1.96*pred_se$se.fit),
                   lci = m_nt$family$linkinv(pred_se$fit - 1.96*pred_se$se.fit))
library(ggplot2)
ggplot(data = pred) +
  geom_ribbon(aes(x = distance, ymax = uci, ymin = lci), fill = 'gray') +
  geom_line(aes(x = distance, y = pred)) +
  xlim(0, 1500) + ylim(0,1) + theme_bw()







#### More models!!!
####
library(dplyr)
library(mgcv)

data <- readRDS('data and imports/rangetest_logit_binary.RDS') %>%
  mutate(array = as.factor(array),
         date_f = as.factor(date),
         cp_phase = ifelse(lubridate::month(date) %in% 5:8, 'p', 'a')) %>%
  rename(noise = `Average noise`,
         tilt = `Tilt angle`,
         bwt = `Average temperature`)


m <- gam(cbind(success, fail) ~ distance + array +
           ti(dt, bs = 'cs', k = 5) +
           ti(noise, bs = 'cs', k = 5) ,
         data = data[data$distance != 0, ],
         family = 'binomial',
         method = 'REML')

m2 <- gam(cbind(success, fail) ~ distance + array +
            ti(dt, noise),
          data = data[data$distance != 0, ],
          family = 'binomial',
          method = 'REML')

m3 <- gam(cbind(success, fail) ~ distance + array +
            ti(dt, bs = 'cs', k = 5) +
            ti(noise, bs = 'cs', k = 5)+
            ti(dt, noise, k = 3),
          data = data[data$distance != 0, ],
          family = 'binomial',
          method = 'REML')


m4 <- gam(cbind(success, fail) ~ distance + array +
            te(dt, noise),
          data = data[data$distance != 0, ],
          family = 'binomial',
          method = 'REML')

# m3 and m4 are the same. m1 is the same as using s(). m2 is maybe just trash

library(itsadug)

fvisgam(m3, view = c('dt', 'noise'), transform = m_nt$family$linkinv,
        cond = list(distance = 600, array = 'MD WEA'))
pvisgam(m,view = c('dt', 'noise'))
ggplot()+
  geom_point(data= data[data$distance == 800, ],
             aes(x = dt, y = noise, color = (success/(success+fail))))+
  scale_color_viridis_c()

plot_smooth(m3, 'noise', transform = m_nt$family$linkinv,
            cond = list(distance = 600, dt = 0, array = 'MD WEA'),
            rug = T, xlab = 'Ambient noise (mV)', ylab = 'Detection Efficiency',
            ylim = c(0, 1))
plot_smooth(m3, 'dt', transform = m_nt$family$linkinv,
            cond = list(distance = 600, noise = 200),
            rug = T, xlab = 'dT (C)', ylab = 'Detection Efficiency',
            ylim = c(0, 1), v0 = 2.2)


pred_data <- data %>%
  distinct(date_f, array, .keep_all = T) %>%
  filter(array == 'MD WEA') %>%
  mutate(distance = 0)


lpm <- predict(m, pred_data, type = 'lpmatrix')
d50 <- (log10(0.5/(1-0.5)) -(lpm[,-2] %*% coef(m)[-2])) / coef(m)[2]


d95 <- (log10(0.95/(1-0.95)) -(lpm[,-2] %*% coef(m)[-2])) / coef(m)[2]


library(ggplot2); library(dplyr); library(mgcv)
data <- readRDS('data and imports/rangetest_no_outliers.RDS')

## Coding different random effects ----
# Random slope
mod <- gam(d50_adj ~ s(average_temperature, bs = 'ts') +
             s(average_temperature, array, bs = 're') +
             s(average_noise, bs = 'ts') +
             s(average_noise, array, bs = 're')+
             s(tilt_angle, bs = 'ts'),
           family = Gamma(),
           method = 'REML',
           data = data, verbosePQL = F)

# Random slope and intercept
mod2 <- gam(d50_adj ~ s(average_temperature, bs = 'ts') +
              s(average_temperature, array, bs = 're') +
              s(average_noise, bs = 'ts') +
              s(average_noise, array, bs = 're')+
              s(tilt_angle, bs = 'ts') +
              s(array, bs = 're'),
            family = Gamma(),
            method = 'REML',
            data = data, verbosePQL = F)

# Random intercept
mod3 <- gam(d50_adj ~ s(average_temperature, bs = 'ts') +
              s(average_noise, bs = 'ts') +
              s(tilt_angle, bs = 'ts') +
              s(array, bs = 're'),
            family = Gamma(),
            method = 'REML',
            data = data, verbosePQL = F)

# Random smooth
mod4 <- gam(d50_adj ~ s(average_temperature, bs = 'ts') +
              s(average_temperature, array, bs = 'fs', m = 1) +
              s(average_noise, bs = 'ts') +
              s(average_noise, array, bs = 'fs', m  =1)+
              s(tilt_angle, bs = 'ts'),
            family = Gamma(),
            method = 'REML',
            data = data, verbosePQL = F)


# Predictions without random effects ----
pred_ts <- function(model){
  p_mod <- predict(model, type = 'response', se = T,
                   # this removes random effects (anything with 'array')
                   exclude = grep('array', row.names(summary(model)$s.table),
                                  value = T))
  p_mod <- cbind(data, p_mod)

  p <- p_mod %>%
    select(date, array, Observed = d50_adj, Estimated = fit, se.fit) %>%
    tidyr::gather(key = 'key', value = 'val', -date, -array, -se.fit)

  ggplot() +
    geom_ribbon(data = p[p$key == 'Estimated',],
                aes(x = date, ymin = val - 2*se.fit, ymax =  val + 2*se.fit),
                fill = 'gray', alpha = 0.5) +
    geom_line(data = p[p$key %in% c('Observed', 'Estimated'),],
              aes(x = date, y = val, linetype = key)) +
    labs(x=NULL, y = 'Estimated D50', linetype = NULL) +
    facet_wrap(~array, nrow = 2) +
    theme_bw()
}

# random slope
pred_ts(mod)
# random slope and intercept
pred_ts(mod2)
# random intercept
pred_ts(mod3)
#random smooth
pred_ts(mod4)


## Marginal effects ----
margplot <- function(predicted, variable, ...){
plot(variable, predicted$fit, type = 'l',
     ylim = c(min(predicted$fit - 2 * predicted$se.fit),
              max(predicted$fit + 2 * predicted$se.fit)),
     ...)
lines(variable, predicted$fit + 2 * predicted$se.fit, lty = 2)
lines(variable, predicted$fit - 2 * predicted$se.fit, lty = 2)
}
# BWT
bwt_inner <- data.frame(average_temperature = seq(0,
                                                max(data$average_temperature),
                                                length = 100),
                      average_noise = median(data$average_noise),
                      tilt_angle = median(data$tilt_angle),
                      array = 'Inner')
bwt_pred_rmranef <- predict(mod2, bwt_inner, type = 'response', se.fit = T,
                    exclude = c('s(average_temperature,array)',
                                's(average_noise,array)',
                                's(array)'))
bwt_pred_inner <- predict(mod2, bwt_inner, type = 'response', se.fit = T)
bwt_wea <- data.frame(average_temperature = seq(0,
                                                  max(data$average_temperature),
                                                  length = 100),
                        average_noise = median(data$average_noise),
                        tilt_angle = median(data$tilt_angle),
                        array = 'MD WEA')
bwt_pred_wea <- predict(mod2, bwt_wea, type = 'response', se.fit = T)

margplot(bwt_pred_rmranef, bwt_inner$average_temperature,
         xlab = 'BWT (°C)', ylab = 'Marginal effect on D50')
rug(data$average_temperature)
lines(bwt_inner$average_temperature, bwt_pred_inner$fit, col = 'red')
lines(bwt_inner$average_temperature,
      bwt_pred_inner$fit + 2*bwt_pred_inner$se.fit, col = 'red', lty = 2)
lines(bwt_inner$average_temperature,
      bwt_pred_inner$fit - 2*bwt_pred_inner$se.fit, col = 'red', lty = 2)
lines(bwt_wea$average_temperature, bwt_pred_wea$fit, col = 'blue')
lines(bwt_wea$average_temperature,
      bwt_pred_wea$fit + 2*bwt_pred_wea$se.fit, col = 'blue', lty = 2)
lines(bwt_wea$average_temperature,
      bwt_pred_wea$fit - 2*bwt_pred_wea$se.fit, col = 'blue', lty = 2)

#Noise
noise_inner <- data.frame(average_temperature = mean(data$average_temperature),
                        average_noise = seq(min(data$average_noise),
                                            max(data$average_noise),
                                            length = 100),
                        tilt_angle = median(data$tilt_angle),
                        array = 'Inner')
noise_pred_rmranef <- predict(mod, noise_inner, type = 'response', se.fit = T,
                            exclude = c('s(average_temperature,array)',
                                        's(average_noise,array)'))
noise_pred_inner <- predict(mod, noise_inner, type = 'response', se.fit = T)
noise_wea <- data.frame(average_temperature = mean(data$average_temperature),
                        average_noise = seq(min(data$average_noise),
                                            max(data$average_noise),
                                            length = 100),
                      tilt_angle = median(data$tilt_angle),
                      array = 'MD WEA')
noise_pred_wea <- predict(mod, noise_wea, type = 'response', se.fit = T)

margplot(noise_pred_rmranef, noise_inner$average_noise,
         xlab = 'Noise (mV)', ylab = 'Marginal effect on D50')
rug(data$average_noise)
lines(noise_inner$average_noise, noise_pred_inner$fit, col = 'red')
lines(noise_inner$average_noise,
      noise_pred_inner$fit + 2*noise_pred_inner$se.fit, col = 'red', lty = 2)
lines(noise_inner$average_noise,
      noise_pred_inner$fit - 2*noise_pred_inner$se.fit, col = 'red', lty = 2)
lines(noise_wea$average_noise, noise_pred_wea$fit, col = 'blue')
lines(noise_wea$average_noise,
      noise_pred_wea$fit + 2*noise_pred_wea$se.fit, col = 'blue', lty = 2)
lines(noise_wea$average_noise,
      noise_pred_wea$fit - 2*noise_pred_wea$se.fit, col = 'blue', lty = 2)

#tilt
tilt_inner <- data.frame(average_temperature = mean(data$average_temperature),
                          average_noise = median(data$average_noise),
                          tilt_angle = seq(min(data$tilt_angle),
                                           max(data$tilt_angle),
                                           length = 100),
                          array = 'Inner')
tilt_pred_rmranef <- predict(mod, tilt_inner, type = 'response', se.fit = T,
                              exclude = c('s(average_temperature,array)',
                                          's(average_noise,array)'))
tilt_pred_inner <- predict(mod, tilt_inner, type = 'response', se.fit = T)
tilt_wea <- data.frame(average_temperature = mean(data$average_temperature),
                       average_noise = median(data$average_noise),
                       tilt_angle = seq(min(data$tilt_angle),
                                        max(data$tilt_angle),
                                        length = 100),
                        array = 'MD WEA')
tilt_pred_wea <- predict(mod, tilt_wea, type = 'response', se.fit = T)

margplot(tilt_pred_rmranef, tilt_inner$tilt_angle,
         xlab = 'Receiver tilt (°)', ylab = 'Marginal effect on D50')
rug(data$tilt_angle)
lines(tilt_inner$tilt_angle, tilt_pred_inner$fit, col = 'red')
lines(tilt_inner$tilt_angle,
      tilt_pred_inner$fit + 2*tilt_pred_inner$se.fit, col = 'red', lty = 2)
lines(tilt_inner$tilt_angle,
      tilt_pred_inner$fit - 2*tilt_pred_inner$se.fit, col = 'red', lty = 2)
lines(tilt_wea$tilt_angle, tilt_pred_wea$fit, col = 'blue')
lines(tilt_wea$tilt_angle,
      tilt_pred_wea$fit + 2*tilt_pred_wea$se.fit, col = 'blue', lty = 2)
lines(tilt_wea$tilt_angle,
      tilt_pred_wea$fit - 2*tilt_pred_wea$se.fit, col = 'blue', lty = 2)





m2 <- gam(d50_adj ~ s(average_temperature, bs = 'ts') +
             s(average_temperature, array, bs = 're') +
             s(average_noise, bs = 'ts') +
             s(average_noise, array, bs = 're')+
             s(tilt_angle, bs = 'ts') +
            s(tilt_angle, array, bs = 're'),
           family = Gamma(),
           method = 'REML',
           data = data, verbosePQL = F, na.action = 'na.fail')

library(MuMIn)

j <- parallel::makeCluster(parallel::detectCores()-2)
clusterEvalQ(j, library(mgcv))
clusterExport(j, "data")

k <- pdredge(m2, j)
parallel::stopCluster(j)
par(mar = c(3,5,6,4))
plot(k, labAsExpr = TRUE)

model.sel(k)






# w/satellite variables
data2 <- data[complete.cases(data[, c(1, 2, 12:16, 18:27)]),]
library(vegan)
pca <- rda(data2[data2$array == 'MD WEA', c(12:16, 18:27, 31),],
           scale = T)
biplot(pca, type = c('text', 'points'))

mod <- gam(d50_adj ~ s(average_temperature, bs = 'ts') +
             s(average_temperature, array, bs = 're') +
             s(average_noise, bs = 'ts') +
             s(average_noise, array, bs = 're')+
             s(tilt_angle, bs = 'ts'),
           family = Gamma(),
           method = 'REML',
           data = data, verbosePQL = F)










## Ranefff specified differently
library(mgcv)
data <- readRDS('data and imports/rangetest_no_outliers.RDS')


mod3 <- gam(d50_adj ~ s(average_temperature, bs = 'ts') +
             s(average_noise, bs = 'ts') +
             s(tilt_angle, bs = 'ts') +
             s(array, bs = 're'),
           family = Gamma(),
           method = 'REML',
           data = data, verbosePQL = F, na.action = 'na.fail')

anova(mod, mod2, test ='Chisq')

# Take the random effects out.
prediction <- predict(mod3, type = 'response', se = T,
                      exclude = c('s(array)'))
prediction <- cbind(data, prediction)


library(MuMIn)

j <- parallel::makeCluster(parallel::detectCores()-2)
parallel::clusterEvalQ(j, library(mgcv))
parallel::clusterExport(j, "data")

k <- pdredge(mod2, j)
parallel::stopCluster(j)
par(mar = c(3,5,6,4))
plot(k, labAsExpr = TRUE)

model.sel(k)

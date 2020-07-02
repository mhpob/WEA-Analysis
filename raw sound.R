library(mgcv)
data <- readRDS('data and imports/rangetest_logit_data.rds')
names(data)[c(12:14)] <- c('noise', 'bwt', 'tilt')
data$array <- as.factor(data$array)
j <- gam(D50 ~  s(dt, by = array),
         data = data, method = 'REML', family = Gamma(link = log))
summary(j)


par(mfrow=c(2,2))
plot(j,residuals=TRUE,pch=19)


k <- gamm(D50 ~ s(bwt) + s(dt), random = list(array = ~1),
          data = data, method = 'REML', family = Gamma(link = log))


absorp <- function(f, t = 10, sal = 35, ph = 8, depth = 0.02){
  # Van Moll et al. 2009; absorption in dB/km
  # t = Temperature (C)
  # sal = salinity (g/kg)
  # depth = depth (km)
  f_1 <- 0.91 * sqrt(sal/35) * exp(t/33)
  boric <- 0.101 * (f_1 * f^2) / (f_1^2 + f^2) * exp((ph - 8) / 0.57)

  f_2 <- 46.6 * exp(t / 18)
  magsul <- 0.56 * (1 + t / 76) * (sal / 35) *
    ((f_2 * f^2) / (f_2^2 + f^2)) * exp(- depth / 4.9)

  A_3 <- ifelse(t <= 20,
                4.937e-4 - 2.59e-5 * t + 9.11e-7 * t^2 - 1.5e-8 * t^3,
                3.964e-4 - 1.146e-5 * t + 1.45e-7 * t^2 - 6.5e-10 * t^3)
  P_3 <- 1 - 3.83e-2 * depth + 4.9e-4 * depth^2
  fresh <- A_3 * P_3 * f^2

  boric + magsul + fresh
}

plot(x = seq(0, 30, 1), y = absorp(f = 69, t = seq(0, 30, 1)),
     xlab = 'Temperature (C)', ylab = 'Absorption (dB/km)', type = 'l')
plot(x = seq(0, 30, 1), y = absorp(f = 0.8, t = seq(0, 30, 1)),
     xlab = 'Temperature (C)', ylab = 'Absorption (dB/km)', type = 'l')

transloss <- function(r, absorp){
  (20 * log10(r)) + (absorp * r)
}

tl_plot <- function(f){
plot(seq(100, 1000, 100), transloss(seq(100, 1000, 100), absorp(f, 30)/1000),
     type = 'l', col = 'red', xlab = 'Radius (m)', ylab = 'Transmission loss (dB)')
lines(seq(100, 1000, 100), transloss(seq(100, 1000, 100), absorp(f, 10)/1000),
      col = 'blue')
lines(seq(100, 1000, 100), transloss(seq(100, 1000, 100), absorp(f, 20)/1000))
lines(seq(100, 1000, 100), 20*log10(seq(100, 1000, 100)), lty = 2)
}

tl_plot(0.8)
tl_plot(69)

k <- expand.grid(r = seq(10, 1500, 10), bwt = seq(0, 27, 0.5))
k$TL <- transloss(k$r, absorp = absorp(69, t = k$bwt)/1000)


ggplot() + geom_tile(data = k, aes(x = r, y = bwt, fill = TL)) +
  geom_contour(data = k, aes(x = r, y = bwt, z = TL),
               col = 'black', binwidth = 5) +
  scale_fill_viridis_c(limits = c(20, 100)) +
  labs(x = 'Distance (m)', y = 'Temperature (C)', fill = 'TL (dB)', color = 'Array') +
  # geom_point(data = data, aes(x = D95, y = sst), size = 3, color = 'red') +
  geom_point(data = data, aes(x = D50, y = bwt, color = array), size = 3) +
  xlim(0, 1500) +
  theme_bw()

data$absorp <- absorp(69, data$bwt)/1000
data$tloss_d50 <- transloss(data$D50, data$absorp)

# Theoretical maximum detectable distance
k <- expand.grid(r = seq(10, 5000, 10), bwt = seq(0, 30, 0.5))
k$TL <- transloss(k$r, absorp = absorp(69, t = k$bwt)/1000)
k <- k[k$TL >= 137.22,]
library(dplyr)
k <- group_by(k, bwt) %>% summarize(min = min(r))
plot(k$min, k$bwt, type = 'l',
     xlab = 'Maximum detectable distance (m)', ylab = 'Temperature (C)')
hist(data[data$bwt > 5 && data$bwt < 15, 'D5'], main = NULL, xlab = 'D5 between 5 and 15 C')

s_vel <- function(t, sal, depth){
  1449.2 +
    4.6 * t - 0.055 * t ^ 2 + 0.00029 * t ^ 3 +
    (1.34 - 0.01 * t) * (sal - 35) +
    0.016 * depth
}

data$s_vel <- ifelse(data$array == 'Inner',
                     s_vel(data$bwt, 32, 17),
                     s_vel(data$bwt, 32, 27))



L <-  32
D <-  65

h <- sqrt((1/8) * (D + L)) * 914.4 #in meters

ggplot() + geom_line(data = data, aes(x = date, y = absorp, color = array))
ggplot() + geom_line(data = data, aes(x = date, y = bwt, color = array))

ggplot() + geom_line(data = data, aes(x = date, y = tloss_d50, color = array))
ggplot() + geom_line(data = data, aes(x = date, y = D50, color = array))

ggplot(data = data, aes(x = absorp, y = D50)) + geom_point(aes(color = array)) +
  geom_smooth( method = 'glm', method.args = list(family = Gamma(link = 'log')))

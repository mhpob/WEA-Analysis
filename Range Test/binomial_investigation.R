library(dplyr)
library(gamm4)

data <- readRDS('data and imports/rangetest_logit_binary.RDS') %>%
  mutate(array = as.factor(array),
         date_f = as.factor(date)) %>%
  rename(noise = `Average noise`,
         tilt = `Tilt angle`,
         bwt = `Average temperature`)



# distance grouped by date ----
m1 <- gamm4(cbind(success, fail) ~ distance,
            random = ~(1 | date_f),
            data = data, family = 'binomial')

# Plot the per-date (random intercept) D50
# Add per-date intercept terms together
head(as.matrix(coef(m1$mer)[[1]]))
plot.ts(-(coef(m1$mer)[[1]][1] + coef(m1$mer)[[1]][2]) / coef(m1$mer)[[1]][3])


# Try random slope/intercept with distance, random array intercept ----
m2 <- gamm4(cbind(success, fail) ~ distance,
            random = ~(distance | date_f) + (1 | array),
            data = data, family = 'binomial')

c.mat <- coef(m2$mer)[[1]]
# Add per-day intercepts
c.mat[, '(Intercept)'] <-  rowSums(c.mat[, grep('Intercept', names(c.mat), value = T)])
# Add per-day distance slopes
c.mat[, 'distance'] <-  rowSums(c.mat[, grep('distance', names(c.mat), value = T)])
# drop unused terms
c.mat <- c.mat[c('(Intercept)', 'distance')]
# Join in array intercepts
c.mat <- list(inner = c.mat,
              wea = c.mat)
c.mat[[1]][, 1] <- c.mat[[1]][,1] + coef(m2$mer)[[2]][1, 1]
c.mat[[2]][, 1] <- c.mat[[1]][,1] + coef(m2$mer)[[2]][2, 1]

# Linear predictor matrix
pred_data <- data %>%
  distinct(date_f, array, .keep_all = T) %>%
  mutate(distance = 1,
         date = as.Date(date_f))
lpm <- predict(m2$gam, pred_data[pred_data$array == 'Inner',],
               type = 'lpmatrix')

# Calculate per-day D50
d50.inn <- NULL
for(i in 1:349){
  d50.inn[i] <- -(sum(lpm[i,-2] * c.mat[[1]][i,-2])) / c.mat[[1]][i,2]
}
# Plot
plot.ts(d50.inn)


lpm <- predict(m2$gam, pred_data[pred_data$array == 'MD WEA',],
               type = 'lpmatrix')

# Calculate per-day D50
d50.wea <- NULL
for(i in 1:349){
  d50.wea[i] <- -(sum(lpm[i,-2] * c.mat[[2]][i,-2])) / c.mat[[2]][i,2]
}
# Plot
plot.ts(d50.wea)





# Try with a smoother and a random slope/intercept with distance, random array intercept ----
m3 <- gamm4(cbind(success, fail) ~ distance +
              s(dt, bs = 'ts'),
            random = ~(distance | date_f/array),
            data = data, family = 'binomial')

# Three sets of coef: one for random date, one for the smoother, one for array
names(coef(m3$mer))

head(coef(m3$mer)[[1]])

# And the overall coefficients are here
coef(m3$gam)


c.mat <- coef(m3$mer)[[1]]

# Add per-day intercepts
c.mat[, '(Intercept)'] <-  rowSums(c.mat[, grep('Intercept', names(c.mat), value = T)])
# Add per-day distance slopes
c.mat[, 'distance'] <-  rowSums(c.mat[, grep('distance', names(c.mat), value = T)])
# drop unused terms
c.mat <- c.mat[c('(Intercept)', 'distance', 's(dt)')]
# Join in smoother coefs
c.mat <- cbind(c.mat[1:2], matrix(rep(coef(m3$gam)[3:11], each = 349), nrow = 349))


# Linear predictor matrix
pred_data <- data %>%
distinct(date_f, array, .keep_all = T) %>%
  mutate(distance = 1,
         date = as.Date(date_f))
lpm <- predict(m3$gam, pred_data[pred_data$array == 'MD WEA'],
               type = 'lpmatrix')


# Calculate per-day D50
d50 <- NULL
for(i in 1:349){
  d50[i] <- -(sum(lpm[i,-2] * c.mat[i,-2])) / c.mat[i,2]
}
# Plot
plot.ts(d50)




# Try with a random slope of distance for date:array combos. No random intercept ----
m4 <- gamm4(cbind(success, fail) ~ distance,
            random = ~(distance - 1| date_f:array),
            data = data, family = 'binomial')

names(coef(m4$mer))
k <- coef(m4$mer)[[1]][grepl('Inner', row.names(coef(m4$mer)[[1]])),]
plot.ts(-k[,2]/(k[,1]+k[,3]))
k <- coef(m4$mer)[[1]][grepl('MD WEA', row.names(coef(m4$mer)[[1]])),]
plot.ts(-k[,2]/(k[,1]+k[,3]))



# Try with a random slope of distance for date:array combos and uncorrelated random intercept ----
m5 <- gamm4(cbind(success, fail) ~ distance,
            random = ~ (1 | date_f:array) + (distance - 1| date_f:array),
            data = data, family = 'binomial')

names(coef(m5$mer)[[1]])
k <- coef(m5$mer)[[1]][grepl('Inner', row.names(coef(m5$mer)[[1]])),]
plot.ts(-(k[,1] + k[, 3]) / (k[, 2]+k[, 4]))
k <- coef(m5$mer)[[1]][grepl('MD WEA', row.names(coef(m5$mer)[[1]])),]
plot.ts(-(k[,1] + k[, 3]) / (k[, 2]+k[, 4]))
# Big issues here.



# Try with a random slope of distance for date:array combos, random array intercept ----
m6 <- gamm4(cbind(success, fail) ~ distance,
            random = ~ (1 | array) + (distance - 1| date_f:array),
            data = data, family = 'binomial')
shell.exec('https://media.giphy.com/media/o2La4Pvf9CdJC/giphy.gif')


names(coef(m6$mer))
k <- coef(m6$mer)[[1]][grepl('Inner', row.names(coef(m6$mer)[[1]])),]
plot.ts(-(k[, 'X(Intercept)'] + coef(m6$mer)[[2]][1, '(Intercept)']) /
          (k[, 'distance'] + k[, 'Xdistance']))
k <- coef(m6$mer)[[1]][grepl('MD WEA', row.names(coef(m6$mer)[[1]])),]
plot.ts(-(k[, 'X(Intercept)'] + coef(m6$mer)[[2]][1, '(Intercept)']) /
          (k[, 'distance'] + k[, 'Xdistance']))




m6.1 <- gamm4(cbind(success, fail) ~ distance +
                s(dt, bs = 'ts') + s(),
              random = ~ (1 | array) + (distance - 1| date_f:array),
              data = data, family = 'binomial')
shell.exec('https://media.giphy.com/media/o2La4Pvf9CdJC/giphy.gif')



# Add per-day distance slopes
c.mat <- ranef(m6.1$mer)[['date_f:array']] + fixef(m6.1$mer)['Xdistance']
c.mat[, '(Intercept)'] <- ifelse(grepl('Inner', row.names(c.mat)),
                                 ranef(m6.1$mer)[[3]]['Inner',] +
                                   fixef(m6.1$mer)['X(Intercept)'],
                                 ranef(m6.1$mer)[[3]]['MD WEA',] +
                                   fixef(m6.1$mer)['X(Intercept)'])

# Join in smoother coefs
c.mat <- cbind(c.mat[, c(2,1)], matrix(rep(coef(m6.1$gam)[grepl('dt', names(coef(m6.1$gam)))],
                                 each = 698), nrow = 698))

c.mat <- cbind(c.mat[, c(2,1)], matrix(rep(t(ranef(m6.1$mer)[[2]]),
                                           each = 698), nrow = 698))



pred_data <- data %>%
  distinct(date_f, array, .keep_all = T) %>%
  filter(array == 'Inner') %>%
  mutate(distance = 250)
lpm <- predict(m6.1$gam, pred_data,
               type = 'lpmatrix')
cm <- c.mat[grepl('Inner', row.names(c.mat)),]


d50 <- NULL
for(i in 1:349){
  d50[i] <- -(sum(lpm[i,-2] * cm[i,-2])) / cm[i,2]
}
# Plot
plot.ts(d50)




### Newest ----
library(dplyr)
library(gamm4)

data <- readRDS('data and imports/rangetest_logit_binary.RDS') %>%
  mutate(array = as.factor(array),
         date_f = as.factor(date),
         date_n = as.numeric(date_f)) %>%
  rename(noise = `Average noise`,
         tilt = `Tilt angle`,
         bwt = `Average temperature`)
mod_data <- data[, names(data) %in%
                   c('date_f', 'array', 'distance', 'success', 'fail', 'noise',
                     'bwt', 'tilt', 'dt', 'wdir', 'wspd', 'wvht', 'dpd', 'apd',
                     'mwd', 'pres')]
mod_data <- mod_data[complete.cases(mod_data),]
email <- gmailr::mime(To = 'obrien@umces.edu',
                     From = 'obrien@umces.edu',
                     Subject = 'Model run finished')


mod <- gamm4(cbind(success, fail) ~ distance +
               s(dt) +
               s(noise) +
               s(tilt) +
               s(wdir)+
               s(wspd) +
               s(wvht) +
               s(dpd) +
               s(apd) +
               s(mwd) +
               s(pres),
             random = ~ (1 | array) + (0 + distance | date_f:array),
             data = mod_data, family = 'binomial')

gmailr::send_message(email)

# library(mgcv)
system.time(
m <- mgcv::gam(cbind(success, fail) ~ distance +
                 s(distance, date_f, array, bs = 're') +
                 s(array, bs = 're') +
                 s(dt) +
                 s(noise) +
                 s(tilt) +
                 s(wdir)+
                 s(wspd) +
                 s(wvht) +
                 s(dpd) +
                 s(apd) +
                 s(mwd) +
                 s(pres),
         data = mod_data, family = 'binomial')
)
system.time(
m2 <- gamm4::gamm4(cbind(success, fail) ~ distance +
                     s(dt) +
                     s(noise) +
                     s(tilt) +
                     s(wdir)+
                     s(wspd) +
                     s(wvht) +
                     s(dpd) +
                     s(apd) +
                     s(mwd) +
                     s(pres),
           random = ~ (0 + distance | date_f:array) + (1 | array),
           data = mod_data, family = 'binomial')
)
gmailr::send_message(email)

saveRDS(m, 'global_gam.rds')
saveRDS(m2, 'global_gamm4.rds')

# options(na.action = "na.fail")
# cl <- parallel::makeCluster(parallel::detectCores(logical = F) - 1)
# parallel::clusterExport(cl, 'mod_data')
# parallel::clusterEvalQ(cl, library(mgcv))
# parallel::clusterEvalQ(cl, options(na.action = "na.fail"))
# all.models <- MuMIn::pdredge(m, cluster = cl, fixed = ~distance)
# parallel::stopcluster(cl)
#
# saveRDS(all.models, 'c:/users/darpa1/desktop/mumin_models.RDS')
models <- readRDS('mumin_models.RDS')

m <- readRDS('global_gam.rds')


library(itsadug)

plot_smooth(m, 'dt')


pred_data <- data %>%
  distinct(date_f, array, .keep_all = T) %>%
  filter(array == 'Inner') %>%
  mutate(distance = 250)
lpm <- predict(m, pred_data,
               type = 'lpmatrix',
               exclude = grep('date_f|array', row.names(summary(m)$s.table),
                              value = T))
cm <- c.mat[grepl('Inner', row.names(c.mat)),]


# OTHER CODE -----


lpm <- predict(m1$gam, pred_data,
               type = 'lpmatrix',
               exclude = grep('date_f', row.names(summary(m1$gam)$s.table),
                              value = T))

 -(lpm[,-2] %*% coef(m1$gam)[-2]) / coef(m1$gam)[2]



ci.pred <- function(n, obj, lpreds){
  # n = number of draws
  # obj = fitted model
  # lpreds = linear predictor matrix of new data

  param_reps <- MASS::mvrnorm(n, coef(obj), vcov(obj))
  estimate_reps <- rep(0, n)

  for (i in 1:n){
    estimate_reps[i] <- -(lpreds[-2] %*% param_reps[i, -2]) / param_reps[i, 2]
  }

  mod_d50 <- -(lpreds[-2] %*% coef(obj)[-2]) / coef(obj)[2]
  sd_estimate <- sqrt(var(estimate_reps))

  out <- c(d50 = mod_d50,
           uci = mod_d50 + 1.96 * sd_estimate,
           lci = mod_d50 - 1.96 * sd_estimate)
  out
}

d50_pred <- t(apply(lpm, 1, function(x) ci.pred(100, m, x)))
pred_data <- cbind(pred_data, d50_pred)

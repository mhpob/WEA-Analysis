# ==================
# MODEL SELECTION
# ==================
# Packages ----
library(data.table); library(mgcv)



# Custom functions ----

## Cross validation
cv <- function(data, model, k, repeats = 1, seed = NULL){

  terms_to_exclude <- grep('site|receiver',
                           row.names(summary(model)$s.table),
                           value = T)

  refit <- function(i, data_shuffle, folds){
    # Segment data by fold using the which() function
    test_ids <- which(folds == i, arr.ind = TRUE)
    test_data <- data_shuffle[test_ids, ]
    train_data <- data_shuffle[-test_ids, ]

    # Train the model.
    CV_mod <- bam(data = train_data,
                  formula = model$formula,
                  family = model$family$family,
                  discrete = T,
                  weights = obs_weights)


    train_data$pred <- predict(CV_mod, type = 'response',
                               exclude = terms_to_exclude)
    # train_data <- train_data[!train_data$distance %in% c(0, 2400),]


    # Test the model
    test_data$pred <- predict(CV_mod, test_data, type = 'response',
                              exclude = c('s(receiver)', 's(noise,site)',
                                          's(dt,site)', 's(bwt,site)'))
    # test_data <- test_data[!test_data$distance %in% c(0, 2400),]


    # Penalty functions
    ## Overall RMSE
    train_overall <- sqrt(mean((train_data$detectability - train_data$pred) ^ 2))
    test_overall <- sqrt(mean((test_data$detectability - test_data$pred) ^ 2))

    ## 250m RMSE
    train_250 <- sqrt(mean((train_data$detectability[train_data$distance == 250] -
                              train_data$pred[train_data$distance == 250]) ^ 2))
    test_250 <- sqrt(mean((test_data$detectability[test_data$distance == 250] -
                             test_data$pred[test_data$distance == 250]) ^ 2))

    ## 550m RMSE
    train_550 <- sqrt(mean((train_data$detectability[train_data$distance == 550] -
                              train_data$pred[train_data$distance == 550]) ^ 2))
    test_550 <- sqrt(mean((test_data$detectability[test_data$distance == 550] -
                             test_data$pred[test_data$distance == 550]) ^ 2))

    ## 800m RMSE
    train_800 <- sqrt(mean((train_data$detectability[train_data$distance == 800] -
                              train_data$pred[train_data$distance == 800]) ^ 2))
    test_800 <- sqrt(mean((test_data$detectability[test_data$distance == 800] -
                             test_data$pred[test_data$distance == 800]) ^ 2))

    c(train_overall = train_overall, test_overall = test_overall,
      train_250 = train_250, test_250 = test_250,
      train_550 = train_550, test_550 = test_550,
      train_800 = train_800, test_800 = test_800)

  }


  if(repeats != 1){
    reps <- function(reps){
      if(!is.null(seed)){
        set.seed(seed * reps)
      }

      data_shuffle <- data[sample(nrow(data)),]
      folds <- cut(seq(1, nrow(data_shuffle)), breaks = k, labels = F)

      folds <- sapply(1:k, refit, data_shuffle, folds)

      data.frame((t(folds)),
                 fold = 1:k,
                 rep = reps)
    }

    cvs <- lapply(1:repeats, reps)
    do.call(rbind, cvs)

  }else{

    data_shuffle <- data[sample(nrow(data)),]
    folds <- cut(seq(1, nrow(data_shuffle)), breaks = k, labels = F)

    folds <- sapply(1:k, refit, data_shuffle, folds)
    data.frame((t(folds)),
               fold = 1:k)

  }

}



# Data ----
data <- fread('range test/manuscript/revisions/additional_file_1.csv')
data[, ':='(

  # Convert character vectors to factors
  site = as.factor(site),
  receiver = as.factor(receiver),

  # "Prior weights on the contribution of the data to the log likelihood."
  #   See ?mgcv::bam
  obs_weights = (success + fail) / mean(success + fail),

  # Note the start of each time series for AR1 modeling
  ts.start = fifelse(date == '2017-12-21' |
                       (date == '2018-08-08' & receiver == 'A' & distance == 800) |
                       (date == '2018-08-08' & receiver == 'B' & distance == 550),
                     T, F))]

## Order data according to each receiver's time series
setorder(data, site, receiver, distance, date)

## Set seed for CV reproducibility; number reflects internal ABIT manuscript number
seed <- 2000031



# Fit global model for rho estimate to fit AR1 model ----
m_global <- bam(detectability ~
                  distance +
                  s(receiver, bs = 're') +
                  s(dt, k = 40, m = 2) +
                  s(dt, site, bs = 'fs', k = 40, m = 2) +
                  s(bwt, k = 40, m = 2) +
                  s(bwt, site, bs = 'fs', k = 40, m = 2) +
                  s(noise, k = 40, m = 2) +
                  s(noise, site, bs = 'fs', k = 40, m = 2) +
                  ti(noise, bwt, k = c(10, 10)) +
                  ti(noise, dt, k = c(10, 10)),
                family = binomial(),
                data = data,
                weights = data$obs_weights,
                discrete = T)

## Guess rho
rho_guess <- acf(resid(m_global), plot = F)$acf[2]



# ~ Distance ----
##  Linear response to distance only (Null model)
m <- bam(detectability ~
           distance +
           s(receiver, bs = 're'),
         family = binomial(),
         data = data,
         weights = data$wgt,
         discrete = T,
         rho = rho_guess,
         AR.start = data$ts.start)

##  Rsq/Dev expl
sapply(summary(m)[c(10,14)], round, 3)

##  5-fold cross-validation
kf_d <- cv(data, m, k = 5, repeats = 5, seed = seed)

##  Summarize the above. Only use RMSE calculated on the test set. If desired,
##    compare RMSE of the test and train fits. If RMSE test >> RMSE train, then
##    the model is over-fit.
colMeans(kf_d[, grepl('test', names(kf_d))]) * 100
apply(kf_d[, grepl('test', names(kf_d))], 2, sd) * 100



# ~ Distance + s(noise) ----
## Nonlinear response to noise
m_n <- bam(detectability ~
             distance +
             s(receiver, bs = 're') +
             s(noise, k = 40, m = 2) +
             s(noise, site, bs = 'fs', k = 40, m = 2),
           family = binomial(),
           data = data,
           weights = data$wgt,
           discrete = T,
           rho = rho_guess,
           AR.start = data$ts.start)


##  Rsq/Dev expl
sapply(summary(m_n)[c(10,14)], round, 3)

##  5-fold cross-validation
kf_n <- cv(data, m_n, k = 5, repeats = 5, seed = seed)
colMeans(kf_n[, grepl('test', names(kf_n))]) * 100
apply(kf_n[, grepl('test', names(kf_n))], 2, sd) * 100



# ~ Distance + s(noise) + s(dt) ----
## Nonlinear response to noise and stratification
m_ndt <- bam(detectability ~
               distance +
               s(receiver, bs = 're') +
               s(noise, k = 40, m = 2) +
               s(noise, site, bs = 'fs', k = 40, m = 2) +
               s(dt, k = 40, m = 2) +
               s(dt, site, bs = 'fs', k = 40, m = 2),
             family = binomial(),
             data = data,
             weights = data$wgt,
             discrete = T,
             rho = rho_guess,
             AR.start = data$ts.start)


##  Rsq/Dev expl
sapply(summary(m_ndt)[c(10,14)], round, 3)

##  5-fold cross-validation
kf_ndt <- cv(data, m_ndt, k = 5, repeats = 5, seed = seed)
colMeans(kf_ndt[, grepl('test', names(kf_ndt))]) * 100
apply(kf_ndt[, grepl('test', names(kf_ndt))], 2, sd) * 100



# ~ Distance + s(noise) + s(dt) + ti(noise, dt)----
## Nonlinear response to noise, modulated by stratification
m_ndti <- bam(detectability ~
                distance +
                s(receiver, bs = 're') +
                s(dt, k = 40, m = 2) +
                s(dt, site, bs = 'fs', k = 40, m = 2) +
                s(noise, k = 40, m = 2) +
                s(noise, site, bs = 'fs', k = 40, m = 2) +
                ti(noise, dt, k = c(10, 10)),
              family = binomial(),
              data = data,
              weights = data$wgt,
              discrete = T,
              rho = rho_guess,
              AR.start = data$ts.start)


##  Rsq/Dev expl
sapply(summary(m_ndti)[c(10,14)], round, 3)

##  5-fold cross-validation
kf_ndti <- cv(data, m_ndti, k = 5, repeats = 5, seed = seed)
colMeans(kf_ndti[, grepl('test', names(kf_ndti))]) * 100
apply(kf_ndti[, grepl('test', names(kf_ndti))], 2, sd) * 100



# ~ Distance + s(dt) ----
## Nonlinear response to stratification
m_dt <- bam(detectability ~
              distance +
              s(receiver, bs = 're') +
              s(dt, k = 40, m = 2) +
              s(dt, site, bs = 'fs', k = 40, m = 2),
            family = binomial(),
            data = data,
            weights = data$wgt,
            discrete = T,
            rho = rho_guess,
            AR.start = data$ts.start)

##  Rsq/Dev expl
sapply(summary(m_dt)[c(10,14)], round, 3)


##  5-fold cross-validation
kf_dt <- cv(data, m_dt, k = 5, repeats = 5, seed = seed)
colMeans(kf_dt[, grepl('test', names(kf_dt))]) * 100
apply(kf_dt[, grepl('test', names(kf_dt))], 2, sd) * 100


# ~ Distance + s(noise) + s(bwt) ----
## Nonlinear response to noise and near-receiver temperature
m_nbt <- bam(detectability ~
               distance +
               s(receiver, bs = 're') +
               s(noise, k = 40, m = 2) +
               s(noise, site, bs = 'fs', k = 40, m = 2) +
               s(bwt, k = 40, m = 2) +
               s(bwt, site, bs = 'fs', k = 40, m = 2),
             family = binomial(),
             data = data,
             weights = data$wgt,
             discrete = T,
             rho = rho_guess,
             AR.start = data$ts.start)

##  Rsq/Dev expl
sapply(summary(m_nbt)[c(10,14)], round, 3)


##  5-fold cross-validation
kf_nbt <- cv(data, m_nbt, k = 5, repeats = 5, seed = seed)
colMeans(kf_nbt[, grepl('test', names(kf_nbt))]) * 100
apply(kf_nbt[, grepl('test', names(kf_nbt))], 2, sd) * 100



# ~ Distance + s(noise) + s(bwt) + ti(noise, bwt)----
## Nonlinear response to noise, modulated by near-receiver temperature
m_nbti <- bam(detectability ~
                distance +
                s(receiver, bs = 're') +
                s(noise, k = 40, m = 2) +
                s(noise, site, bs = 'fs', k = 40, m = 2) +
                s(bwt, k = 40, m = 2) +
                s(bwt, site, bs = 'fs', k = 40, m = 2) +
                ti(noise, bwt, k = c(10, 10)),
              family = binomial(),
              data = data,
              weights = data$wgt,
              discrete = T,
              rho = rho_guess,
              AR.start = data$ts.start)


##  Rsq/Dev expl
sapply(summary(m_nbti)[c(10,14)], round, 3)

##  5-fold cross-validation
kf_nbti <- cv(data, m_nbti, k = 5, repeats = 5, seed = seed)
colMeans(kf_nbti[, grepl('test', names(kf_nbti))]) * 100
apply(kf_nbti[, grepl('test', names(kf_nbti))], 2, sd) * 100



# ~ Distance + s(bwt) ----
## Nonlinear response to near-receiver temperature
m_bt <- bam(detectability ~
              distance +
              s(receiver, bs = 're') +
              s(bwt, k = 40, m = 2) +
              s(bwt, site, bs = 'fs', k = 40, m = 2),
            family = binomial(),
            data = data,
            weights = data$wgt,
            discrete = T,
            rho = rho_guess,
            AR.start = data$ts.start)

##  Rsq/Dev expl
sapply(summary(m_bt)[c(10,14)], round, 3)


##  5-fold cross-validation
kf_bt <- cv(data, m_bt, k = 5, repeats = 5, seed = seed)
colMeans(kf_bt[, grepl('test', names(kf_bt))]) * 100
apply(kf_bt[, grepl('test', names(kf_bt))], 2, sd) * 100



# AIC comparison ---
ic <- data.frame(
  model = c('m', 'm_n', 'm_ndt', 'm_ndti', 'm_dt', 'm_nbt',
            'm_nbti', 'm_bt'),
  AIC = sapply(list(m, m_n, m_ndt, m_ndti,
                    m_dt, m_nbt, m_nbti, m_bt),
               AIC)
)

ic$dAIC <- ic$AIC - min(ic$AIC)
ic$wAIC <- exp(-0.5 * ic$dAIC)
ic$wAIC <- ic$wAIC / sum(ic$wAIC)
ic <- ic[order(ic$dAIC),]
ic





# =================
# D50 CALCULATION
#==================

# Packages ----
library(data.table); library(changepoint); library(mcp); library(mgcv)


# Import data ----
data <- fread('range test/manuscript/revisions/additional_file_1.csv')
data[, ':='(

  # Convert character vectors to factors
  site = as.factor(site),
  receiver = as.factor(receiver),

  # "Prior weights on the contribution of the data to the log likelihood."
  #   See ?mgcv::bam
  obs_weights = (success + fail) / mean(success + fail),

  # Note the start of each time series for AR1 modeling
  ts.start = fifelse(date == '2017-12-21' |
                       (date == '2018-08-08' & receiver == 'A' & distance == 800) |
                       (date == '2018-08-08' & receiver == 'B' & distance == 550),
                     T, F))]

## Order data according to each receiver's time series
setorder(data, site, receiver, distance, date)



# Run winning model ----
##  see above for selection
mod <- bam(detectability ~
             distance +
             s(receiver, bs = 're') +
             s(noise, k = 40, m = 2) +
             s(noise, site, bs = 'fs', k = 40, m = 2) +
             s(dt, k = 40, m = 2) +
             s(dt, site, bs = 'fs', k = 40, m = 2),
           family = binomial(),
           data = data,
           weights = data$wgt,
           discrete = T,
           rho = 0.5136916,
           AR.start = data$ts.start)



# Posterior simulation ----
## Use the posterior distribution to simulate smooths
### Calculate variance-covariance matrix
vcm <- vcov(mod, unconditional = T)


### Simulate parameters (x10000)
sims <- rmvn(10000, mu = coef(mod), V = vcm)


## Create new data
new_data <- data[distance == 800,]
new_data <- new_data[, .(noise = mean(noise),
                         dt = mean(dt),
                         distance = 800,
                         receiver = 'A'), # dummy receiver for predict.bam
                     by = c('date', 'site')]



# Calculate the linear predictor matrix ----
lpm <- predict(mod, new_data, type = 'lpmatrix',
               exclude = c('s(receiver)', 's(dt,site)', 's(noise,site)'))



# D50 calculation ----
# log10(50% / (100% - 50%)) = intercept +
#                             distance_coef * D50 +
#                             squiggly_coefs * squiggly_bits
# 0 = intercept + distance_coef * D50 + squiggly_coefs * squiggly_bits
# D50 = -(intercept + squiggly_coefs * squiggly_bits) / distance_coef
#
# Note that the non-coefficient terms come from the linear predictor matrix
#
# "2" is the coefficient of distance
d50_sim <- data.table(
  (log10(0.5 / (1 - 0.5)) -
     (lpm[, -2] %*% t(sims)[-2,])) /
    t(sims)[2,]
)


# Caclulate median and 95% CIs
d50 <- data.table(d50_med = apply(d50_sim, 1, median),
                  d50_lci = apply(d50_sim, 1, quantile, 0.025),
                  d50_uci = apply(d50_sim, 1, quantile, 0.975))

d50 <- data.table(new_data, d50)



# Calculate change points ----
##  At this time (2020-08), mcp does not offer automatic identification of number of
##    change points, so using changepoint package first, then mcp in order to
##    get an idea of the error around the estimate.
##    Note: need JAGS and rjags installed to run mcp.

set.seed(2000031)

d50[, ':='(n_date = as.numeric(as.Date(date)) - as.numeric(min(as.Date(date))),
           d50_lci = fifelse(d50_lci < 0, 0, d50_lci))]
setorder(d50, site, n_date)



## Use changepoint package to find number of change points ----
cp_d50ns <- cpt.meanvar(d50[site == 'nearshore']$d50_med,
                        method = 'PELT', penalty = 'CROPS',
                        pen.value = c(15, 500))

plot(cp_d50ns, diagnostic = T)
abline(v = 2, col = 'blue')
### Two change points found.

cp_d50ms <- cpt.meanvar(d50[site == 'midshelf']$d50_med,
                        method = 'PELT', penalty = 'CROPS',
                        pen.value = c(15, 500))

plot(cp_d50ms, diagnostic = T)
abline(v = 2, col = 'blue')
### Two change points found.



# Use mcp to model change points ----
## Fit nearshore model
cp_ns_fit <-  mcp(
  ## A model with 3 intercepts and abrupt changes in between
  model = list(
    d50_med ~ 1,
    ~ 1,
    ~ 1
  ),
  par_x = 'n_date',
  ## Dirichlet priors to push the change points away from each other
  prior = list(
    cp_1 = "dirichlet(2)",
    cp_2 = "dirichlet(2)"
  ),
  ## Initialize intercept estimates to help the model converge
  inits = list(int_1 = 600,
               int_2 = 800,
               int_3 = 600),
  data = d50[site == 'nearshore'],
  cores = 3,
  chains = 3,
  iter = 10000
)
summary(cp_ns_fit)


## Fit mid-shelf model
cp_ms_fit <-  mcp(
  ## A model with 3 intercepts and abrupt changes in between
  model = list(
    d50_med ~ 1,
    ~ 1,
    ~ 1
  ),
  par_x = 'n_date',
  ## Dirichlet priors to push the change points away from each other
  prior = list(
    cp_1 = "dirichlet(2)",
    cp_2 = "dirichlet(2)"
  ),
  ## Initialize intercept estimates to help the model converge
  inits = list(int_1 = 600,
               int_2 = 800,
               int_3 = 600),
  data = d50[site == 'midshelf'],
  cores = 3,
  chains = 3,
  iter = 10000
)
summary(cp_ms_fit)

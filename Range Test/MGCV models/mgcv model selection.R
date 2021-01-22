# Packages ----
library(data.table); library(mgcv)



# Functions ----
cv <- function(data, model, k, repeats = 1, seed = NULL){

  terms_to_exclude <- grep('array|station',
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
                  weights = wgt)


    train_data$pred <- predict(CV_mod, type = 'response',
                               exclude = terms_to_exclude)
    # train_data <- train_data[!train_data$distance %in% c(0, 2400),]


    # Test the model
    test_data$pred <- predict(CV_mod, test_data, type = 'response',
                              exclude = c('s(station)', 's(average_noise,array)',
                                          's(dt,array)', 's(average_temperature,array)'))
    # test_data <- test_data[!test_data$distance %in% c(0, 2400),]


    # Penalty functions
    ## Overall RMSE
    train_overall <- sqrt(mean((train_data$freq - train_data$pred) ^ 2))
    test_overall <- sqrt(mean((test_data$freq - test_data$pred) ^ 2))

    ## 250m RMSE
    train_250 <- sqrt(mean((train_data$freq[train_data$distance == 250] -
                            train_data$pred[train_data$distance == 250]) ^ 2))
    test_250 <- sqrt(mean((test_data$freq[test_data$distance == 250] -
                           test_data$pred[test_data$distance == 250]) ^ 2))

    ## 550m RMSE
    train_550 <- sqrt(mean((train_data$freq[train_data$distance == 550] -
                              train_data$pred[train_data$distance == 550]) ^ 2))
    test_550 <- sqrt(mean((test_data$freq[test_data$distance == 550] -
                             test_data$pred[test_data$distance == 550]) ^ 2))

    ## 800m RMSE
    train_800 <- sqrt(mean((train_data$freq[train_data$distance == 800] -
                              train_data$pred[train_data$distance == 800]) ^ 2))
    test_800 <- sqrt(mean((test_data$freq[test_data$distance == 800] -
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



qbn_hack <- function (link = "logit") {
  linktemp <- substitute(link)
  if (!is.character(linktemp))
    linktemp <- deparse(linktemp)
  okLinks <- c("logit", "probit", "cloglog",
               "cauchit", "log")
  family <- "quasibinomial"
  if (linktemp %in% okLinks)
    stats <- make.link(linktemp)
  else if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  }
  else {
    if (inherits(link, "link-glm")) {
      stats <- link
      if (!is.null(stats$name))
        linktemp <- stats$name
    }
    else {
      stop(gettextf("link \"%s\" not available for %s family; available links are %s",
                    linktemp, family, paste(sQuote(okLinks), collapse = ", ")),
           domain = NA)
    }
  }
  aic <- function(y, n, mu, wt, dev) {
    m <- if (any(n > 1))
      n
    else wt
    -2 * sum(ifelse(m > 0, (wt/m), 0) * dbinom(round(m *
                                                       y), round(m), mu, log = TRUE))
  }
  structure(list(family = family, link = linktemp, linkfun = stats$linkfun,
                 linkinv = stats$linkinv, variance = function(mu) mu *
                   (1 - mu), dev.resids = function(y, mu, wt) .Call(stats:::C_binomial_dev_resids,
                                                                    y, mu, wt), aic = aic,
                 mu.eta = stats$mu.eta, initialize = stats:::binomInitialize(family),
                 validmu = function(mu) all(is.finite(mu)) && all(0 <
                                                                    mu & mu < 1), valideta = stats$valideta), class = "family")
}



# Data ----
data <- data.table(
  readRDS('data and imports/rangetest_median_sitespec.RDS')
)
data <- data[distance > 0 & distance < 2400]

setnames(data, gsub(' ', '_', tolower(names(data))))

data[, ':='(array = as.factor(gsub(' ', '', array)),
            station = as.factor(station),
            freq = success / (success + fail),
            wgt = (success + fail) / mean(success + fail),
            date = as.factor(date),
            ts.start = fifelse(date == '2017-12-21' |
                                 (date == '2018-08-08' & station == 'AN3' & distance == 800) |
                                 (date == '2018-08-08' & station == 'AN3_250' & distance == 550),
                               T, F))]

setorder(data, station, distance, date)


# set seed for CV reproducibility; number reflects internal ABIT ms number
seed <- 2000031

# Check global model for evidence of over/underdispersion
m_gq <- bam(freq ~
              distance +
              s(station, bs = 're') +
              s(dt, k = 40, m = 2) +
              s(dt, array, bs = 'fs', k = 40, m = 2) +
              s(average_temperature, k = 40, m = 2) +
              s(average_temperature, array, bs = 'fs', k = 40, m = 2) +
              s(average_noise, k = 40, m = 2) +
              s(average_noise, array, bs = 'fs', k = 40, m = 2) +
              ti(average_noise, average_temperature, k = c(10, 10)) +
              ti(average_noise, dt, k = c(10, 10)),
            family = qbn_hack(),
            data = data,
            weights = data$wgt,
            discrete = T)

summary(m_gq)$dispersion
# 0.1623854
# Evidence of underdispersion. See if AR term helps

# Temporal autocorrelation ----
rho_guess <- acf(resid(m_gq), plot = F)$acf[2]
m_gq_ar <- bam(freq ~
                 distance +
                 s(station, bs = 're') +
                 s(dt, k = 40, m = 2) +
                 s(dt, array, bs = 'fs', k = 40, m = 2) +
                 s(average_temperature, k = 40, m = 2) +
                 s(average_temperature, array, bs = 'fs', k = 40, m = 2) +
                 s(average_noise, k = 40, m = 2) +
                 s(average_noise, array, bs = 'fs', k = 40, m = 2) +
                 ti(average_noise, average_temperature, k = c(10, 10)) +
                 ti(average_noise, dt, k = c(10, 10)),
               family = qbn_hack(),
               data = data,
               weights = data$wgt,
               discrete = T,
               rho = rho_guess,
               AR.start = data$ts.start)

summary(m_gq_ar)$dispersion
# Doesn't really, but autocorrelation is apparent so it's likely needed


#Going back to binomial family
m_g <- bam(freq ~
              distance +
              s(station, bs = 're') +
              s(dt, k = 40, m = 2) +
              s(dt, array, bs = 'fs', k = 40, m = 2) +
              s(average_temperature, k = 40, m = 2) +
              s(average_temperature, array, bs = 'fs', k = 40, m = 2) +
              s(average_noise, k = 40, m = 2) +
              s(average_noise, array, bs = 'fs', k = 40, m = 2) +
              ti(average_noise, average_temperature, k = c(10, 10)) +
              ti(average_noise, dt, k = c(10, 10)),
            family = binomial(),
            data = data,
            weights = data$wgt,
            discrete = T)

rho_guess <- acf(resid(m_g), plot = F)$acf[2]

# ~ Distance ----
##  Linear response to distance only (Null model)
m <- bam(freq ~
           distance +
           s(station, bs = 're'),
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
colMeans(kf_d[, grepl('test', names(kf_d))]) * 100
apply(kf_d[, grepl('test', names(kf_d))], 2, sd) * 100


# ~ Distance + s(noise) ----
## Nonlinear response to noise
m_n <- bam(freq ~
             distance +
             s(station, bs = 're') +
             s(average_noise, k = 40, m = 2) +
             s(average_noise, array, bs = 'fs', k = 40, m = 2),
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
m_ndt <- bam(freq ~
               distance +
               s(station, bs = 're') +
               s(average_noise, k = 40, m = 2) +
               s(average_noise, array, bs = 'fs', k = 40, m = 2) +
               s(dt, k = 40, m = 2) +
               s(dt, array, bs = 'fs', k = 40, m = 2),
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



# ~ Distance + te(noise, dt)----
## Nonlinear response to noise, modulated by stratification
m_ndti <- bam(freq ~
                distance +
                s(station, bs = 're') +
                s(dt, k = 40, m = 2) +
                s(dt, array, bs = 'fs', k = 40, m = 2) +
                s(average_noise, k = 40, m = 2) +
                s(average_noise, array, bs = 'fs', k = 40, m = 2) +
                ti(average_noise, dt, k = c(10, 10)),
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
m_dt <- bam(freq ~
              distance +
              s(station, bs = 're') +
              s(dt, k = 40, m = 2) +
              s(dt, array, bs = 'fs', k = 40, m = 2),
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
m_nbt <- bam(freq ~
               distance +
               s(station, bs = 're') +
               s(average_noise, k = 40, m = 2) +
               s(average_noise, array, bs = 'fs', k = 40, m = 2) +
               s(average_temperature, k = 40, m = 2) +
               s(average_temperature, array, bs = 'fs', k = 40, m = 2),
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



# ~ Distance + te(noise, bwt)----
## Nonlinear response to noise, modulated by near-receiver temperature
m_nbti <- bam(freq ~
                distance +
                s(station, bs = 're') +
                s(average_noise, k = 40, m = 2) +
                s(average_noise, array, bs = 'fs', k = 40, m = 2) +
                s(average_temperature, k = 40, m = 2) +
                s(average_temperature, array, bs = 'fs', k = 40, m = 2) +
                ti(average_noise, average_temperature, k = c(10, 10)),
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
m_bt <- bam(freq ~
              distance +
              s(station, bs = 're') +
              s(average_temperature, k = 40, m = 2) +
              s(average_temperature, array, bs = 'fs', k = 40, m = 2),
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



# QAIC comparison ---
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

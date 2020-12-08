# Packages ----
library(mgcv)



# Functions ----
cv <- function(data, model, k, repeats = 1, seed = NULL){

  refit <- function(i, data_shuffle, folds){
    # Segment data by fold using the which() function
    test_ids <- which(folds == i, arr.ind = TRUE)
    test_data <- data_shuffle[test_ids, ]
    train_data <- data_shuffle[-test_ids, ]

    # Train the model.
    if('gam' %in% class(model)){
      CV_mod <- gam(data = train_data,
                    formula = model$formula,
                    family = model$family$family,
                    method = model$method,
                    control = model$control)
    } else{
      CV_mod <- glm(data = train_data,
                    formula = model$formula,
                    family = model$family$family,
                    control = model$control)
    }


    train_data$pred <- predict(CV_mod, type = 'response',
                               exclude = 'array')
    # train_data <- train_data[!train_data$distance %in% c(0, 2400),]


    # Test the model
    test_data$pred <- predict(CV_mod, test_data, type = 'response',
                              exclude = 'array')
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

QAIC <- function(model, c_hat){
  k <- length(model$coefficients) + 1      # add 1 for the estimate of c_hat
  (model$deviance / c_hat) + 2 * k
}



# Data ----
data <- readRDS('data and imports/rangetest_logit_binary_pt0.RDS')
names(data) <- gsub(' ', '_', tolower(names(data)))
data$array <- as.factor(gsub(' ', '', data$array))
data$freq <- data$success / (data$success + data$fail)
data <- data[data$distance > 0 & data$distance < 2400,]



# Overdispersion ----
## Fit global model it see if data is overdispersed.
## Use quasibinomial family since it returns the dispersion
mod_global <- gam(cbind(success, fail) ~
                         distance + array +
                         s(dt) +
                         s(average_temperature) +
                         s(average_noise) +
                         ti(average_noise, average_temperature) +
                         ti(average_noise, dt),
                       family = quasibinomial(),
                       data = data,
                       method = 'REML',
                       control = gam.control(nthreads = 4))

## c-hat:
c_hat <- summary(mod_global)$dispersion
# [1] 2.498263

## c-hat estimated to be >1, so overdispersed. Moving forward with quasibinomial
##  family and QAIC as outlined in functions section



# ~ Distance ----
##  Linear response to distance only (Null model)
##  There are some convergence issues when fitting as a GAM, so fitting this as
##    a GLM.
mod_d <- glm(freq ~
               distance + array,
             family = quasibinomial(),
             data = data)

##  Calculate adj-R^2 and dev expl to compare with GAM fits
summary_func <- function(model){
  w <- model$prior.weights
  nobs <- nrow(model$model)
  mean.y <- sum(w * model$y) / sum(w)
  residual.df <- length(model$y) - sum(model$edf)

  data.frame(
    adjR2 = 1 - var(w * (as.numeric(model$y) - model$fitted.values)) *
      (nobs - 1) / (var(w * (as.numeric(model$y) - mean.y)) * model$df.residual),
    dev.expl = (model$null.deviance - model$deviance) / model$null.deviance,
    qaic = QAIC(mod_d, c_hat)
  )
}

summary_func(mod_d)
##  5-fold cross-validation
kf_d <- cv(data, mod_d, k = 5, repeats = 5)
colMeans(kf_d[, grepl('test', names(kf_d))]) * 100
apply(kf_d[, grepl('test', names(kf_d))], 2, sd) * 100


# ~ Distance + s(noise) ----
## Nonlinear response to noise
mod_dn <- gam(freq ~
                distance + array +
                s(average_noise, k = 12),
              family = quasibinomial(),
              data = data,
              method = 'REML',
              control = gam.control(nthreads = 4))


##  5-fold cross-validation
kf_dn <- cv(data, mod_dn, k = 5, repeats = 5)
colMeans(kf_dn[, grepl('test', names(kf_dn))]) * 100
apply(kf_dn[, grepl('test', names(kf_dn))], 2, sd) * 100



# ~ Distance + s(noise) + s(dt) ----
## Nonlinear response to noise and stratification
mod_dndt <- gam(freq ~
                  distance + array +
                  s(average_noise) +
                  s(dt),
                family = quasibinomial(),
                data = data,
                method = 'REML',
                control = gam.control(nthreads = 4))


##  5-fold cross-validation
kf_dndt <- cv(data, mod_dndt, k = 5, repeats = 5)
colMeans(kf_dndt[, grepl('test', names(kf_dndt))]) * 100
apply(kf_dndt[, grepl('test', names(kf_dndt))], 2, sd) * 100



# ~ Distance + te(noise, dt)----
## Nonlinear response to noise, modulated by stratification
mod_dndt_int <- gam(freq ~
                      distance + array +
                      s(average_noise) +
                      s(dt) +
                      ti(average_noise, dt),
                    family = quasibinomial(),
                    data = data,
                    method = 'REML',
                    control = gam.control(nthreads = 4))


##  5-fold cross-validation
kf_dndt_int <- cv(data, mod_dndt_int, k = 5, repeats = 5)
colMeans(kf_dndt_int[, grepl('test', names(kf_dndt_int))]) * 100
apply(kf_dndt_int[, grepl('test', names(kf_dndt_int))], 2, sd) * 100



# ~ Distance + s(dt) ----
## Nonlinear response to stratification
mod_ddt <- gam(freq ~
                 distance + array +
                 s(dt),
               family = quasibinomial(),
               data = data,
               method = 'REML',
               control = gam.control(nthreads = 4))


##  5-fold cross-validation
kf_ddt <- cv(data, mod_ddt, k = 5, repeats = 5)
colMeans(kf_ddt[, grepl('test', names(kf_ddt))]) * 100
apply(kf_ddt[, grepl('test', names(kf_ddt))], 2, sd) * 100


# ~ Distance + s(noise) + s(bwt) ----
## Nonlinear response to noise and near-receiver temperature
mod_dnbt <- gam(freq ~
                  distance + array +
                  s(average_noise) +
                  s(average_temperature),
                family = quasibinomial(),
                data = data,
                method = 'REML',
                control = gam.control(nthreads = 4))


##  5-fold cross-validation
kf_dnbt <- cv(data, mod_dnbt, k = 5, repeats = 5)
colMeans(kf_dnbt[, grepl('test', names(kf_dnbt))]) * 100
apply(kf_dnbt[, grepl('test', names(kf_dnbt))], 2, sd) * 100



# ~ Distance + te(noise, bwt)----
## Nonlinear response to noise, modulated by near-receiver temperature
mod_dnbt_int <- gam(freq ~
                      distance + array +
                      s(average_noise) +
                      s(average_temperature) +
                      ti(average_noise, average_temperature),
                    family = quasibinomial(),
                    data = data,
                    method = 'REML',
                    control = gam.control(nthreads = 4))


##  5-fold cross-validation
kf_dnbt_int <- cv(data, mod_dnbt_int, k = 5, repeats = 5)
colMeans(kf_dnbt_int[, grepl('test', names(kf_dnbt_int))]) * 100
apply(kf_dnbt_int[, grepl('test', names(kf_dnbt_int))], 2, sd) * 100



# ~ Distance + s(bwt) ----
## Nonlinear response to near-receiver temperature
mod_dbt <- gam(freq ~
                 distance + array +
                 s(average_temperature),
               family = quasibinomial(),
               data = data,
               method = 'REML',
               control = gam.control(nthreads = 4))


##  5-fold cross-validation
kf_dbt <- cv(data, mod_dbt, k = 5, repeats = 5)
colMeans(kf_dbt[, grepl('test', names(kf_dbt))]) * 100
apply(kf_dbt[, grepl('test', names(kf_dbt))], 2, sd) * 100



# QAIC comparison ---
ic <- data.frame(
  model = c('mod_d', 'mod_dn', 'mod_dndt', 'mod_dndt_int', 'mod_ddt', 'mod_dnbt',
            'mod_dnbt_int', 'mod_dbt'),
  QAIC = sapply(list(mod_d, mod_dn, mod_dndt, mod_dndt_int,
                     mod_ddt, mod_dnbt, mod_dnbt_int, mod_dbt),
                QAIC, c_hat = c_hat)
)

ic$dQAIC <- ic$QAIC - min(ic$QAIC)
ic$wQAIC <- exp(-0.5 * ic$dQAIC)
ic$wQAIC <- ic$wQAIC / sum(ic$wQAIC)
ic <- ic[order(ic$dQAIC),]
ic



mod_global <- gam(freq ~
                    distance + array +
                    s(dt) +
                    s(average_temperature) +
                    s(average_noise) +
                    ti(average_noise, average_temperature) +
                    ti(average_noise, dt),
                  family = quasibinomial(),
                  data = data,
                  method = 'REML',
                  control = gam.control(nthreads = 4))

c_hat <- summary(mod_global2)$dispersion
QAIC <- function(model.candidate, c_hat) {
  k <- length(model.candidate$coefficients) + 1      # add 1 for the estimate of c_hat
  (model.candidate$deviance/c_hat) + 2 * k
}
QAIC(mod)
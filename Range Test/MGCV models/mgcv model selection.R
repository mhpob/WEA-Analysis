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
    CV_mod <- gam(data = train_data,
                  formula = model$formula,
                  family = model$family$family)


    # Test the model
    pred <- predict(CV_mod, test_data, type = "response")

    # Penalty functions
    # RMSE of test > RMSE of train => OVER FITTING of the data.
    # RMSE of test < RMSE of train => UNDER FITTING of the data.
    rmse_train <- sqrt(mean((CV_mod$y - CV_mod$fitted.values) ^ 2))
    rmse_test <- sqrt(mean((
      (test_data$success / (test_data$success + test_data$fail)) - pred) ^ 2))

    c(rmse_train = rmse_train, rmse_test = rmse_test)
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
data <- readRDS('data and imports/rangetest_logit_binary_pt0.RDS')
names(data) <- gsub(' ', '_', tolower(names(data)))
data$array <- as.factor(gsub(' ', '', data$array))



# ~ Distance ----
##  Linear response to distance only (Null model)
mod_d <- gam(cbind(success, fail) ~
               distance +
               s(array, bs = 're'),
             family = binomial(),
             data = data,
             method = 'REML',
             control = gam.control(nthreads = 4))


##  5-fold cross-validation
kf_d <- cv(data, mod_d, k = 5, repeats = 5)
mean(kf_d['rmse_test',])



# ~ Distance + s(noise) ----
## Nonlinear response to noise
mod_dn <- gam(cbind(success, fail) ~
                distance +
                s(average_noise, k = 6, bs = 'tp') +
                s(array, bs = 're'),
              family = binomial(),
              data = data,
              method = 'REML',
              control = gam.control(nthreads = 4))


##  5-fold cross-validation
kf_dn <- cv(data, mod_dn, k = 5, repeats = 5)
mean(kf_dn['rmse_test',])



# ~ Distance + s(noise) + s(dt) ----
## Nonlinear response to noise and stratification
mod_dndt <- gam(cbind(success, fail) ~
                  distance +
                  s(average_noise, k = 6, bs = 'tp') +
                  s(dt, k = 6, bs = 'tp') +
                  s(array, bs = 're'),
                family = binomial(),
                data = data,
                method = 'REML',
                control = gam.control(nthreads = 4))


##  5-fold cross-validation
kf_dndt <- cv(data, mod_dndt, k = 5, repeats = 5)
mean(kf_dndt['rmse_test',])



# ~ Distance + s(noise) + s(dt) + te(noise, dt)----
## Nonlinear response to noise, modulated by stratification
mod_dndt_int <- gam(cbind(success, fail) ~
                      distance +
                      s(average_noise, k = 6, bs = 'tp') +
                      s(dt, k = 6, bs = 'tp') +
                      te(average_noise, dt, bs = 'tp') +
                      s(array, bs = 're'),
                    family = binomial(),
                    data = data,
                    method = 'REML',
                    control = gam.control(nthreads = 4))


##  5-fold cross-validation
kf_dndt_int <- cv(data, mod_dndt_int, k = 5, repeats = 5)
mean(kf_dndt_int['rmse_test',])



# ~ Distance + s(dt) ----
## Nonlinear response to stratification
mod_ddt <- gam(cbind(success, fail) ~
                 distance +
                 s(dt, k = 6, bs = 'tp') +
                 s(array, bs = 're'),
               family = binomial(),
               data = data,
               method = 'REML',
               control = gam.control(nthreads = 4))


##  5-fold cross-validation
kf_ddt <- cv(data, mod_ddt, k = 5, repeats = 5)
mean(kf_ddt['rmse_test',])



# ~ Distance + s(noise) + s(bwt) ----
## Nonlinear response to noise and near-receiver temperature
mod_dnbt <- gam(cbind(success, fail) ~
                  distance +
                  s(average_noise, k = 6, bs = 'tp') +
                  s(average_temperature, k = 6, bs = 'tp') +
                  s(array, bs = 're'),
                family = binomial(),
                data = data,
                method = 'REML',
                control = gam.control(nthreads = 4))


##  5-fold cross-validation
kf_dnbt <- cv(data, mod_dnbt, k = 5, repeats = 5)
mean(kf_dnbt['rmse_test',])



# ~ Distance + s(noise) + s(bwt) + te(noise, bwt)----
## Nonlinear response to noise, modulated by near-receiver temperature
mod_dnbt_int <- gam(cbind(success, fail) ~
                      distance +
                      s(average_noise, k = 6, bs = 'tp') +
                      s(average_temperature, k = 6, bs = 'tp') +
                      te(average_noise, average_temperature, bs = 'tp') +
                      s(array, bs = 're'),
                    family = binomial(),
                    data = data,
                    method = 'REML',
                    control = gam.control(nthreads = 4))


##  5-fold cross-validation
kf_dnbt_int <- cv(data, mod_dnbt_int, k = 5, repeats = 5)
mean(kf_dnbt_int['rmse_test',])



# ~ Distance + s(bwt) ----
## Nonlinear response to near-receiver temperature
mod_dbt <- gam(cbind(success, fail) ~
                 distance +
                 s(average_temperature, k = 6, bs = 'tp') +
                 s(array, bs = 're'),
               family = binomial(),
               data = data,
               method = 'REML',
               control = gam.control(nthreads = 4))


##  5-fold cross-validation
kf_dbt <- cv(data, mod_dbt, k = 5, repeats = 5)
mean(kf_dbt['rmse_test',])

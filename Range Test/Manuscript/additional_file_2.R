# ==================
# MODEL SELECTION
# ==================
# Packages ----
library(mgcv)



# Custom functions ----

## Cross validation
cv <- function(data, model, k, repeats = 1, seed = NULL){

  # Create refitting function
  refit <- function(i, data_shuffle, folds){
    # Segment data by fold
    test_ids <- which(folds == i, arr.ind = TRUE)
    test_data <- data_shuffle[test_ids, ]
    train_data <- data_shuffle[-test_ids, ]

    # Train the model
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

    train_data$pred <- predict(CV_mod,
                               unconditional = T, type = 'response',
                               exclude = 'sitenearshore')


    # Test the model
    test_data$pred <- predict(CV_mod, test_data,
                              unconditional = T, type = 'response',
                              exclude = 'sitenearshore')


    # Errors
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

  # Run repeated folds, if requested.
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

## Calculate adj-R^2 and deviance explained from GLM to compare with GAM fits
summary_func <- function(glm_model){
  w <- glm_model$prior.weights
  nobs <- nrow(glm_model$model)
  mean.y <- sum(w * glm_model$y) / sum(w)
  residual.df <- length(glm_model$y) - sum(glm_model$edf)

  data.frame(
    adjR2 = 1 - var(w * (as.numeric(glm_model$y) - glm_model$fitted.values)) *
      (nobs - 1) / (var(w * (as.numeric(glm_model$y) - mean.y)) * glm_model$df.residual),
    dev.expl = (glm_model$null.deviance - glm_model$deviance) / glm_model$null.deviance
  )
}

# Data ----
data <- read.csv('additional_file_1/stratification_range_test_data.csv')


# ~ Distance ----
##  Linear response to distance only (Null model)
### Fitting as a GLM
mod_null <- glm(detectability ~
               distance + site,
             family = quasibinomial,
             data = data)

##  5-fold cross-validation
kf_null <- cv(data, mod_null, k = 5, repeats = 5, seed = 8675309)

##  Summarize the above. Only use RMSE calculated on the test set. If desired,
##    compare RMSE of the test and train fits. If RMSE test >> RMSE train, then
##    the model is over-fit.
colMeans(kf_null[, grepl('test', names(kf_null))]) * 100
apply(kf_null[, grepl('test', names(kf_null))], 2, sd) * 100
summary_func(mod_null)


## ~ Distance + s(noise) ----
### Nonlinear response to noise
mod_n <- gam(detectability ~
                distance + site +
                s(noise),
              family = quasibinomial,
              data = data,
              method = 'REML')


###  5-fold cross-validation
kf_n <- cv(data, mod_n, k = 5, repeats = 5, seed = 8675309)
colMeans(kf_n[, grepl('test', names(kf_n))]) * 100
apply(kf_n[, grepl('test', names(kf_n))], 2, sd) * 100



## ~ Distance + s(noise) + s(dt) ----
### Nonlinear response to noise and stratification
mod_ndt <- gam(detectability ~
                  distance + site +
                  s(noise) +
                  s(dt),
                family = quasibinomial,
                data = data,
                method = 'REML')


###  5-fold cross-validation
kf_ndt <- cv(data, mod_ndt, k = 5, repeats = 5, seed = 8675309)
colMeans(kf_ndt[, grepl('test', names(kf_ndt))]) * 100
apply(kf_ndt[, grepl('test', names(kf_ndt))], 2, sd) * 100



## ~ Distance + te(noise, dt)----
### Nonlinear response to noise, modulated by stratification
mod_ndt_int <- gam(detectability ~
                      distance + site +
                      te(noise, dt),
                    family = quasibinomial,
                    data = data,
                    method = 'REML')


###  5-fold cross-validation
kf_ndt_int <- cv(data, mod_ndt_int, k = 5, repeats = 5, seed = 8675309)
colMeans(kf_ndt_int[, grepl('test', names(kf_ndt_int))]) * 100
apply(kf_ndt_int[, grepl('test', names(kf_ndt_int))], 2, sd) * 100



## ~ Distance + s(dt) ----
### Nonlinear response to stratification
mod_dt <- gam(detectability ~
                 distance + site +
                 s(dt),
               family = quasibinomial,
               data = data,
               method = 'REML')


###  5-fold cross-validation
kf_dt <- cv(data, mod_dt, k = 5, repeats = 5, seed = 8675309)
colMeans(kf_dt[, grepl('test', names(kf_dt))]) * 100
apply(kf_dt[, grepl('test', names(kf_dt))], 2, sd) * 100


## ~ Distance + s(noise) + s(bwt) ----
### Nonlinear response to noise and near-receiver temperature
mod_nbt <- gam(detectability ~
                  distance + site +
                  s(noise) +
                  s(bwt),
                family = quasibinomial,
                data = data,
                method = 'REML')


###  5-fold cross-validation
kf_nbt <- cv(data, mod_nbt, k = 5, repeats = 5, seed = 8675309)
colMeans(kf_nbt[, grepl('test', names(kf_nbt))]) * 100
apply(kf_nbt[, grepl('test', names(kf_nbt))], 2, sd) * 100



## ~ Distance + te(noise, bwt)----
### Nonlinear response to noise, modulated by near-receiver temperature
mod_nbt_int <- gam(detectability ~
                      distance + site +
                      te(noise, bwt),
                    family = quasibinomial,
                    data = data,
                    method = 'REML')


###  5-fold cross-validation
kf_nbt_int <- cv(data, mod_nbt_int, k = 5, repeats = 5, seed = 8675309)
colMeans(kf_nbt_int[, grepl('test', names(kf_nbt_int))]) * 100
apply(kf_nbt_int[, grepl('test', names(kf_nbt_int))], 2, sd) * 100



## ~ Distance + s(bwt) ----
### Nonlinear response to near-receiver temperature
mod_bt <- gam(detectability ~
                 distance + site +
                 s(bwt),
               family = quasibinomial,
               data = data,
               method = 'REML')


###  5-fold cross-validation
kf_bt <- cv(data, mod_bt, k = 5, repeats = 5, seed = 8675309)
colMeans(kf_bt[, grepl('test', names(kf_bt))]) * 100
apply(kf_bt[, grepl('test', names(kf_bt))], 2, sd) * 100





# =================
# D50 CALCULATION
#==================

# Packages ----
library(changepoint); library(mcp); library(mgcv)



# Custom functions ----
## MVN random deviates
rmvn <- function(n, mu, sig){
  L <- mroot(sig)
  m <- ncol(L)
  t(mu + L %*% matrix(rnorm(m*n), m, n))
}



# Import data ----
data <- read.csv('additional_file_1/stratification_range_test_data.csv')
data$date <- as.Date(data$date)



# Run winning model ----
##  see above for selection
mod_int <- gam(detectability ~
                 distance + site +
                 te(noise, dt),
               family = quasibinomial,
               data = data,
               method = 'REML')



# Posterior simulation ----
## Use the posterior distribution to simulate smooths
### Calculate variance-covariance matrix
vcm <- vcov(mod_int, unconditional = T)


### Simulate parameters (x10000)
sims <- rmvn(10000, mu = coef(mod_int), sig = vcm)


## Create new data
new_data <- data[data$distance == 800,]
new_data <- new_data[!duplicated(new_data[, c('date', 'site')]),]
new_data <- split(new_data, new_data$site)



# Calculate the linear predictor matrix ----
lpm <- lapply(new_data,
              function(.) predict(mod_int, newdata = .,
                                  type = 'lpmatrix', unconditional = T)
)



# D50 calculation ----
# log10(50% / (100% - 50%)) = intercept +
#                             distance_coef * D50 + site_coef * site +
#                             squiggly_coefs * squiggly_bits
# 0 = intercept + distance_coef * D50 + site_coef * site + squiggly_coefs * squiggly_bits
# D50 = -(intercept + site_coef * site + squiggly_coefs * squiggly_bits) / distance_coef
#
# Note that the non-coefficient terms come from the linear predictor matrix
#
# "2" is the coefficient of distance
d50_sim <- lapply(lpm,
                  function(.){
                    (log10(0.5 / (1 - 0.5)) -
                       (.[, -2] %*% t(sims)[-2,])) /
                      t(sims)[2,]
                  }
)


# Pool site-specific simulations
d50_sim <- do.call(cbind, d50_sim)


# Caclulate median and 95% CIs
d50 <- list(d50_med = apply(d50_sim, 1, median),
            d50_lci = apply(d50_sim, 1, quantile, 0.025),
            d50_uci = apply(d50_sim, 1, quantile, 0.975))
d50 <- as.data.frame(d50)


# Add back dates
d50 <- data.frame(date = unique(data$date), d50)



# Calculate change points ----
##  At this time (2020-08), mcp does not offer automatic identification of number of
##    change points, so using changepoint package first, then mcp in order to
##    get an idea of the error around the estimate.
##    Note: need JAGS and rjags installed to run mcp.

d50$n_date <- as.numeric(d50$date) - as.numeric(min(d50$date))



## Use changepoint package to find number of change points ----
cp_d50 <- cpt.meanvar(d50$d50_med,
                      method = 'PELT', penalty = 'CROPS',
                      pen.value = c(15, 500))

plot(cp_d50, diagnostic = T)
abline(v = 2, col = 'blue')
### Two change points found.



# Use mcp to model change points ----
## Fit model
cp_fit <-  mcp(
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
  inits = list(int_1 = 500,
               int_2 = 800,
               int_3 = 500),
  data = d50,
  cores = 3,
  iter = 10000
)
summary(cp_fit)
plot(cp_fit)

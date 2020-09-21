data <- readRDS('p:/obrien/biotelemetry/md wea habitat/wea-analysis/data and imports/rangetest_logit_binary_pt0.RDS')
data$n_date <- as.factor(as.numeric(data$date) - as.numeric(min(data$date)))
data$sc_dist <- scale(data$distance)

options(mc.cores = parallel::detectCores() - 2)

library(tidybayes); library(dplyr); library(brms)

mod_d <- brm(success | trials(success + fail) ~
              sc_dist +
              (1|array),
            data = data,
            family = 'binomial',
            iter = 7500
            )

# saveRDS(mod_d, 'range test/brms models/mod_d.rds')

# mod_d <- readRDS('range test/brms models/mod_d.rds')

library(future)
options(mc.cores = parallel::detectCores() - 2)
plan(multisession)
rstan_options(auto_write = TRUE)
kf_test <- kfold(mod_d, chains = 2, future = T)
saveRDS(kf_test, 'range test/brms models/kf_d.rds')

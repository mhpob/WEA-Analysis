data <- readRDS('p:/obrien/biotelemetry/md wea habitat/wea-analysis/data and imports/rangetest_logit_binary_pt0.RDS')
data$n_date <- as.factor(as.numeric(data$date) - as.numeric(min(data$date)))
data$sc_dist <- scale(data$distance)
names(data) <- gsub(' ', '_', tolower(names(data)))

library(mgcv); library(dplyr); library(brms)

options(mc.cores = parallel::detectCores() - 2)
mod_dn <- brm(success | trials(success + fail) ~
                distance +
                s(average_noise, bs = 'ts') +
                (1|array),
              data = data,
              family = 'binomial',
              chains = 4,
              iter = 5000,
              prior = set_prior('student_t(3, 0, 2.5)', 'sd', 'Intercept', 'array'))
# saveRDS(mod_dn, 'range test/brms models/mod_dist_noise.rds')
mod_dn <- readRDS('range test/brms models/mod_dist_noise.rds')

library(future)
options(mc.cores = parallel::detectCores() - 2)
plan(multiprocess)
kf_dn <- kfold(mod_dn, chains = 2, future = T,
               save_fits = T)
# saveRDS(kf_dn, 'range test/brms models/kf_dist_noise.rds')
kf_dn <- readRDS('range test/brms models/kf_dist_noise.rds')

# ***
# PLEASE CONTACT MIKE OBRIEN (OBRIEN@UMCES.EDU//267-970-1973) IF YOU NEED TO
# KILL THIS PROGRAM.
# ***
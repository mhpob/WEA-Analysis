data <- readRDS('data and imports/rangetest_logit_binary.RDS')
data$n_date <- as.factor(as.numeric(data$date) - as.numeric(min(data$date)))
data$sc_dist <- scale(data$distance)


library(brms)

j <- brm(success | trials(success + fail) ~ sc_dist + (sc_dist|n_date) + (1|array),
         data = data,
         family = 'binomial',
         chains = 3,
         iter = 7000,
         cores = 10,
         file = 'two_array_d50_brms_model')

plot(j) # look for caterpillars, not waves
pp_check(j) #look for overlay

plot.ts((-coef(j)$n_date[, 1, 'Intercept'] / coef(j)$n_date[, 1, 'sc_dist'] *
          sd(data$distance)) + mean(data$distance))
lines((-coef(j)$n_date[, 3, 'Intercept'] / coef(j)$n_date[, 3, 'sc_dist'] *
         sd(data$distance)) + mean(data$distance))
lines((-coef(j)$n_date[, 4, 'Intercept'] / coef(j)$n_date[, 4, 'sc_dist'] *
         sd(data$distance)) + mean(data$distance))


k <- coef(j)$n_date[, , 'Intercept']
k <- data.frame('int' = coef(j)$n_date[, 1, 'Intercept'],
                'slope' = coef(j)$n_date[, 1, 'sc_dist'],
                'r_int' = ranef(j)$n_date[, 1, 'Intercept'],
                'r_slope' = ranef(j)$n_date[, 1, 'sc_dist'],
                'r_slope_lb' = ranef(j)$n_date[, 3, 'sc_dist'],
                'r_slope_ub' = ranef(j)$n_date[, 4, 'sc_dist'])
k$d50 <- (-k$int/k$slope) * sd(data$distance) + mean(data$distance)
View(k)


# we dont want slope to get too positive (>2); set sd prior that truncates at higher values


names(j)




set_prior('normal(0,10)', class = "b", nlpar = "sc_dist", lb = 0)







m <- brm(success | trials(success + fail) ~ sc_dist + (sc_dist|n_date),
         data = data,
         family = 'binomial',
         chains = 3,
         iter = 7000,
         cores = 10,
         file = 'two_array_norandom_d50_brms_model')


plot(m)
pp_check(m)

plot.ts((-coef(m)$n_date[, 1, 'Intercept'] / coef(m)$n_date[, 1, 'sc_dist'] *
           sd(data$distance)) + mean(data$distance))
lines((-coef(m)$n_date[, 3, 'Intercept'] / coef(m)$n_date[, 3, 'sc_dist'] *
         sd(data$distance)) + mean(data$distance),
      col = 'gray')
lines((-coef(m)$n_date[, 4, 'Intercept'] / coef(m)$n_date[, 4, 'sc_dist'] *
         sd(data$distance)) + mean(data$distance),
      col = 'gray')


library(tidybayes)
j <- spread_draws(m, r_n_date[, Intercept])

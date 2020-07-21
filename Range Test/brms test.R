data <- readRDS('p:/obrien/biotelemetry/md wea habitat/wea-analysis/data and imports/rangetest_logit_binary_pt0.RDS')
data$n_date <- as.factor(as.numeric(data$date) - as.numeric(min(data$date)))
data$sc_dist <- scale(data$distance)


library(brms)

j <- brm(success | trials(success + fail) ~ sc_dist + (sc_dist|n_date) + (1|array),
         data = data,
         family = 'binomial',
         chains = 6,
         iter = 5000,
         cores = 6,
         control = list(adapt_delta = 0.98),
         file = 'two_array_d50_brms_model_pt0_adaptd')

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












mod <- readRDS('p:/obrien/biotelemetry/md wea habitat/wea-analysis/two_array_d50_brms_model.rds')
r_eff <- ranef(mod)

library(tidybayes)
get_variables(mod)
k <- mod %>%
  spread_draws(r_n_date[day, term], r_array[array,], b_Intercept, b_sc_dist) %>%
  tidyr::pivot_wider(names_from = term, values_from = r_n_date)
k <- k %>%
  mutate(d50 = (-(b_Intercept + Intercept + r_array) / (b_sc_dist + sc_dist)) *
           sd(data$distance) + mean(data$distance))

l <- k %>%
  summarize(mean = mean(d50),
            median = median(d50),
            pct_25 = quantile(d50 , 0.25),
            pct_95 = quantile(d50, 0.95))


library(ggplot2)
ggplot() +
  geom_ribbon(data = l, aes(x = day+ as.Date('2017-12-21'), ymin = pct_25, ymax = pct_95), fill = 'gray') +
  geom_line(data = l, aes(x = day+ as.Date('2017-12-21'), y = mean)) +
  geom_hline(yintercept = 800) +
  labs(x = NULL, y = 'D50') +
  scale_x_date(date_breaks = 'month', date_labels = '%m-%Y') +
  theme_bw()

ggplot() +
  geom_histogram(data = l, aes(y = median, fill = array), position = 'dodge')




mod_noranarray <- readRDS('p:/obrien/biotelemetry/md wea habitat/wea-analysis/two_array_norandom_d50_brms_model.rds')

d50_draws <- mod_noranarray %>%
  spread_draws(r_n_date[day, term], b_Intercept, b_sc_dist) %>%
  tidyr::pivot_wider(names_from = term, values_from = r_n_date) %>%
  mutate(d50 = (-(b_Intercept + Intercept) / (b_sc_dist + sc_dist)) *
           sd(data$distance) + mean(data$distance)) %>%
  summarize(mean = mean(d50),
            median = median(d50),
            pct_25 = quantile(d50 , 0.25),
            pct_95 = quantile(d50, 0.95),
            max = max(d50),
            min = min(d50))


library(ggplot2)
ggplot(data = d50_draws) +
  geom_ribbon(aes(x = day, ymin = min, ymax = max), fill = 'gray') +
  geom_line(aes(x = day, y = median)) +
  geom_hline(yintercept = 800)

# Packages ----
library(ggplot2); library(data.table); library(mgcv)



# Custom functions ----
## MVN random deviates
rmvn <- function(n, mu, sig){
  L <- mroot(sig)
  m <- ncol(L)
  t(mu + L %*% matrix(rnorm(m*n), m, n))
}



# Import data ----
data <- readRDS('data and imports/rangetest_logit_binary_pt0.RDS')
names(data) <- gsub(' ', '_', tolower(names(data)))
data$array <- as.factor(gsub(' ', '', data$array))
data <- data[data$distance > 0,]
data$freq <- data$success / (data$success + data$fail)



# Import winning model ----
##  see "mgcv model selection.R" for selection

mod_int <- gam(cbind(success, fail) ~
                 distance +
                 te(average_noise, dt) +
                 s(array, bs = 're'),
               family = binomial(),
               data = data,
               method = 'REML',
               control = gam.control(nthreads = 4))



# Posterior simulation ----
## Most code/notes copied from offset_calculation.R
## Code for posterior simulation copied from:
#   https://fromthebottomoftheheap.net/2014/06/16/simultaneous-confidence-intervals-for-derivatives/
#   https://fromthebottomoftheheap.net/2016/12/15/simultaneous-interval-revisited/

## Create new data
new_data <- data.table(data)[distance == 800,]
new_data <- unique(new_data, by = c('date', 'array'))



## Use the posterior distribution to simulate smooths
### Calculate unconditional variance-covariance matrix
vcm <- vcov(mod_int, unconditional = T)


### Simulate parameters (x10000)
sims <- rmvn(10000, mu = coef(mod_int), sig = vcm)



# Calculate the linear predictor matrix ----
# I don't actually know what I'm talking about and this is probably an incorrect
# way to look at it, but I view this as the number that is multiplied by the
# fitted model coefficients.
#
# The coefs of  the linear terms are always multiplied by 1 since that slope
# (the response) is constant no matter where the point is that we're trying to
# predict. The coefs of the curvy terms are multiplied by a whole bunch of funky
# numbers since the response is different depending on the value of the point we're
# trying to fit.
#
# unconditional = T corrects covariance matrix so that it's not conditional on the
# smooth being correct

lpm <- predict(mod_int, new_data, type = 'lpmatrix')

# Note that since we didn't specify knots in the tensor product, mgcv defaulted
#   to (5 ^ 2) - 1 knots.

colnames(lpm)



# D50 calculation ----
# log10(50% / (100% - 50%)) = intercept +
#                             distance_coef * D50 +
#                             squiggly_coefs * squiggly_bits
# **ALGEBRAAAAaaaaa**
# 0 = intercept + distance_coef * D50 + squiggly_coefs * squiggly_bits
# D50 = -(intercept + squiggly_coefs * squiggly_bits) / distance_coef
#
# Note that the non-coefficient terms come from the linear predictor matrix
#
# -c(2, 27, 28) removes the distance coefficients and influence of the random effect
d50_sim <- data.table(
  (log10(0.5 / (1 - 0.5)) -
              (lpm[, -c(2, 27, 28)] %*% t(sims)[-c(2, 27, 28),])) /
  t(sims)[2,]
)

d50 <- list(d50_med = apply(d50_sim, 1, median),
            d50_lci = apply(d50_sim, 1, quantile, 0.025),
            d50_uci = apply(d50_sim, 1, quantile, 0.975))
d50 <- as.data.frame(d50)

d50 <- data.frame(new_data, d50)



# Vis ----
ggplot(data = d50) +
  geom_ribbon(aes(x = date, ymin = d50_lci, ymax = d50_uci, fill = array),
              alpha = 0.5) +
  geom_line(aes(x = date, y = d50_med, color = array)) +
  labs(x = NULL, y = 'Distance at 50% detection probability (m)') +
  ylim(0, 1100)+
  theme_bw()


# D50, but aggregating across arrays ----
new_inn <- new_data[array == 'Inner']
new_wea <- new_data[array == 'MDWEA']

lpm_inn <- predict(mod_int, new_inn, type = 'lpmatrix', unconditional = T)
lpm_wea <- predict(mod_int, new_wea, type = 'lpmatrix', unconditional = T)


d50_sim_inn <- data.table(
  (log10(0.5 / (1 - 0.5)) -
     (lpm_inn[, -c(2, 27, 28)] %*% t(sims)[-c(2, 27, 28),])) /
    t(sims)[2,]
)
d50_sim_wea <- data.table(
  (log10(0.5 / (1 - 0.5)) -
     (lpm_wea[, -c(2, 27, 28)] %*% t(sims)[-c(2, 27, 28),])) /
    t(sims)[2,]
)

d50_sim <- cbind(d50_sim_inn, d50_sim_wea)

d50 <- list(d50_med = apply(d50_sim, 1, median),
            d50_lci = apply(d50_sim, 1, quantile, 0.025),
            d50_uci = apply(d50_sim, 1, quantile, 0.975))
d50 <- as.data.frame(d50)

d50 <- data.frame(date = unique(new_data$date), d50)



# Vis ----
# d50_agg_plot <-
  ggplot(data = d50) +
  geom_ribbon(aes(x = date, ymin = d50_lci, ymax = d50_uci),
              fill = 'gray', size = 1) +
  geom_line(aes(x = date, y = d50_med)) +
  labs(x = NULL, y = 'Distance at 50% detection probability (m)') +
  scale_y_continuous(limits=c(0, 1100), expand = c(0, 0)) +
  scale_x_date(date_breaks = 'month', date_labels = '%b-%y', expand = c(0, 0)) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))



# D50, but including influence of random effects ----
lpm_inn <- predict(mod_int, new_inn, type = 'lpmatrix', unconditional = T)
lpm_wea <- predict(mod_int, new_wea, type = 'lpmatrix', unconditional = T)


d50_sim_inn <- data.table(
  (log10(0.5 / (1 - 0.5)) -
     (lpm_inn[, -2] %*% t(sims)[-2,])) /
    t(sims)[2,]
)
d50_sim_wea <- data.table(
  (log10(0.5 / (1 - 0.5)) -
     (lpm_wea[, -2] %*% t(sims)[-2,])) /
    t(sims)[2,]
)

d50_inn <- apply(d50_sim_inn, 1, median)
d50_wea <- apply(d50_sim_wea, 1, median)

range(d50_inn)
range(d50_wea)

d50 <- data.frame(d50 = c(d50_inn, d50_wea),
                  array = c(rep('inn', times = length(d50_inn)),
                            rep('wea', times = length(d50_wea))),
                  date = rep(unique(new_data$date), times = 2))

# fwrite(d50, 'range test/mgcv models/d50_20200803.csv')




# Calculate change points ----
##  At this time, mcp does not offer automatic identification of number of
##    change points, so using changepoint package first, then mcp in order to
##    get an idea of the error around the estimate.
##    Note: need JAGS and rjags installed to run mcp.
library(changepoint); library(mcp); library(data.table)


d50 <- fread('range test/mgcv models/d50_20200803.csv')
d50 <- d50[, n_date := as.numeric(date) - as.numeric(min(date))]



## Use changepoint package to find number of change points ----
cp_d50wea <- cpt.meanvar(d50[array == 'wea',]$d50,
                         method = 'PELT', penalty = 'CROPS',
                         pen.value = c(15, 500))

plot(cp_d50wea, diagnostic = T)
abline(v = 2, col = 'blue')
### Two change points found.


cp_d50inn <- cpt.meanvar(d50[array == 'inn',]$d50,
                         method = 'PELT', penalty = 'CROPS',
                         pen.value = c(15, 500))

plot(cp_d50inn, diagnostic = T)
abline(v = 2, col = 'blue')
### Two change points found



# Use mcp to model change points ----
## A model with 3 intercepts and abrupt changes in between
inn_model = list(
  d50 ~ 1,
  ~ 1,
  ~ 1
)

## Dirichlet priors to push the change points away from each other
prior = list(
  cp_1 = "dirichlet(2)",
  cp_2 = "dirichlet(2)"
)


## Fit model
inn_fit <-  mcp(inn_model, data = d50[array == 'inn'], prior = prior,
                par_x = 'n_date', iter = 10000, chains = 3, cores = 3,
                adapt = 9000)
plot(inn_fit)


wea_model = list(
  d50 ~ 1,   # constant rate
  ~ 1,
  ~ 1
)

prior = list(
  cp_1 = "dirichlet(2)",
  cp_2 = "dirichlet(2)"
)

wea_fit <- mcp(wea_model, data = d50[array == 'wea'], prior = prior,
               par_x = 'n_date', iter = 10000, adapt = 5000, chains = 3, cores = 3)
wea_fit
plot(wea_fit)



# Stack density onto previous D50 vis ----


if (which_y == "ct" && geom_data != FALSE) {
  limits_y = c(min(fit$data[, fit$pars$y]),
               max(fit$data[, fit$pars$y]))
} else if (any(q_predict != FALSE)) {
  limits_y = c(min(samples_expanded$.predicted),
               max(samples_expanded$.predicted))
} else if (as.character(yvar) %in% names(samples_expanded)) {
  limits_y = c(min(dplyr::pull(samples_expanded, as.character(yvar))),
               max(dplyr::pull(samples_expanded, as.character(yvar))))
} else {
  stop("Failed to draw change point density for this plot. Please raise an error on GitHub.")
}



d50_agg_plot + geom_cp_density(wea_fit,c(0,1100)) +
  ggplot2::coord_cartesian(
    ylim = c(limits_y[1], NA),  # Remove density flat line from view
    xlim = c(min(wea_fit$data[, wea_fit$pars$x]), max(wea_fit$data[, wea_fit$pars$x]))  # Very broad varying change point posteriors can expand beyond observed range. TO DO
  )


# From internals of plot.mcp (not exported)
geom_cp_density = function(fit, facet_by = NULL, limits_y) {
  dens_scale = 0.2  # Proportion of plot height
  dens_cut = 0.05  # How much to move density down. 5% is ggplot default. Move a bit further.

  # Get varying and population change point parameter names
    varying = NULL
    population = wea_fit$.other$ST$cp_name[-1]



  # Get samples in long format

  samples = tidybayes::spread_draws(fit$mcmc_prior, population = population, varying = varying, absolute = TRUE)
  samples = tidyr::pivot_longer(samples, cols = tidyselect::starts_with("cp_"), names_to = "cp_name", values_to = "value")

  # Make the geom!
  ggplot2::stat_density(aes(
    x = value,
    y = ..scaled.. * diff(limits_y) * dens_scale +  # Scale to proportion of view
      0 -  # Put on x-axis
      diff(limits_y) * dens_cut,  # Move a bit further down to remove zero-density line from view.
    group = paste0(.chain, cp_name),  # Apply scaling for each chain X cp_i combo
    color = .chain
  ),
  data = samples,
  position = "identity",
  geom = "line",
  show.legend = FALSE
  )
}


k <- do.call(rbind, wea_fit$mcmc_post)
k <- data.table(k)
k <- melt(k[,1:2], measure.vars = c('cp_1', 'cp_2'), variable.name = 'cp',value.name = 'date')
k <- k[, date := as.Date(date + as.numeric(min(d50$date)))]

d50_agg_plot + ggplot2::stat_density(aes(
  x = date,
  y = ..scaled.. * diff(c(0,1100)) * dens_scale +  # Scale to proportion of view
    0 -  # Put on x-axis
    diff(c(0,1100)) * dens_cut,  # Move a bit further down to remove zero-density line from view.
  group = cp
),data = k,
position = "identity",
geom = "line",
show.legend = FALSE)

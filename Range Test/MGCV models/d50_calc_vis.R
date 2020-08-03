# Packages ----
library(ggplot2); library(patchwork); library(data.table)
library(changepoint); library(mcp); library(mgcv)



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


# fwrite(d50, 'range test/mgcv models/d50_20200803.csv')


# Vis ----
d50_agg_plot <-
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

d50_ran <- data.frame(d50 = c(d50_inn, d50_wea),
                      array = c(rep('inn', times = length(d50_inn)),
                                rep('wea', times = length(d50_wea))),
                      date = rep(unique(new_data$date), times = 2))






# Calculate change points ----
##  At this time, mcp does not offer automatic identification of number of
##    change points, so using changepoint package first, then mcp in order to
##    get an idea of the error around the estimate.
##    Note: need JAGS and rjags installed to run mcp.

d50 <- fread('range test/mgcv models/d50_20200803.csv')
d50 <- d50[, n_date := as.numeric(date) - as.numeric(min(date))]



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
  inits = list(int_1 = 550,
               int_2 = 800,
               int_3 = 600),
  data = d50,
  cores = 3,
  iter = 10000
)
summary(cp_fit)
plot(cp_fit)



# Stack density onto previous D50 vis ----
##  Code modified from geom_cp_density() within internals of plot.mcp().
##  geom_cp_density is not exported, check Github.

## Aggregate posterior samples
posts <- do.call(rbind, cp_fit$mcmc_post)
posts <- data.table(posts)
posts <- melt(posts[, c('cp_1', 'cp_2')],
              measure.vars = c('cp_1', 'cp_2'),
              variable.name = 'cp',
              value.name = 'date')
posts <- posts[, date := as.Date(round(date) + as.numeric(min(d50$date)))]


## Grab y axis limits from the time series plot
base_lims <- ggplot_build(d50_agg_plot)$layout$panel_scales_y[[1]]$limits


## Plot
TS <- d50_agg_plot +
  ggplot2::stat_density(aes(x = date,
                            # scale to 20% y axis range
                            y = ..scaled.. * diff(base_lims) * 0.2 - 1e-13,
                            group = cp),
                        data = posts, linetype = 'dashed', size = 1,
                        position = "identity",
                        geom = "line",
                        show.legend = FALSE) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))



# Create and append histogram ----
hist <- ggplot() +
  geom_histogram(data= d50_ran, aes(y = d50, fill = array),
                 binwidth = 50,
                 position = 'dodge',
                 show.legend = F) +
  scale_fill_viridis_d() +
  scale_y_continuous(limits = c(0, 1100), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 60), expand = c(0,0)) +
  labs(x = 'Count', y = NULL) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        panel.grid.minor.x = element_blank())


TS + hist + plot_layout(widths = c(4, 1))

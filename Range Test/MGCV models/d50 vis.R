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

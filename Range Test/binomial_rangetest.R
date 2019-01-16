library(dplyr)
library(mgcv)

data <- readRDS('data and imports/rangetest_logit_binary.RDS') %>%
  mutate(array = as.factor(array),
         date_f = as.factor(date))

mod <- gam(cbind(success, fail) ~ distance + s(dt, bs = "ts") +
            s(array, bs = 're') + s(date_f, bs = 're'),
         data = data, family = 'binomial', method = 'REML')
summary(mod)
anova(mod)

j <- predict(mod, data.frame(distance = seq(1, 1000, 10),
                           dt = 0,
                           array = 'MD WEA'),
             type = 'response', exclude = c('s(date_f)', 's(array)'))
j_se <- predict(mod, data.frame(distance = seq(1, 1000, 10),
                              dt = 0,
                              array = 'MD WEA'),
             type = 'link', exclude = c('s(date_f)', 's(array)'), se = T)
jj <- predict(mod, data.frame(distance = seq(1, 1000, 10),
                            dt = 5,
                            array = 'MD WEA'),
              type = 'response')
jjj <- predict(mod, data.frame(distance = seq(1, 1000, 10),
                            dt = 1,
                            array = 'MD WEA'),
              type = 'response', )
plot(j ~ seq(1, 1000, 10), type = 'l')
lines(mod$family$linkinv(j_se$fit + 1.96*j_se$se.fit) ~ seq(1, 1000, 10))
lines(mod$family$linkinv(j_se$fit - 1.96*j_se$se.fit) ~ seq(1, 1000, 10))
lines(jj ~ seq(1, 1000, 10), col = 'blue')
lines(jjj ~ seq(1, 1000, 10), col = 'red')


j_se2 <- predict(mod, data.frame(distance = d50,
                              dt = 0),
                type = 'link', se = T)
plot(seq(1, 1000, 10) ~ j_se$fit, type = 'l')
lines(seq(1, 1000, 10) ~ j_se$fit + 1.96 * j_se$se.fit, lty = 'dashed')
lines(seq(1, 1000, 10) ~ j_se$fit - 1.96 * j_se$se.fit, lty = 'dashed')


# sum of linear predictor matrix without the distance term multiplied by their coefficients,
# divided by the distance coef
lpm <- predict(mod, data.frame(distance = 500,
                             dt = 5,
                             date_f = '2018-12-02',
                             array = 'MD WEA'),
               type = 'lpmatrix', exclude = c('s(date_f)', 's(array)'))
d50 <- -sum(lpm[1, -2] * coef(mod)[-2]) / coef(mod)[2]



mod2 <- gamm(cbind(success, fail) ~ distance + s(dt, bs = "ts"),
             random = list(date_f = ~1),
             data = data, family = 'binomial', method = 'REML')

lpm <- predict(mod2$gam, data.frame(distance = 1,
                               dt = 10),
               type = 'lpmatrix', exclude = c('s(date_f)'))
d50 <- -(lpm[, -2] %*% coef(mod2$gam)[-2]) / coef(mod2$gam)[2]



rmvn <- function(n,mu,sig) { ## MVN random deviates
  L <- mroot(sig)
  m <- ncol(L)
  t(mu + L %*% matrix(rnorm(m * n), m, n))
}

br <- MASS::mvrnorm(10000, coef(mod2$gam), vcov(mod2$gam)) ## 1000 replicate param. vectors
res <- rep(0, 10000)

for (i in 1:10000){
  res[i] <- -(lpm[, -2] %*% br[i, -2]) / br[i, 2] ## example non-linear function
}






# SE through bootstrap ----
## Bootstrapping
library(boot)

# Function to be bootstrapped
boot_fun <- function(x, indices, form){
  d <- x[indices,]
  form <- as.formula(form)
  fit <- mgcv::gam(form, family = 'binomial', data = d, method = 'REML')
  coeffs <- coef(fit)

  lpm <- predict(fit, data.frame(distance = seq(1, 1000, 10),
                               dt = 0),
                 type = 'lpmatrix')

  # https://stackoverflow.com/a/35589393/7496818
  # solve the following for dist, holding all terms but dist at a constant:
  # logit(p) ~ intercept + (GLM/GAM solved terms) + dist_coef * dist
  # logit(0.5) ~ intercept + (GLM/GAM solved terms) + dist_coef * D50
  # 0 ~ intercept + (GLM/GAM solved terms) + dist_coef * D50
  # -(intercept + (GLM/GAM solved terms)) ~ dist_coef * D50
  # D50 ~ -(intercept + (GLM/GAM solved terms)) / dist_coef
  #
  # In words:
  # negative sum of linear predictor matrix without the distance term,
  # multiplied by their coefficients,
  # then divided by the distance coef

  d50 <- setNames(-sum(lpm[1, -2] * coef(fit)[-2]) / coef(fit)[2], 'd50')

  # Returns the distance at 50% detection
  d50
}

# Apply bootstrap
l_e_boot <- function(data, form){
  ncore <- as.integer(Sys.getenv('NUMBER_OF_PROCESSORS')) - 2

  boot(data = data, statistic = boot_fun, R = 2500, form = form,
       parallel = 'snow', ncpus = ncore)
}

test <- l_e_boot(p, pct_det ~ as.numeric(dist) + s(dt, bs = "ts"))

# BCa (adjusted bootstrap percentile) 95% CIs
boot.ci(test, type = 'basic') # BCa won't work unless R is large enough (thousands)


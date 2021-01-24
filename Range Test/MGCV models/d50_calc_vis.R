# Packages ----
library(ggplot2); library(patchwork); library(data.table) ;library(mgcv)



# Import data ----
data <- data.table(
  readRDS('data and imports/rangetest_median_sitespec.RDS')
)
data <- data[distance > 0 & distance < 2400]

setnames(data, gsub(' ', '_', tolower(names(data))))

data[, ':='(array = as.factor(gsub(' ', '', array)),
            station = as.factor(station),
            freq = success / (success + fail),
            wgt = (success + fail) / mean(success + fail),
            date = as.factor(date),
            ts.start = fifelse(date == '2017-12-21' |
                                 (date == '2018-08-08' & station == 'AN3' & distance == 800) |
                                 (date == '2018-08-08' & station == 'AN3_250' & distance == 550),
                               T, F))]

setorder(data, station, distance, date)


# Import winning model ----
##  see "mgcv model selection.R" for selection
m_ndt <- bam(freq ~ distance +
               s(station, bs = 're') +
               s(average_noise, k = 40, m = 2) +
               s(dt, k = 40, m = 2) +
               s(average_noise, array, k = 40, m = 2, bs = 'fs') +
               s(dt, array, k = 40, m = 2, bs = 'fs'),
             family = binomial(),
             data = data,
             weights = data$wgt,
             discrete = T,
             rho = 0.5,
             AR.start = data$ts.start)



# Posterior simulation ----
## Most code/notes copied from offset_calculation.R
## Code for posterior simulation copied from:
#   https://fromthebottomoftheheap.net/2014/06/16/simultaneous-confidence-intervals-for-derivatives/
#   https://fromthebottomoftheheap.net/2016/12/15/simultaneous-interval-revisited/

## Create new data
new_data <- data[distance == 800,]
new_data <- new_data[, .(average_noise = mean(average_noise),
                         average_temperature = mean(average_temperature),
                         dt = mean(dt),
                         freq = mean(freq),
                         distance = 800,
                         station = 'IS2'),
                     by = c('date', 'array')]
setorder(new_data, date, array)

# new_data <- unique(new_data, by = c('date', 'array'))



## Use the posterior distribution to simulate smooths
### Calculate unconditional variance-covariance matrix
vcm <- vcov(m_ndt, unconditional = T)


### Simulate parameters (x10000)
sims <- rmvn(10000, mu = coef(m_ndt), V = vcm)


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

lpm <- predict(m_ndt, new_data, type = 'lpmatrix',
               exclude = c('s(station)', 's(dt,array)', 's(average_noise,array)'))

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
              (lpm[, -2] %*% t(sims)[-2,])) /
  t(sims)[2,]
)

d50 <- data.table(d50_med = apply(d50_sim, 1, median),
                  d50_lci = apply(d50_sim, 1, quantile, 0.025),
                  d50_uci = apply(d50_sim, 1, quantile, 0.975))

d50 <- data.table(new_data, d50)



# Vis ----
ggplot(data = d50) +
  geom_ribbon(aes(x = as.Date(date), ymin = d50_lci, ymax = d50_uci),
              alpha = 0.5) +
  geom_line(aes(x = as.Date(date), y = d50_med)) +
  facet_wrap(~array, ncol = 1)+
  labs(x = NULL, y = 'Distance at 50% detection probability (m)') +
  coord_cartesian(ylim=c(0, 1250), expand = F) +
  # ylim(0, 1100)+
  theme_bw()


# d50 <- fwrite(d50, 'range test/mgcv models/d50_20210120.csv')



# Calculate change points ----
##  At this time, mcp does not offer automatic identification of number of
##    change points, so using changepoint package first, then mcp in order to
##    get an idea of the error around the estimate.
##    Note: need JAGS and rjags installed to run mcp.
library(ggplot2); library(changepoint); library(mcp); library(data.table)

set.seed(2000031)

d50 <- fread('range test/mgcv models/d50_20210120.csv')
d50 <- d50[, ':='(n_date = as.numeric(as.Date(date)) - as.numeric(min(as.Date(date))),
                  array = fifelse(array == 'Inner', 'Nearshore', 'Mid-shelf'),
                  d50_lci = fifelse(d50_lci < 0, 0, d50_lci))]



## Use changepoint package to find number of change points ----
cp_d50ns <- cpt.meanvar(d50[array == 'Nearshore']$d50_med,
                         method = 'PELT', penalty = 'CROPS',
                         pen.value = c(15, 500))

plot(cp_d50ns, diagnostic = T)
abline(v = 2, col = 'blue')
### Two change points found.

cp_d50ms <- cpt.meanvar(d50[array == 'Mid-shelf']$d50_med,
                         method = 'PELT', penalty = 'CROPS',
                         pen.value = c(15, 500))

plot(cp_d50ms, diagnostic = T)
abline(v = 2, col = 'blue')
### Two change points found.



# Use mcp to model change points ----
## Fit Nearshore model
cp_ns_fit <-  mcp(
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
  inits = list(int_1 = 600,
               int_2 = 800,
               int_3 = 600),
  data = d50[array == 'Nearshore'],
  cores = 3,
  chains = 3,
  iter = 10000
)
summary(cp_ns_fit)
# plot(cp_fit)


## Fit Inner model
cp_ms_fit <-  mcp(
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
  inits = list(int_1 = 600,
               int_2 = 800,
               int_3 = 600),
  data = d50[array == 'Mid-shelf'],
  cores = 3,
  chains = 3,
  iter = 10000
)
summary(cp_ms_fit)


# Stack density onto previous D50 vis ----
##  Code modified from geom_cp_density() within internals of plot.mcp().
##  geom_cp_density is not exported, check Github.

## Aggregate posterior samples
posts_ns <- do.call(rbind, cp_ns_fit$mcmc_post)
posts_ns <- data.table(posts_ns)
posts_ns <- melt(posts_ns[, c('cp_1', 'cp_2')],
              measure.vars = c('cp_1', 'cp_2'),
              variable.name = 'cp',
              value.name = 'date')
posts_ns <- posts_ns[, ':='(date = as.Date(round(date) + as.numeric(min(d50$date))),
                      array = 'Nearshore')]


posts_ms <- do.call(rbind, cp_ms_fit$mcmc_post)
posts_ms <- data.table(posts_ms)
posts_ms <- melt(posts_ms[, c('cp_1', 'cp_2')],
                  measure.vars = c('cp_1', 'cp_2'),
                  variable.name = 'cp',
                  value.name = 'date')
posts_ms <- posts_ms[, ':='(date = as.Date(round(date) + as.numeric(min(d50$date))),
                              array = 'Mid-shelf')]

posts <- rbind(
  posts_ns,
  posts_ms
)

## Misc. stats
### Median and credible interval
posts[, .(median = median(date),
          lower = as.Date(quantile(as.numeric(date), 0.025)),
          upper = as.Date(quantile(as.numeric(date), 0.975))),
      by = c('array', 'cp')]


# Data for wind rug ----
ndbc_dl <- function(url){
  buoy_header <- fread(url,
                       nrows = 2, header = F)
  buoy_header <- buoy_header[, lapply(.SD, function(.) gsub('#', '', .))]
  buoy_header <- buoy_header[, lapply(.SD, function(.) gsub('/', '_', .))]
  buoy_header <- buoy_header[, lapply(.SD, function(.) paste(., collapse = '_'))]

  buoy_data <- fread(url,
                     skip = 2,
                     col.names = as.character(buoy_header))
  buoy_data <- buoy_data[WSPD_m_s != 99,]


  buoy_data <- buoy_data[, date := as.Date(paste(YY_yr, MM_mo, DD_dy, sep = '-'))]
}


coast_2017 <- ndbc_dl('https://www.ndbc.noaa.gov/view_text_file.php?filename=44009h2017.txt.gz&dir=data/historical/stdmet/')
coast_2018 <- ndbc_dl('https://www.ndbc.noaa.gov/view_text_file.php?filename=44009h2018.txt.gz&dir=data/historical/stdmet/')

coast <- rbind(coast_2017, coast_2018)
coast <- coast[date >= '2017-12-21' & date <= '2018-12-04']


coast <- coast[, .(max = max(WSPD_m_s)),
               by = date]



# No records after 2018-11-15; supplement with Ocean City Inlet records
ocmd <- ndbc_dl('https://www.ndbc.noaa.gov/view_text_file.php?filename=ocim2h2018.txt.gz&dir=data/historical/stdmet/')

ocmd <- ocmd[date > '2018-11-15' & date <= '2018-12-04']

# Coastal buoy are hourly averages, so do the same here
ocmd <- ocmd[, .(wspd = mean(WSPD_m_s),
                 date = as.Date(paste(YY_yr, MM_mo, DD_dy, sep = '-'))),
             by = c('YY_yr', 'MM_mo', 'DD_dy', 'hh_hr')]

ocmd <- ocmd[, .(max = max(wspd)), by = date]

coast <- rbind(coast, ocmd)

coast[, array := 'Nearshore']



## Plot
d50_plot <-
  ggplot(data = d50) +

  # D50 time series
  geom_ribbon(aes(x = date, ymin = d50_lci, ymax = d50_uci), fill = 'gray63') +
  geom_line(aes(x = date, y = d50_med, color = array), size = 0.3) +
  labs(x = NULL, y = 'Distance at 50% detectability (m)') +
  facet_wrap(~array, ncol = 1, strip.position = 'right') +
  scale_color_manual(values = c('#0072B2', '#D55E00')) +

  # Label receiver tending
  geom_rug(
    data = data.table(
      date = as.Date(c('2018-04-11', '2018-08-08',  '2018-08-09')),
      array = 'Nearshore'
    ),
    aes(x = date), length = unit(0.1, "npc"),
    size = 0.11, color = 'red', linetype = 'dashed',
    inherit.aes = T,
    show.legend = F,
    sides = 't', outside = T) +

  # Wind rug
  geom_rug(data = coast[max >= 17.2],
           aes(x = date), length = unit(0.07, "npc"),
           size = 0.1,
           inherit.aes = T,
           show.legend = F,
           sides = 't', outside = T) +
  geom_rug(data = coast[max >= 10.8 & max < 17.2],
           aes(x = date), length = unit(0.03, "npc"),
           size = 0.1,
           inherit.aes = T,
           show.legend = F,
           sides = 't', outside = T) +

  coord_cartesian(ylim = c(0, 1250), expand = F, clip = 'off') +
  scale_x_date(date_breaks = 'month', date_labels = '%b', expand = c(0, 0)) +
  theme_bw() +
  theme(legend.position = 'none')

## Grab y axis limits from the time series plot
base_lims <-
  ## Un-comment if using scale_y_... above
  # ggplot_build(d50_agg_plot)$layout$panel_scales_y[[1]]$limits

  ## Un-comment if using coord_cartesian above
  ggplot_build(d50_plot)$plot$coordinates$limits$y

compiled_plot <-
  d50_plot +
  ggplot2::stat_density(aes(x = date,
                            # scale to 20% y axis range
                            y = ..scaled.. * diff(base_lims) * 0.2 - 1e-13,
                            group = cp),
                        data = posts, color = 'brown', size = 0.3,
                        position = "identity",
                        geom = "line",
                        show.legend = FALSE) +
  theme(axis.text.y = element_text(size = 5, hjust = 1, angle = 30,margin = margin(r = 0.5)),
        axis.text.x = element_text(size = 5, margin = margin(t = 0)),
        axis.title.y = element_text(size = 6, margin = margin(0, 0, -3, 0)),
        plot.margin = unit(c(0.6, 0.05, 0, 0.1), 'mm'),
        panel.border = element_rect(size = 0.3),
        panel.grid.minor.y = element_blank(),
        panel.grid.major = element_line(size = 0.3),
        panel.grid.minor.x = element_line(size = 0.1),
        axis.ticks = element_line(size = 0.1),
        axis.ticks.length = unit(0.5, 'mm'),
        strip.background = element_blank(),
        strip.text = element_text(size = 6, margin = margin(0, 0, 0.1, 0)))


library(ragg)

agg_tiff("range test/manuscript/revisions/figures/Figure4.tif",
     width = 85, height = 40, units = 'mm', compression = 'lzw', res = 600)

compiled_plot

dev.off()

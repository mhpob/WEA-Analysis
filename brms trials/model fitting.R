library(data.table); library(brms)

basswl <- fread('p:/obrien/bassweek.csv')
basswl[, freq_p := freq / 7]
basswl[, Origin := factor(Origin)]
basswl[, transmitter := factor(transmitter)]
basswl[, ar_start := fifelse(time == 1, T, F)]
basswl[, trials := 7]
basswl[, freqL1_trials := freqL1/trials]



# SBwL1_zib <- brm(
#   freq | trials(7) ~ s(SST,bs='tp',k=11)+
#     s(WE, bs='cp',k=27)+
#     s(CHLA, bs= 'tp',k=6)+
#     s(freqL1,bs='tp',k=6)+
#     Origin+
#     s(transmitter, bs= 're'),
#
#   family = zero_inflated_binomial(),
#   data = basswl,
#   iter = 500,
#   cores = 4,
#   chains = 4,
#
#   file = 'SBwL1_zib'
#
# )
#
# SBwL1_bin <- brm(
#   freq | trials(7) ~ s(SST,bs='tp',k=11)+
#     s(WE, bs='cp',k=27)+
#     s(CHLA, bs= 'tp',k=6)+
#     s(freqL1,bs='tp',k=6)+
#     Origin+
#     s(transmitter, bs= 're'),
#
#   family = binomial(),
#   data = basswl,
#   iter = 500,
#   cores = 4,
#   chains = 4,
#
#   file = 'SBwL1_bin'
#
# )

SBwL1_zib <- readRDS('SbwL1_zib.rds')
SBwL1_bin <- readRDS('SbwL1_bin.rds')

plot(conditional_smooths(SBwL1_zib))

hypothesis(
  SBwL1_zib,
  c('Intercept - OriginHudsonRiverNY = 0',
    'Intercept - OriginKennebecRiverME = 0',
    'Intercept - OriginMassachusetts = 0',
    'Intercept - OriginPotomacRiverMD = 0',
    'OriginHudsonRiverNY - OriginKennebecRiverME = 0',
    'OriginHudsonRiverNY - OriginMassachusetts = 0',
    'OriginHudsonRiverNY - OriginPotomacRiverMD = 0',
    'OriginKennebecRiverME - OriginMassachusetts = 0',
    'OriginKennebecRiverME - OriginPotomacRiverMD = 0',
    'OriginMassachusetts - OriginPotomacRiverMD = 0'
  )

)



# job::job({
#   SBwL1_zib <- brm(
#     freq | trials(7) ~
#       s(SST,bs='tp',k=11) +
#       s(WE, bs='cc',k=27) +
#       s(CHLA, bs= 'tp',k = 6) +
#       ar(gr = transmitter) +
#       Origin+
#       (1|transmitter),
#
#     family = zero_inflated_binomial(),
#     data = basswl,
#     iter = 500,
#     cores = 4,
#     chains = 4,
#
#     file = 'SBwL1_zib_group_ar'
#   )
# })
#
# prior_zip <- c(set_prior('beta(2.5,2.6)', class = 'zi'),
#                )
#
# job::job({
#   SBwL1_zip <- brm(
#     freq ~
#       s(SST,bs='tp',k=11) +
#       s(WE, bs='cc',k=27) +
#       s(CHLA, bs= 'tp',k = 6) +
#       s(freqL1, bs = 'tp', k = 6) +
#       Origin+
#       (1|transmitter),
#
#     family = zero_inflated_poisson(),
#     data = basswl,
#     cores = 4,
#     chains = 4,
#     file = 'SBwL1_zip_group'
#   )
# })


job::job({
  m_fullprior_binom <- brm(
    freq | trials(7) ~ Origin +
      s(SST,bs='tp',k = 11) +
      s(WE, bs='cc',k = 27) +
      s(CHLA, bs= 'tp',k = 6) +
      ar(gr = transmitter) +
      (1|transmitter),

    family = binomial(),
    data = basswl,
    cores = 4,
    chains = 4,
    thin = 5,

    prior = c(
      set_prior('normal(0.6,0.5)', class = 'ar'),
      set_prior('normal(0,20)', class = 'Intercept'),
      set_prior('student_t(5,0,1)', class = 'b', coef = 'OriginHudsonRiverNY'),
      set_prior('student_t(5,0,1)', class = 'b', coef = 'OriginKennebecRiverME'),
      set_prior('student_t(5,0,1)', class = 'b', coef = 'OriginMassachusetts'),
      set_prior('student_t(5,0,1)', class = 'b', coef = 'OriginPotomacRiverMD'),
      set_prior('normal(0,20)', class = 'b', coef = 'sCHLA_1'),
      set_prior('normal(0,20)', class = 'b', coef = 'sSST_1')
    ),
    sample_prior = 'yes',

    file = 'brms trials/bin_full_prior.rds'
  )
})




SBwL1_zib_prior <- readRDS('brms trials/SbwL1_zib_group_ar_prior.rds')
SBwL1_zip <- readRDS('SbwL1_zip_group.rds')




job::job({
  SBwL1_bin <- brm(
    freq | trials(7) ~
      s(SST,bs='tp',k = 11) +
      s(WE, bs='cc',k = 27) +
      s(CHLA, bs= 'tp',k = 6) +
      ar(gr = transmitter) +
      Origin+
      (1|transmitter),

    family = binomial(),
    data = basswl,
    iter=500,
    cores = 4,
    chains = 4,

    prior = c(
      set_prior('normal(0,0.75)', class = 'ar'),
      set_prior('normal(0,100)', class = 'b')
    ),

    file = 'SBwL1_bin_ar.rds'

  )
})

SBwL1_bin <- readRDS('SBwL1_bin.rds')

qbn_hack <- function (link = "logit") {
  linktemp <- substitute(link)
  if (!is.character(linktemp))
    linktemp <- deparse(linktemp)
  okLinks <- c("logit", "probit", "cloglog",
               "cauchit", "log")
  family <- "quasibinomial"
  if (linktemp %in% okLinks)
    stats <- make.link(linktemp)
  else if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  }
  else {
    if (inherits(link, "link-glm")) {
      stats <- link
      if (!is.null(stats$name))
        linktemp <- stats$name
    }
    else {
      stop(gettextf("link \"%s\" not available for %s family; available links are %s",
                    linktemp, family, paste(sQuote(okLinks), collapse = ", ")),
           domain = NA)
    }
  }
  aic <- function(y, n, mu, wt, dev) {
    m <- if (any(n > 1))
      n
    else wt
    -2 * sum(ifelse(m > 0, (wt/m), 0) * dbinom(round(m *
                                                       y), round(m), mu, log = TRUE))
  }
  structure(list(family = family, link = linktemp, linkfun = stats$linkfun,
                 linkinv = stats$linkinv, variance = function(mu) mu *
                   (1 - mu), dev.resids = function(y, mu, wt) .Call(stats:::C_binomial_dev_resids,
                                                                    y, mu, wt), aic = aic,
                 mu.eta = stats$mu.eta, initialize = stats:::binomInitialize(family),
                 validmu = function(mu) all(is.finite(mu)) && all(0 <
                                                                    mu & mu < 1), valideta = stats$valideta), class = "family")
}

job::job({
  search <- sapply(seq(0, 0.999, 0.001), function(x){
    AIC(mgcv::bam(cbind(freq, trials) ~ s(SST,bs='tp',k = 11) +
                    s(WE, bs='cc', k = 28) +
                    s(CHLA, bs= 'tp', k = 6) +
                    Origin+
                    s(transmitter, bs = 're'),
                  rho = x, AR.start = ar_start,
                  family = qbn_hack(), discrete = T,
                  data=basswl))
  })})
# winner is rho = 0.865

job::job({m3 <- mgcv::bam(cbind(freq, trials) ~ s(SST,bs='ts',k = 11) +
                 s(WE, bs='cc',k = 28) +
                 s(CHLA, bs= 'ts',k = 6) +
                 Origin+
                   # s(freqL1_trials, bs = 'ts', k = 7)+
                 s(transmitter, bs = 're'),
               # rho = 0.865, AR.start = basswl$ar_start,
               family = binomial, discrete = T,
               data=basswl)})


SBLSBwL1_zip <- bam(freq ~ s(SST,bs='tp',k=11)+
                   s(WE, bs='cp',k=27)+
                   s(CHLA, bs= 'tp',k=6)+
                   s(freqL1,bs='tp',k=6)+
                   Origin+
                   s(transmitter, bs= 're'), discrete = T,
                 family=ziP,
                 data=basswl)

appraise(SBwL1)

# what I'm using for HGAM:

# SST
get_prior(freq | trials(7) ~
            s(SST, bs='cr', k = 11, m = 2) +
            s(SST, by = Origin, k=11, m = 1, bs='cs') +
            ar(gr = transmitter) +
            (1|Origin),
          family = zero_inflated_binomial(),
          data = basswl)


job::job({
  h_sst_cubic <- brm(
    freq | trials(7) ~
      s(SST, bs='cr', k = 11, m = 2) +
      s(SST, by = Origin, k=11, m = 1, bs='cs') +
      ar(gr = transmitter) +
      (1|Origin),

    family = zero_inflated_binomial(),
    data = basswl,
    iter = 500,
    warmup = 250,
    cores = 4,
    chains = 4,

    prior = c(
      set_prior('normal(0.5,0.4)', class = 'ar'),
      set_prior('normal(0,50)', class = 'b'),
      set_prior('normal(0,50)', class = 'sds')
    ),

    file = 'brms trials/HGAM_SST_cubic_lowiter.rds'

  )
})








job::job({
  SBh1L_ar2 <- sapply(seq(0, 0.99, 0.01), function(x){
    AIC(mgcv::bam(cbind(freq, trials) ~ s(SST,bs='tp',k=11, m=2)+
                    s(SST, by=Origin, k=11, m=1, bs='ts')+
                    s(Origin, bs='re') +
                    s(transmitter, bs = 're'),
                  rho = x, AR.start = ar_start,
                  family=quasibinomial, discrete = T,
                  data=basswl))
  })})
kkk <- mgcv::bam(cbind(freq, trials) ~ s(SST,bs='tp',k=11, m=2)+
                   s(SST, by=Origin, k=11, m=1, bs='ts')+
                   s(Origin, bs='re') +
                   s(transmitter, bs = 're'),
                 # s(freqL1,k=8),
                 rho = 0.86, AR.start = ar_start,
                 family=quasibinomial, discrete = T,
                 data=basswl)


# Week
job::job({
  h_week_cubic <- brm(
    freq | trials(7) ~
      s(WE,bs='cc',k=27, m=2)+
      s(WE, by=Origin, k=27, m=1, bs='cc') +
      ar(gr = transmitter) +
      (1|Origin),

    family = zero_inflated_binomial(),
    data = basswl,
    iter = 500,
    warmup = 250,
    cores = 4,
    chains = 4,

    prior = c(
      set_prior('normal(0.5,0.4)', class = 'ar'),
      set_prior('normal(0,10)', class = 'sds')
    ),

    file = 'brms trials/HGAM_week_cubic_lowiter.rds'

  )
})

SBh2L <- gam(freq ~ s(WE,bs='cc',k=27, m=2)+
               s(WE, by=Origin, k=27, m=1, bs='cc')+
               s(Origin, bs='re',k=6)+
               s(freqL1,k=8),
             family=ziP,
             data=basswl)

# CHL-A
job::job({
  h_chla_cubic <- brm(
    freq | trials(7) ~
      s(CHLA,bs='cr',k=6, m=2)+
      s(CHLA, by=Origin, k=6, m=1, bs='cr') +
      ar(gr = transmitter) +
      (1|Origin),

    family = zero_inflated_binomial(),
    data = basswl,
    iter = 500,
    warmup = 250,
    cores = 4,
    chains = 4,

    prior = c(
      set_prior('normal(0.5,0.4)', class = 'ar'),
      set_prior('normal(0,10)', class = 'sds')
    ),

    file = 'brms trials/HGAM_chla_cubic_lowiter.rds'

  )
})

SBh3L <- gam(freq ~ s(CHLA,bs='tp',k=6, m=2)+
               s(CHLA, by=Origin, k=6, m=1, bs='tp')+
               s(Origin, bs='re',k=6)+
               s(freqL1,k=8),
             family=ziP,
             data=basswl)
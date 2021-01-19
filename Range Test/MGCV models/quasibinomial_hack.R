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
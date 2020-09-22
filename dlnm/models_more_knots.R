library(dplyr); library(dlnm); library(mgcv)


sturgd <- read.csv("dlnm/data/sturg_w_lags.csv")


## For visualization and cross-validation later, we need to have all the data
##    saved in a list rather than a data.frame + 2 matrices (Q & L)
##    -MOB
sturgd_list <- as.list(sturgd[, 1:12])
sturgd_list[['Qs2']] <- as.matrix(sturgd[c(13:42)])
sturgd_list[['Ls2']] <- matrix(0:(ncol(sturgd_list$Qs2) - 1),
                               nrow(sturgd_list$Qs2),
                               ncol(sturgd_list$Qs2),
                               byrow = TRUE)


# Log transform
sturgd_list$CHLA <- log(sturgd_list$CHLA)

# Site needs to be a factor or mgcv throws an error
#   Note that since R v4.0, the default of read.csv has changed from
#   stringsAsFactors = T to stringsAsFactors = F -- this might break some of your
#   old code!! If you see the following error when running a model ...
#   "Error in names(dat) <- object$term : 'names' attribute [1] must be the same length as the vector [0]"
#   ... chances are that you need to change something from a character to a factor.
sturgd_list$Site <- as.factor(sturgd_list$Site)


# Pulling out all of the formulas to fit the models in parallel
forms <- c("freq~s(Qs2,Ls2, bs='cb', k = c(10, 14)) +
             s(Site, bs = 're') +
             s(Year, bs = 're') +
             s(CHLA, bs='tp', k = 9) +
             s(DOY, bs='cs', k = 26) +
             s(Depth, bs='tp', k = 6) +
             offset(log(d50.s))",
           "freq~s(Qs2,Ls2, bs='cb',k = c(10, 14)) +
             s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(CHLA,bs='tp', k = 9)+
             s(DOY,bs='cs', k = 26)+
             offset(log(d50.s))",
           "freq~s(Qs2,Ls2,bs='cb', k = c(10, 14)) +
             s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(CHLA,bs='tp', k = 9)+
             s(Depth,bs='tp', k = 6)+
             offset(log(d50.s))",
           "freq~s(Qs2,Ls2,bs='cb',k=c(10, 14)) +
             s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(DOY,bs='cs',k=26)+
             s(Depth,bs='tp', k = 6)+
             offset(log(d50.s))",
           "freq~s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(CHLA,bs='tp', k = 9)+
             s(DOY,bs='cs',k=26)+
             s(Depth,bs='tp', k = 6)+
             offset(log(d50.s))",
           "freq~s(Qs2,Ls2,bs='cb', k=c(10, 14)) +
           s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(CHLA,bs='tp', k = 9)+
             te(DOY,Depth,bs=c('cp', 'tp'),k=c(26,6))+
             offset(log(d50.s))",
           "freq~s(Qs2,Ls2,bs='cb', k = c(10, 14)) +
           s(Site, bs = 're')+
             s(Year, bs = 're')+
             te(DOY,Depth,bs=c('cp', 'tp'),k=c(26, 6))+
             offset(log(d50.s))",
           "freq~s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(CHLA,bs='tp', k = 9)+
             te(DOY,Depth,bs=c('cp', 'tp'),k=c(26, 6))+
             offset(log(d50.s))",
           "freq~s(Site, bs = 're')+
             s(Year, bs = 're')+
             te(DOY,Depth,bs=c('cp', 'tp'),k=c(26, 6))+
             offset(log(d50.s))",
           "freq~s(Qs2,Ls2,bs='cb', k = c(10, 14)) +
           s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(Depth,bs='tp', k = 6)+
              te(DOY,CHLA,bs=c('cp', 'tp'),k= c(26, 9))+
              offset(log(d50.s))",
           "freq~s(Qs2,Ls2,bs='cb',k = c(10, 14)) +
           s(Site, bs = 're')+
              s(Year, bs = 're')+
              te(DOY,CHLA,bs=c('cp', 'tp'),k=c(26, 9))+
              offset(log(d50.s))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(Depth,bs='tp', k = 6)+
              te(DOY,CHLA,bs=c('cp', 'tp'),k=c(26, 9))+
              offset(log(d50.s))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              te(DOY,CHLA,bs=c('cp', 'tp'),k=c(26, 9))+
              offset(log(d50.s))",
           "freq~s(Qs2,Ls2,bs='cb',k = c(10, 14)) +
           s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs='cs',k=26)+
              te(Depth,CHLA,bs=c('tp', 'tp'),k=c(6, 9))+
              offset(log(d50.s))",
           "freq~s(Qs2,Ls2,bs='cb', k = c(10, 14)) +
           s(Site, bs = 're')+
              s(Year, bs = 're')+
              te(Depth,CHLA,bs=c('tp', 'tp'),k=c(6, 9))+
              offset(log(d50.s))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs='cs',k=26)+
              te(Depth,CHLA,bs=c('tp', 'tp'),k=c(6, 9))+
              offset(log(d50.s))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              te(Depth,CHLA,bs=c('tp', 'tp'),k=c(6, 9))+
              offset(log(d50.s))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(CHLA,bs='tp',k=9)+
              s(Depth,bs='tp', k=6)+
              offset(log(d50.s))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs='cs',k=26)+
              s(Depth,bs='tp', k =6)+
              offset(log(d50.s))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs='cs',k=26)+
              s(CHLA,bs='tp',k=9)+
              offset(log(d50.s))",
           "freq~s(Qs2,Ls2,bs='cb', k = c(10, 14)) +
           s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(CHLA,bs='tp', k=9)+
              offset(log(d50.s))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(CHLA,bs='tp', k=9)+
              offset(log(d50.s))",
           "freq~s(Qs2,Ls2,bs='cb', k = c(10, 14)) +
           s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs='cs',k=26)+
              offset(log(d50.s))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs='cs',k=26)+
              offset(log(d50.s))",
           "freq~s(Qs2,Ls2,bs='cb',k = c(10, 14)) +
           s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(Depth,bs='tp', k=6)+
              offset(log(d50.s))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(Depth,bs='tp', k=6)+
              offset(log(d50.s))",
           "freq ~ s(Qs2, Ls2, bs='cb', k = c(10, 14)) +
              s(Site, bs = 're')+
              s(Year, bs = 're')+
              offset(log(d50.s))")


forms <- gsub('\n', ' ', forms)


# Set up cluster
library(parallel)
cl <- detectCores(logical = F) - 1
cl <- makeCluster(cl)
clusterEvalQ(cl, list(library(dlnm), library(mgcv)))
clusterExport(cl, 'sturgd_list')


# Fit models
library(pbapply)
as_models <- pblapply(forms,
                      function(.){
                        gam(formula = as.formula(.),
                            data = sturgd_list,
                            family = ziP,
                            method = 'REML')
                      },
                      cl = cl)

# Close cluster
stopCluster(cl)

names(as_models) <- paste0('SM', 1:27)

# saveRDS(as_models, 'dlnm/sturg_lag_models.RDS')
# as_models <- readRDS('dlnm/sturg_lag_models.RDS')

summary(as_models$SM1)

s_ic <- data.frame(AIC = sapply(as_models, AIC))
s_ic$dAIC <- s_ic$AIC - min(s_ic$AIC)
s_ic <- s_ic[order(s_ic$dAIC),]
s_ic$wAIC <- exp(-0.5 * s_ic$dAIC)
s_ic$wAIC <- s_ic$wAIC / sum(s_ic$wAIC)
s_ic$cwAIC <- cumsum(s_ic$wAIC)
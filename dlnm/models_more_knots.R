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
sturgd_list$Year <- as.factor(sturgd_list$Year)

# Pulling out all of the formulas to fit the models in parallel
forms <- c("freq~s(Qs2,Ls2, bs='cb', k = c(11, 15)) +
             s(Site, bs = 're') +
             s(Year, bs = 're') +
             s(CHLA, bs='tp', k = 10) +
             s(DOY, bs='cc', k = 27) +
             s(Depth, bs='tp', k = 15) +
             offset(log(d50.s))",
           "freq~s(Qs2,Ls2, bs='cb',k = c(11, 15)) +
             s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(CHLA,bs='tp', k = 10)+
             s(DOY,bs='cc', k = 27)+
             offset(log(d50.s))",
           "freq~s(Qs2,Ls2,bs='cb', k = c(11, 15)) +
             s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(CHLA,bs='tp', k = 10)+
             s(Depth,bs='tp', k = 15)+
             offset(log(d50.s))",
           "freq~s(Qs2,Ls2,bs='cb',k=c(11, 15)) +
             s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(DOY,bs='cc',k=27)+
             s(Depth,bs='tp', k = 15)+
             offset(log(d50.s))",
           "freq~s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(CHLA,bs='tp', k = 10)+
             s(DOY,bs='cc',k=27)+
             s(Depth,bs='tp', k = 15)+
             offset(log(d50.s))",
           "freq~s(Qs2,Ls2,bs='cb', k=c(11, 15)) +
           s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(CHLA,bs='tp', k = 10)+
             te(DOY,Depth,bs=c('cc', 'tp'),k=c(27,15))+
             offset(log(d50.s))",
           "freq~s(Qs2,Ls2,bs='cb', k = c(11, 15)) +
           s(Site, bs = 're')+
             s(Year, bs = 're')+
             te(DOY,Depth,bs=c('cc', 'tp'),k=c(27, 15))+
             offset(log(d50.s))",
           "freq~s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(CHLA,bs='tp', k = 10)+
             te(DOY,Depth,bs=c('cc', 'tp'),k=c(27, 15))+
             offset(log(d50.s))",
           "freq~s(Site, bs = 're')+
             s(Year, bs = 're')+
             te(DOY,Depth,bs=c('cc', 'tp'),k=c(27, 15))+
             offset(log(d50.s))",
           "freq~s(Qs2,Ls2,bs='cb', k = c(11, 15)) +
           s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(Depth,bs='tp', k = 15)+
              te(DOY,CHLA,bs=c('cc', 'tp'),k= c(27, 10))+
              offset(log(d50.s))",
           "freq~s(Qs2,Ls2,bs='cb',k = c(11, 15)) +
           s(Site, bs = 're')+
              s(Year, bs = 're')+
              te(DOY,CHLA,bs=c('cc', 'tp'),k=c(27, 10))+
              offset(log(d50.s))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(Depth,bs='tp', k = 15)+
              te(DOY,CHLA,bs=c('cc', 'tp'),k=c(27, 10))+
              offset(log(d50.s))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              te(DOY,CHLA,bs=c('cc', 'tp'),k=c(27, 10))+
              offset(log(d50.s))",
           "freq~s(Qs2,Ls2,bs='cb',k = c(11, 15)) +
           s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs='cc',k=27)+
              te(Depth,CHLA,bs=c('tp', 'tp'),k=c(15, 10))+
              offset(log(d50.s))",
           "freq~s(Qs2,Ls2,bs='cb', k = c(11, 15)) +
           s(Site, bs = 're')+
              s(Year, bs = 're')+
              te(Depth,CHLA,bs=c('tp', 'tp'),k=c(15, 10))+
              offset(log(d50.s))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs='cc',k=27)+
              te(Depth,CHLA,bs=c('tp', 'tp'),k=c(15, 10))+
              offset(log(d50.s))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              te(Depth,CHLA,bs=c('tp', 'tp'),k=c(15, 10))+
              offset(log(d50.s))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(CHLA,bs='tp',k=10)+
              s(Depth,bs='tp', k=15)+
              offset(log(d50.s))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs='cc',k=27)+
              s(Depth,bs='tp', k =15)+
              offset(log(d50.s))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs='cc',k=27)+
              s(CHLA,bs='tp',k=10)+
              offset(log(d50.s))",
           "freq~s(Qs2,Ls2,bs='cb', k = c(11, 15)) +
           s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(CHLA,bs='tp', k=10)+
              offset(log(d50.s))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(CHLA,bs='tp', k=10)+
              offset(log(d50.s))",
           "freq~s(Qs2,Ls2,bs='cb', k = c(11, 15)) +
           s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs='cc',k=27)+
              offset(log(d50.s))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs='cc',k=27)+
              offset(log(d50.s))",
           "freq~s(Qs2,Ls2,bs='cb',k = c(11, 15)) +
           s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(Depth,bs='tp', k=15)+
              offset(log(d50.s))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(Depth,bs='tp', k=15)+
              offset(log(d50.s))",
           "freq ~ s(Qs2, Ls2, bs='cb', k = c(11, 15)) +
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

# saveRDS(as_models, 'dlnm/sturg_lag_models_moreknots.RDS')
# as_models <- readRDS('dlnm/sturg_lag_models_moreknots.RDS')

summary(as_models$SM1)

s_ic <- data.frame(AIC = sapply(as_models, AIC))
s_ic$dAIC <- s_ic$AIC - min(s_ic$AIC)
s_ic <- s_ic[order(s_ic$dAIC),]
s_ic$wAIC <- exp(-0.5 * s_ic$dAIC)
s_ic$wAIC <- s_ic$wAIC / sum(s_ic$wAIC)
s_ic$cwAIC <- cumsum(s_ic$wAIC)




bassd2 <- read.csv('dlnm/data/bass_w_lags.csv')

# Converted to list, as for sturgeon, above   -MOB
bassd_list <- as.list(bassd2[, 1:12])
bassd_list[['Qb2']] <- as.matrix(bassd2[c(13:42)])
bassd_list[['Lb2']] <- matrix(0:(ncol(bassd_list$Qb2) - 1),
                              nrow(bassd_list$Qb2),
                              ncol(bassd_list$Qb2),
                              byrow = TRUE)


# Log transform
bassd_list$CHLA <- log(bassd_list$CHLA)

# Site and Year to factor, as for sturgeon, above    -MOB
bassd_list$Site <- as.factor(bassd_list$Site)
bassd_list$Year <- as.factor(bassd_list$Year)

forms <- c("freq~s(Qb2,Lb2, bs='cb', k = c(11, 15)) +
             s(Site, bs = 're') +
             s(Year, bs = 're') +
             s(CHLA, bs='tp', k = 10) +
             s(DOY, bs='cc', k = 27) +
             s(Depth, bs='tp', k = 15) +
             offset(log(d50.b))",
           "freq~s(Qb2,Lb2, bs='cb',k = c(11, 15)) +
             s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(CHLA,bs='tp', k = 10)+
             s(DOY,bs='cc', k = 27)+
             offset(log(d50.b))",
           "freq~s(Qb2,Lb2,bs='cb', k = c(11, 15)) +
             s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(CHLA,bs='tp', k = 10)+
             s(Depth,bs='tp', k = 15)+
             offset(log(d50.b))",
           "freq~s(Qb2,Lb2,bs='cb',k=c(11, 15)) +
             s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(DOY,bs='cc',k=27)+
             s(Depth,bs='tp', k = 15)+
             offset(log(d50.b))",
           "freq~s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(CHLA,bs='tp', k = 10)+
             s(DOY,bs='cc',k=27)+
             s(Depth,bs='tp', k = 15)+
             offset(log(d50.b))",
           "freq~s(Qb2,Lb2,bs='cb', k=c(11, 15)) +
           s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(CHLA,bs='tp', k = 10)+
             te(DOY,Depth,bs=c('cc', 'tp'),k=c(27,15))+
             offset(log(d50.b))",
           "freq~s(Qb2,Lb2,bs='cb', k = c(11, 15)) +
           s(Site, bs = 're')+
             s(Year, bs = 're')+
             te(DOY,Depth,bs=c('cc', 'tp'),k=c(27, 15))+
             offset(log(d50.b))",
           "freq~s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(CHLA,bs='tp', k = 10)+
             te(DOY,Depth,bs=c('cc', 'tp'),k=c(27, 15))+
             offset(log(d50.b))",
           "freq~s(Site, bs = 're')+
             s(Year, bs = 're')+
             te(DOY,Depth,bs=c('cc', 'tp'),k=c(27, 15))+
             offset(log(d50.b))",
           "freq~s(Qb2,Lb2,bs='cb', k = c(11, 15)) +
           s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(Depth,bs='tp', k = 15)+
              te(DOY,CHLA,bs=c('cc', 'tp'),k= c(27, 10))+
              offset(log(d50.b))",
           "freq~s(Qb2,Lb2,bs='cb',k = c(11, 15)) +
           s(Site, bs = 're')+
              s(Year, bs = 're')+
              te(DOY,CHLA,bs=c('cc', 'tp'),k=c(27, 10))+
              offset(log(d50.b))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(Depth,bs='tp', k = 15)+
              te(DOY,CHLA,bs=c('cc', 'tp'),k=c(27, 10))+
              offset(log(d50.b))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              te(DOY,CHLA,bs=c('cc', 'tp'),k=c(27, 10))+
              offset(log(d50.b))",
           "freq~s(Qb2,Lb2,bs='cb',k = c(11, 15)) +
           s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs='cc',k=27)+
              te(Depth,CHLA,bs=c('tp', 'tp'),k=c(15, 10))+
              offset(log(d50.b))",
           "freq~s(Qb2,Lb2,bs='cb', k = c(11, 15)) +
           s(Site, bs = 're')+
              s(Year, bs = 're')+
              te(Depth,CHLA,bs=c('tp', 'tp'),k=c(15, 10))+
              offset(log(d50.b))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs='cc',k=27)+
              te(Depth,CHLA,bs=c('tp', 'tp'),k=c(15, 10))+
              offset(log(d50.b))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              te(Depth,CHLA,bs=c('tp', 'tp'),k=c(15, 10))+
              offset(log(d50.b))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(CHLA,bs='tp',k=10)+
              s(Depth,bs='tp', k=15)+
              offset(log(d50.b))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs='cc',k=27)+
              s(Depth,bs='tp', k =15)+
              offset(log(d50.b))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs='cc',k=27)+
              s(CHLA,bs='tp',k=10)+
              offset(log(d50.b))",
           "freq~s(Qb2,Lb2,bs='cb', k = c(11, 15)) +
           s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(CHLA,bs='tp', k=10)+
              offset(log(d50.b))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(CHLA,bs='tp', k=10)+
              offset(log(d50.b))",
           "freq~s(Qb2,Lb2,bs='cb', k = c(11, 15)) +
           s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs='cc',k=27)+
              offset(log(d50.b))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs='cc',k=27)+
              offset(log(d50.b))",
           "freq~s(Qb2,Lb2,bs='cb',k = c(11, 15)) +
           s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(Depth,bs='tp', k=15)+
              offset(log(d50.b))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(Depth,bs='tp', k=15)+
              offset(log(d50.b))",
           "freq ~ s(Qb2, Lb2, bs='cb', k = c(11, 15)) +
              s(Site, bs = 're')+
              s(Year, bs = 're')+
              offset(log(d50.b))")

# remove the carriage returns
forms <- gsub('\n', ' ', forms)


library(parallel)
cl <- detectCores(logical = F) - 1
cl <- makeCluster(cl)
clusterEvalQ(cl, list(library(dlnm), library(mgcv)))
clusterExport(cl, 'bassd_list')

library(pbapply)
sb_models <- pblapply(forms,
                      function(.){
                        # One of the models doesn't want to converge, so using
                        # tryCatach to make it save the error and move on rather
                        # than killing the whole process.     -MOB
                        tryCatch({
                          gam(formula = as.formula(.),
                              data = bassd_list,
                              family = ziP,
                              method = 'REML')
                        }, error = function(e) return(paste('Error in', as.formula(.))))
                      },
                      cl = cl)
stopCluster(cl)

names(sb_models) <-  paste0('BM', 1:27)

# saveRDS(sb_models, 'dlnm/sb_lag_models_moreknots.RDS')
# sb_models <- readRDS('dlnm/sb_lag_models_moreknots.RDS')

summary(sb_models$BM1)

b_ic <- data.frame(AIC = sapply(sb_models, AIC))
b_ic$dAIC <- b_ic$AIC - min(b_ic$AIC)
b_ic <- b_ic[order(b_ic$dAIC),]
b_ic$wAIC <- exp(-0.5 * b_ic$dAIC)
b_ic$wAIC <- b_ic$wAIC / sum(b_ic$wAIC)
b_ic$cwAIC <- cumsum(b_ic$wAIC)

## Added line to load packages      -MOB
library(dplyr); library(dlnm); library(mgcv)

## I edited the column name of "Day" in sturg_full.csv directly before loading in,
##    as the column had wound up called "Day.x".    -MOB
# sturgn<-read.csv("dlnm/data/Sturg_full.csv",header=T)
# sturgn <- subset(sturgn,select=-c(1))
#
# sturgn[,2]<-as.POSIXct(sturgn[,2],format="%Y-%m-%d")
#
# sturgd<- aggregate(sturgn$freq,
#                    by = list(Site=sturgn$Site,
#                              Day=sturgn$Day,
#                              DOY=sturgn$DOY,
#                              CHLA=sturgn$CHLA,
#                              SST=sturgn$SST,
#                              Depth=sturgn$Depth,
#                              d50.s=sturgn$d50.s,
#                              time=sturgn$time,
#                              Season=sturgn$Season,
#                              Seasonb=sturgn$Seasonb,
#                              Year=sturgn$Year),FUN= "sum")
#
# colnames(sturgd)[12] <- "freq"
# sturgd$Year<-as.factor(sturgd$Year)
#
# Lags<-read.csv("dlnm/data/AllLags.csv",header=T)
#
# ## Also removing SST column -- it's redundant after the join. -- MOB
# Lags <- subset(Lags,select=-c(1, 4))
#
# Lags[,2]<-as.POSIXct(Lags[,2],format="%Y-%m-%d")
#
# sturgda<-left_join(sturgd, Lags, by=c("Day","Site"))
# sturgd2<-sturgda[!duplicated(sturgda),]
# #
# sturgd2 <- sturgd2[with(sturgd2,order(time,Site)),]
# write.csv(sturgd2, 'dlnm/sturg_w_lags.csv', row.names = F)

sturgd2 <- read.csv("dlnm/data/sturg_w_lags.csv")


## For visualization and cross-validation later, we need to have all the data
##    saved in a list rather than a data.frame + 2 matrices (Q & L)
##    -MOB
sturgd_list <- as.list(sturgd2[, 1:12])
sturgd_list[['Qs2']] <- as.matrix(sturgd2[c(13:42)])
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
forms <- c("freq~s(Qs2,Ls2,bs='cb',k=c(6, 6)) +s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(CHLA,bs='tp', k = 6)+
             s(DOY,bs='cs',k=6)+
             s(Depth,bs='tp', k = 6)+
             offset(log(d50.s))",
           "freq~s(Qs2,Ls2,bs='cb',k=c(6, 6))+s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(CHLA,bs='tp', k = 6)+
             s(DOY,bs='cs',k=6)+
             offset(log(d50.s))",
           "freq~s(Qs2,Ls2,bs='cb',k=c(6, 6))+s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(CHLA,bs='tp', k = 6)+
             s(Depth,bs='tp', k = 6)+
             offset(log(d50.s))",
           "freq~s(Qs2,Ls2,bs='cb',k=c(6, 6))+s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(DOY,bs='cs',k=6)+
             s(Depth,bs='tp', k = 6)+
             offset(log(d50.s))",
           "freq~s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(CHLA,bs='tp', k = 6)+
             s(DOY,bs='cs',k=6)+
             s(Depth,bs='tp', k = 6)+
             offset(log(d50.s))",
           "freq~s(Qs2,Ls2,bs='cb',k=c(6, 6))+s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(CHLA,bs='tp', k = 6)+
             te(DOY,Depth,bs=c('cp', 'tp'),k=c(6,6))+
             offset(log(d50.s))",
           "freq~s(Qs2,Ls2,bs='cb',k=c(6, 6))+s(Site, bs = 're')+
             s(Year, bs = 're')+
             te(DOY,Depth,bs=c('cp', 'tp'),k=c(6, 6))+
             offset(log(d50.s))",
           "freq~s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(CHLA,bs='tp', k = 6)+
             te(DOY,Depth,bs=c('cp', 'tp'),k=c(6, 6))+
             offset(log(d50.s))",
           "freq~s(Site, bs = 're')+
             s(Year, bs = 're')+
             te(DOY,Depth,bs=c('cp', 'tp'),k=c(6, 6))+
             offset(log(d50.s))",
           "freq~s(Qs2,Ls2,bs='cb',k=c(6, 6))+s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(Depth,bs='tp', k = 6)+
              te(DOY,CHLA,bs=c('cp', 'tp'),k= c(6, 6))+
              offset(log(d50.s))",
           "freq~s(Qs2,Ls2,bs='cb',k=c(6, 6))+s(Site, bs = 're')+
              s(Year, bs = 're')+
              te(DOY,CHLA,bs=c('cp', 'tp'),k=c(6, 6))+
              offset(log(d50.s))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(Depth,bs='tp', k = 6)+
              te(DOY,CHLA,bs=c('cp', 'tp'),k=c(6, 6))+
              offset(log(d50.s))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              te(DOY,CHLA,bs=c('cp', 'tp'),k=c(6, 6))+
              offset(log(d50.s))",
           "freq~s(Qs2,Ls2,bs='cb',k=c(6, 6))+s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs='cs',k=6)+
              te(Depth,CHLA,bs=c('tp', 'tp'),k=c(6, 6))+
              offset(log(d50.s))",
           "freq~s(Qs2,Ls2,bs='cb',k=c(6, 6))+s(Site, bs = 're')+
              s(Year, bs = 're')+
              te(Depth,CHLA,bs=c('tp', 'tp'),k=c(6, 6))+
              offset(log(d50.s))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs='cs',k=6)+
              te(Depth,CHLA,bs=c('tp', 'tp'),k=c(6, 6))+
              offset(log(d50.s))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              te(Depth,CHLA,bs=c('tp', 'tp'),k=c(6, 6))+
              offset(log(d50.s))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(CHLA,bs='tp',k=6)+
              s(Depth,bs='tp', k=6)+
              offset(log(d50.s))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs='cs',k=6)+
              s(Depth,bs='tp', k =6)+
              offset(log(d50.s))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs='cs',k=6)+
              s(CHLA,bs='tp',k=6)+
              offset(log(d50.s))",
           "freq~s(Qs2,Ls2,bs='cb',k=c(6, 6))+s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(CHLA,bs='tp', k=6)+
              offset(log(d50.s))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(CHLA,bs='tp', k=6)+
              offset(log(d50.s))",
           "freq~s(Qs2,Ls2,bs='cb',k=c(6, 6))+s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs='cs',k=6)+
              offset(log(d50.s))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs='cs',k=6)+
              offset(log(d50.s))",
           "freq~s(Qs2,Ls2,bs='cb',k = c(6, 6))+s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(Depth,bs='tp', k=6)+
              offset(log(d50.s))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(Depth,bs='tp', k=6)+
              offset(log(d50.s))",
           "freq~s(Qs2,Ls2,bs='cb',k = c(6, 6))+s(Site, bs = 're')+
              s(Year, bs = 're')+
              offset(log(d50.s))")

# remove the carriage returns (I was too lazy to do anything but copy/paste)
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

## SM6 and SM7 seem to be somewhat-comparabe candidates.



## Now for striped bass
# bassn<-read.csv("dlnm/data/Bass_full.csv",header=T)
# bassn <- subset(bassn,select=-c(1))
#
# bassn[,2]<-as.POSIXct(bassn[,2],format="%Y-%m-%d")
#
# bassd<- aggregate(bassn$freq,
#                   by = list(Site=bassn$Site,
#                             Day=bassn$Day,
#                             DOY=bassn$DOY,
#                             CHLA=bassn$CHLA,
#                             SST=bassn$SST,
#                             Depth=bassn$Depth,
#                             d50.b=bassn$d50.b,
#                             time=bassn$time,
#                             Season=bassn$Season,
#                             Seasonb=bassn$Seasonb,
#                             Year=bassn$Year),FUN= "sum")
#
# colnames(bassd)[12] <- "freq"
# bassd$Year<-as.factor(bassd$Year)
#
# bassda<-left_join(bassd, Lags, by=c("Day","Site"))
# bassd2<-bassda[!duplicated(bassda),]
# #
# bassd2 <- bassd2[with(bassd2,order(time,Site)),]
# write.csv(bassd2, 'dlnm/data/bass_w_lags.csv', row.names = F)

bassd2 <- read.csv('dlnm/data/bass_w_lags.csv')

# Converted to list, as for sturgeon, above   -MOB
bassd_list <- as.list(bassd2[, 1:12])
bassd_list[['Qb2']] <- as.matrix(bassd2[c(13:42)])
bassd_list[['Lb2']] <- matrix(0:(ncol(bassd_list$Qb2) - 1),
                              nrow(bassd_list$Qb2),
                              ncol(bassd_list$Qb2),
                              byrow = TRUE)


ts.start <- proc.time()
# Log transform
bassd_list$CHLA <- log(bassd_list$CHLA)

# Site to factor, as for sturgeon, above    -MOB
bassd_list$Site <- as.factor(bassd_list$Site)


forms <- c("freq~s(Qb2,Lb2,bs='cb',k = c(6, 6))+s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(CHLA,bs='tp', k=6)+
             s(DOY,bs='cs',k=6)+
             s(Depth,bs='tp', k=6)+
             offset(log(d50.b))",
           "freq~s(Qb2,Lb2,bs='cb',k = c(6, 6))+s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(CHLA,bs='tp', k=6)+
             s(DOY,bs='cs',k=6)+
             offset(log(d50.b))",
           "freq~s(Qb2,Lb2,bs='cb',k = c(6, 6))+s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(CHLA,bs='tp', k=6)+
             s(Depth,bs='tp', k=6)+
             offset(log(d50.b))",
           "freq~s(Qb2,Lb2,bs='cb',k = c(6, 6))+s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(DOY,bs='cs',k=6)+
             s(Depth,bs='tp', k=6)+
             offset(log(d50.b))",
           "freq~s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(CHLA,bs='tp', k=6)+
             s(DOY,bs='cs',k=6)+
             s(Depth,bs='tp', k=6)+
             offset(log(d50.b))",
           "freq~s(Qb2,Lb2,bs='cb',k = c(6, 6))+s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(CHLA,bs='tp', k=6)+
             te(DOY,Depth,bs=c('cp', 'tp'),k=c(6, 6))+
             offset(log(d50.b))",
           "freq~s(Qb2,Lb2,bs='cb',k = c(6, 6))+s(Site, bs = 're')+
             s(Year, bs = 're')+
             te(DOY,Depth,bs=c('cp', 'tp'),k=c(6, 6))+
             offset(log(d50.b))",
           "freq~s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(CHLA,bs='tp', k=6)+
             te(DOY,Depth,bs=c('cp', 'tp'),k=c(6, 6))+
             offset(log(d50.b))",
           "freq~s(Site, bs = 're')+
             s(Year, bs = 're')+
             te(DOY,Depth,bs=c('cp', 'tp'),k=c(6, 6))+
             offset(log(d50.b))",
           "freq~s(Qb2,Lb2,bs='cb',k = c(6, 6))+s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(Depth,bs='tp', k=6)+
              te(DOY,CHLA,bs=c('cp', 'tp'),k=c(6, 6))+
              offset(log(d50.b))",
           "freq~s(Qb2,Lb2,bs='cb',k = c(6, 6))+s(Site, bs = 're')+
              s(Year, bs = 're')+
              te(DOY,CHLA,bs=c('cp', 'tp'),k=c(6, 6))+
              offset(log(d50.b))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(Depth,bs='tp', k=6)+
              te(DOY,CHLA,bs=c('cp', 'tp'),k=c(6, 6))+
              offset(log(d50.b))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              te(DOY,CHLA,bs=c('cp', 'tp'),k=c(6, 6))+
              offset(log(d50.b))",
           "freq~s(Qb2,Lb2,bs='cb',k = c(6, 6))+s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs='cs',k=6)+
              te(Depth,CHLA,bs=c('tp', 'tp'),k=6)+
              offset(log(d50.b))",
           "freq~s(Qb2,Lb2,bs='cb',k = c(6, 6))+s(Site, bs = 're')+
              s(Year, bs = 're')+
              te(Depth,CHLA,bs=c('tp', 'tp'),k=c(6, 6))+
              offset(log(d50.b))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs='cs',k=6)+
              te(Depth,CHLA,bs=c('tp', 'tp'),k=c(6, 6))+
              offset(log(d50.b))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              te(Depth,CHLA,bs=c('tp', 'tp'),k=c(6, 6))+
              offset(log(d50.b))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(CHLA,bs='tp',k=6)+
              s(Depth,bs='tp', k=6)+
              offset(log(d50.b))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs='cs',k=6)+
              s(Depth,bs='tp', k=6)+
              offset(log(d50.b))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(CHLA, bs = 'tp', k = 6) +
              s(DOY, bs = 'cs', k = 6) +
              offset(log(d50.b))",
           "freq~s(Qb2,Lb2,bs='cb',k = c(6, 6))+s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(CHLA,bs='tp', k=6)+
              offset(log(d50.b))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(CHLA,bs='tp', k=6)+
              offset(log(d50.b))",
           "freq~s(Qb2,Lb2,bs='cb',k = c(6, 6))+s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs='cs',k=6)+
              offset(log(d50.b))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs='cs',k=6)+
              offset(log(d50.b))",
           "freq~s(Qb2,Lb2,bs='cb',k = c(6, 6))+s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(Depth,bs='tp', k=6)+
              offset(log(d50.b))",
           "freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(Depth,bs='tp', k=6)+
              offset(log(d50.b))",
           "freq~s(Qb2,Lb2,bs='cb',k = c(6, 6))+s(Site, bs = 're')+
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

# saveRDS(sb_models, 'dlnm/sb_lag_models.RDS')

summary(sb_models$BM1)

b_ic <- data.frame(AIC = sapply(sb_models, AIC))
b_ic$dAIC <- b_ic$AIC - min(b_ic$AIC)
b_ic <- b_ic[order(b_ic$dAIC),]
b_ic$wAIC <- exp(-0.5 * b_ic$dAIC)
b_ic$wAIC <- b_ic$wAIC / sum(b_ic$wAIC)
b_ic$cwAIC <- cumsum(b_ic$wAIC)

## BM6 is clear winner

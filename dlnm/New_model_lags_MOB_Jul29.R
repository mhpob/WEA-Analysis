## Added line to load packages      -MOB 20200729
library(dplyr); library(dlnm); library(mgcv)

## I edited the column name of "Day" in sturg_full.csv directly before loading in,
##    as the column had wound up called "Day.x".    -MOB 20200729
sturgn<-read.csv("dlnm/Sturg_full.csv",header=T)
sturgn <- subset(sturgn,select=-c(1))

sturgn[,2]<-as.POSIXct(sturgn[,2],format="%Y-%m-%d")

sturgd<- aggregate(sturgn$freq,
                   by = list(Site=sturgn$Site,
                             Day=sturgn$Day,
                             DOY=sturgn$DOY,
                             CHLA=sturgn$CHLA,
                             SST=sturgn$SST,
                             Depth=sturgn$Depth,
                             d50.s=sturgn$d50.s,
                             time=sturgn$time,
                             Season=sturgn$Season,
                             Seasonb=sturgn$Seasonb,
                             Year=sturgn$Year),FUN= "sum")

colnames(sturgd)[12] <- "freq"
sturgd$Year<-as.factor(sturgd$Year)

Lags<-read.csv("dlnm/AllLags.csv",header=T)

## Also removing SST column -- it's redundant after the join. -- MOB 20200729
Lags <- subset(Lags,select=-c(1, 4))

Lags[,2]<-as.POSIXct(Lags[,2],format="%Y-%m-%d")

sturgda<-left_join(sturgd, Lags, by=c("Day","Site"))
sturgd2<-sturgda[!duplicated(sturgda),]
#
sturgd2 <- sturgd2[with(sturgd2,order(time,Site)),]


## For visualization and cross-validation later, we need to have all the data
##    saved in a list rather than a data.frame + 2 matrices (Q & L)
##    -MOB 20200729
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

SM1 <- gam(freq~s(Qs2,Ls2,bs="cb",k=10)+s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(CHLA,bs='tp')+
             s(DOY,bs="cs",k=6)+
             s(Depth,bs='tp')+
             offset(log(d50.s)),
           data=sturgd2,family=ziP,method="REML")

SM2 <- gam(freq~s(Qs2,Ls2,bs="cb",k=10)+s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(CHLA,bs='tp')+
             s(DOY,bs="cs",k=6)+
             offset(log(d50.s)),
           data=sturgd2,family=ziP,method="REML")

SM3 <- gam(freq~s(Qs2,Ls2,bs="cb",k=10)+s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(CHLA,bs='tp')+
             s(Depth,bs='tp')+
             offset(log(d50.s)),
           data=sturgd2,family=ziP,method="REML")

SM4 <- gam(freq~s(Qs2,Ls2,bs="cb",k=10)+s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(DOY,bs="cs",k=6)+
             s(Depth,bs='tp')+
             offset(log(d50.s)),
           data=sturgd2,family=ziP,method="REML")

SM5 <- gam(freq~s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(CHLA,bs='tp')+
             s(DOY,bs="cs",k=6)+
             s(Depth,bs='tp')+
             offset(log(d50.s)),
           data=sturgd2,family=ziP,method="REML")

SM6 <- gam(freq~s(Qs2,Ls2,bs="cb",k=10)+s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(CHLA,bs='tp')+
             t2(DOY,Depth,bs=c("cp", "tp"),k=6)+
             offset(log(d50.s)),
           data=sturgd2,family=ziP,method="REML")

SM7 <- gam(freq~s(Qs2,Ls2,bs="cb",k=10)+s(Site, bs = 're')+
             s(Year, bs = 're')+
             t2(DOY,Depth,bs=c("cp", "tp"),k=6)+
             offset(log(d50.s)),
           data=sturgd2,family=ziP,method="REML")

SM8 <- gam(freq~s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(CHLA,bs='tp')+
             t2(DOY,Depth,bs=c("cp", "tp"),k=6)+
             offset(log(d50.s)),
           data=sturgd2,family=ziP,method="REML")

SM9 <- gam(freq~s(Site, bs = 're')+
             s(Year, bs = 're')+
             t2(DOY,Depth,bs=c("cp", "tp"),k=6)+
             offset(log(d50.s)),
           data=sturgd2,family=ziP,method="REML")

SM10 <- gam(freq~s(Qs2,Ls2,bs="cb",k=10)+s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(Depth,bs='tp')+
              t2(DOY,CHLA,bs=c("cp", "tp"),k=6)+
              offset(log(d50.s)),
            data=sturgd2,family=ziP,method="REML")

SM11 <- gam(freq~s(Qs2,Ls2,bs="cb",k=10)+s(Site, bs = 're')+
              s(Year, bs = 're')+
              t2(DOY,CHLA,bs=c("cp", "tp"),k=6)+
              offset(log(d50.s)),
            data=sturgd2,family=ziP,method="REML")

SM12 <- gam(freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(Depth,bs='tp')+
              t2(DOY,CHLA,bs=c("cp", "tp"),k=6)+
              offset(log(d50.s)),
            data=sturgd2,family=ziP,method="REML")

SM13 <- gam(freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              t2(DOY,CHLA,bs=c("cp", "tp"),k=6)+
              offset(log(d50.s)),
            data=sturgd2,family=ziP,method="REML")

SM14 <- gam(freq~s(Qs2,Ls2,bs="cb",k=10)+s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs="cs",k=6)+
              t2(Depth,CHLA,bs=c("tp", "tp"),k=6)+
              offset(log(d50.s)),
            data=sturgd2,family=ziP,method="REML")

SM15 <- gam(freq~s(Qs2,Ls2,bs="cb",k=10)+s(Site, bs = 're')+
              s(Year, bs = 're')+
              t2(Depth,CHLA,bs=c("tp", "tp"),k=6)+
              offset(log(d50.s)),
            data=sturgd2,family=ziP,method="REML")

SM16 <- gam(freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs="cs",k=6)+
              t2(Depth,CHLA,bs=c("tp", "tp"),k=6)+
              offset(log(d50.s)),
            data=sturgd2,family=ziP,method="REML")

SM17 <- gam(freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              t2(Depth,CHLA,bs=c("tp", "tp"),k=6)+
              offset(log(d50.s)),
            data=sturgd2,family=ziP,method="REML")

SM18 <- gam(freq~s(Qs2,Ls2,bs="cb",k=10)+s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs="cs",k=6)+
              s(Depth,bs='tp')+
              offset(log(d50.s)),
            data=sturgd2,family=ziP,method="REML")

SM19 <- gam(freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs="cs",k=6)+
              s(Depth,bs='tp')+
              offset(log(d50.s)),
            data=sturgd2,family=ziP,method="REML")

SM20 <- gam(freq~s(Qs2,Ls2,bs="cb",k=10)+s(Site, bs = 're')+
              s(Year, bs = 're')+
              t2(DOY,Depth,bs=c("cp", "tp"),k=6)+
              offset(log(d50.s)),
            data=sturgd2,family=ziP,method="REML")

SM21 <- gam(freq~s(Qs2,Ls2,bs="cb",k=10)+s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(CHLA,bs='tp')+
              offset(log(d50.s)),
            data=sturgd2,family=ziP,method="REML")

SM22 <- gam(freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(CHLA,bs='tp')+
              offset(log(d50.s)),
            data=sturgd2,family=ziP,method="REML")

SM23 <- gam(freq~s(Qs2,Ls2,bs="cb",k=10)+s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs="cs",k=6)+
              offset(log(d50.s)),
            data=sturgd2,family=ziP,method="REML")

SM24 <- gam(freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs="cs",k=6)+
              offset(log(d50.s)),
            data=sturgd2,family=ziP,method="REML")

SM25 <- gam(freq~s(Qs2,Ls2,bs="cb",k=10)+s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(Depth,bs='tp')+
              offset(log(d50.s)),
            data=sturgd2,family=ziP,method="REML")

SM26 <- gam(freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(Depth,bs='tp')+
              offset(log(d50.s)),
            data=sturgd2,family=ziP,method="REML")


summary(SM1)
AIC(SM1,SM2,SM3,SM4,SM5,SM6,SM7,SM8,SM9,SM10,
    SM11,SM12,SM13,SM14,SM15,SM16,SM17,SM18,SM19,SM20,
    SM21,SM22,SM23,SM24,SM25,SM26)

summary(SM6)
summary(SM7)
summary(SM10)
summary(SM11)
summary(SM1)

plot(SM6)

## Now for striped bass
bassn<-read.csv("dlnm/Bass_full.csv",header=T)
bassn <- subset(bassn,select=-c(1))

bassn[,2]<-as.POSIXct(bassn[,2],format="%Y-%m-%d")

bassd<- aggregate(bassn$freq,
                   by = list(Site=bassn$Site,
                             Day=bassn$Day,
                             DOY=bassn$DOY,
                             CHLA=bassn$CHLA,
                             SST=bassn$SST,
                             Depth=bassn$Depth,
                             d50.b=bassn$d50.b,
                             time=bassn$time,
                             Season=bassn$Season,
                             Seasonb=bassn$Seasonb,
                             Year=bassn$Year),FUN= "sum")

colnames(bassd)[12] <- "freq"
bassd$Year<-as.factor(bassd$Year)

bassda<-left_join(bassd, Lags, by=c("Day","Site"))
bassd2<-bassda[!duplicated(bassda),]
#
bassd2 <- bassd2[with(bassd2,order(time,Site)),]

# Converted to list, as for sturgeon, above   -MOB 20200729
bassd_list <- as.list(bassd2[, 1:12])
bassd_list[['Qb2']] <- as.matrix(bassd2[c(13:42)])
bassd_list[['Lb2']] <- matrix(0:(ncol(bassd_list$Qb2) - 1),
                               nrow(bassd_list$Qb2),
                               ncol(bassd_list$Qb2),
                               byrow = TRUE)


ts.start <- proc.time()
# Log transform
bassd_list$CHLA <- log(bassd_list$CHLA)

# Site to factor, as for sturgeon, above    -MOB 20200729
bassd_list$Site <- as.factor(bassd_list$Site)

BM1 <- gam(freq~s(Qb2,Lb2,bs="cb",k=10)+s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(CHLA,bs='tp')+
             s(DOY,bs="cs",k=6)+
             s(Depth,bs='tp')+
             offset(log(d50.b)),
           data=bassd_list,family=ziP,method="REML")

BM2 <- gam(freq~s(Qb2,Lb2,bs="cb",k=10)+s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(CHLA,bs='tp')+
             s(DOY,bs="cs",k=6)+
             offset(log(d50.b)),
           data=bassd_list,family=ziP,method="REML")

BM3 <- gam(freq~s(Qb2,Lb2,bs="cb",k=10)+s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(CHLA,bs='tp')+
             s(Depth,bs='tp')+
             offset(log(d50.b)),
           data=bassd_list,family=ziP,method="REML")

BM4 <- gam(freq~s(Qb2,Lb2,bs="cb",k=10)+s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(DOY,bs="cs",k=6)+
             s(Depth,bs='tp')+
             offset(log(d50.b)),
           data=bassd_list,family=ziP,method="REML")

BM5 <- gam(freq~s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(CHLA,bs='tp')+
             s(DOY,bs="cs",k=6)+
             s(Depth,bs='tp')+
             offset(log(d50.b)),
           data=bassd_list,family=ziP,method="REML")

BM6 <- gam(freq~s(Qb2,Lb2,bs="cb",k=10)+s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(CHLA,bs='tp')+
             t2(DOY,Depth,bs=c("cp", "tp"),k=6)+
             offset(log(d50.b)),
           data=bassd_list,family=ziP,method="REML")

BM7 <- gam(freq~s(Qb2,Lb2,bs="cb",k=10)+s(Site, bs = 're')+
             s(Year, bs = 're')+
             t2(DOY,Depth,bs=c("cp", "tp"),k=6)+
             offset(log(d50.b)),
           data=bassd_list,family=ziP,method="REML")

BM8 <- gam(freq~s(Site, bs = 're')+
             s(Year, bs = 're')+
             s(CHLA,bs='tp')+
             t2(DOY,Depth,bs=c("cp", "tp"),k=6)+
             offset(log(d50.b)),
           data=bassd_list,family=ziP,method="REML")

BM9 <- gam(freq~s(Site, bs = 're')+
             s(Year, bs = 're')+
             t2(DOY,Depth,bs=c("cp", "tp"),k=6)+
             offset(log(d50.b)),
           data=bassd_list,family=ziP,method="REML")

BM10 <- gam(freq~s(Qb2,Lb2,bs="cb",k=10)+s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(Depth,bs='tp')+
              t2(DOY,CHLA,bs=c("cp", "tp"),k=6)+
              offset(log(d50.b)),
            data=bassd_list,family=ziP,method="REML")

BM11 <- gam(freq~s(Qb2,Lb2,bs="cb",k=10)+s(Site, bs = 're')+
              s(Year, bs = 're')+
              t2(DOY,CHLA,bs=c("cp", "tp"),k=6)+
              offset(log(d50.b)),
            data=bassd_list,family=ziP,method="REML")

BM12 <- gam(freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(Depth,bs='tp')+
              t2(DOY,CHLA,bs=c("cp", "tp"),k=6)+
              offset(log(d50.b)),
            data=bassd_list,family=ziP,method="REML")

BM13 <- gam(freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              t2(DOY,CHLA,bs=c("cp", "tp"),k=6)+
              offset(log(d50.b)),
            data=bassd_list,family=ziP,method="REML")

BM14 <- gam(freq~s(Qb2,Lb2,bs="cb",k=10)+s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs="cs",k=6)+
              t2(Depth,CHLA,bs=c("tp", "tp"),k=6)+
              offset(log(d50.b)),
            data=bassd_list,family=ziP,method="REML")

BM15 <- gam(freq~s(Qb2,Lb2,bs="cb",k=10)+s(Site, bs = 're')+
              s(Year, bs = 're')+
              t2(Depth,CHLA,bs=c("tp", "tp"),k=6)+
              offset(log(d50.b)),
            data=bassd_list,family=ziP,method="REML")

BM16 <- gam(freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs="cs",k=6)+
              t2(Depth,CHLA,bs=c("tp", "tp"),k=6)+
              offset(log(d50.b)),
            data=bassd_list,family=ziP,method="REML")

BM17 <- gam(freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              t2(Depth,CHLA,bs=c("tp", "tp"),k=6)+
              offset(log(d50.b)),
            data=bassd_list,family=ziP,method="REML")

BM18 <- gam(freq~s(Qb2,Lb2,bs="cb",k=10)+s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs="cs",k=6)+
              s(Depth,bs='tp')+
              offset(log(d50.b)),
            data=bassd_list,family=ziP,method="REML")

BM19 <- gam(freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs="cs",k=6)+
              s(Depth,bs='tp')+
              offset(log(d50.b)),
            data=bassd_list,family=ziP,method="REML")

BM20 <- gam(freq~s(Qb2,Lb2,bs="cb",k=10)+s(Site, bs = 're')+
              s(Year, bs = 're')+
              t2(DOY,Depth,bs=c("cp", "tp"),k=6)+
              offset(log(d50.b)),
            data=bassd_list,family=ziP,method="REML")

BM21 <- gam(freq~s(Qb2,Lb2,bs="cb",k=10)+s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(CHLA,bs='tp')+
              offset(log(d50.b)),
            data=bassd_list,family=ziP,method="REML")

BM22 <- gam(freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(CHLA,bs='tp')+
              offset(log(d50.b)),
            data=bassd_list,family=ziP,method="REML")

BM23 <- gam(freq~s(Qb2,Lb2,bs="cb",k=10)+s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs="cs",k=6)+
              offset(log(d50.b)),
            data=bassd_list,family=ziP,method="REML")

BM24 <- gam(freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(DOY,bs="cs",k=6)+
              offset(log(d50.b)),
            data=bassd_list,family=ziP,method="REML")

BM25 <- gam(freq~s(Qb2,Lb2,bs="cb",k=10)+s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(Depth,bs='tp')+
              offset(log(d50.b)),
            data=bassd_list,family=ziP,method="REML")

BM26 <- gam(freq~s(Site, bs = 're')+
              s(Year, bs = 're')+
              s(Depth,bs='tp')+
              offset(log(d50.b)),
            data=bassd_list,family=ziP,method="REML")

AIC(BM1,BM2,BM3,BM4,BM5,BM6,BM7,BM8,BM9,BM10,
    BM11,BM12,BM13,BM14,BM15,BM16,BM17,BM18,BM19,BM20,
    BM21,BM22,BM23,BM24,BM25,BM26)

summary(BM6)
summary(BM7)
summary(BM11)
summary(BM10)
summary(BM14)

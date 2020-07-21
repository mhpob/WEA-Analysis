library(dlnm) ; library(mgcv); library(dplyr)

sturgo<-read.csv("dlnm/sturgoSSTlags.csv",header=T)

# Drop rownames.
sturgo <- sturgo[, -1]

sturgo[,2]<-as.POSIXct(sturgo[,2], format="%Y-%m-%d")

sturgd<- aggregate(sturgo$freq,
                   by = list(Site=sturgo$Site,
                             Day=sturgo$Day,
                             DOY=sturgo$DOY,
                             CHLA=sturgo$CHLA,
                             SST=sturgo$SST,
                             Depth=sturgo$Depth,
                             d50.s=sturgo$d50.s,
                             time=sturgo$time,
                             Season=sturgo$Season,
                             Seasonb=sturgo$Seasonb,
                             Year=sturgo$Year),FUN= "sum")

colnames(sturgd)[12] <- "freq"

Lags <- sturgo[c(1,2,28:57)]

sturgda<-left_join(sturgd, Lags, by=c("Day","Site"))
sturgd2<-sturgda[!duplicated(sturgda),]

sturgd2 <- sturgd2[with(sturgd2,order(time,Site)),]

# Log transform
sturgd2$CHLA<-(log(sturgd2$CHLA))
sturgd2$Site <- as.factor(sturgd2$Site)
sturgd2$Year <- as.factor(sturgd2$Year)



####################
## New RMSE, with suggestions from Dong
set.seed(12345)
sturgdata<-sturgd2[sample(nrow(sturgd2)),]
# Create 5 equally size folds
folds <- cut(seq(1,nrow(sturgdata)),breaks=5,labels=FALSE)
# Perform 5 fold cross validation
storage <- vector("list",5) ## DL added

pb <- txtProgressBar(min = 0, max = length(unique(folds)), style = 3)
for(i in 1:5){
  # Segment data by fold using the which() function
  StestIndexes <- which(folds==i,arr.ind=TRUE)
  StestData <- sturgdata[StestIndexes, ]
  StrainData <- sturgdata[-StestIndexes, ]
  ## THE DLNM WILL BE TRAINED HERW USING BTRAIN DATA
  ## AND THE TRAINED MODEL WITH USED TO PREDICT ON BTEST DATA
  ## RECORD THE PREDICTION AND COMPARE WITH THE FREUQNCYS
  ## THIS FOLLOWING IS AN EXAMPLE USING YOUR CODE
  ## BUT REPLACE WITH DLNM

  create_input_data <- function(model_data){
    # Lag matrix
    Q <- as.matrix(model_data[c(13:42)])
    L <- matrix(0:(ncol(Q) - 1), nrow(Q), ncol(Q), byrow = TRUE)

    # Convert predictors to elements in a list
    model_data <- as.list(model_data)

    # Add Q and L to that list
    model_data[['Q']] <- Q
    model_data[['L']] <- L

    model_data
  }

  # Train the model.
  train <- create_input_data(StrainData)
  Sfinal <- gam(freq ~ s(Q, L, bs="cb", k = 10) +
                  s(Site, bs = 're') +
                  s(Year, bs = 're') +
                  s(CHLA, bs='tp') +
                  t2(DOY, Depth, bs=c("cp", "tp"), k = 6) +
                  offset(log(d50.s)),
                data = train, family = ziP, method = "REML")


  # Test the model
  test <- create_input_data(StestData)

  StestData$T1 = predict(Sfinal, test, type = "response")
  storage[[i]] <- StestData ## DL added

  setTxtProgressBar(pb, i)
}
close(pb)


rmse <- lapply(storage, function(.) sqrt(mean((.$freq - .$T1) ^ 2, na.rm = T)))
mean(unlist(rmse))
sd(unlist(rmse))


# CVfinal <- do.call(rbind,storage) ## DL added
CVRMSE=with(CVfinal,sqrt(mean((freq - T1)^2))) ## DL added

StestData$T1= predict(Sfinal,StestData,type = "response")

SRMSE=sqrt(mean((StestData$freq - StestData$T1)^2))
# Calculated values= 0.8055, 0.8517, 0.9373, 0.8346, 0.7492

mean(c(0.8055, 0.8517, 0.9373, 0.8346, 0.7492))
sd(c(0.8055, 0.8517, 0.9373, 0.8346, 0.7492))

BAE=mean((StestData$freq - StestData$T1))
# Calculated values= 0.0069, -0.0232, 0.0366, -0.0017, -0.0036

mean(c(0.0069, -0.0232, 0.0366, -0.0017, -0.0036))
sd(c(0.0069, -0.0232, 0.0366, -0.0017, -0.0036))

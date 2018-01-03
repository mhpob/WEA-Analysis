library(ggplot2); library(reshape2)

times <- read.csv(
  'p:/obrien/biotelemetry/md wea habitat/data/rec_deploydates.csv')
times$d.st <- lubridate::ymd(times$d.st)
times$d.end <- lubridate::ymd(times$d.end)

tend <- read.csv(
  'p:/obrien/biotelemetry/md wea habitat/data/rec_tenddates.csv')
tend <- melt(tend, id = c('Station'), value.name = 'd.t', na.rm = T)
tend$d.t <- lubridate::ymd(tend$d.t)
tend <- merge(tend, times)
tend <- unique(tend)


j <- list(NULL)
for(i in 1:dim(times)[1]){
  k <- data.frame(Station = times[i, 'Station'],
                  order = times[i, 'ord'],
                  date = seq(times[i, 'd.st'], times[i, 'd.end'], by = 'day'))
  j[[i]] <- k
}
times <- do.call(rbind.data.frame, j)


labs <- c('IN1', 'IN2', 'IS1', 'IS2', 'IS2_250', 'IS2_800', 'AN1', 'AN2', 'AN3',
          'AN3_250', 'AN3_800', 'AN4', 'AM1', 'AM2', 'AM3', 'AM4', 'AS1', 'AS2',
          'AS3', 'AS4', 'ON1', 'ON2', 'OS1', 'OS2')



ggplot() + geom_raster(data = times,
                       aes(x = date, y = as.factor(order)),
                       fill = 'darkgreen') +
  geom_segment(data = tend,
               aes(x = d.t, y = ord - 0.45,
                   xend = d.t, yend = ord + 0.45),
               size = 2) +
  scale_y_discrete(labels = rev(labs)) +
  theme_bw() +
  labs(x = 'Date', y = 'Site')

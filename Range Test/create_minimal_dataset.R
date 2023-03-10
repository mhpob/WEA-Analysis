library(data.table)

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


setnames(data,
         c('array', 'station', 'average_noise', 'average_temperature', 'freq'),
         c('site', 'receiver', 'noise', 'bwt', 'detectability'))


data <- data[, .(date, site, receiver, distance, success, fail, detectability, noise, bwt, sst, dt)]

data[, ':='(site = fifelse(site == 'MDWEA', 'midshelf', 'nearshore'),
            receiver = fcase(grepl('^...$', receiver), 'A',
                             grepl('_2', receiver), 'B',
                             grepl('_8', receiver), 'C'))]

fwrite(data, 'range test/manuscript/revisions/additional_file_1.csv')

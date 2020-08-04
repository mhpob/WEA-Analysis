# Packages ----
library(ggplot); library(data.table)


# Import data ----
data <- setDT(
  readRDS('data and imports/rangetest_logit_binary_pt0.RDS')
)
names(data) <- gsub(' ', '_', tolower(names(data)))
data <- unique(data, by = c('date', 'array'))



# Visual inspection of cold pool ramp-up ----
## Meant to be quick and dirty.
cp_start <- data[date %between% c('2018-04-10', '2018-05-01')]

ggplot() +
  geom_line(data = cp_start, aes(x = date, y = dt, color = array)) +
  scale_x_date(date_breaks = 'day') +
  scale_y_continuous(breaks = -2:20) +
  theme(axis.text.x = element_text(angle = 45))

## Ramp-up starts on April 15


# Visual inspection and stats for cold pool formation ----
## Meant to be quick and dirty.
cp_form <- data[date %between% c('2018-04-28', '2018-05-07')]

ggplot() +
  geom_line(data = cp_start, aes(x = date, y = dt, color = array)) +
  scale_x_date(date_breaks = 'day') +
  scale_y_continuous(breaks = -2:20) +
  theme(axis.text.x = element_text(angle = 45))

## Formation occurs from May 2 - May 5

## Deg. per day
cp_form[date %between% c('2018-05-02', '2018-05-05'),
        (max(dt) - min(dt)) / (.N - 1), by = 'array']



# Visual inspection and stats for cold pool collapse ----
cp_collapse <- data[date %between% c('2018-09-01', '2018-09-14')]

ggplot() +
  geom_line(data = cp_collapse, aes(x = date, y = dt, color = array)) +
  scale_x_date(date_breaks = 'day') +
  scale_y_continuous(breaks = -2:20) +
  theme(axis.text.x = element_text(angle = 45))

## Collapse occurs from Sept 7 - Sept 10

## Deg. per day
### Inner collapses Sep 7 - 9
cp_collapse[date %between% c('2018-09-07', '2018-09-09') & array == 'Inner',
            (min(dt) - max(dt)) / (.N - 1)]

### MDWEA collapses Sep 8 - 10
cp_collapse[date %between% c('2018-09-08', '2018-09-10') & array == 'MD WEA',
            (min(dt) - max(dt)) / (.N - 1)]



# Range of CP stratification ----
data[date %between% c('2018-05-05', '2018-09-07'),
     range(dt),
     by = 'array']



# Visual inspection of storm-driven decreases ----
## Storm in May
cp_storm <- data[date %between% c('2018-05-01', '2018-05-31')]

ggplot() +
  geom_line(data = cp_storm, aes(x = date, y = dt, color = array)) +
  scale_x_date(date_breaks = 'day') +
  scale_y_continuous(breaks = -2:20) +
  theme(axis.text.x = element_text(angle = 45))

### Partial destratification due to storm from May 16-19
### Deg. per day
cp_storm[date %between% c('2018-05-16', '2018-05-19'),
         (min(dt) - max(dt)) / (.N - 1),
         by = 'array']


## Storm in June
cp_storm <- data[date %between% c('2018-06-01', '2018-06-30')]

ggplot() +
  geom_line(data = cp_storm, aes(x = date, y = dt, color = array)) +
  scale_x_date(date_breaks = 'day') +
  scale_y_continuous(breaks = -2:20) +
  theme(axis.text.x = element_text(angle = 45))

### Partial destratification due to storm from June 2-4
### Deg. per day
cp_storm[date %between% c('2018-06-02', '2018-06-04'),
         (min(dt) - max(dt)) / (.N - 1),
         by = 'array']


## Storm in July
cp_storm <- data[date %between% c('2018-07-01', '2018-07-31')]

ggplot() +
  geom_line(data = cp_storm, aes(x = date, y = dt, color = array)) +
  scale_x_date(date_breaks = 'day') +
  scale_y_continuous(breaks = -2:20) +
  theme(axis.text.x = element_text(angle = 45))

### Partial destratification due to storm from July 5-9
### Deg. per day
cp_storm[date %between% c('2018-07-05', '2018-07-09'),
         (min(dt) - max(dt)) / (.N - 1),
         by = 'array']


## Storm in August
cp_storm <- data[date %between% c('2018-08-01', '2018-08-31')]

ggplot() +
  geom_line(data = cp_storm, aes(x = date, y = dt, color = array)) +
  scale_x_date(date_breaks = 'day') +
  scale_y_continuous(breaks = -2:20) +
  theme(axis.text.x = element_text(angle = 45))

### Partial destratification due to storm from August 19-24
### Deg. per day
cp_storm[date %between% c('2018-08-19', '2018-08-24'),
         (min(dt) - max(dt)) / (.N - 1),
         by = 'array']



# Noise stats and tests ----
data <- data[, ':='(cp = ifelse(date %between% c('2018-05-02', '2018-09-09'),
                                'present',
                                'absent'),
                    mask = ifelse(average_noise >= 300, 'T', 'F'))]


## ANOVA b/w noise in CP+ and CP- periods
summary(lm(average_noise ~ cp, data = data))


## Chi-square test b/w proportion of masked (>300 mV) hours in CP+ and CP- periods
chisq.test(xtabs(~ cp + mask, data = data), correct = F)

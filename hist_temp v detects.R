library(ggplot2); library(dplyr)
# '%>%', filter(), mutate(), group_by(), and summarize() are all in the dplyr package

# Read events (temperatures) ----
WT <- readRDS("P:/Wiernicki/Striped Bass/BOEM Wind Energy Area/rec_events.rds")

WT <- WT %>%
  # Filter data to fewer rows in order to save computing time
  filter(Description == 'Average temperature') %>%
  mutate(
    # Define the dates as POSIX
    Date.Time = lubridate::ymd_hms(Date.Time),
    # Pull out date information
    date = lubridate::date(Date.Time),
    # Define the data as numeric
    Data = as.numeric(Data),
    season = ifelse(lubridate::month(Date.Time) %in% 4:9, 'SprSum',
                   'AutWin')) %>%
  # Group by date to prepare for aggregation
  group_by(date, season) %>%
  # Aggregate the temperature data by date
  summarize(BWT = mean(Data))

WT$bins <- cut(WT$BWT,
               breaks = seq(floor(min(WT$BWT)), ceiling(max(WT$BWT)), 1))
# Save levels for later
bin_levels <- levels(WT$bins)


# Load detection data ----
data <- TelemetryR::vemsort('p:/obrien/biotelemetry/detections/offshore md/fish migration')
load(file = 'p:/obrien/randomr/ACTactive.rda')

data <- data %>%
  # Join ACT database to detection data
  left_join(ACTactive, by = c('transmitter' = 'Tag.ID.Code.Standard')) %>%
  mutate(season = ifelse(lubridate::month(date.local) %in% 4:9, 'SprSum',
                         'AutWin'),
    # Make sure SB, AS, and WS are named the same thing. Use GREP pattern matching.
    Common.Name = ifelse(grepl('striped', Common.Name, ignore.case = T),
                         'Striped bass',
                         ifelse(grepl('c stur|^stur', Common.Name, ignore.case = T),
                                'Atlantic sturgeon',
                         ifelse(grepl('white shark',  Common.Name, ignore.case = T),
                                'White shark',
                                Common.Name))),
    # Pull out date information
    date = lubridate::date(date.utc)) %>%
  # Keep only target species
  filter(Common.Name %in% c('Striped bass', 'Atlantic sturgeon', 'White shark')) %>%
  # Group by date/season to prepare for aggregation
  group_by(season, date, Common.Name) %>%
  # Aggregate detections per species by date
  summarize(n_fish = n_distinct(transmitter)) %>%
  # Join corresponding temperature data
  full_join(WT)

# Aggregate both temperature (count) and detections (sum) by WT bin ----
# Treat temperature first by removing redundant dates...
hold_BWT <- data_daily[!duplicated(data_daily$date),
                       names(data_daily) %in% c('bins', 'BWT')]
# ...and then aggregating
hold_BWT <- aggregate(BWT ~ bins, FUN = length, data = hold_BWT)

# Treat detection data second by summing detections. Reshape2 can handle this
library(reshape2)
hold_detects <- dcast(data = data_daily, bins ~ Common.Name,
                      fun.aggregate = sum,
                      value.var = 'n_fish')

data_daily <- merge(hold_detects[,!names(hold_detects) == 'NA'], hold_BWT)

# Scale by total number
data_daily[, 2:5] <- sweep(data_daily[, 2:5], 2, colSums(data_daily[, 2:5]), `/`)

# Re-order the factor levels
data_daily$bins <- factor(data_daily$bins, levels = bin_levels)

# Melt the data frame for use in ggplot
data_daily <- melt(data_daily, id = 'bins', variable.name = 'type')

# All-data plotting ----
ggplot() + geom_col(data = filter(data_daily, grepl('BWT|Atl', type)),
                    aes(x = bins, y = value, fill = type),
                    position = 'dodge') +
  labs(title = 'Atlantic sturgeon', fill = NULL,
       x = 'Temperature Bin', y = 'Relative frequency') +
  scale_fill_manual(labels = c('Detections', 'Bottom\nTemperature'),
                    values = c('red', 'blue')) +
  ylim(0, 0.24) +
  theme_bw() +
  theme(legend.position = c(0.9, 0.9))

ggplot() + geom_col(data = filter(data_daily, grepl('BWT|Str', type)),
                    aes(x = bins, y = value, fill = type),
                    position = 'dodge') +
  labs(title = 'Striped bass', fill = NULL,
       x = 'Temperature Bin', y = 'Relative frequency') +
  scale_fill_manual(labels = c('Detections', 'Bottom\nTemperature'),
                    values = c('red', 'blue')) +
  ylim(0, 0.24) +
  theme_bw() +
  theme(legend.position = c(0.9, 0.9))

ggplot() + geom_col(data = filter(data_daily, grepl('BWT|Wh', type)),
                    aes(x = bins, y = value, fill = type),
                    position = 'dodge') +
  labs(title = 'White shark', fill = NULL,
       x = 'Temperature Bin', y = 'Relative frequency') +
  scale_fill_manual(labels = c('Detections', 'Bottom\nTemperature'),
                    values = c('red', 'blue')) +
  ylim(0, 0.24) +
  theme_bw() +
  theme(legend.position = c(0.9, 0.9))





library(dplyr)

library(ggplot2)
ggplot() +
  geom_line(data = data, aes(x = date, y = d50_adj), lwd = 1) +
  facet_wrap(~array, nrow = 2)+
  labs(x=NULL, y = 'Adjusted Modeled D50 (m)')+
  theme_bw()

ggplot() +
  geom_line(data = j, aes(x = date, y = value), lwd = 1) +
  facet_wrap(type ~ array, scales = 'free') +
  labs(x=NULL, y = 'D50 (m)')+
  theme_bw()

j <- data %>%
  select(Modeled = d50, `Outliers Removed` = d50_adj, array, date) %>%
  tidyr::gather(key = 'type', value ='value',
                -array, -date, Modeled, `Outliers Removed`) %>%
  mutate(key = paste(array, type, sep = ', '))

ggplot() +
  geom_line(data = j, aes(x = date, y = value), lwd = 1) +
  facet_wrap(~key, scales = 'free', dir = 'v') +
  labs(x=NULL, y = 'D50 (m)')+
  theme_bw()

data %>%
  group_by(array) %>%
  summarize(min(d95_adj),
            max(d95_adj))

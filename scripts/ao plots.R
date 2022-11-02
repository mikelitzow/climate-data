# plot AO time series

library(tidyverse)

dat <- read.csv("./downloads/AO.csv")

plot <- dat %>%
  group_by(year) %>%
  summarise(ao = mean(ao))

ggplot(plot, aes(year, ao)) +
  geom_line() +
  geom_point()

dat <- dat %>%
  mutate(dec.yr = year + (month - 0.5) / 12)

ggplot(dat, aes(dec.yr, ao)) +
  geom_line() +
  geom_point()

arrange(dat, desc(ao))

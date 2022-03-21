# plot annual EBS / GOA SST for outreach talk

library(tidyverse)

# load summarized ersst data

dat <- read.csv("./output data/regional_north_pacific_ersst_time_series.csv")

head(dat)

dat <- dat %>%
  filter(region %in% c("Eastern_Bering_Sea", "Gulf_of_Alaska"),
         year >= 1950) %>%
  mutate(temp = (annual.unsmoothed*9/5)+32) %>%
  select(region, year, temp)

ggplot(dat, aes(year, temp)) +
  geom_line() +
  geom_point() +
  facet_wrap(~region, ncol = 1, scales = "free_y") +
  ylab("Sea surface temperature (°F)") +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  scale_x_continuous(breaks = seq(1950,2020,10)) 

ggsave("./figs/fahrenheit_EBS_GOA_sst.png", width = 4, height = 5, units = 'in')

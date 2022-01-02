library(tidyverse)

theme_set(theme_bw())

# colorblind palette
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# download MASIE daily ice extent data

download.file("https://masie_web.apps.nsidc.org/pub/DATASETS/NOAA/G02186/masie_4km_allyears_extent_sqkm.csv", "./downloads/MAISIE.csv")

dat <- read.csv("./downloads/MAISIE.csv", skip = 1)

head(dat)

# Bering Sea only

keep <- grep("Bering", names(dat))

dat <- dat %>%
  select(1, keep) 

names(dat)[2] <- "Bering"

dat$year <- floor(dat$yyyyddd / 1000)
dat$day <- dat$yyyyddd - dat$year*1000

# now establish winter year corresponding to January

dat$winter <- if_else(dat$day < 244, dat$year, dat$year+1)

# now get day of winter (Sept. 1 = winter day 1)

winters <- unique(dat$winter)

win.temp <- data.frame()

for(i in 1:length(winters)) {
  # i <- 2
  
  temp <- dat %>%
    filter(winter == winters[i])
  
  temp$winter.day <- 1:nrow(temp)
  
  # bump up winter 2006 to account for missing data from 2005
  if(i == 1){temp$winter.day = temp$winter.day + 121}
  
  win.temp <- rbind(win.temp, temp)
    
}


# put back into main data

win.temp <- win.temp %>%
  select(yyyyddd, winter.day)

dat <- left_join(dat, win.temp)

dat$winter <- as.factor(dat$winter)

# get 2006:2013 mean
climatology <- dat %>%
  filter (winter %in% 2006:2013) %>%
  group_by(winter.day) %>%
  summarise(mean = mean(Bering))

clim.dat <- data.frame(Bering = climatology$mean,
                       winter = "mean", 
                       winter.day = climatology$winter.day)

dat <- dat %>%
  select(Bering, winter, winter.day)


dat <- rbind(dat, clim.dat)

dat$period <- case_when(
  
  dat$winter %in% 2006:2021 ~ "2006-2021",
  dat$winter == "mean" ~ "2006-2013 mean",
  dat$winter == 2022 ~ "2022"
)
  
# set line width for plot
dat$line.width <- if_else(dat$winter %in% 2006:2021, 0.3,
                          if_else(dat$winter == "mean", 0.8, 1))

ggplot(dat, aes(winter.day, Bering/1e4, group = winter, color = period)) +
  geom_line(lwd = dat$line.width) +
  scale_color_manual(values = c("black", "grey", cb[6])) +
  labs(x = "Winter day (September 1 = day 1)",
       y = expression(Ice~extent~(10^4~km^2)),
       title = "Bering Sea ice extent",
       subtitle = "MAISIE data, NSIDC") +
  theme(legend.title = element_blank(),
        legend.position = c(0.2, 0.8))

ggsave("./figs/maisie.png", width = 5, height = 4, units = 'in')

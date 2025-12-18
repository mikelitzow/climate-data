library(tidyverse)
theme_set(theme_bw())
library(chron)
library(fields)
library(maps)
library(mapdata)
library(oce)
library(ncdf4)

# colorblind palette
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# identify most recent month/date for desired data
last_month <- 11
last_year <- 2025

# define URL
URL <- 
  "https://coastwatch.pfeg.noaa.gov/erddap/griddap/nceiErsstv5.nc?sst%5B(1920-01-01):1:(2025-11-15T00:00:00Z)%5D%5B(0.0):1:(0.0)%5D%5B(20):1:(66)%5D%5B(110):1:(250)%5D"

# URL <-
  # paste("https://apdrc.soest.hawaii.edu/erddap/griddap/hawaii_soest_31a3_72d5_401e.nc?sst[(1920-01-15):1:(",
        # last_year, "-", last_month, "-15)][(20):1:(66)][(110):1:(250)]", sep = "")

# download ERSSTv5 data
# download.file(URL, "./downloads/ERSST.nc")

# load data
nc <- nc_open("./downloads/nceiErsstv5_17a9_06de_e3f8.nc")
# nc <- nc_open("./downloads/hawaii_soest_31a3_72d5_401e_4e48_f10a_e9b0.nc")

# nc <- nc_open("./downloads/hawaii_soest_31a3_72d5_401e_4a7c_bdaf_7af6.nc")

# process

ncvar_get(nc, "time")   # seconds since 1-1-1970
raw <- ncvar_get(nc, "time")
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))
m <- months(d)
yr <- years(d)

x <- ncvar_get(nc, "longitude")
y <- ncvar_get(nc, "latitude")

SST <- ncvar_get(nc, "sst", verbose = F)

SST <- aperm(SST, 3:1)

SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))
lon <- rep(x, each = length(y))
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# plot to check
temp.mean <- colMeans(SST, na.rm=T)
z <- t(matrix(temp.mean,length(y)))
image.plot(x,y,z, col=oceColorsPalette(64), xlab = "°E Longitude", ylab = "°N Latitude",
           legend.lab = "°C")
contour(x, y, z, add=T)
map('world2Hires',c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China', 'Mongolia'), fill=T,add=T, lwd=1, col="lightyellow3")

box()

# identify polygons for each area


# EBS/NBS
# 54-66 deg. N, 188-202 deg. E
ebs.x <- c(183, 183, 203, 203, 191) 
ebs.y <- c(53, 65, 65, 57.5, 53)

# GOA polygon
goa.x <- c(201, 201, 205, 208, 225, 231, 201)
goa.y <- c(55, 56.5, 59, 61, 61, 55, 55)

# define cells within each polygon and plot to check

#first, ebs
ebs.sst <- as.data.frame(SST)

xp <- cbind(ebs.x, ebs.y)
loc=cbind(lon, lat)
check <- in.poly(loc, xp=xp)

ebs.sst[,!check] <- NA

# plot to check
temp.mean <- colMeans(ebs.sst, na.rm=T)
z <- t(matrix(temp.mean,length(y)))
image.plot(x,y,z, col=oceColorsPalette(64), xlab = "", ylab = "", xlim = c(170, 240), ylim = c(40, 66))
contour(x, y, z, add=T)
map('world2Hires',c('Canada', 'usa', 'USSR'), fill=T,add=T, lwd=1, col="lightyellow3")

# now, goa

goa.sst <- as.data.frame(SST)

xp <- cbind(goa.x, goa.y)
loc=cbind(lon, lat)
check <- in.poly(loc, xp=xp)

goa.sst[,!check] <- NA

# plot to check
temp.mean <- colMeans(goa.sst, na.rm=T)
z <- t(matrix(temp.mean,length(y)))
image.plot(x,y,z, col=oceColorsPalette(64), xlab = "", ylab = "", xlim = c(170, 240), ylim = c(40, 66))
contour(x, y, z, add=T)
map('world2Hires',c('Canada', 'usa', 'USSR'), fill=T,add=T, lwd=1, col="lightyellow3")

# create vector of latitudes to weight mean sst by cell area

cell.weight <- sqrt(cos(lat*pi/180))

# plot to check
hist(cell.weight, breaks = 50)
unique(cell.weight) # looks right

# create a weighted mean function
weighted.cell.mean <- function(x) weighted.mean(x, w = cell.weight, na.rm = T)

# create a function to compute monthly mean for pre-1950 data
climatology.mean <- function(x) tapply(x, m[yr < 1950], mean) 

# create data frame to catch anomaly time series
ersst_monthly_anomalies <- data.frame()

## start with N Pac ------------------

# first, calculate monthly anomalies for pre-1950 data
climatology <- SST[yr < 1950,]

mu <- apply(climatology, 2, climatology.mean)	# compute monthly means for each time series (cell)
mu_all <- mu[rep(1:12, length(m)/12),]  # replicate means matrix for each year at each location

# add incomplete year if needed
if(last_month < 12){
  mu_all <- rbind(mu_all, mu[1:last_month,])}

n_pac_anomaly <- SST - mu_all

n_pac_monthly_anomaly <- apply(n_pac_anomaly, 1, weighted.cell.mean)

ersst_monthly_anomalies <- rbind(ersst_monthly_anomalies,
                                 data_frame(region = "North Pacific",
                                            year = yr,
                                            month = m,
                                            dec.year = as.numeric(as.character(yr)) +
                                              (as.numeric(m)-0.5) / 12,
                                            anomaly = n_pac_monthly_anomaly))

## add EBS ------------------

# first, calculate monthly anomalies for pre-1950 data
ebs.climatology <- ebs.sst[yr < 1950,]

ebs.mu <- apply(ebs.climatology, 2, climatology.mean)	# compute monthly means for each time series (cell)
ebs.mu_all <- ebs.mu[rep(1:12, length(m)/12),]  # replicate means matrix for each year at each location

# add incomplete year if needed
if(last_month < 12){
  ebs.mu_all <- rbind(ebs.mu_all, ebs.mu[1:last_month,])}

ebs_anomaly <- ebs.sst - ebs.mu_all

ebs_monthly_anomaly <- apply(ebs_anomaly, 1, weighted.cell.mean)

ersst_monthly_anomalies <- rbind(ersst_monthly_anomalies,
                                 data_frame(region = "Bering Sea",
                                            year = yr,
                                            month = m,
                                            dec.year = as.numeric(as.character(yr)) +
                                              (as.numeric(m)-0.5) / 12,
                                            anomaly = ebs_monthly_anomaly))

## add GOA ------------------

# first, calculate monthly anomalies for pre-1950 data
goa.climatology <- goa.sst[yr < 1950,]

goa.mu <- apply(goa.climatology, 2, climatology.mean)	# compute monthly means for each time series (cell)
goa.mu_all <- goa.mu[rep(1:12, length(m)/12),]  # replicate means matrix for each year at each location

# add incomplete year if needed
if(last_month < 12){
  goa.mu_all <- rbind(goa.mu_all, goa.mu[1:last_month,])}

goa_anomaly <- goa.sst - goa.mu_all

goa_monthly_anomaly <- apply(goa_anomaly, 1, weighted.cell.mean)

ersst_monthly_anomalies <- rbind(ersst_monthly_anomalies,
                                 data_frame(region = "Gulf of Alaska",
                                            year = yr,
                                            month = m,
                                            dec.year = as.numeric(as.character(yr)) +
                                              (as.numeric(m)-0.5) / 12,
                                            anomaly = goa_monthly_anomaly))

# and plot
ggplot(ersst_monthly_anomalies, aes(dec.year, anomaly)) +
  geom_line(color = "dark grey") +
  geom_smooth(method = "gam", se = F) +
  facet_wrap(~region, scales = "free_y", ncol = 1) +
  geom_hline(yintercept = 0) + 
  scale_x_continuous(breaks = seq(1920, 2020, 10)) +
  theme(axis.title.x = element_blank()) + 
  ylab("Anomaly (°C, 1920-1949 climatology)")


# monthly ar1 values for 1950-present
ar1 <- ersst_monthly_anomalies %>%
  filter(year >= 1950) %>%
  mutate(era = as.factor(
           case_when(year %in% 1950:1988 ~ "1950-1988",
                     year %in% 1989:1999 ~ "1989-1999",
                     year %in% 2000:2025 ~ "2000-2025"))) %>%
  group_by(region, era) %>%
  reframe(ar1 = acf(anomaly, lag.max = 60, plot = F)$acf) %>%
  mutate(lag = rep(0:60, 9)) 


# and plot
ggplot(ar1, aes(lag, ar1, color = era)) +
  geom_line() +
  geom_point() +
  facet_wrap(~region, ncol = 1) +
  geom_hline(yintercept = 0, lty = 2) +
  xlab("Lag (months)") +
  scale_x_continuous(breaks = c(0, 12, 24, 36, 48, 60))

# pretty wild pattern! 
# detrend and see if the pattern remains

detr.ar1 <- ersst_monthly_anomalies %>%
  filter(year >= 1950) %>%
  group_by(region) %>%
  reframe(time_step = 1:n(),
    detr.anom = resid(lm(anomaly ~ dec.year)))

# plot to check I'm coding correctly
ggplot(detr.ar1, aes(time_step, detr.anom)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0) +
  facet_wrap(~region, scales = "free")

detr.ar1 <- detr.ar1 %>%
  mutate(year = floor(1950+(time_step-0.5)/12),
    era = as.factor(
    case_when(year %in% 1950:1988 ~ "1950-1988",
              year %in% 1989:1999 ~ "1989-1999",
              year %in% 2000:2025 ~ "2000-2025"))) %>%
  group_by(region, era) %>%
  reframe(ar1 = acf(detr.anom, lag.max = 60, plot = F)$acf) %>%
  mutate(lag = rep(0:60, 9)) 


# and plot
ggplot(detr.ar1, aes(lag, ar1, color = era)) +
  geom_line() +
  geom_point() +
  facet_wrap(~region, ncol = 1) +
  geom_hline(yintercept = 0, lty = 2) +
  xlab("Lag (months)") +
  ylab("ar1 (detrended anomalies)") +
  scale_x_continuous(breaks = c(0, 12, 24, 36, 48, 60))

## calculate annual anomalies
ersst_annual_anomalies <- ersst_monthly_anomalies %>%
  dplyr::group_by(region, year) %>%
  dplyr::summarise(annual_anomaly = mean(anomaly, na.rm = T)) %>%
  mutate(year = as.numeric(as.character(year)))

# and plot
ggplot(filter(ersst_annual_anomalies), aes(year, annual_anomaly)) +
  geom_line() +
  geom_point() +
  geom_smooth(method = "gam", se = F) +
  facet_wrap(~region, scales = "free_y", ncol = 1) +
  geom_hline(yintercept = 0, lty = 2) + 
  scale_x_continuous(breaks = seq(1920, 2020, 10)) +
  theme(axis.title.x = element_blank()) + 
  ylab("Anomaly (°C, 1920-1949 climatology)")



# now plot annual ar1 values by era, for non-detrended and detrended data

# annual ar1 values for 1950-present
annual.ar1 <- ersst_annual_anomalies %>%
  filter(year >= 1950) %>%
  mutate(era = as.factor(
    case_when(year %in% 1950:1988 ~ "1950-1988",
              year %in% 1989:1999 ~ "1989-1999",
              year %in% 2000:2025 ~ "2000-2025"))) %>%
  group_by(region, era) %>%
  reframe(ar1 = acf(annual_anomaly, lag.max = 5, plot = F)$acf) %>%
  mutate(lag = rep(0:5, 9)) 


# and plot
ggplot(annual.ar1, aes(lag, ar1, color = era)) +
  geom_line() +
  geom_point() +
  facet_wrap(~region, ncol = 1) +
  geom_hline(yintercept = 0, lty = 2) +
  xlab("Lag (years)") 

# detrend 

detr.annual.ar1 <- ersst_annual_anomalies %>%
  filter(year >= 1950) %>%
  group_by(region) %>%
  reframe(time_step = 1:n(),
          detr.anom = resid(lm(annual_anomaly ~ year)))

# plot to check I'm coding correctly
ggplot(detr.annual.ar1, aes(time_step, detr.anom)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0) +
  facet_grid(~region)

detr.annual.ar1 <- detr.annual.ar1 %>%  
  mutate(year = time_step + 1949,
          era = as.factor(
  case_when(year %in% 1950:1988 ~ "1950-1988",
            year %in% 1989:1999 ~ "1989-1999",
            year %in% 2000:2025 ~ "2000-2025"))) %>%
  group_by(region, era) %>%
  reframe(ar1 = acf(detr.anom, lag.max = 5, plot = F)$acf) %>%
  mutate(lag = rep(0:5, 9)) 


# and plot
ggplot(detr.annual.ar1, aes(lag, ar1, color = era)) +
  geom_line() +
  geom_point() +
  facet_wrap(~region, ncol = 1) +
  geom_hline(yintercept = 0, lty = 2) +
  xlab("Lag (years)") +
  ylab("ar1 (detrended anomalies)") 


# now plot winter (NDJFM) mean anomalies

## calculate annual anomalies
ersst_winter_anomalies <- ersst_monthly_anomalies %>%
  filter(month %in% c("Nov", "Dec", "Jan", "Feb", "Mar")) %>%
  mutate(year = as.numeric(as.character(year)), 
         winter.year = if_else(month %in% c("Nov", "Dec"), year+1, year)) %>%
  dplyr::group_by(region, winter.year) %>%
  dplyr::summarise(winter_anomaly = mean(anomaly, na.rm = T)) 

# and plot
ggplot(filter(ersst_winter_anomalies, winter.year < 2026), aes(winter.year, winter_anomaly)) +
  geom_line() +
  geom_point() +
  geom_smooth(method = "gam", se = F) +
  facet_wrap(~region, scales = "free_y", ncol = 1) +
  geom_hline(yintercept = 0, lty = 2) + 
  scale_x_continuous(breaks = seq(1920, 2020, 10)) +
  theme(axis.title.x = element_blank()) + 
  ylab("Anomaly (°C, 1920-1949 climatology)")

# now plot winter ar1 values by era, for non-detrended and detrended data

# winter ar1 values for 1950-present
winter.ar1 <- ersst_winter_anomalies %>%
  filter(winter.year %in% 1950:2025) %>%
  mutate(era = as.factor(
    case_when(winter.year %in% 1950:1988 ~ "1950-1988",
              winter.year %in% 1989:1999 ~ "1989-1999",
              winter.year %in% 2000:2025 ~ "2000-2025"))) %>%
  group_by(region, era) %>%
  reframe(ar1 = acf(winter_anomaly, lag.max = 5, plot = F)$acf) %>%
  mutate(lag = rep(0:5, 9)) 


# and plot
ggplot(winter.ar1, aes(lag, ar1, color = era)) +
  geom_line() +
  geom_point() +
  facet_wrap(~region, ncol = 1) +
  geom_hline(yintercept = 0, lty = 2) +
  xlab("Lag (years)") +
  ylab("Winter ar1")

# detrend 

detr.winter.ar1 <- ersst_winter_anomalies %>%
  filter(winter.year %in% 1950:2025) %>%
  group_by(region) %>%
  reframe(time_step = 1:n(),
          detr.anom = resid(lm(winter_anomaly ~ winter.year)))

# plot to check I'm coding correctly
ggplot(detr.winter.ar1, aes(time_step, detr.anom)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0) +
  facet_grid(~region)

detr.winter.ar1 <- detr.winter.ar1 %>%  
  mutate(winter.year = time_step + 1949,
         era = as.factor(
           case_when(winter.year %in% 1950:1988 ~ "1950-1988",
                     winter.year %in% 1989:1999 ~ "1989-1999",
                     winter.year %in% 2000:2025 ~ "2000-2025"))) %>%
  group_by(region, era) %>%
  reframe(ar1 = acf(detr.anom, lag.max = 5, plot = F)$acf) %>%
  mutate(lag = rep(0:5, 9)) 


# and plot
ggplot(detr.winter.ar1, aes(lag, ar1, color = era)) +
  geom_line() +
  geom_point() +
  facet_wrap(~region, ncol = 1) +
  geom_hline(yintercept = 0, lty = 2) +
  xlab("Lag (years)") +
  ylab("ar1 (detrended anomalies)") 


# now plot winter (NDJFM) mean anomalies











# get rate of change for EBS during 2000-2023
mod <- lm(annual_anomaly ~ year, data = filter(ersst_annual_anomalies, year %in% 2000:2023))
summary(mod)


## annual EBS means for full time series

nc <- nc_open("./downloads/hawaii_soest_31a3_72d5_401e_4a52_c192_7222.nc")

# process

ncvar_get(nc, "time")   # seconds since 1-1-1970
raw <- ncvar_get(nc, "time")
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))
m <- months(d)
yr <- years(d)

x <- ncvar_get(nc, "longitude")
y <- ncvar_get(nc, "latitude")

SST <- ncvar_get(nc, "sst", verbose = F)

SST <- aperm(SST, 3:1)

SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))
lon <- rep(x, each = length(y))
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# EBS/NBS
# 54-66 deg. N, 188-202 deg. E
ebs.x <- c(183, 183, 203, 203, 191) 
ebs.y <- c(53, 65, 65, 57.5, 53)


# define cells within each polygon and plot to check

ebs.sst <- as.data.frame(SST)

xp <- cbind(ebs.x, ebs.y)
loc=cbind(lon, lat)
check <- in.poly(loc, xp=xp)

ebs.sst[,!check] <- NA


# create vector of latitudes to weight mean sst by cell area

cell.weight <- sqrt(cos(lat*pi/180))

# plot to check
hist(cell.weight, breaks = 50)
unique(cell.weight) # looks right

# create a weighted mean function
weighted.cell.mean <- function(x) weighted.mean(x, w = cell.weight, na.rm = T)

weighted.mean.monthly.ebs.sst <- apply(ebs.sst, 1, weighted.cell.mean)

annual.mean <- tapply(weighted.mean.monthly.ebs.sst, yr, mean)

plot_temp <- data.frame(year = 1854:2023,
                        temp = annual.mean)

ggplot(plot_temp, aes(year, temp)) +
  geom_point() +
  geom_line() +
  geom_smooth(se = F) +
  theme(axis.title.x = element_blank()) + 
  scale_x_continuous(breaks = seq(1850, 2020, by = 10), minor_breaks = NULL) +
  scale_y_continuous(breaks = seq(2, 6, by = 0.5), minor_breaks = NULL) +
  labs(y = "Annual mean sea surface temperature (°C)",
       title = "Bering Sea temperature, 1854-2023",
       subtitle = "Data source = thermometer readings, collated by NOAA")

ggsave("./figs/bering_sst_1854_2023.png", width = 7.5, height = 4, units = 'in')


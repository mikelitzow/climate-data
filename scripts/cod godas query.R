library(ncdf4)
library(zoo)
library(gplots)
library(dplyr)
library(maps)
library(mapdata)
library(chron)
library(fields)
library(ggplot2)
library(oce)
library(tidyr)
cb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
theme_set(theme_bw())
# 55-60 N, 204-216 E
# URL <- "http://apdrc.soest.hawaii.edu/erddap/griddap/hawaii_soest_d346_28ac_fccf.nc?potdsl[(1980-01-01):1:(2019-08-01T00:00:00Z)][(105):1:(105)][(55.03699999999998):1:(60.03199999999998)][(204.5):1:(216.5)]"
URL <- "http://apdrc.soest.hawaii.edu/erddap/griddap/hawaii_soest_d346_28ac_fccf.nc?potdsl[(1980-01-01):1:(2019-12-01T00:00:00Z)][(5):1:(205)][(55.03699999999998):1:(60.03199999999998)][(204.5):1:(216.5)]"
# retrieving the data from here:
# https://psl.noaa.gov/data/gridded/data.godas.html
# last accessed 10/3/20
# ##############################################################################################
# # examine larger area
# nc <- nc_open("downloads/test2.nc")
# nc
# 
# # extract dates
# raw <- ncvar_get(nc, "time") # seconds since 1-1-1970
# d <- dates(raw, origin = c(1,1,1800))
# 
# x <- ncvar_get(nc, "lon")     # view longitudes (degrees East)
# y <- ncvar_get(nc, "lat")     # view latitudes
# lat <- rep(y, length(x))   # Vector of latitudes
# lon <- rep(x, each = length(y))   # Vector of longitudes
# # check depth
# ncvar_get(nc, "level")
# 
# # will query the data one depth at a time!
# temp.105 <- ncvar_get(nc, "pottmp") # ºK
# dim(temp.105)
# 
# # process!
# temp.105 <- aperm(temp.105, 3:1)  # First, reverse order of dimensions ("transpose" array)
# 
# temp.105 <- matrix(temp.105, nrow=dim(temp.105)[1], ncol=prod(dim(temp.105)[2:3]))  # Change to matrix
# 
# dimnames(temp.105) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))
# 
# # convert to C
# temp.105 <- temp.105-273.15
# 
# # check
# temp.105.mean <- colMeans(temp.105, na.rm=T)
# z <- t(matrix(temp.105.mean,length(y)))  
# image.plot(x,y,z, col=tim.colors(64))
# contour(x, y, z, add=T)  
# map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
# ################################################################################################################
# process
nc <- nc_open("downloads/X69.161.17.52.281.12.11.48.nc")
nc
# extract dates
raw <- ncvar_get(nc, "time") # seconds since 1-1-1970
d <- dates(raw, origin = c(1,1,1800))
x <- ncvar_get(nc, "lon")     # view longitudes (degrees East)
y <- ncvar_get(nc, "lat")     # view latitudes
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes
# check depth
ncvar_get(nc, "level")
# will query the data one depth at a time!
temp1.15 <- ncvar_get(nc, "pottmp") # ºK
dim(temp.15)
# process!
temp1.15 <- aperm(temp1.15, 3:1)  # First, reverse order of dimensions ("transpose" array)
temp1.15 <- matrix(temp1.15, nrow=dim(temp1.15)[1], ncol=prod(dim(temp1.15)[2:3]))  # Change to matrix
dimnames(temp1.15) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))
# now the second half of the TS
nc <- nc_open("downloads/X69.161.17.52.281.12.12.36.nc")
nc
# extract dates
raw <- ncvar_get(nc, "time") # seconds since 1-1-1970
d2 <- dates(raw, origin = c(1,1,1800))
# check depth
ncvar_get(nc, "level")
# will query the data one depth at a time!
temp2.15 <- ncvar_get(nc, "pottmp") # ºK
# process!
temp2.15 <- aperm(temp2.15, 3:1)  # First, reverse order of dimensions ("transpose" array)
temp2.15 <- matrix(temp2.15, nrow=dim(temp2.15)[1], ncol=prod(dim(temp2.15)[2:3]))  # Change to matrix
dimnames(temp2.15) <- list(as.character(d2), paste("N", lat, "E", lon, sep=""))
temp.15 <- rbind(temp1.15, temp2.15)
# convert to C
temp.15 <- temp.15-273.15
# drop some of the abyssal area
drop <- lat < y[6] & lon == 207.5
temp.15[,drop] <- NA
drop <- lat < y[5] & lon %in% c(206.5, 205.5)
temp.15[,drop] <- NA
drop <- lat < y[4] & lon == 204.5
temp.15[,drop] <- NA
drop <- lat < y[3] & lon %in% c(202.5, 203.5)
temp.15[,drop] <- NA
drop <- lat < y[2] & lon == 201.5
temp.15[,drop] <- NA
drop <- lat > y[3] & lon == 199.5
temp.15[,drop] <- NA
# check
temp.15.mean <- colMeans(temp.15, na.rm=T)
z <- t(matrix(temp.15.mean,length(y)))  
image.plot(x,y,z, col=tim.colors(64))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)

# and save the plot
png("figs/godas temperature query area.png", 5, 4, units="in", res=300)
temp.mean <- colMeans(temp.15, na.rm=T)
z <- t(matrix(temp.mean,length(y)))  
image.plot(x,y,z, col=oceColorsPalette(64), xlim=c(198,209), ylim=c(54,59), ylab="", xlab="")
# contour(x, y, z, add=T)  
map('world2Hires',c('Canada', 'usa'), fill=T,add=T, lwd=1, col="lightyellow3")
dev.off()

# now 105m
nc <- nc_open("downloads/X69.161.17.52.281.12.12.8.nc")
nc
# extract dates
raw <- ncvar_get(nc, "time") # seconds since 1-1-1970
d <- dates(raw, origin = c(1,1,1800))
# check depth
ncvar_get(nc, "level")
# will query the data one depth at a time!
temp1.105 <- ncvar_get(nc, "pottmp") # ºK
# process!
temp1.105 <- aperm(temp1.105, 3:1)  # First, reverse order of dimensions ("transpose" array)
temp1.105 <- matrix(temp1.105, nrow=dim(temp1.105)[1], ncol=prod(dim(temp1.105)[2:3]))  # Change to matrix
dimnames(temp1.105) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))
# now the second half of the TS
nc <- nc_open("downloads/X69.161.17.52.281.12.12.52.nc")
nc
# extract dates
raw <- ncvar_get(nc, "time") # seconds since 1-1-1970
d2 <- dates(raw, origin = c(1,1,1800))
# check depth
ncvar_get(nc, "level")
# will query the data one depth at a time!
temp2.105 <- ncvar_get(nc, "pottmp") # ºK
# process!
temp2.105 <- aperm(temp2.105, 3:1)  # First, reverse order of dimensions ("transpose" array)
temp2.105 <- matrix(temp2.105, nrow=dim(temp2.105)[1], ncol=prod(dim(temp2.105)[2:3]))  # Change to matrix
dimnames(temp2.105) <- list(as.character(d2), paste("N", lat, "E", lon, sep=""))
temp.105 <- rbind(temp1.105, temp2.105)
# convert to C
temp.105 <- temp.105-273.15
# drop some of the abyssal area
drop <- lat < y[6] & lon == 207.5
temp.105[,drop] <- NA
drop <- lat < y[5] & lon %in% c(206.5, 205.5)
temp.105[,drop] <- NA
drop <- lat < y[4] & lon == 204.5
temp.105[,drop] <- NA
drop <- lat < y[3] & lon %in% c(202.5, 203.5)
temp.105[,drop] <- NA
drop <- lat < y[2] & lon == 201.5
temp.105[,drop] <- NA
drop <- lat > y[3] & lon == 199.5
temp.105[,drop] <- NA
# check
temp.105.mean <- colMeans(temp.105, na.rm=T)
z <- t(matrix(temp.105.mean,length(y)))  
image.plot(x,y,z, col=tim.colors(64))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)

#################
# now the monthly anomalies for both!
d <- c(d, d2)
m <- months(d)
yr <- years(d)
f <- function(x) tapply(x, m[yr %in% 1980:2009], mean)
ff <- function(x) tapply(x, m[yr %in% 1980:2009], sd)
mu.15 <- apply(temp.15[yr %in% 1980:2009,], 2, f)	# Compute monthly means for each time series (location)
mu.15 <- mu.15[rep(1:12, floor(length(d)/12)),] 
xtra <- 12*((length(d)/12)-floor(length(d)/12))
mu.15 <- rbind(mu.15, mu.15[1:xtra,])
sd.15 <- apply(temp.15[yr %in% 1980:2009,], 2, ff)	# Compute monthly sd for each time series (location)
sd.15 <- sd.15[rep(1:12, floor(length(d)/12)),] 
sd.15 <- rbind(sd.15, sd.15[1:xtra,])
anom.15 <- rowMeans((temp.15 - mu.15)/sd.15, na.rm=T)   # Compute matrix of anomalies!
raw.15 <- rowMeans(temp.15, na.rm=T)

# now 105m
mu.105 <- apply(temp.105[yr %in% 1980:2009,], 2, f)	# Compute monthly means for each time series (location)
mu.105 <- mu.105[rep(1:12, floor(length(d)/12)),] 
xtra <- 12*((length(d)/12)-floor(length(d)/12))
mu.105 <- rbind(mu.105, mu.105[1:xtra,])
sd.105 <- apply(temp.105[yr %in% 1980:2009,], 2, ff)	# Compute monthly sd for each time series (location)
sd.105 <- sd.105[rep(1:12, floor(length(d)/12)),] 
sd.105 <- rbind(sd.105, sd.105[1:xtra,])
anom.105 <- rowMeans((temp.105 - mu.105)/sd.105, na.rm=T)   # Compute matrix of anomalies!
raw.105 <- rowMeans(temp.105, na.rm=T)

yr <- as.numeric(as.character(yr))
dec.yr <- yr + (as.numeric(m)-0.5)/12 
plot.dat <- data.frame(year=dec.yr,
                       anom.15=anom.15,
                       anom.105=anom.105)
plot.dat <- plot.dat %>%
  pivot_longer(cols = -year)

ggplot(plot.dat, aes(year, value, color=name)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values=cb[c(2,4)])

# summarize Jan-April 105m
jan.apr.m <- m[m %in% c("Jan", "Feb", "Mar", "Apr")]
jan.apr.yr <- yr[m %in% c("Jan", "Feb", "Mar", "Apr")]


jan.apr.105 <- anom.105[m %in% c("Jan", "Feb", "Mar", "Apr")]
jan.apr.105 <- tapply(jan.apr.105, jan.apr.yr, mean)

raw.jan.apr.105 <- raw.105[m %in% c("Jan", "Feb", "Mar", "Apr")]
raw.jan.apr.105 <- tapply(raw.jan.apr.105, jan.apr.yr, mean)


# and Apr-Jun 15m
apr.jun.m <- m[m %in% c("Apr", "May", "Jun")]
apr.jun.yr <- yr[m %in% c("Apr", "May", "Jun")]


apr.jun.15 <- anom.15[m %in% c("Apr", "May", "Jun")]
apr.jun.15 <- tapply(apr.jun.15, apr.jun.yr, mean)

raw.apr.jun.15 <- raw.15[m %in% c("Apr", "May", "Jun")]
raw.apr.jun.15 <- tapply(raw.apr.jun.15, apr.jun.yr, mean)


godas.anom <- data.frame(year=1980:2020,
                         apr.jun.15=apr.jun.15,
                         jan.apr.105=jan.apr.105)

plot.anom <- godas.anom %>%
  pivot_longer(cols=-year)

ggplot(plot.anom, aes(year, value, color=name)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values=cb[c(2,4)])

save.dat <- godas.anom

save.dat$mean.anom <- apply(save.dat[,2:3], 1, mean)

plot.dat <- save.dat %>%
  pivot_longer(cols = -year)

ggplot(plot.dat, aes(year, value, color=name)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values=cb[c(2:4)])

write.csv(save.dat, "output data/godas anomalies.csv")

## and raw
godas.raw <- data.frame(year=1980:2020,
                         apr.jun.15=raw.apr.jun.15,
                         jan.apr.105=raw.jan.apr.105)

plot.raw <- godas.raw %>%
  pivot_longer(cols=-year)

ggplot(plot.raw, aes(year, value, color=name)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values=cb[c(2,4)])

save.dat <- godas.raw

plot.dat <- save.dat %>%
  pivot_longer(cols = -year)

ggplot(plot.dat, aes(year, value, color=name)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values=cb[c(2:4)])

write.csv(save.dat, "output data/godas raw.csv")


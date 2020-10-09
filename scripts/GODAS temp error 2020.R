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

# 55-60 N, 204-216 E

# URL <- "http://apdrc.soest.hawaii.edu/erddap/griddap/hawaii_soest_d346_28ac_fccf.nc?potdsl[(1980-01-01):1:(2019-08-01T00:00:00Z)][(105):1:(105)][(55.03699999999998):1:(60.03199999999998)][(204.5):1:(216.5)]"
URL <- "http://apdrc.soest.hawaii.edu/erddap/griddap/hawaii_soest_d346_28ac_fccf.nc?potdsl[(1980-01-01):1:(2020-08-01T00:00:00Z)][(105):1:(105)][(55.03699999999998):1:(60.03199999999998)][(204.5):1:(216.5)]"

# retrieving the data from here:
# https://psl.noaa.gov/data/gridded/data.godas.html
# last accessed 10/3/20

nc <- nc_open("downloads/hawaii_soest_d346_28ac_fccf_f083_0d6a_0ab1.nc")
temp.105 <- ncvar_get(nc, "potdsl") # ºK
temp.105[,,488]
############################
# process
nc <- nc_open("downloads/hawaii_soest_d346_28ac_fccf_f083_0d6a_0ab1.nc")
nc

# extract dates
raw <- ncvar_get(nc, "time") # seconds since 1-1-1970
h <- raw/(60*60*24)
d <- dates(h, origin = c(1,1,1970))

x <- ncvar_get(nc, "longitude")     # view longitudes (degrees East)
y <- ncvar_get(nc, "latitude")     # view latitudes
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes
# check depth
ncvar_get(nc, "LEV")

# will query the data one depth at a time!
temp.105 <- ncvar_get(nc, "potdsl") # ºK
dim(temp.105)

temp.105[,,488]

# process!
temp.105 <- aperm(temp.105, 3:1)  # First, reverse order of dimensions ("transpose" array)

temp.105 <- matrix(temp.105, nrow=dim(temp.105)[1], ncol=prod(dim(temp.105)[2:3]))  # Change to matrix

dimnames(temp.105) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# convert to C
temp.105 <- temp.105-273.15

# check
temp.105.mean <- colMeans(temp.105, na.rm=T)
z <- t(matrix(temp.105.mean,length(y)))  
image.plot(x,y,z, col=tim.colors(64), xlim=c(195,230), ylim=c(53,62))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)

hist(temp.105)

drop <- temp.105 < -5.0e+33
temp.105[drop] <- NA


png("figs/godas temperature error query area .png", 6, 4, units="in", res=300)
temp.mean <- colMeans(temp.105, na.rm=T)
z <- t(matrix(temp.mean,length(y)))  
image.plot(x,y,z, col=oceColorsPalette(64), xlim=c(195,230), ylim=c(53,62))
contour(x, y, z, add=T)  
map('world2Hires',c('Canada', 'usa'), fill=T,xlim=c(130,250), ylim=c(20,66),add=T, lwd=1, col="lightyellow3")
dev.off()

# good 'nuff for now
dec.yr <- as.numeric(as.character(years(d))) + (as.numeric(months(d))-0.5)/12
temp.105ts <- rowMeans(temp.105, na.rm = T)
plot.temp <- data.frame(year=dec.yr, temp=temp.105ts) 

ggplot(plot.temp, aes(year, temp)) +
  theme_bw() +
  geom_line() +
  geom_point()

ggsave("figs/godas 105m pottemp time series.png", width=6, height=4, units='in')

View(temp.105)
keep <- !is.na(colMeans(temp.105))

length <- nrow(temp.105)
share <- temp.105[((length-20):length),keep]

write.csv(share, "downloads/godas.temp.csv")

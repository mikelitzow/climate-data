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


# script for calculating GOA sst anomalies wrt 1951-1980
# download.file("https://coastwatch.pfeg.noaa.gov/erddap/griddap/nceiErsstv5.nc?sst[(1950-01-01):1:(2021-8-01T00:00:00Z)][(0.0):1:(0.0)][(54):1:(62)][(200):1:(226)]", "~temp")

# paste into browser for windows!


# load and process SST data
# nc <- nc_open("~temp")

nc <- nc_open("downloads/nceiErsstv5_1fbf_1ae0_53f2.nc")

# extract dates

ncvar_get(nc, "time")   # seconds since 1-1-1970
raw <- ncvar_get(nc, "time")
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))

# extract study area
# 54-62 deg. N, 200-226 deg. E
x <- ncvar_get(nc, "longitude")
y <- ncvar_get(nc, "latitude")

SST <- ncvar_get(nc, "sst", verbose = F)

# Change data from a 3-D array to a matrix of monthly data by grid point:
# First, reverse order of dimensions ("transpose" array)
SST <- aperm(SST, 3:1)  

# Change to matrix with column for each grid point, rows for monthly means
SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   
lon <- rep(x, each = length(y))   
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# plot to check

# need to drop Bristol Bay cells
BB <- c("N58E200", "N58E202", "N56E200")
SST[,BB] <- NA

# and check
temp.mean <- colMeans(SST, na.rm=T)
z <- t(matrix(temp.mean,length(y)))  
image.plot(x,y,z, col=oceColorsPalette(64), xlim=c(195,230), ylim=c(53,62))
contour(x, y, z, add=T)  
map('world2Hires',c('Canada', 'usa'), fill=T,xlim=c(130,250), ylim=c(20,66),add=T, lwd=1, col="lightyellow3")


###################################################
# WGOA for cod/pollock paper

drop <- lon > 212 | lat < 56
wSST <- SST
wSST[,drop] <- NA


# and check - save spatial domain map
png("./figs/wgoa.ersstv5.domain.annual.means.png", 6, 4, units = 'in', res = 300)
temp.mean <- colMeans(wSST, na.rm=T)
z <- t(matrix(temp.mean,length(y)))  
image.plot(x,y,z, col=oceColorsPalette(64), xlim=c(195,230), ylim=c(53,62))
contour(x, y, z, add=T)  
map('world2Hires',c('Canada', 'usa'), fill=T,xlim=c(130,250), ylim=c(20,66),add=T, lwd=1, col="lightyellow3")
dev.off()

# calculate monthly anomaly

wSST <- rowMeans(wSST, na.rm = T)

# and Apr-Jun
yr <- as.numeric(as.character(years(d)))
m <- months(d)
apr.jun.m <- m[m %in% c("Apr", "May", "Jun")]
apr.jun.yr <- yr[m %in% c("Apr", "May", "Jun")]

apr.jun.wSST <- wSST[m %in% c("Apr", "May", "Jun")]
apr.jun.wSST <- tapply(apr.jun.wSST, apr.jun.yr, mean)

xprt <- data.frame(year=1950:2021,
                   apr.jun.wSST=apr.jun.wSST)

ggplot(xprt, aes(year, apr.jun.wSST)) +
  geom_line() +
  geom_point()

write.csv(xprt, "output data/wgoa.apr.jun.sst.csv", row.names = F)

 # now JFM western GOA
jan.mar.m <- m[m %in% c("Jan", "Feb", "Mar")]
jan.mar.yr <- yr[m %in% c("Jan", "Feb", "Mar")]

jan.mar.wSST <- wSST[m %in% c("Jan", "Feb", "Mar")]
jan.mar.wSST <- tapply(jan.mar.wSST, jan.mar.yr, mean)

xprt <- data.frame(year=1950:2021,
                   jan.feb.mar.SST=jan.mar.wSST)

ggplot(xprt, aes(year, jan.feb.mar.SST)) +
  geom_line() +
  geom_point()

write.csv(xprt, "output data/wgoa.jan.feb.mar.sst.csv", row.names = F)

# f <- function(x) tapply(x, m, mean)
# mu <- apply(SST, 2, f)	# Compute monthly means for each time series (location)
# 
# 
# mu <- mu[rep(1:12, floor(length(d)/12)),] 
# 
# # add trailing months
#  
# xtra <- 12*((length(d)/12)-floor(length(d)/12))
# 
# mu <- rbind(mu, mu[1:xtra,])
# 
# anom <- SST - mu # Compute matrix of anomalies!

# and winter for GOA-wide
SST <- rowMeans(SST, na.rm = T)
 

# set up winter year
win.yr <- yr
win.yr[m %in% c("Nov", "Dec")] <- win.yr[m %in% c("Nov", "Dec")]+1

sst.win <- SST[m %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]
win.yr <- win.yr[m %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]

win.sst <- tapply(sst.win, win.yr, mean)

plot(1951:2020, win.sst[2:(length(win.sst))], type="o")


xprt <- data.frame(year=1951:2020,
                   ndjfm.sst=win.sst[2:(length(win.sst))])

write.csv(xprt, "output data/goa.winter.sst.csv", row.names = F)





png("GOA SST anom 1950-2019.png", 6, 4, units="in", res=300)
par(las=1, mar=c(3,4,1,1), mgp=c(2,0.8,0))
plot(dec.yr, anom, type="l", xlab="", ylab = "Anomaly (ºC)", lwd=0.5)
abline(h=0)
#abline(h=seq(-3,3,0.5), lty=3, lwd=0.5)
#abline(v=seq(1920,2020,5), lty=3, lwd=0.5)
lines(dec.yr, rollmean(anom,36,fill=NA), col="red", lwd=1.5)
dev.off()

sst.sd <- rollapply(anom, 120, sd, align="right", fill=NA)

par(las=1)
plot(dec.yr, sst.sd, type="l", xlab="", ylab = "10-yr running SD (ºC)", col="grey", xlim=c(1920,max(dec.yr)), ylim=range(sst.sd[yr >= 1920]))
abline(h=mean(sst.sd[yr >= 1920]))

par(las=1)
plot(dec.yr, anom, type="l", xlab="", ylab = "Anomaly wrt 1951-1980 (ºC)", col="grey", xlim=c(1951,max(dec.yr)), ylim=range(anom[yr >= 1951]))
abline(h=0)
abline(h=seq(-3,3,0.5), lty=3, lwd=0.5)
abline(v=seq(1920,2020,5), lty=3, lwd=0.5)
lines(dec.yr, rollmean(anom,25,fill=NA), col="red")

annual <- tapply(anom, yr, mean)
plot(names(annual), annual, type="o", pch=19, xlim=c(1951,2018))

# get means in Fº for presentation
short.yr <- yr[yr < 2019]

short.SST <- SST[m %in% c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct"),]

short.annual.mean <- tapply(rowMeans(short.SST, na.rm = T), yr[m %in% c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct")], mean)

# change to F

short.annual.mean.F <- (short.annual.mean*9/5) + 32

# now estimate 2019 F value from annual mean
annual <- tapply(rowMeans(SST, na.rm = T), yr, mean)
annual.F <- (annual*9/5) + 32

mod <- lm(annual.F[names(annual.F) < 2019] ~ short.annual.mean.F[names(short.annual.mean.F) < 2019])
summary(mod)

est <- annual.F[70]*coef(mod)[2] + coef(mod)[1]

plot.annual.F <- data.frame(year=1950:2019,
                            temp=c(annual.F[1:69], est))

ggplot(plot.annual.F, aes(year, temp)) +
  theme_bw() +
  geom_line() + 
  geom_point() +
  geom_hline(yintercept = mean(annual.F[names(annual.F) %in% 1981:2010]), lty=2) +
  ggtitle("Gulf of Alaska annual sea surface temperature") +
  theme(axis.title.x = element_blank()) +
  scale_x_continuous(breaks=seq(1950, 2020, 10)) +
  ylab("º Fahrenheit")

ggsave("GOA annual SST fahrenheit.png", height = 3, width = 5, units="in")

ann.anom <- data.frame(year=as.numeric(names(annual)), anomaly=annual)

short.anom <- filter(ann.anom, year >=1990)

png("GOA SST anom 1990-2018.png", 6, 4, units="in", res=300)
ggplot(short.anom, aes(year, anomaly)) + geom_hline(yintercept = 0, col="dark grey") + geom_line() + geom_point() + 
  xlab("") + ylab("ºC") + ggtitle("GOA SST anomaly: difference from 1981-2010 mean", subtitle= "Data: ERSST v5")
dev.off()


arrange(short.anom, desc(anomaly))

win.anom <- anom[m %in% 1:3]
win.yr <- yr[m %in% 1:3]

win.annual <- tapply(win.anom, win.yr, mean)
plot(names(win.annual), win.annual, type="o", pch=19, xlim=c(1951,2018), ylim=c(-1.3,2))

# plot as ºfahrenheit

# get anomaly for 1951:1980

yr <- as.numeric(as.character(years(d)))

m <- months(d[yr %in% 1951:1980])

SST.f <- SST*1.8 +32

f <- function(x) tapply(x, m, mean)
mu <- apply(SST.f[yr %in% 1951:1980,], 2, f)	# Compute monthly means for each time series (location)
mu <- mu[rep(1:12, floor(length(d)/12)),] 

xtra <- 12*((length(d)/12)-floor(length(d)/12))

mu <- rbind(mu, mu[1:xtra,])

anom.f <- rowMeans(SST.f - mu, na.rm=T)   # Compute matrix of anomalies!

m <- as.numeric(months(d))

dec.yr <- yr + (m-0.5)/12


par(las=1)

annual.f <- tapply(anom.f, yr, mean)
png("GOA annual anomaly F.png", 6.5, 4, units="in", res=300)
par(las=1, mar=c(3,4,1,1), cex.axis=1.3)
plot(names(annual.f), annual.f, type="o", pch=19, xlim=c(1951,2018), ylim=c(-2.3,3.5),
     ylab="Difference from 1951-1980 mean (ºF)", xlab="")
abline(h=0, lty=2)
dev.off()

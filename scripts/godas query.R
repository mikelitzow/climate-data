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
temp.mean <- colMeans(temp.105, na.rm=T)
z <- t(matrix(temp.mean,length(y)))  
image.plot(x,y,z, col=oceColorsPalette(64), xlim=c(200,215), ylim=c(53,62))
contour(x, y, z, add=T)  
map('world2Hires',c('Canada', 'usa'), fill=T,xlim=c(130,250), ylim=c(20,66),add=T, lwd=1, col="lightyellow3")
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

# now 105m
mu.105 <- apply(temp.105[yr %in% 1980:2009,], 2, f)	# Compute monthly means for each time series (location)
mu.105 <- mu.105[rep(1:12, floor(length(d)/12)),] 

xtra <- 12*((length(d)/12)-floor(length(d)/12))

mu.105 <- rbind(mu.105, mu.105[1:xtra,])

sd.105 <- apply(temp.105[yr %in% 1980:2009,], 2, ff)	# Compute monthly sd for each time series (location)
sd.105 <- sd.105[rep(1:12, floor(length(d)/12)),] 

sd.105 <- rbind(sd.105, sd.105[1:xtra,])


anom.105 <- rowMeans((temp.105 - mu.105)/sd.105, na.rm=T)   # Compute matrix of anomalies!

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



##

# summarize Jan-April 105m
jan.apr.m <- m[m %in% c("Jan", "Feb", "Mar", "Apr")]
jan.apr.yr <- yr[m %in% c("Jan", "Feb", "Mar", "Apr")]


jan.apr.105 <- anom.105[m %in% c("Jan", "Feb", "Mar", "Apr")]
jan.apr.105 <- tapply(jan.apr.105, jan.apr.yr, mean)

# and Apr-Jun 15m
apr.jun.m <- m[m %in% c("Apr", "May", "Jun")]
apr.jun.yr <- yr[m %in% c("Apr", "May", "Jun")]


apr.jun.15 <- anom.15[m %in% c("Apr", "May", "Jun")]
apr.jun.15 <- tapply(apr.jun.15, apr.jun.yr, mean)

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

#################
# old stuff below!
m2 <- months(d2)
yr2 <- years(d2)

temp.105.combined <- rbind(temp2.105, temp.105)




jan.apr.m <- m[m %in% c("Jan", "Feb", "Mar", "Apr")]
jan.apr.yr <- yr[m %in% c("Jan", "Feb", "Mar", "Apr")]
jan.apr.temp <- temp.105[m %in% c("Jan", "Feb", "Mar", "Apr")]

jan.apr.temp <- tapply(jan.apr.temp, jan.apr.yr, mean)

plot.spawn.temp <- data.frame(year=as.numeric(names(jan.apr.temp)), temp=jan.apr.temp)

ggplot(plot.spawn.temp, aes(year, temp)) +
  theme_bw() +
  geom_line() +
  geom_point() +
  ggtitle("GODAS 105m temp, annual Jan-Apr means")

temp.data <- plot.spawn.temp
names(temp.data)[2] <- "godas.105m.Jan.Apr"

###########
# now 205m
nc <- nc_open("downloads/X161.55.246.38.277.19.46.13.nc")
nc

# extract dates
raw <- ncvar_get(nc, "time") # seconds since 1-1-1970
d <- dates(raw, origin = c(1,1,1800))

ncvar_get(nc, "level")

# will query the data one depth at a time!
temp.205 <- ncvar_get(nc, "pottmp") # ºK
dim(temp.205)

# process!
temp.205 <- aperm(temp.205, 3:1)  # First, reverse order of dimensions ("transpose" array)

temp.205 <- matrix(temp.205, nrow=dim(temp.205)[1], ncol=prod(dim(temp.205)[2:3]))  # Change to matrix

dimnames(temp.205) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# convert to C
temp.205 <- temp.205-273.15

temp.205[mask.drop] <- NA

# check
temp.205.mean <- colMeans(temp.205, na.rm=T)
z <- t(matrix(temp.205.mean,length(y)))  
image.plot(x,y,z, col=tim.colors(64), xlim=c(195,230), ylim=c(53,62))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)


temp.205 <- rowMeans(temp.205, na.rm = T)
plot.temp <- data.frame(year=dec.yr, temp=temp.205) 

ggplot(plot.temp, aes(year, temp)) +
  theme_bw() +
  geom_line() +
  geom_point()


jan.apr.temp <- temp.205[m %in% c("Jan", "Feb", "Mar", "Apr")]

jan.apr.temp <- tapply(jan.apr.temp, jan.apr.yr, mean)

plot.spawn.temp <- data.frame(year=as.numeric(names(jan.apr.temp)), temp=jan.apr.temp)

ggplot(plot.spawn.temp, aes(year, temp)) +
  theme_bw() +
  geom_line() +
  geom_point() +
  ggtitle("GODAS 105m temp, annual Jan-Apr means")

temp.data$godas.205m.Jan.Apr <- jan.apr.temp

# now 15m
nc <- nc_open("downloads/X161.55.246.38.277.19.46.37.nc")
nc

ncvar_get(nc, "level")

# will query the data one depth at a time!
temp.15 <- ncvar_get(nc, "pottmp") # ºK
dim(temp.15)

# process!
temp.15 <- aperm(temp.15, 3:1)  # First, reverse order of dimensions ("transpose" array)

temp.15 <- matrix(temp.15, nrow=dim(temp.15)[1], ncol=prod(dim(temp.15)[2:3]))  # Change to matrix

dimnames(temp.15) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# convert to C
temp.15 <- temp.15-273.15

temp.15[mask.drop] <- NA

# check
temp.15.mean <- colMeans(temp.15, na.rm=T)
z <- t(matrix(temp.15.mean,length(y)))  
image.plot(x,y,z, col=tim.colors(64), xlim=c(195,230), ylim=c(53,62))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)

temp.15 <- rowMeans(temp.15, na.rm = T)
plot.temp <- data.frame(year=dec.yr, temp=temp.15) 

ggplot(plot.temp, aes(year, temp)) +
  theme_bw() +
  geom_line() +
  geom_point()


apr.jun.m <- m[m %in% c("Apr", "May", "Jun")]
apr.jun.yr <- yr[m %in% c("Apr", "May", "Jun")]
apr.jun.temp <- temp.15[m %in% c("Apr", "May", "Jun")]

apr.jun.temp <- tapply(apr.jun.temp, apr.jun.yr, mean)

plot.spawn.temp <- data.frame(year=as.numeric(names(apr.jun.temp )), temp=apr.jun.temp)

ggplot(plot.spawn.temp, aes(year, temp)) +
  theme_bw() +
  geom_line() +
  geom_point() +
  ggtitle("GODAS 15m temp, annual Apr-Jun means")

temp.data$godas.15m.Apr.Jun <- apr.jun.temp

# add ERSST GOA-wide for Apr-Jun
wgoa.sst <- read.csv("output data/wgoa.apr.jun.sst.csv")

temp.data$ersst.Apr.Jun <- wgoa.sst$apr.jun.wSST[match(temp.data$year, wgoa.sst$year)]

plot.dat <- temp.data %>%
  pivot_longer(cols=-year)

ggplot(plot.dat, aes(year, value, color=name)) +
  geom_point() +
  geom_line()

# so godas 15m functionally the same as ersst!

# now compare with GAK1!
gak <- read.csv("downloads/gak1.csv")
gak <- gak[-(1:2),]
gak$date <- lubridate::date_decimal(as.numeric(as.character(gak$Year)))
gak <- gak %>%
  select(-Year)
gak$year <- lubridate::year(gak$date)
gak$month <- lubridate::month(gak$date)
gak$Depth <- as.numeric(as.character(gak$Depth))
gak$Temp <- as.numeric(as.character(gak$Temp))

gak.jan.apr <- gak %>%
  filter(month %in% 1:4)

tapply(gak.jan.apr$Temp, list(gak.jan.apr$year, gak.jan.apr$Depth), mean, na.rm=T)

ff <- function(x) sum(!is.na(x))

tapply(gak.jan.apr$Temp, list(gak.jan.apr$year, gak.jan.apr$Depth), ff)

gak.100 <- gak.jan.apr %>%
  filter(Depth==100) %>%
  group_by(year) %>%
  summarise(mu=mean(Temp))

gak.200 <- gak.jan.apr %>%
  filter(Depth==200) %>%
  group_by(year) %>%
  summarise(mu=mean(Temp))

temp.data$gak.100.Jan.Apr <- gak.100$mu[match(temp.data$year , gak.100$year)]
temp.data$gak.200.Jan.Apr <- gak.200$mu[match(temp.data$year , gak.200$year)]

plot.dat <- temp.data %>%
  pivot_longer(cols=-year)

ggplot(plot.dat, aes(year, value, color=name)) +
  geom_point() +
  geom_line()

# plot as a multi-panel!
temp.surf <- data.frame(year=2006:2020,
                        depth="Surface Apr-Jun",
                        observed=temp.data$ersst.Apr.Jun,
                        GODAS=temp.data$godas.15m.Apr.Jun)

temp.100 <- data.frame(year=2006:2020,
                        depth="100m Jan-Apr",
                        observed=temp.data$gak.100.Jan.Apr,
                        GODAS=temp.data$godas.105m.Jan.Apr)

temp.200 <- data.frame(year=2006:2020,
                       depth="200m Jan-Apr",
                       observed=temp.data$gak.200.Jan.Apr,
                       GODAS=temp.data$godas.205m.Jan.Apr)

new.plot <- rbind(temp.surf, temp.100, temp.200)

new.plot <- new.plot %>%
  pivot_longer(cols=c(observed, GODAS))

ggplot(new.plot, aes(year, value, color=name)) +
         geom_point() +
         geom_line() +
         facet_wrap(~depth) +
        ylab("°C")

ggsave("figs/godas ersst gak1.png", width=6, height=3, units='in')


# and fit a DFA model
library(MARSS)

# change to matrix
# dfa.dat <- t(as.matrix(temp.data[,2:7]))
# dfa.dat <- t(as.matrix(temp.data[,c(2,4,6)])) # restricting to surface and 100m
dfa.dat <- t(as.matrix(temp.data[,c(2,4)]))

# changing convergence criterion to ensure convergence
cntl.list = list(minit=200, maxit=20000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)

# set up forms of R matrices
levels.R = c("diagonal and equal",
             "diagonal and unequal",
             "equalvarcov",
             "unconstrained")
model.data = data.frame()

# fit models & store results
for(R in levels.R) {
  for(m in 1) {  # allowing up to 3 trends
    dfa.model = list(A="zero", R=R, m=m)
    kemz = MARSS(dfa.dat, model=dfa.model, control=cntl.list,
                 form="dfa", z.score=TRUE)
    model.data = rbind(model.data,
                       data.frame(R=R,
                                  m=m,
                                  logLik=kemz$logLik,
                                  K=kemz$num.params,
                                  AICc=kemz$AICc,
                                  stringsAsFactors=FALSE))
    assign(paste("kemz", m, R, sep="."), kemz)
  } # end m loop
} # end R loop

# calculate delta-AICc scores, sort in descending order, and compare
model.data$dAICc <- model.data$AICc-min(model.data$AICc)
model.data <- model.data %>%
  arrange(dAICc)
model.data

# plot
model.list = list(A="zero", m=1, R="diagonal and equal")
mod = MARSS(dfa.dat, model=model.list, z.score=TRUE, form="dfa", control=cntl.list)

# get CI and plot loadings...
modCI <- MARSSparamCIs(mod)
modCI

plot.CI <- data.frame(names=rownames(dfa.dat),
                      mean=modCI$par$Z,
                      upCI=modCI$par.upCI$Z,
                           lowCI=modCI$par.lowCI$Z)

plot.CI. <- arrange(plot.CI, mean)
dodge <- position_dodge(width=0.9)



ggplot(plot.CI, aes(x=names, y=mean)) + 
  geom_bar(position="dodge", stat="identity", fill=cb[2]) +
  geom_errorbar(aes(ymax=upCI, ymin=lowCI), position=dodge, width=0.5) +
  ylab("") + xlab("") + theme_bw() + theme(axis.text.x  = element_text(angle=45, hjust=1, size=12)) +
  geom_hline(yintercept = 0) + ggtitle("1989-2012 community") 

# and plot

plot.dat <- data.frame(year=2006:2020,
                       godas.15m.Apr.Jun=scale(temp.data$godas.15m.Apr.Jun),
                       godas.105m.Jan.Apr=scale(temp.data$godas.105m.Jan.Apr), 
                       trend=as.vector(mod$states))

plot.dat <- plot.dat %>%
  pivot_longer(cols = -year) %>%
  arrange(desc(name))

plot.dat$trendUCI = c(as.vector(mod$states+1.96*mod$states.se), rep(NA, 30))
plot.dat$trendLCI = c(as.vector(mod$states-1.96*mod$states.se), rep(NA, 30))

ggplot(plot.dat, aes(year, value, color=name)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values=cb[2:4]) + 
  geom_ribbon(aes(ymin=trendLCI, ymax=trendUCI), alpha=0.2)

# this plot is problematic - and proabably not needed!!

# finally! compare with winter ERSST we were using before
win.sst <- read.csv("output data/goa.winter.sst.csv")

compare <- data.frame(year=2006:2020,
                      winter.ersst=win.sst$ndjfm.sst[win.sst$year %in% 2006:2020],
                      dfa.trend=as.vector(mod$states))

ggplot(compare, aes(winter.ersst, dfa.trend)) +
  geom_point()

comp2 <- compare %>%
  pivot_longer(cols= -year)

ggplot(comp2, aes(year, value, color=name)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values=cb[c(2,4)])

write.csv(compare, "output data/winter sst and godas dfa trend.csv")

# might be better to just take the average anomaly of the two!
plot.dat <- data.frame(year=2006:2020,
                       godas.15m.Apr.Jun=scale(temp.data$godas.15m.Apr.Jun),
                       godas.105m.Jan.Apr=scale(temp.data$godas.105m.Jan.Apr), 
                       trend=as.vector(mod$states))

plot.dat$mean.anomaly=apply(plot.dat[,2:3], 1, mean)

ggplot(plot.dat, aes(trend, mean.anomaly)) +
  geom_point()

plot2 <- plot.dat %>%
  select(year, trend, mean.anomaly) %>%
  pivot_longer(cols=-year)

ggplot(plot2, aes(year, value, color=name)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values=cb[c(2,4)])

write.csv(plot.dat, "output data/winter sst and godas dfa trend.csv")


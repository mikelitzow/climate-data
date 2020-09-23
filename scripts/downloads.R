library(ncdf4)
library(maps)
library(maptools)
library(mapdata)
library(fields)



# load ERSSTv5 data for the N. Pacific
# 20º-70ºN, 120º-250ºE, 1854-present

# identify latest year and month needed
year <- 2019
month <- "07"

URL <- paste("https://coastwatch.pfeg.noaa.gov/erddap/griddap/nceiErsstv5.nc?sst[(1854-01-01):1:(", year, "-", month, "-01T00:00:00Z)][(0.0):1:(0.0)][(20):1:(70)][(120):1:(250)]", sep="")

download.file(URL, "data/North.Pacific.ersst")

################
# download GODAS
################

# identify latest year and month needed
year <- 2019
month <- "0b"
query <- c("93a6_4907_ab65.nc?dbssbitl",
           "13c5_96ef_e443.nc?dbssbmxl",
           "2ee3_0bfa_a8d6.nc?sshgsfc",
           "5cd4_e4be_d3dd.nc?uflxsfc",
           "5644_bc32_c3fd.nc?vflxsfc")

variable <- c("dbssbitl", "dbssbmxl", "sshgsfc", "uflxsfc", "vflxsfc")

for(i in 1:length(query)){
URL <- paste("http://apdrc.soest.hawaii.edu/erddap/griddap/hawaii_soest_", query[i], "[(1980-01-01):1:(", year, "-",
             month, "-01T00:00:00Z)][(40.05199999999999):1:(64.36099999999999)][(160.5):1:(240.5)]", sep="")

download.file(URL, paste("data/North.Pacific.godas.", variable[i], sep=""))
}

# and test these downloads!
test <- nc_open("data/North.Pacific.godas.sshgsfc")
test

x <- ncvar_get(test, "longitude")
y <- ncvar_get(test, "latitude")
ssh <- ncvar_get(test, "sshgsfc", verbose = F)
# Change data from a 3-D array to a matrix of monthly data by grid point:
# First, reverse order of dimensions ("transpose" array)
ssh <- aperm(ssh, 3:1)  

# Change to matrix with column for each grid point, rows for monthly means
ssh <- matrix(ssh, nrow=dim(ssh)[1], ncol=prod(dim(ssh)[2:3]))  


z <- colMeans(ssh, na.rm=T)
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
#image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image(x,y,z, col=tim.colors(64), xlab = "", ylab = "", yaxt="n", xaxt="n")

contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 

################
# download NCEP/NCAR
################

# identify latest year and month needed
year <- 2019
month <- "05"
query <- c("e025_be03_a4bb.nc?vflx",
           "803d_41d9_b553.nc?uflx",
           "f19d_3925_d70b.nc?slp")

variable <- c("vflx", "uflx", "slp")

for(i in 1:length(query)){
  URL <- paste("http://apdrc.soest.hawaii.edu/erddap/griddap/hawaii_soest_", query[i], "[(1948-01-01):1:(", year, "-",
               month, "-01T00:00:00Z)][(19.99970054626):1:(69.52169799805)][(120):1:(249.375)]", sep="")
  
  download.file(URL, paste("data/North.Pacific.NCEP.NCAR.", variable[i], sep=""))
}

# and test
test <- nc_open("data/North.Pacific.NCEP.NCAR.uflx")
test

x <- ncvar_get(test, "longitude")
y <- ncvar_get(test, "latitude")
z <- ncvar_get(test, "uflx", verbose = F)
# Change data from a 3-D array to a matrix of monthly data by grid point:
# First, reverse order of dimensions ("transpose" array)
z <- aperm(z, 3:1)  

# Change to matrix with column for each grid point, rows for monthly means
z <- matrix(z, nrow=dim(z)[1], ncol=prod(dim(z)[2:3]))  

z <- colMeans(z, na.rm=T)
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
#image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image(x,y,z, col=tim.colors(64), xlab = "", ylab = "", yaxt="n", xaxt="n")

contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires', add=T, lwd=1)

################
# download SODA 2.24
################

URL <- "http://apdrc.soest.hawaii.edu/erddap/griddap/hawaii_soest_f29d_0f00_8b34.nc?taux[(1948-01-15):1:(2010-12-15T00:00:00Z)][(20.25):1:(70.25)][(120.25):1:(250.25)],tauy[(1948-01-15):1:(2010-12-15T00:00:00Z)][(20.25):1:(70.25)][(120.25):1:(250.25)],ssh[(1948-01-15):1:(2010-12-15T00:00:00Z)][(20.25):1:(70.25)][(120.25):1:(250.25)]"
download.file(URL, "data/North.Pacific.SODA.2.2.4")

#####
# ECO-FOCI M2 mooring

URL <- "https://ferret.pmel.noaa.gov/pmel/erddap/tabledap/EcoFOCI_Bering_Sea_timeseries_data.nc?id%2Clongitude%2Clon360%2Clatitude%2Cdepth%2Ctime%2Cocean_velocity_u_eastward_true_north%2Cocean_velocity_v_northward_true_north%2Cocean_current_speed_component%2Cocean_current_direction_component_true_north%2Cocean_temperature_1%2Cocean_depth%2Cocean_chlorophyll_a_concentration_factoryCal%2Cocean_chlorophyll_fluorescence_raw%2Cocean_chlorophyll_fluorescence_raw_standard_deviation%2Cocean_practical_salinity_1&longitude%3E=-164.0663&latitude%3E=56.8633&depth%3E=9.0&time%3E=2010-10-05T03%3A00%3A00Z&time%3C=2017-04-16T11%3A20%3A00Z"
download.file(URL, "~temp")

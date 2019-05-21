library(ncdf4)
library(zoo)
library(gplots)
library(dplyr)
library(maps)
library(mapdata)
library(chron)
library(fields)
library(tidyr)
library(nlme)
library(pracma)
library(ggplot2)
library(car)
library(FactoMineR)
library(ggpubr)
library(mgcv)


# using monthly NCEP/NCAR!
# these data are in the google drive data folder
nc.slp <- nc_open("/Users/MikeLitzow 1/Documents/R/climate scripts/ncep_ncar_monthly_slp.nc")

# now process SLP data - first, extract dates
raw <- ncvar_get(nc.slp, "TIME")  # seconds since 1-1-1970
# h <- raw/(24*60*60)
d <- dates(raw, origin = c(1,1,0001))
yr <- years(d)

x <- ncvar_get(nc.slp, "LON53_101")
y <- ncvar_get(nc.slp, "LAT45_69")

SLP <- ncvar_get(nc.slp, "SLP", verbose = F)
# Change data from a 3-D array to a matrix of monthly data by grid point:
# First, reverse order of dimensions ("transpose" array)
SLP <- aperm(SLP, 3:1)  

# Change to matrix with column for each grid point, rows for monthly means
SLP <- matrix(SLP, nrow=dim(SLP)[1], ncol=prod(dim(SLP)[2:3]))  

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   
lon <- rep(x, each = length(y))   
dimnames(SLP) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# plot to check
z <- colMeans(SLP)   # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image(x,y,z, col=tim.colors(64), xlab = "", ylab = "", yaxt="n", xaxt="n")
contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,80),add=T, lwd=1)

# load pdo
names <- read.table("~pdo", skip=30, nrows=1, as.is = T)
pdo <- read.table("~pdo", skip=31, nrows=119, fill=T, col.names = names)
pdo$YEAR <- 1900:(1899+nrow(pdo)) # drop asterisks!
pdo <- pdo %>%
  gather(month, value, -YEAR) %>%
  arrange(YEAR)

# load npgo
# download.file("http://www.oces.us/npgo/npgo.php", "~npgo")
npgo <- read.table("~npgo", skip=10, nrows=828, fill=T, col.names = c("Year", "month", "value"))


# now, AL sd
# subset SLP to AL area...45-55 N, 192.5-207.5 E
latAL <- lat >= 46.25 & lat <=56.25
lonAL <- lon >= 186.25 & lon <=206

SLP.AL <- SLP
SLP.AL[,!latAL] <- NA
SLP.AL[,!lonAL] <- NA

# check
z <- colMeans(SLP.AL)   # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image(x,y,z, col=tim.colors(64), xlab = "", ylab = "", yaxt="n", xaxt="n")
contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,80),add=T, lwd=1)

# remove seasonal signal
m <- months(d)
yr <- years(d)

f <- function(x) tapply(x, m, mean)
mu <- apply(SLP.AL, 2, f)	# Compute monthly means for each time series (location)

mu <- mu[rep(1:12, floor(length(d)/12)),] 
xtra <- 12*((length(d)/12)-floor(length(d)/12))
mu <- rbind(mu, mu[1:xtra,])

SLP.ALanom <- SLP.AL - mu   # Compute matrix of anomalies - dropping year and month!

# get average anomaly across the area
SLP.ALanom <- rowMeans(SLP.ALanom, na.rm=T)

# smooth with 11-m rolling mean
SLP.sm <- rollmean(SLP.ALanom,11,fill=NA)

dec.yr <- as.numeric(as.character(yr))+(as.numeric(m)-0.5)/12

plot(dec.yr, SLP.sm, type="l")
plot(dec.yr, SLP.sd, type="l")

# put together time series plots
plot.dat <- data.frame(dec.yr=dec.yr, year=as.numeric(as.character(yr)), month=as.numeric(m), AL.sd=SLP.sd, PDO.NPGO.cor=NA, SLP.PDO.NS=NA, SLP.NPGO=NA)

# and get full-field SD values by era to plot

# recalculate anomalies for the entire era
mu <- apply(SLP, 2, f)	# Compute monthly means for each time series (location)

mu <- mu[rep(1:12, floor(length(d)/12)),] 
xtra <- 12*((length(d)/12)-floor(length(d)/12))
mu <- rbind(mu, mu[1:xtra,])

SLPanom <- SLP - mu   # Compute matrix of anomalies - dropping year and month!

ff <- function(x) rollmean(x,11,fill=NA)

# smooth with 11-m rolling mean
SLP.sm <- apply(SLPanom, 2, ff)

SLPsd1 <- apply(SLPanom[yr<=1988,], 2, sd, na.rm=T)
SLPsd2 <- apply(SLPanom[yr %in% 1989:2013,], 2, sd, na.rm=T)
SLPsd.diff <- SLPsd2-SLPsd1



# and attempt some plots that depict the declining independent predictive skill of the NPGO
# reload SST
nc <- nc_open("~updated.sst")

# get lat/long
x.t <- ncvar_get(nc, "longitude")
y.t <- ncvar_get(nc, "latitude")
lat.t <- rep(y.t, length(x.t))   # Vector of latitudes
lon.t <- rep(x.t, each = length(y.t))   # Vector of longitudes

# assign dates
raw <- ncvar_get(nc, "time") # seconds since January 1, 1970
h <- raw/(24*60*60)
d.t <- dates(h, origin = c(1,1,1970))

# year for processing later
m <- as.numeric(months(d.t))
yr <- as.numeric(as.character(years(d.t)))

# get required sst data
SST <- ncvar_get(nc, "sst")

# need to change to matrix for easier use
SST <- aperm(SST, 3:1) # transpose array

SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  # Change to matrix

# plot to check
z <- colMeans(SST)   # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y.t)))  # Convert vector to matrix and transpose for plotting
image(x.t,y.t,z, col=tim.colors(64), xlab = "", ylab = "", yaxt="n", xaxt="n")
contour(x.t,y.t,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)

# set names
dimnames(SST) <- list(as.character(d.t), paste("N", lat.t, "E", lon.t, sep=""))


f <- function(x) tapply(x, m[yr %in% 1951:1980], mean)  # function to compute monthly means for a single time series

mu <- apply(SST[yr %in% 1951:1980,], 2, f)	# Compute monthly means for each time series (location)

mu <- mu[rep(1:12, floor(length(d.t)/12)),] 

xtra <- 12*((length(d.t)/12)-floor(length(d.t)/12))

mu <- rbind(mu, mu[1:xtra,])

SST.anom <- SST-mu

# and remove land
land <- is.na(colMeans(SST.anom))
SST.anom <- SST.anom[,!land]

# separate into 1950:1988 and 1989:2013
SST.anom1 <- SST.anom[yr <= 1988,]
SST.anom2 <- SST.anom[yr %in% 1989:2013,]

# and separate pdo/npgo to same periods
pdo1 <- pdo$value[pdo$YEAR %in% 1950:1988]
pdo2 <- pdo$value[pdo$YEAR %in% 1989:2013]

npgo1 <- npgo$value[npgo$Year %in% 1950:1988]
npgo2 <- npgo$value[npgo$Year %in% 1989:2013]


# look at NPGO regression on SSTa
npgo.sst.regr1 <- npgo.sst.regr2 <- NA

for(i in 1:ncol(SST.anom1)){
  # i <- 1
  mod <- lm(SST.anom1[,i] ~ npgo1)
  npgo.sst.regr1[i] <- summary(mod)$coefficients[2,1]
  
  mod <- lm(SST.anom2[,i] ~ npgo2)
  npgo.sst.regr2[i] <- summary(mod)$coefficients[2,1]
}

# plot

lim <- range(npgo.sst.regr1, npgo.sst.regr2)

par(mfrow=c(2,2), mar=c(0.5, 0.5,2,2))
z <- rep(NA, ncol(SST))

z[!land] <- npgo.sst.regr1

z <- t(matrix(z, length(y.t)))  # Convert vector to matrix and transpose for plotting
image.plot(x.t,y.t,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x.t,y.t,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("NPGO regr. SSTa 1950-1988", cex=0.8)

z <- rep(NA, ncol(SST))
z[!land] <- npgo.sst.regr2
z <- t(matrix(z, length(y.t)))  # Convert vector to matrix and transpose for plotting
image.plot(x.t,y.t,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))


contour(x.t,y.t,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("NPGO regr. EOF2 1989-2102", cex=0.8)

z <- rep(NA, ncol(SST))
# z[!land] <- npgo.pattern1  # replace elements NOT corresponding to land with loadings!
# z[!land] <- npgo.eof.r1  # replace elements NOT corresponding to land with loadings!
z[!land] <- npgo.eof2.regr1  # replace elements NOT corresponding to land with loadings!

z <- t(matrix(z, length(y.t)))  # Convert vector to matrix and transpose for plotting
image.plot(x.t,y.t,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))


contour(x.t,y.t,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("NPGO regr. EOF2 1950-1988", cex=0.8)

z <- rep(NA, ncol(SST))
# z[!land] <- npgo.pattern2  # replace elements NOT corresponding to land with loadings!
# z[!land] <- npgo.eof.r2
z[!land] <- npgo.eof2.regr2
z <- t(matrix(z, length(y.t)))  # Convert vector to matrix and transpose for plotting
image.plot(x.t,y.t,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))


contour(x.t,y.t,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("NPGO regr. EOF2 1989-2102", cex=0.8)

###
# and SLP-PDO
# reload SLP to make things easy...
nc.slp <- nc_open("/Users/MikeLitzow 1/Documents/R/climate scripts/64AA9803F7345DFE8991A731014FFB01_ferret_listing.nc")

# now process SLP data - first, extract dates
raw <- ncvar_get(nc.slp, "TIME")  # seconds since 1-1-1970
# h <- raw/(24*60*60)
d <- dates(raw, origin = c(1,1,0001))
yr <- years(d)


x <- ncvar_get(nc.slp, "LON53_101")
y <- ncvar_get(nc.slp, "LAT45_69")

# save to plot below
x.slp <- x
y.slp <- y

SLP <- ncvar_get(nc.slp, "SLP", verbose = F)
# Change data from a 3-D array to a matrix of monthly data by grid point:
# First, reverse order of dimensions ("transpose" array)
SLP <- aperm(SLP, 3:1)  

# Change to matrix with column for each grid point, rows for monthly means
SLP <- matrix(SLP, nrow=dim(SLP)[1], ncol=prod(dim(SLP)[2:3]))  

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   
lon <- rep(x, each = length(y))   
dimnames(SLP) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# plot to check
z <- colMeans(SLP)   # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image(x,y,z, col=tim.colors(64), xlab = "", ylab = "", yaxt="n", xaxt="n")
contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,80),add=T, lwd=1)

# load pdo
names <- read.table("~pdo", skip=30, nrows=1, as.is = T)
pdo <- read.table("~pdo", skip=31, nrows=119, fill=T, col.names = names)
pdo$YEAR <- 1900:(1899+nrow(pdo)) # drop asterisks!
pdo <- pdo %>%
  gather(month, value, -YEAR) %>%
  arrange(YEAR)

# load npgo
download.file("http://www.oces.us/npgo/npgo.php", "~npgo")
npgo <- read.table("~npgo", skip=10, nrows=828, fill=T, col.names = c("Year", "month", "value"))

# smooth SLP with three month rolling mean
f <- function(x) rollmean(x,3,fill=NA)

SLPsm <- apply(SLP, 2, f)

rownames(SLPsm)

# limit to 1950:2017
SLPsm <- SLPsm[yr %in% 1950:2017,]

# limit PDO and NPGO to March 1950 - Feb 2018
# also expanding to split out 1950-1988 and 1989-2013
pdo[c(603,1418),]
pdo[c(603,1068),]
pdo[c(1069,1368),]

pdoTS <- pdo$value[603:1418]
pdoTS1 <- pdo$value[603:1068]
pdoTS2 <- pdo$value[1069:1368]

npgo[c(3,818),]
npgo[c(3,468),]
npgo[c(469,768),]

npgoTS <- npgo$value[3:818]
npgoTS1 <- npgo$value[3:468]
npgoTS2 <- npgo$value[469:768]

# and two era time series
rownames(SLPsm)[c(1,length(pdoTS1))]
rownames(SLPsm)[c(length(pdoTS1),length(pdoTS1)+length(pdoTS2))]
SLPsm1 <- SLPsm[1:length(pdoTS1),]
SLPsm2 <- SLPsm[(1+length(pdoTS1)):(length(pdoTS1)+length(pdoTS2)),]

# separate regressions in each era!
pdo.regr1 <- pdo.regr2 <- npgo.regr1 <- npgo.regr2 <- NA

for(i in 1:ncol(SLPsm)){
  #  i <- 1
  mod <- lm(SLPsm1[,i] ~ pdoTS1)
  pdo.regr1[i] <- summary(mod)$coef[2,1]
  
  mod <- lm(SLPsm2[,i] ~ pdoTS2)
  pdo.regr2[i] <- summary(mod)$coef[2,1]
  
  mod <- lm(SLPsm1[,i] ~ npgoTS1)
  npgo.regr1[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(SLPsm2[,i] ~ npgoTS2)
  npgo.regr2[i] <- summary(mod)$coef[2,1] 
}

# and plot
# set up color schemes
new.col <- my.col <- tim.colors(64)
grays <- c("gray98", "gray97", "gray96", "gray95", "gray94", "gray93", "gray92", "gray91", "gray90", "gray89", "gray88",
           "gray87", "gray86", "gray85", "gray84", "gray83", "gray82", "gray81", "gray80", "gray79", "gray78", "gray77",
           "gray76", "gray75", "gray74", "gray73", "gray72", "gray71", "gray71", "gray71", "gray70", "gray69", "gray68")


my.col[1:33] <- grays
my.col[22:43] <- c(grays[11:1], grays[1:11])

new.col[27:36] <- c(grays[5:1], grays[1:5])

lim <- range(pdo.regr1, pdo.regr2, npgo.regr1, npgo.regr2)


# setup the layout
mt.cex <- 1.1
l.mar <- 3
l.cex <- 0.8
l.l <- 0.2
tc.l <- -0.2

par(mar=c(1.5,1.5,1.5,1),  tcl=tc.l, mgp=c(1.5,0.3,0), las=1, mfrow=c(2,2), cex.axis=0.8, cex.lab=0.8, oma=c(0,0,0,0.2))

# PDO first era
z <- pdo.regr1  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
#image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
      xlab = "", ylab = "")

contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("PDO forcing pattern 1950-1988", cex=0.8)

xx1 <- c(183.5, 183.5, 204, 204, 183.5)
yy1 <- c(46.5, 54, 54, 46.5, 46.5)

xx2 <- c(196, 196, 216.5, 216.5, 196)
yy2 <- c(39, 46.5, 46.5, 39, 39)

lines(xx1, yy1, lwd=1.5, col="magenta")
lines(xx2, yy2, lwd=1.5, col="magenta")

# PDO second era
z <- pdo.regr2  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
#image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
      xlab = "", ylab = "")

contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("PDO forcing pattern 1989-2013", cex=0.8)

lines(xx1, yy1, lwd=1.5, col="magenta")
lines(xx2, yy2, lwd=1.5, col="magenta")


# npgo first era
z <- npgo.regr1  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
# image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#            xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
      xlab = "", ylab = "")

contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("NPGO forcing pattern 1950-1988", cex=0.8)

xx3 <- c(194, 194, 199, 211, 216.5, 216.5, 194)
yy3 <- c(51.5, 54, 56.5, 61.5, 61.5, 51.5, 51.5)

lines(xx3, yy3, lwd=1.5, col="magenta")

# npgo second era
z <- npgo.regr2  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
# image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#            xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
      xlab = "", ylab = "")

contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("NPGO forcing pattern 1989-2013", cex=0.8)

xx3 <- c(194, 194, 199, 211, 216.5, 216.5, 194)
yy3 <- c(51.5, 54, 56.5, 61.5, 61.5, 51.5, 51.5)

lines(xx3, yy3, lwd=1.5, col="magenta")

###########################
# get rolling regression coefficients for each square

PDOkeep1 <- PDOkeep2 <- NPGOkeep <- NA
pdoSLP1 <- pdoSLP2 <- npgoSLP <- SLPsm

for(i in 1:length(lat)){
  # i <- 1
  PDOkeep1[i] <- inpolygon(lon[i], lat[i], xx1, yy1)
  PDOkeep2[i] <- inpolygon(lon[i], lat[i], xx2, yy2)
  NPGOkeep[i] <- inpolygon(lon[i], lat[i], xx3, yy3)
}

pdoSLP1[,!PDOkeep1] <- NA
pdoSLP2[,!PDOkeep2] <- NA
npgoSLP[,!NPGOkeep] <- NA


# plot to check
# first PDO1
z <- colMeans(pdoSLP1, na.rm=T)  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
#image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image(x,y,z, col=tim.colors(64), ylim=c(20,68),
      xlab = "", ylab = "")

contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
lines(xx1, yy1, lwd=1.5, col="magenta")

# PDO2 
z <- colMeans(pdoSLP2, na.rm=T)  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
#image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image(x,y,z, col=tim.colors(64), ylim=c(20,68),
      xlab = "", ylab = "")

contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
lines(xx2, yy2, lwd=1.5, col="magenta")

# now NPGO
z <- colMeans(npgoSLP, na.rm=T)  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
#image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image(x,y,z, col=tim.colors(64), ylim=c(20,68),
      xlab = "", ylab = "")

contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
lines(xx3, yy3, lwd=1.5, col="magenta")

# all look good

# now try 21-yr rolling regressions!
pdo.regr1 <- pdo.regr2 <- npgo.regr <- NA

# get monthly means of each
pdoSLP1 <- rowMeans(pdoSLP1, na.rm=T)
pdoSLP2 <- rowMeans(pdoSLP2, na.rm=T)
npgoSLP <- rowMeans(npgoSLP, na.rm=T)

names(pdoSLP1)[1:(12*21)]
# will use 253 month (21 yr + 1 mo) rolling windows

roll <- data.frame(dec.yr=as.numeric(as.character(years(names(pdoSLP1))))+(as.numeric(months(names(pdoSLP1)))-0.5)/12, pdo.regr1=NA, pdo.regr2=NA, npgo.regr=NA)

# and the data for regressions
dat <- data.frame(year=as.numeric(as.character(years(names(pdoSLP1))))+(as.numeric(months(names(pdoSLP1)))-0.5)/12,
                  pdoSLP1=pdoSLP1, pdoSLP2=pdoSLP2, npgoSLP=npgoSLP, pdoTS=pdoTS, npgoTS=npgoTS, pdo.regr1=NA, pdo.regr2=NA, npgo.regr=NA) # so the correct lags are built in


for(i in 127:(nrow(dat)-126)){
  # i <- 127
  temp <- dat[(i-126):(i+126),]
  
  mod <- lm(temp$pdoSLP1 ~ temp$pdoTS)
  dat$pdo.regr1[i] <- summary(mod)$coef[2,1]
  
  mod <- lm(temp$pdoSLP2 ~ temp$pdoTS)
  dat$pdo.regr2[i] <- summary(mod)$coef[2,1]
  
  mod <- lm(temp$npgoSLP ~ temp$npgoTS)
  dat$npgo.regr[i] <- summary(mod)$coef[2,1]
  
}

dat$pdoNS <- dat$pdo.regr1/dat$pdo.regr2
  
plot.dat$SLP.PDO.NS  <- dat$pdoNS[match(plot.dat$dec.yr, dat$year)]
plot.dat$SLP.NPGO  <- dat$npgo.regr[match(plot.dat$dec.yr, dat$year)]

plot.dat <- dat %>%
  select(year, pdo.regr1, pdo.regr2, npgo.regr) %>%
  gather(key, value, -year)

ggplot(plot.dat, aes(year, value, color=key)) +
  theme_linedraw() +
  geom_line() +
  ylim(c(-0.5, -2.8)) + 
  xlim(c(1962,2005))

dat$pdo.ratio <- ifelse(is.na(dat$pdo.regr1),NA, dat$pdo.regr1/dat$pdo.regr2)

plot.dat <- dat %>%
  select(year, pdo.ratio, npgo.regr) %>%
  gather(key, value, -year)

ggplot(plot.dat, aes(year, value)) +
  theme_linedraw() +
  geom_line() +
  facet_wrap(~key, scales="free") + 
  xlim(c(1962,2005))







#######
# combined plot
png("reduced basin scale combined plot.png", 6,6, units="in", res=300)

# setup the layout
mt.cex <- 1.1
l.mar <- 3
l.cex <- 0.8
l.l <- 0.2
tc.l <- -0.2

cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
my.col <- tim.colors(64)
grays <- c("gray98", "gray98", "gray97", "gray97", "gray96", "gray96", "gray95", "gray95", "gray94", "gray94", "gray93",
           "gray93", "gray92", "gray92", "gray91", "gray91", "gray90", "gray90", "gray89", "gray89", "gray88", "gray88",
           "gray87", "gray87", "gray86", "gray86", "gray85", "gray85", "gray84", "gray84", "gray83", "gray83", "gray82")


my.col[1:33] <- grays
# my.col[22:43] <- c(grays[11:1], grays)
# new.col[27:36] <- c(grays[5:1], grays[1:5])

xlim <- c(min(dec.yr), max(dec.yr.t))

par(mar=c(1.25,1.25,1.25,1),  tcl=tc.l, mgp=c(1.5,0.3,0), las=1, mfrow=c(4,3), cex.axis=0.8, cex.lab=0.8, oma=c(0,2,2,0))


###
# now AL SD

lim <- range(SLPsd1, SLPsd2)

xx <- c(186.25, 186.25, 206, 206, 186.25)
yy <- c(46.25, 56.25, 56.25, 46.25, 46.25)


z <- SLPsd1   # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image.plot(x,y,z, col=my.col, xlab = "", ylab = "", zlim=lim, ylim=c(20,70), yaxt="n", xaxt="n",
           legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,80),add=T, lwd=1)
lines(xx,yy, lwd=2, col="magenta")
mtext("a", adj=0.05, line=-1.4, cex=1)

z <- SLPsd2   # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image.plot(x,y,z, col=my.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=lim, ylim=c(20,70), yaxt="n", xaxt="n", 
           legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,80),add=T, lwd=1)
lines(xx,yy, lwd=2, col="magenta")
mtext("b", adj=0.05, line=-1.4, cex=1)

par(mar=c(1.5,3,1.5,1))

plot(dec.yr, SLP.sd, type="l", xlab="", ylab="Standard dev. (Pa)", xlim=xlim, col=cb[6])
abline(h=mean(SLP.sd, na.rm=T))
abline(v=1989, lty=2)
mtext("c", adj=0.05, line=-1.4, cex=1)

par(mar=c(1.25,1.25,1.25,1))

# now atmospheric forcing of PDO/NPGO
lim <- range(pdo.regr1, pdo.regr2, npgo.regr1, npgo.regr2)
# PDO first era
z <- pdo.regr1  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y.slp)))  # Convert vector to matrix and transpose for plotting
#image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image.plot(x.slp,y.slp,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
           xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x.slp,y.slp,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("d", adj=0.05, line=-1.4, cex=1)

xx1 <- c(183.5, 183.5, 204, 204, 183.5)
yy1 <- c(46.5, 54, 54, 46.5, 46.5)

xx2 <- c(196, 196, 216.5, 216.5, 196)
yy2 <- c(39, 46.5, 46.5, 39, 39)

lines(xx1, yy1, lwd=1.5, col="magenta")
lines(xx2, yy2, lwd=1.5, col="magenta")

# PDO second era
z <- pdo.regr2  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y.slp)))  # Convert vector to matrix and transpose for plotting
#image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image.plot(x.slp,y.slp,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
           xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x.slp,y.slp,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("e", adj=0.05, line=-1.4, cex=1)

lines(xx1, yy1, lwd=1.5, col="magenta")
lines(xx2, yy2, lwd=1.5, col="magenta")

par(mar=c(1.5,3,1.5,1))

plot(dat$year, dat$pdoNS, type="l", xlab="", ylab="North:South forcing ratio", xlim=xlim, col=cb[6])
abline(h=mean(dat$pdoNS, na.rm=T))
abline(v=1989, lty=2)
mtext("f", adj=0.05, line=-1.4, cex=1)

par(mar=c(1.25,1.25,1.25,1))

# npgo first era
z <- npgo.regr1  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y.slp)))  # Convert vector to matrix and transpose for plotting
# image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#            xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image.plot(x.slp,y.slp,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
           xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x.slp,y.slp,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("g", adj=0.05, line=-1.4, cex=1)

xx3 <- c(194, 194, 199, 211, 216.5, 216.5, 194)
yy3 <- c(51.5, 54, 56.5, 61.5, 61.5, 51.5, 51.5)

lines(xx3, yy3, lwd=1.5, col="magenta")

# npgo second era
z <- npgo.regr2  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y.slp)))  # Convert vector to matrix and transpose for plotting
# image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#            xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image.plot(x.slp,y.slp,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
           xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x.slp,y.slp,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("h", adj=0.05, line=-1.4, cex=1)

xx3 <- c(194, 194, 199, 211, 216.5, 216.5, 194)
yy3 <- c(51.5, 54, 56.5, 61.5, 61.5, 51.5, 51.5)

lines(xx3, yy3, lwd=1.5, col="magenta")

par(mar=c(1.5,3,1.5,1))

plot(dat$year, dat$npgo.regr, type="l", xlab="", ylab="Regression coef (Pa)", xlim=xlim, col=cb[6],
     ylim=c(max(dat$npgo.regr, na.rm=T), min(dat$npgo.regr, na.rm=T)))
abline(h=mean(dat$npgo.regr, na.rm=T))
abline(v=1989, lty=2)
mtext("i", adj=0.05, line=-1.4, cex=1)

par(mar=c(1.25,1.25,1.25,1))


# finally, the NPGO-SST plots
lim <- range(npgo.sst.regr1, npgo.sst.regr2)
z <- rep(NA, ncol(SST))
# z[!land] <- npgo.pattern1  # replace elements NOT corresponding to land with loadings!
# z[!land] <- npgo.eof.r1  # replace elements NOT corresponding to land with loadings!
z[!land] <- npgo.sst.regr1  # replace elements NOT corresponding to land with loadings!

z <- t(matrix(z, length(y.t)))  # Convert vector to matrix and transpose for plotting
image.plot(x.t,y.t,z, col=new.col, zlim=c(lim[1], -lim[1]),
           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))


contour(x.t,y.t,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("j", adj=0.05, line=-1.4, cex=1)

z <- rep(NA, ncol(SST))
# z[!land] <- npgo.pattern2  # replace elements NOT corresponding to land with loadings!
# z[!land] <- npgo.eof.r2
z[!land] <- npgo.sst.regr2
z <- t(matrix(z, length(y.t)))  # Convert vector to matrix and transpose for plotting
image.plot(x.t,y.t,z, col=new.col, zlim=c(lim[1], -lim[1]),
           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))


contour(x.t,y.t,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("k", adj=0.05, line=-1.4, cex=1)

par(mar=c(1.5,3,1.5,1))

plot(plot.dat$dec.yr, plot.dat$PDO.NPGO.cor, type="l", xlab="", ylab="PDO-NPGO correlation", xlim=xlim, col=cb[6])
abline(h=0)
abline(v=1989, lty=2)
mtext("l", adj=0.05, line=-1.4, cex=1)
dev.off()


########
# older version
png("older larger basin scale combined plot.png", 7,10, units="in", res=300)

# setup the layout
mt.cex <- 1.1
l.mar <- 3
l.cex <- 0.8
l.l <- 0.2
tc.l <- -0.2

cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
 my.col <- tim.colors(64)
grays <- c("gray98", "gray98", "gray97", "gray97", "gray96", "gray96", "gray95", "gray95", "gray94", "gray94", "gray93",
           "gray93", "gray92", "gray92", "gray91", "gray91", "gray90", "gray90", "gray89", "gray89", "gray88", "gray88",
           "gray87", "gray87", "gray86", "gray86", "gray85", "gray85", "gray84", "gray84", "gray83", "gray83", "gray82")


my.col[1:33] <- grays
# my.col[22:43] <- c(grays[11:1], grays)
# new.col[27:36] <- c(grays[5:1], grays[1:5])

xlim <- c(min(dec.yr), max(dec.yr.t))

par(mar=c(1.5,1.5,1.5,1),  tcl=tc.l, mgp=c(1.5,0.3,0), las=1, mfrow=c(5,3), cex.axis=0.8, cex.lab=0.8, oma=c(0,2,2,0))

lim <- range(SSTanom1, SSTanom2, na.rm=T)

z <- SSTanom1   # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y.t)))  # Convert vector to matrix and transpose for plotting
image.plot(x.t,y.t,z, col=tim.colors(64), xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim[2], lim[2]), 
           legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x.t,y.t,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,80),add=T, lwd=1)
mtext("a", adj=0.05, line=-1.4, cex=1.1)

z <- SSTanom2   # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y.t)))  # Convert vector to matrix and transpose for plotting
image.plot(x.t,y.t,z, col=tim.colors(64), xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim[2], lim[2]),
           legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x.t,y.t,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,80),add=T, lwd=1)
mtext("b", adj=0.05, line=-1.4, cex=1.1)

par(mar=c(1.5,3,1.5,1))
plot(dec.yr.t, SST.anomTS, type="l", xlab="", ylab="ºC wrt 1951-1980", xlim=xlim, col=cb[6])
abline(h=0)
abline(v=1989, lty=2)
mtext("c", adj=0.05, line=-1.4, cex=1.1)

par(mar=c(1.5,1.5,1.5,1))

###
# now AL SD

lim <- range(SLPsd1, SLPsd2)

xx <- c(186.25, 186.25, 206, 206, 186.25)
yy <- c(46.25, 56.25, 56.25, 46.25, 46.25)


z <- SLPsd1   # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image.plot(x,y,z, col=my.col, xlab = "", ylab = "", zlim=lim, ylim=c(20,70), yaxt="n", xaxt="n",
           legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,80),add=T, lwd=1)
lines(xx,yy, lwd=2, col="magenta")
mtext("d", adj=0.05, line=-1.4, cex=1.1)

z <- SLPsd2   # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image.plot(x,y,z, col=my.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=lim, ylim=c(20,70), yaxt="n", xaxt="n", 
           legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,80),add=T, lwd=1)
lines(xx,yy, lwd=2, col="magenta")
mtext("e", adj=0.05, line=-1.4, cex=1.1)

par(mar=c(1.5,3,1.5,1))

plot(dec.yr, SLP.sd, type="l", xlab="", ylab="Standard dev. (Pa)", xlim=xlim, col=cb[6])
abline(h=mean(SLP.sd, na.rm=T))
abline(v=1989, lty=2)
mtext("f", adj=0.05, line=-1.4, cex=1.1)

par(mar=c(1.5,1.5,1.5,1))

# now atmospheric forcing of PDO/NPGO
lim <- range(pdo.regr1, pdo.regr2, npgo.regr1, npgo.regr2)
# PDO first era
z <- pdo.regr1  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y.slp)))  # Convert vector to matrix and transpose for plotting
#image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image.plot(x.slp,y.slp,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
      xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x.slp,y.slp,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("g", adj=0.05, line=-1.4, cex=1.1)

xx1 <- c(183.5, 183.5, 204, 204, 183.5)
yy1 <- c(46.5, 54, 54, 46.5, 46.5)

xx2 <- c(196, 196, 216.5, 216.5, 196)
yy2 <- c(39, 46.5, 46.5, 39, 39)

lines(xx1, yy1, lwd=1.5, col="magenta")
lines(xx2, yy2, lwd=1.5, col="magenta")

# PDO second era
z <- pdo.regr2  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y.slp)))  # Convert vector to matrix and transpose for plotting
#image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image.plot(x.slp,y.slp,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
           xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x.slp,y.slp,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("h", adj=0.05, line=-1.4, cex=1.1)

lines(xx1, yy1, lwd=1.5, col="magenta")
lines(xx2, yy2, lwd=1.5, col="magenta")

par(mar=c(1.5,3,1.5,1))

plot(dat$year, dat$pdoNS, type="l", xlab="", ylab="North:South forcing ratio", xlim=xlim, col=cb[6])
abline(h=mean(dat$pdoNS, na.rm=T))
abline(v=1989, lty=2)
mtext("i", adj=0.05, line=-1.4, cex=1.1)

par(mar=c(1.5,1.5,1.5,1))

# npgo first era
z <- npgo.regr1  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y.slp)))  # Convert vector to matrix and transpose for plotting
# image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#            xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image.plot(x.slp,y.slp,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
           xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x.slp,y.slp,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("j", adj=0.05, line=-1.4, cex=1.1)

xx3 <- c(194, 194, 199, 211, 216.5, 216.5, 194)
yy3 <- c(51.5, 54, 56.5, 61.5, 61.5, 51.5, 51.5)

lines(xx3, yy3, lwd=1.5, col="magenta")

# npgo second era
z <- npgo.regr2  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y.slp)))  # Convert vector to matrix and transpose for plotting
# image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#            xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image.plot(x.slp,y.slp,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
           xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x.slp,y.slp,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("k", adj=0.05, line=-1.4, cex=1.1)

xx3 <- c(194, 194, 199, 211, 216.5, 216.5, 194)
yy3 <- c(51.5, 54, 56.5, 61.5, 61.5, 51.5, 51.5)

lines(xx3, yy3, lwd=1.5, col="magenta")

par(mar=c(1.5,3,1.5,1))

plot(dat$year, dat$npgo.regr, type="l", xlab="", ylab="Regression coef (Pa)", xlim=xlim, col=cb[6],
     ylim=c(max(dat$npgo.regr, na.rm=T), min(dat$npgo.regr, na.rm=T)))
abline(h=mean(dat$npgo.regr, na.rm=T))
abline(v=1989, lty=2)
mtext("l", adj=0.05, line=-1.4, cex=1.1)

par(mar=c(1.5,1.5,1.5,1))


# finally, the NPGO-EOF2 plots
lim <- range(npgo.eof2.regr1, npgo.eof2.regr2)
z <- rep(NA, ncol(SST))
# z[!land] <- npgo.pattern1  # replace elements NOT corresponding to land with loadings!
# z[!land] <- npgo.eof.r1  # replace elements NOT corresponding to land with loadings!
z[!land] <- npgo.eof2.regr1  # replace elements NOT corresponding to land with loadings!

z <- t(matrix(z, length(y.t)))  # Convert vector to matrix and transpose for plotting
image.plot(x.t,y.t,z, col=new.col, zlim=c(lim[1], -lim[1]),
           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))


contour(x.t,y.t,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("m", adj=0.05, line=-1.4, cex=1.1)

z <- rep(NA, ncol(SST))
# z[!land] <- npgo.pattern2  # replace elements NOT corresponding to land with loadings!
# z[!land] <- npgo.eof.r2
z[!land] <- npgo.eof2.regr2
z <- t(matrix(z, length(y.t)))  # Convert vector to matrix and transpose for plotting
image.plot(x.t,y.t,z, col=new.col, zlim=c(lim[1], -lim[1]),
           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))


contour(x.t,y.t,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("n", adj=0.05, line=-1.4, cex=1.1)

par(mar=c(1.5,3,1.5,1))

plot(plot.dat$dec.yr, plot.dat$PDO.NPGO.cor, type="l", xlab="", ylab="PDO-NPGO correlation", xlim=xlim, col=cb[6])
abline(h=0)
abline(v=1989, lty=2)
mtext("o", adj=0.05, line=-1.4, cex=1.1)
dev.off()


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

# download ncep/ncar slp and wind data for N. Pacific: 20-67.5 deg. N, 132.5-250 deg. E
# breaking into individual files for each variable as the files are so large!
download.file("https://coastwatch.pfeg.noaa.gov/erddap/griddap/esrlNcepRe.nc?slp[(1948-01-01):1:(2018-10-15T00:00:00Z)][(67.5):1:(20)][(132.5):1:(250)]", "~slptemp")
download.file("https://coastwatch.pfeg.noaa.gov/erddap/griddap/esrlNcepRe.nc?uwnd[(1948-01-01):1:(2018-10-15T00:00:00Z)][(67.5):1:(20)][(132.5):1:(250)]", "~uwindtemp")
download.file("https://coastwatch.pfeg.noaa.gov/erddap/griddap/esrlNcepRe.nc?vwnd[(1948-01-01):1:(2018-10-15T00:00:00Z)][(67.5):1:(20)][(132.5):1:(250)]", "~vwindtemp")
nc.slp <- nc_open("~slptemp")
nc.uwnd <- nc_open("~uwindtemp")
nc.vwnd <- nc_open("~vwindtemp")

# load pdo and npgo
download.file("http://jisao.washington.edu/pdo/PDO.latest", "~pdo")
names <- read.table("~pdo", skip=32, nrows=1, as.is = T)
pdo <- read.table("~pdo", skip=33, nrows=119, fill=T, col.names = names)
pdo$YEAR <- 1900:(1899+nrow(pdo)) # drop asterisks!
pdo <- pdo %>%
  gather(month, value, -YEAR) %>%
  arrange(YEAR)


download.file("http://www.oces.us/npgo/npgo.php", "~npgo")
npgo <- read.table("~npgo", skip=10, nrows=818, fill=T, col.names = c("Year", "month", "value"))

# start with regressions on SLP
# process SLP data - first, extract dates
raw <- ncvar_get(nc.slp, "time")  # seconds since 1-1-1970
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))

# extract study area
# 54-62 deg. N, 200-226 deg. E
x <- ncvar_get(nc.slp, "longitude")
y <- ncvar_get(nc.slp, "latitude")

SLP <- ncvar_get(nc.slp, "slp", verbose = F)
# Change data from a 3-D array to a matrix of monthly data by grid point:
# First, reverse order of dimensions ("transpose" array)
SLP <- aperm(SLP, 3:1)  

# Reverse order of latitudes to be increasing for convenience (in later plotting)
SLP <- SLP[,20:1,]  
# Also reverse corresponding vector of latitudes
y <- rev(y)  
# Change to matrix with column for each grid point, rows for monthly means
SLP <- matrix(SLP, nrow=dim(SLP)[1], ncol=prod(dim(SLP)[2:3]))  

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   
lon <- rep(x, each = length(y))   
dimnames(SLP) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# now we need to get monthly means!
m.y <- paste(years(d), as.numeric(months(d)), sep="-") # make a month-year factor from the dates

f <- function(x) tapply(x, m.y, mean)
SLP.m <- as.data.frame(apply(SLP, 2, f))

vv <- matrix(unlist(strsplit(as.character(rownames(SLP.m)), "-")),ncol=2, byrow = T)

SLP.m$year <- as.numeric(vv[,1])
SLP.m$month <- as.numeric(vv[,2])

SLP.m <- SLP.m %>%
  arrange(month) %>%
  arrange(year)

# remove seasonal signal
# and set up vector of winter years (identify winters by the year corresponding to Jan.)
m <- SLP.m$month
yr <- SLP.m$year
win.yr <- yr
win.yr[m %in% c(11, 12)] <- win.yr[m %in% c(11, 12)] +1

f <- function(x) tapply(x, m, mean)
mu <- apply(SLP.m, 2, f)	# Compute monthly means for each time series (location)

mu <- mu[rep(1:12, floor(length(d)/12)),] 
xtra <- 12*((length(d)/12)-floor(length(d)/12))
mu <- rbind(mu, mu[1:xtra,])

SLP.anom <- SLP.m[,1:960] - mu   # Compute matrix of anomalies - dropping year and month!

# restrict to relevant months
p.win <- c(11,12,1) # months for SLP data
SLP.anom <- SLP.anom[SLP.m$month %in% p.win,]
win.yr <- win.yr[SLP.m$month %in% p.win]

# clean up
rownames(SLP.anom) <- 1:nrow(SLP.anom)

# restrict PDO to relevant months
t.win <- c("FEB", "MAR", "APR")
pdo <- pdo[pdo$month %in% t.win,]
rownames(pdo) <- 1:nrow(pdo)

r1 <- r2 <- r3 <- NA # vectors to catch regression coefficients


pdo.FMA <- tapply(pdo$value, pdo$YEAR, mean) # mean values for winter year corresponding to Jan.

ff <- function(x) tapply(x, win.yr, mean)

SLP.NDJ <- apply(SLP.anom, 2, ff) # mean values for winter year corresponding to Jan. Note that 1949 is incomplete and will not be used!
# (1 mo and 2 mo, respectively!)

for(j in 1:ncol(SLP.anom)){
 # j <- 1
  
  # subset the data for only the cell of interest and set up the early and late eras (pre/post 1988/89)
  temp <- data.frame(slp=SLP.NDJ[2:71, j], pdo=pdo.FMA[50:119], era=c(rep(1, 40), rep(2, 25), rep(3,5)))

  # now fit era-specific regressions to plot coefficients
  mod <- lm(slp ~ pdo, data=temp[temp$era==1,])
  r1[j] <- summary(mod)$coefficients[2,1]
  
  mod <- lm(slp ~ pdo, data=temp[temp$era==2,])
  r2[j] <- summary(mod)$coefficients[2,1]
  
  mod <- lm(slp ~ pdo, data=temp[temp$era==3,])
  r3[j] <- summary(mod)$coefficients[2,1]
  
}

# convert Pa to hPa
r1 <- r1/100
r2 <- r2/100
r3 <- r3/100

# now the same regression approach for the NPGO
# restrict NPGO to relevant months
t.win <- c(2,3,4)
npgo <- npgo[npgo$month %in% t.win,]
rownames(npgo) <- 1:nrow(npgo)

n1 <- n2 <- n3 <- NA # vectors to catch regression coefficients


npgo.FMA <- tapply(npgo$value, npgo$Year, mean) # mean values for winter year corresponding to Jan.

for(j in 1:ncol(SLP.anom)){
  #j <- 1
  
  # subset the data for only the cell of interest and set up the early and late eras (pre/post 1988/89)
  temp <- data.frame(slp=SLP.NDJ[3:71, j], npgo=npgo.FMA, era=c(rep(1, 39), rep(2, 25), rep(3,5)))
  
  # now fit era-specific regressions to plot coefficients
  mod <- lm(slp ~ npgo, data=temp[temp$era==1,])
  n1[j] <- summary(mod)$coefficients[2,1]
  
  mod <- lm(slp ~ npgo, data=temp[temp$era==2,])
  n2[j] <- summary(mod)$coefficients[2,1]
  
  mod <- lm(slp ~ npgo, data=temp[temp$era==3,])
  n3[j] <- summary(mod)$coefficients[2,1]
  
}

# convert Pa to hPa
n1 <- n1/100
n2 <- n2/100
n3 <- n3/100

# set up color schemes
new.col <- my.col <- tim.colors(64)
grays <- c("gray98", "gray97", "gray96", "gray95", "gray94", "gray93", "gray92", "gray91", "gray90", "gray89", "gray88")

my.col[22:43] <- c(grays[11:1], grays)
new.col[27:36] <- c(grays[5:1], grays[1:5])

png("SLP-index regression by era.png", 8, 5, units="in", res=300)

# setup the layout
mt.cex <- 1.1
l.mar <- 3
l.cex <- 0.8
l.l <- 0.2
tc.l <- -0.2

par(mar=c(0.5,0.5,1.5,1),  tcl=tc.l, mgp=c(1.5,0.3,0), las=1, mfrow=c(2,3), cex.axis=0.8, cex.lab=0.8, oma=c(0,0,0,0.2))

# set the limit for plotting 
lim <- range(r1, r2, r3)

z <- r1   
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting

image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("SLP-PDO 1949-1988", cex=0.8)

z <- r2  
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting

image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("SLP-PDO 1989-2013", cex=0.8)

z <- r3  
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting

image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("SLP-PDO 2014-2018", cex=0.8)

# and npgo plots
# using NEW limits for z!

lim <- range(n1, n2, n3)
z <- n1   
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting

image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("SLP-NPGO 1950-1988", cex=0.8)

z <- n2  
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting

image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("SLP-NPGO 1989-2013", cex=0.8)

z <- n3  
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting

image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("SLP-NPGO 2014-2018", cex=0.8)

dev.off()

#####################
# now u-wind!

u.w <- ncvar_get(nc.uwnd, "uwnd", verbose = F)
# Change data from a 3-D array to a matrix of monthly data by grid point:
# First, reverse order of dimensions ("transpose" array)
UW <- aperm(u.w, 3:1)  

# Reverse order of latitudes to be increasing for convenience (in later plotting)
UW <- UW[,20:1,]  

# Change to matrix with column for each grid point, rows for monthly means
UW <- matrix(UW, nrow=dim(UW)[1], ncol=prod(dim(UW)[2:3]))  

# Keep track of corresponding latitudes and longitudes of each column:
dimnames(UW) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# convert to daily wind stress!
UW <- 0.0013*1.22*UW^2
# now we need to get monthly means!
m.y <- paste(years(d), as.numeric(months(d)), sep="-") # make a month-year factor from the dates

f <- function(x) tapply(x, m.y, mean)
UW.m <- as.data.frame(apply(UW, 2, f))

vv <- matrix(unlist(strsplit(as.character(rownames(UW.m)), "-")),ncol=2, byrow = T)

UW.m$year <- as.numeric(vv[,1])
UW.m$month <- as.numeric(vv[,2])

UW.m <- UW.m %>%
  arrange(month) %>%
  arrange(year)

# remove seasonal signal
# and set up vector of winter years (identify winters by the year corresponding to Jan.)
m <- UW.m$month
yr <- UW.m$year
win.yr <- yr
win.yr[m %in% c(11, 12)] <- win.yr[m %in% c(11, 12)] +1

f <- function(x) tapply(x, m, mean)
mu <- apply(UW.m, 2, f)	# Compute monthly means for each time series (location)

mu <- mu[rep(1:12, floor(length(d)/12)),] 
xtra <- 12*((length(d)/12)-floor(length(d)/12))
mu <- rbind(mu, mu[1:xtra,])

UW.anom <- UW.m[,1:960] - mu   # Compute matrix of anomalies - dropping year and month!

# restrict to relevant months
p.win <- c(11,12,1) # months for UW data
UW.anom <- UW.anom[UW.m$month %in% p.win,]
win.yr <- win.yr[UW.m$month %in% p.win]

# clean up
rownames(UW.anom) <- 1:nrow(UW.anom)

# and regressions
r1 <- r2 <- r3 <- p.vals <- NA # vectors to catch regression coefficients

UW.NDJ <- apply(UW.anom, 2, ff) # mean values for winter year corresponding to Jan. Note that 1949 is incomplete and will not be used!
# (1 mo and 2 mo, respectively!)

for(j in 1:ncol(UW.anom)){
  # j <- 1
  
  # subset the data for only the cell of interest and set up the early and late eras (pre/post 1988/89)
  temp <- data.frame(uw=UW.NDJ[2:71, j], pdo=pdo.FMA[50:119], era=c(rep(1, 40), rep(2, 25), rep(3,5)))
  
  # now fit era-specific regressions to plot coefficients
  mod <- lm(uw ~ pdo, data=temp[temp$era==1,])
  r1[j] <- summary(mod)$coefficients[2,1]
  
  mod <- lm(uw ~ pdo, data=temp[temp$era==2,])
  r2[j] <- summary(mod)$coefficients[2,1]
  
  mod <- lm(uw ~ pdo, data=temp[temp$era==3,])
  r3[j] <- summary(mod)$coefficients[2,1]
  
  # and add p-vals
  mod <- gls(uw ~ pdo*era, data=temp, correlation = corAR1())
  p.vals[j] <- summary(mod)$tTable[4,4]
}


# now the same regression approach for the NPGO

n1 <- n2 <- n3 <- n.vals <- NA # vectors to catch regression coefficients

for(j in 1:ncol(SLP.anom)){
  #j <- 1
  
  # subset the data for only the cell of interest and set up the early and late eras (pre/post 1988/89)
  temp <- data.frame(uw=UW.NDJ[3:71, j], npgo=npgo.FMA, era=c(rep(1, 39), rep(2, 25), rep(3,5)))
  
  # now fit era-specific regressions to plot coefficients
  mod <- lm(uw ~ npgo, data=temp[temp$era==1,])
  n1[j] <- summary(mod)$coefficients[2,1]
  
  mod <- lm(uw ~ npgo, data=temp[temp$era==2,])
  n2[j] <- summary(mod)$coefficients[2,1]
  
  mod <- lm(uw ~ npgo, data=temp[temp$era==3,])
  n3[j] <- summary(mod)$coefficients[2,1]
  
  # and add p-vals
  mod <- gls(uw ~ npgo*era, data=temp, correlation = corAR1())
  n.vals[j] <- summary(mod)$tTable[4,4]
  
}


png("U wind stress-index regression by era.png", 10, 5, units="in", res=300)

# setup the layout
mt.cex <- 1.1
l.mar <- 3
l.cex <- 0.8
l.l <- 0.2
tc.l <- -0.2

par(mar=c(0.5,0.5,1.5,1),  tcl=tc.l, mgp=c(1.5,0.3,0), las=1, mfrow=c(2,4), cex.axis=0.8, cex.lab=0.8, oma=c(0,0,0,0.2))

# set the limit for plotting 
lim <- range(r1, r2, r3)

z <- r1   
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting

image.plot(x,y,z, col=new.col, zlim=c(-lim[2], lim[2]), xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("U wind-PDO 1949-1988", cex=0.8)

z <- r2  
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting

image.plot(x,y,z, col=new.col, zlim=c(-lim[2], lim[2]), xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("U wind-PDO 1989-2013", cex=0.8)

z <- r3  
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting

image.plot(x,y,z, col=new.col, zlim=c(-lim[2], lim[2]), xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("U wind-PDO 2014-2018", cex=0.8)

# now p values
z <- r1 
z <- t(matrix(z,length(y)))  
image(x,y,z, col=new.col, zlim=c(999,9999), ylab="", xlab="", yaxt="n", xaxt="n")

z <- p.vals
z <- t(matrix(z,length(y))) 
contour(x, y, z, add=T, drawlabels = F, levels = seq(0.05, 0, length.out = 1000), col="#56B4E9", lwd=2.5) 


map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)

mtext("PDO: Era effect P < 0.05", cex=0.8)


# and npgo plots
# using same limits for z!

z <- n1   
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting

image.plot(x,y,z, col=new.col, zlim=c(-lim[2], lim[2]), xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("U wind-NPGO 1950-1988", cex=0.8)

z <- n2  
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting

image.plot(x,y,z, col=new.col, zlim=c(-lim[2], lim[2]), xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("U wind-NPGO 1989-2013", cex=0.8)

z <- n3  
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting

image.plot(x,y,z, col=new.col, zlim=c(-lim[2], lim[2]), xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("U wind-NPGO 2014-2018", cex=0.8)

# now p values
z <- r1 
z <- t(matrix(z,length(y)))  
image(x,y,z, col=new.col, zlim=c(999,9999), ylab="", xlab="", yaxt="n", xaxt="n")

z <- n.vals
z <- t(matrix(z,length(y))) 
contour(x, y, z, add=T, drawlabels = F, levels = seq(0.05, 0, length.out = 1000), col="#56B4E9", lwd=2.5) 

map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)

mtext("NPGO: Era effect P < 0.05", cex=0.8)

dev.off()

#############
# and finally, V-wind
v.w <- ncvar_get(nc.vwnd, "vwnd", verbose = F)
# Change data from a 3-D array to a matrix of monthly data by grid point:
# First, reverse order of dimensions ("transpose" array)
VW <- aperm(v.w, 3:1)  

# Reverse order of latitudes to be increasing for convenience (in later plotting)
VW <- VW[,20:1,]  

# Change to matrix with column for each grid point, rows for monthly means
VW <- matrix(VW, nrow=dim(VW)[1], ncol=prod(dim(VW)[2:3]))  

# Keep track of corresponding latitudes and longitudes of each column:
dimnames(VW) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# convert to daily wind stress!
VW <- 0.0013*1.22*VW^2
# now we need to get monthly means!
m.y <- paste(years(d), as.numeric(months(d)), sep="-") # make a month-year factor from the dates

f <- function(x) tapply(x, m.y, mean)
VW.m <- as.data.frame(apply(VW, 2, f))

vv <- matrix(unlist(strsplit(as.character(rownames(VW.m)), "-")),ncol=2, byrow = T)

VW.m$year <- as.numeric(vv[,1])
VW.m$month <- as.numeric(vv[,2])

VW.m <- VW.m %>%
  arrange(month) %>%
  arrange(year)

# remove seasonal signal
# and set up vector of winter years (identify winters by the year corresponding to Jan.)
m <- VW.m$month
yr <- VW.m$year
win.yr <- yr
win.yr[m %in% c(11, 12)] <- win.yr[m %in% c(11, 12)] +1

f <- function(x) tapply(x, m, mean)
mu <- apply(VW.m, 2, f)	# Compute monthly means for each time series (location)

mu <- mu[rep(1:12, floor(length(d)/12)),] 
xtra <- 12*((length(d)/12)-floor(length(d)/12))
mu <- rbind(mu, mu[1:xtra,])

VW.anom <- VW.m[,1:960] - mu   # Compute matrix of anomalies - dropping year and month!

# restrict to relevant months
p.win <- c(11,12,1) # months for VW data
VW.anom <- VW.anom[VW.m$month %in% p.win,]
win.yr <- win.yr[VW.m$month %in% p.win]

# clean up
rownames(VW.anom) <- 1:nrow(VW.anom)

# and regressions
r1 <- r2 <- r3 <- p.vals <- NA # vectors to catch regression coefficients

VW.NDJ <- apply(VW.anom, 2, ff) # mean values for winter year corresponding to Jan. Note that 1949 is incomplete and will not be used!
# (1 mo and 2 mo, respectively!)

for(j in 1:ncol(VW.anom)){
  # j <- 1
  
  # subset the data for only the cell of interest and set up the early and late eras (pre/post 1988/89)
  temp <- data.frame(vw=VW.NDJ[2:71, j], pdo=pdo.FMA[50:119], era=c(rep(1, 40), rep(2, 25), rep(3,5)))
  
  # now fit era-specific regressions to plot coefficients
  mod <- lm(vw ~ pdo, data=temp[temp$era==1,])
  r1[j] <- summary(mod)$coefficients[2,1]
  
  mod <- lm(vw ~ pdo, data=temp[temp$era==2,])
  r2[j] <- summary(mod)$coefficients[2,1]
  
  mod <- lm(vw ~ pdo, data=temp[temp$era==3,])
  r3[j] <- summary(mod)$coefficients[2,1]
  
  # and add p-vals
  mod <- gls(vw ~ pdo*era, data=temp, correlation = corAR1())
  p.vals[j] <- summary(mod)$tTable[4,4]
}


# now the same regression approach for the NPGO

n1 <- n2 <- n3 <- n.vals <- NA # vectors to catch regression coefficients

for(j in 1:ncol(SLP.anom)){
  #j <- 1
  
  # subset the data for only the cell of interest and set up the early and late eras (pre/post 1988/89)
  temp <- data.frame(vw=VW.NDJ[3:71, j], npgo=npgo.FMA, era=c(rep(1, 39), rep(2, 25), rep(3,5)))
  
  # now fit era-specific regressions to plot coefficients
  mod <- lm(vw ~ npgo, data=temp[temp$era==1,])
  n1[j] <- summary(mod)$coefficients[2,1]
  
  mod <- lm(vw ~ npgo, data=temp[temp$era==2,])
  n2[j] <- summary(mod)$coefficients[2,1]
  
  mod <- lm(vw ~ npgo, data=temp[temp$era==3,])
  n3[j] <- summary(mod)$coefficients[2,1]
  
  # and add p-vals
  mod <- gls(vw ~ npgo*era, data=temp, correlation = corAR1())
  n.vals[j] <- summary(mod)$tTable[4,4]
  
}


png("V wind stress-index regression by era.png", 10, 5, units="in", res=300)

# setup the layout
mt.cex <- 1.1
l.mar <- 3
l.cex <- 0.8
l.l <- 0.2
tc.l <- -0.2

par(mar=c(0.5,0.5,1.5,1),  tcl=tc.l, mgp=c(1.5,0.3,0), las=1, mfrow=c(2,4), cex.axis=0.8, cex.lab=0.8, oma=c(0,0,0,0.2))

# set the limit for plotting 
lim <- range(r1, r2, r3)

z <- r1   
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting

image.plot(x,y,z, col=new.col, zlim=c(-lim[2], lim[2]), xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("V wind-PDO 1949-1988", cex=0.8)

z <- r2  
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting

image.plot(x,y,z, col=new.col, zlim=c(-lim[2], lim[2]), xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("V wind-PDO 1989-2013", cex=0.8)

z <- r3  
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting

image.plot(x,y,z, col=new.col, zlim=c(-lim[2], lim[2]), xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("V wind-PDO 2014-2018", cex=0.8)

# now p values
z <- r1 
z <- t(matrix(z,length(y)))  
image(x,y,z, col=new.col, zlim=c(999,9999), ylab="", xlab="", yaxt="n", xaxt="n")

z <- p.vals
z <- t(matrix(z,length(y))) 
contour(x, y, z, add=T, drawlabels = F, levels = seq(0.05, 0, length.out = 1000), col="#56B4E9", lwd=2.5) 


map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)

mtext("PDO: Era effect P < 0.05", cex=0.8)


# and npgo plots
# using same limits for z!

z <- n1   
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting

image.plot(x,y,z, col=new.col, zlim=c(-lim[2], lim[2]), xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("V wind-NPGO 1950-1988", cex=0.8)

z <- n2  
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting

image.plot(x,y,z, col=new.col, zlim=c(-lim[2], lim[2]), xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("V wind-NPGO 1989-2013", cex=0.8)

z <- n3  
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting

image.plot(x,y,z, col=new.col, zlim=c(-lim[2], lim[2]), xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("V wind-NPGO 2014-2018", cex=0.8)

# now p values
z <- r1 
z <- t(matrix(z,length(y)))  
image(x,y,z, col=new.col, zlim=c(999,9999), ylab="", xlab="", yaxt="n", xaxt="n")

z <- n.vals
z <- t(matrix(z,length(y))) 
contour(x, y, z, add=T, drawlabels = F, levels = seq(0.05, 0, length.out = 1000), col="#56B4E9", lwd=2.5) 

map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)

mtext("NPGO: Era effect P < 0.05", cex=0.8)

dev.off()
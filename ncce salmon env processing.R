library(tidyverse)
library(chron)
library(fields)
library(maps)
library(mapdata)

bifur <- read.csv("bifurcation-index.csv", row.names = 1, skip=15)

head(bifur)

# upwelling
upwell.dat<-read.csv("upwelling_PFEL.csv")
mod.upwell.dat<-arrange(upwell.dat,year,lat_n, lon_w)

mod.upwell.dat$prev.nov = NA
mod.upwell.dat$prev.dec = NA

start.row<-dim(mod.upwell.dat[mod.upwell.dat$year==1946,])[1]+1
end.row<-nrow(mod.upwell.dat)-dim(mod.upwell.dat[mod.upwell.dat$year==1946,])[1]

mod.upwell.dat$prev.nov[start.row:nrow(mod.upwell.dat)] = mod.upwell.dat$nov[1:end.row]
mod.upwell.dat$prev.dec[start.row:nrow(mod.upwell.dat)] = mod.upwell.dat$dec[1:end.row]

# Then, calculate seasonal averages for upwelling

melt.uw.dat<-melt(mod.upwell.dat, id.vars=c("lat_n","lon_w","year"), variable.name = "month",
                  value.name = "uw")

#For seabirds both winter and spring ocean conditions are imporant so I'm including months
#Dec-March, but we might want to apply winter (Dec-Feb) and Spring (Mar-May) separately
#in the analysis

win.mo<-c("prev.nov","prev.dec","jan","feb", "mar")
spr.sum.mo<-c("apr","may", "jun", "jul")



mean.uw<-melt.uw.dat %>%
  group_by(year) %>%
  summarise(uw.win.51 = mean(uw[lat_n == 51 & month %in% win.mo],na.rm=TRUE),
            uw.spr.sum.51 = mean(uw[lat_n == 51 & month %in% spr.sum.mo],na.rm=TRUE),
            uw.win.48 = mean(uw[lat_n == 48 & month %in% win.mo],na.rm=TRUE),
            uw.spr.sum.48 = mean(uw[lat_n == 48 & month %in% spr.sum.mo],na.rm=TRUE))


# now ERSST
unloadNamespace("tidyverse")
library(ncdf4)

download.file("https://coastwatch.pfeg.noaa.gov/erddap/griddap/nceiErsstv5.nc?sst[(1949-01-01):1:(2019-03-01T00:00:00Z)][(0.0):1:(0.0)][(46):1:(54)][(230):1:(236)]", "~ncc.sst")

nc <- nc_open("~ncc.sst")

# 42-50 deg. N, 234-242 deg. E
x <- ncvar_get(nc, "longitude")
y <- ncvar_get(nc, "latitude")
x; y 


# get required sst data
SST <- ncvar_get(nc, "sst")

# extract dates
# seconds since 1-1-1970
raw <- ncvar_get(nc, "time")
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))
m <- months(d)
yr <- as.numeric(as.character(years(d)))

# need to change to matrix for easier use
SST <- aperm(SST, 3:1) # transpose array

SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  # Change to matrix
dim(SST)  # Matrix with column for each grid point, rows for monthly means

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes

# set names
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# plot mean temperature pattern to check
SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=tim.colors(64), xlim=c(225,238), ylim=c(44,57))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)

# blank out the ones we don't want

blank <- c("N50E230", "N48E230", "N48E232", 
           "N46E230", "N46E232", "N46E234", "N54E230")

SST[,blank] <- NA

# plot mean temperature pattern to check
SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=tim.colors(64), xlim=c(225,238), ylim=c(44,57))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
# looks good!

# now remove seasonal signal and scale 
f <- function(x) tapply(x, m, mean)
mu <- apply(SST, 2, f)	# Compute monthly means for each time series (location)
mu <- mu[rep(1:12, floor(length(d)/12)),]  # Replicate means matrix for each year at each location 

#now need to account for trailing months (i.e., fractions of a year that are available with updated data)
add <- length(d)-12*floor(length(d)/12)

# add in additional required months
mu <- rbind(mu, mu[1:add,])
# check
identical(nrow(mu), nrow(SST)) #true!

# now calculate anomalies...
SST.anom <- SST - mu

# check with some plots!
par(las=1)
ncc.sst <- ts(rowMeans(SST.anom, na.rm=T), start=c(1950,1), frequency=12)
plot(ncc.sst, type="l", col="grey", lwd=0.8, xlab="", ylab="SST anomaly")
sm <- ts(rollmean(rowMeans(SST.anom, na.rm=T), 13, fill=NA), start=c(1950,1), frequency=12)
lines(sm, col="red", lwd=1.5)
abline(h=0)
mtext("NCC SST & 13-mo rolling mean")

SST.anom <- rowMeans(SST.anom, na.rm=T)

# set up winter year and get spring/summer and winter means
win.yr <- ifelse(m %in% c("Nov", "Dec"), yr+1, yr)

temp <- SST.anom[m %in% c("Apr", "May", "Jun", "Jul")]
temp.yr <- yr[m %in% c("Apr", "May", "Jun", "Jul")]

spr.sum.SST <- tapply(temp, temp.yr, mean, na.rm=T)

# winter
temp <- SST.anom[m %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]
temp.yr <- yr[m %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]

win.SST <- tapply(temp, temp.yr, mean, na.rm=T)

ncc.dat <- data.frame(year=1950:2013, 
                      win.sst=win.SST[names(win.SST) %in% 1950:2013], 
                      spr.sum.sst=spr.sum.SST[names(spr.sum.SST) %in% 1950:2013],
                      bifur=bifur$bi_stnd[match(c(1950:2013), rownames(bifur))],
                      win.UW.51=mean.uw$uw.win.51[match(c(1950:2013), mean.uw$year)],
                      spr.sum.UW.51=mean.uw$uw.spr.sum.51[match(c(1950:2013), mean.uw$year)],
                      win.UW.48=mean.uw$uw.win.48[match(c(1950:2013), mean.uw$year)],
                      spr.sum.UW.48=mean.uw$uw.spr.sum.48[match(c(1950:2013), mean.uw$year)])

write.csv(ncc.dat, "ncc.env.dat.csv")

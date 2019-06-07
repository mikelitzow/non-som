library(reshape2)
library(ncdf4)
library(MuMIn)
library(zoo)
library(nlme)
library(dplyr)
library(chron)
library(fields)
library(maps)
library(mapdata)

# load the pdo index
pdo <- read.csv("pdo.csv")

# the PDO file that I used requires some processing
# drop winter mean and change pdo to long form data

pdo <- pdo[,-2]
pdo <- melt(pdo, id.vars = "YEAR")
pdo <- pdo[order(pdo$YEAR),]

# and get FMA mean values
pdo <- filter(pdo, variable %in% c("FEB", "MAR", "APR"))
pdo.FMA <- tapply(pdo$value, pdo$YEAR, mean) 

# load the SODA wind stress data
nc <- nc_open("/Users/MikeLitzow 1/Documents/R/pdo-npgo paper/hawaii_3e19_7ccd_16ff_ad5d_5cb7_0e6f.nc")

# view dates (middle of month):
raw <- ncvar_get(nc, "time")
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970)) # jan 1949 - dec 2010

# first, eastward wind stress!
# get all the data - they have already been subsetted by date and area in my version
tauX <- ncvar_get(nc, "taux") 

x <- ncvar_get(nc, "longitude")     # view longitudes (degrees East)
y <- ncvar_get(nc, "latitude")     # view latitudes

# process!
tauX <- aperm(tauX, 3:1)  # First, reverse order of dimensions ("transpose" array)

tauX <- matrix(tauX, nrow=dim(tauX)[1], ncol=prod(dim(tauX)[2:3]))  # Change to matrix

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes
dimnames(tauX) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

m1 <- months(d)
y1 <- years(d)
dec.yr1 <- as.numeric(as.character(y1)) + (as.numeric(m1)-0.5)/12

# and define the seasons for analysis
win <- c("Nov", "Dec", "Jan") # using NDJ as wind period to relate to FMA PDO

# define winter years
win.y1 <- as.numeric(as.character(y1))
win.y1[m1 %in% c("Nov", "Dec")] <- win.y1[m1 %in% c("Nov", "Dec")] + 1

# restrict to our selected winter months
tauX <- tauX[m1 %in% win,]

# restrict the indexing vector of winter years
win.y1 <- win.y1[m1 %in% win]

# and get annual means of these winter values
ff <- function(x) tapply(x, win.y1, mean)

tauX <- apply(tauX, 2, ff)

# now regress on the PDO for 1950:1988 and 1989:2010

# get rid of NAs for regression
land <- is.na(colMeans(tauX))  # Logical vector that's true over land!

# For analysis, we only use the columns of the matrix with non-missing values:
tauX <- tauX[,!land] 

regr.early.X <- regr.late.X <- NA # vectors for regression coefficients in both eras
X.pvals <- NA # object to catch p values

for(j in 1:ncol(tauX)){
  
  # subset for cell of interest
  temp <- data.frame(tauX=tauX[2:62, j], pdo=pdo.FMA[51:111], era=c(rep("early", 39), rep("late", 22)))
  mod <- gls(tauX ~ pdo*era, data=temp, corAR1()) # again, autocorrelated residuals allowed
  regr.early.X[j] <- summary(mod)$tTable[2,1]
  regr.late.X[j] <- regr.early.X[j] + summary(mod)$tTable[4,1]
  X.pvals[j] <- summary(mod)$tTable[4,4]
}

# And now the northward wind stress.
# load the data 
tauY <- ncvar_get(nc, "tauy") # get all the data!

# process!
tauY <- aperm(tauY, 3:1)  # First, reverse order of dimensions ("transpose" array)

tauY <- matrix(tauY, nrow=dim(tauY)[1], ncol=prod(dim(tauY)[2:3]))  # Change to matrix

dimnames(tauY) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# re-define winter years
win.y1 <- as.numeric(as.character(y1))
win.y1[m1 %in% c("Nov", "Dec")] <- win.y1[m1 %in% c("Nov", "Dec")] + 1

# restrict to our selected winter months
tauY <- tauY[m1 %in% win,]

# restrict the indexing vector of winter years
win.y1 <- win.y1[m1 %in% win]

# and get annual means of these winter values
tauY <- apply(tauY, 2, ff)

# now regress on the PDO for 1950:1988 and 1989:2010

# For analysis, we only use the columns of the matrix with non-missing values:
tauY <- tauY[,!land] 

regr.early.Y <- regr.late.Y <- NA # vectors for regression coefficients in both eras
Y.pvals <- NA # object to catch p values

for(j in 1:ncol(tauY)){
  
  # again subset by cell
  temp <- data.frame(tauY=tauY[2:62, j], pdo=pdo.FMA[51:111], era=c(rep("early", 39), rep("late", 22)))
  mod <- gls(tauY ~ pdo*era, data=temp, corAR1()) 
  regr.early.Y[j] <- summary(mod)$tTable[2,1]
  regr.late.Y[j] <- regr.early.Y[j] + summary(mod)$tTable[4,1]
  Y.pvals[j] <- summary(mod)$tTable[4,4]
}


# Now plot the combined regression coefficients

# combine the regression coefficients for the two directions
regr.early.XY <- sqrt(regr.early.X^2 + regr.early.Y^2)
regr.late.XY <- sqrt(regr.late.X^2 + regr.late.Y^2)

# set up the color scheme
new.col <- tim.colors(64)
grays <- c("gray98", "gray97", "gray96", "gray95", "gray94", "gray93", "gray92", "gray91", "gray90", "gray89", "gray88")
new.col[27:36] <- c(grays[5:1], grays[1:5])

# setup the layout
l.mar <- 3
l.cex <- 0.8
l.l <- 0.2
tc.l <- -0.2

par(mar=c(1.5,2.5,1,0.5),  tcl=tc.l, mfrow=c(1,2), oma=c(0,0,0,0.2))

zlim <- range(regr.early.XY, regr.late.XY)

z <- rep(NA, ncol(tauY))
z[!land] <- regr.early.XY 
z <- t(matrix(z,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image.plot(x,y,z, col=new.col, zlim=c(-zlim[2],zlim[2]), ylab="", xlab="", yaxt="n", xaxt="n",legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("Wind stress-PDO 1950-1988", cex=0.8)

z <- rep(NA, ncol(tauY))
z[!land] <- regr.late.XY
z <- t(matrix(z,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image.plot(x,y,z, col=new.col, zlim=c(-zlim[2],zlim[2]), ylab="", xlab="", yaxt="n", xaxt="n",legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("Wind stress-PDO 1989-2012", cex=0.8)

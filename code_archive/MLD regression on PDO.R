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
library(tidyr)

# load and process the PDO
names <- read.table("~pdo", skip=30, nrows=1, as.is = T)
pdo <- read.table("~pdo", skip=31, nrows=119, fill=T, col.names = names)
pdo$YEAR <- 1900:(1899+nrow(pdo)) # drop asterisks!
pdo <- pdo %>%
  gather(month, value, -YEAR) %>%
  arrange(YEAR)

# and get FMA mean values
pdo <- filter(pdo, month %in% c("JAN", "FEB", "MAR"))
pdo.FMA <- tapply(pdo$value, pdo$YEAR, mean) 

# load the MLD data
nc <- nc_open("MLD_ORAS3.nc")

# view dates (middle of month):
d <- ncvar_get(nc, "TIME")
d <- dates(d, origin = c(1,1,0001)) 

# get all the data - they have already been subsetted by date and area in my version
mld <- ncvar_get(nc, "MLD") 

x <- ncvar_get(nc, "LON122_168")     # view longitudes (degrees East)
y <- ncvar_get(nc, "LAT157_175")     # view latitudes

# process!
mld <- aperm(mld, 3:1)  # First, reverse order of dimensions ("transpose" array)

mld <- matrix(mld, nrow=dim(mld)[1], ncol=prod(dim(mld)[2:3]))  # Change to matrix

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes
dimnames(mld) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

m1 <- months(d)
y1 <- years(d)
dec.yr1 <- as.numeric(as.character(y1)) + (as.numeric(m1)-0.5)/12

# and define the seasons for analysis
use <- c("Jan", "Feb", "Mar") # using NDJ as wind period to relate to FMA PDO

# restrict to our selected months
mld <- mld[m1 %in% use,]

# restrict the indexing vector of winter years
y1 <- y1[m1 %in% use]

# and get annual means of these winter values
ff <- function(x) tapply(x, y1, mean)

mld <- apply(mld, 2, ff)

# now regress on the PDO for 1959:1988 and 1989:2011

# get rid of NAs for regression
land <- is.na(colMeans(mld))  # Logical vector that's true over land!


# For analysis, we only use the columns of the matrix with non-missing values:
mld <- mld[,!land] 

regr.early <- regr.late <- NA # vectors for regression coefficients in both eras
X.pvals <- NA # object to catch p values

mld <- log(mld) # log transform...

for(j in 1:ncol(mld)){
 # j <- 1
  # subset for cell of interest
  temp <- data.frame(mld=mld[2:53,j], pdo=pdo.FMA[61:112], era=c(rep("early", 29), rep("late", 23)))
  mod <- gls(mld ~ pdo*era, data=temp, corAR1()) # again, autocorrelated residuals allowed
  regr.early[j] <- summary(mod)$tTable[2,1]
  regr.late[j] <- regr.early[j] + summary(mod)$tTable[4,1]
  # X.pvals[j] <- summary(mod)$tTable[4,4]
}


# Now plot 

# set up the color scheme
new.col <- tim.colors(64)
grays <- c("gray98", "gray97", "gray96", "gray95", "gray94", "gray93", "gray92", "gray91", "gray90", "gray89", "gray88")
new.col[27:36] <- c(grays[5:1], grays[1:5])

# setup the layout
l.mar <- 3
l.cex <- 0.8
l.l <- 0.2
tc.l <- -0.2

# two panel layout
png("mld vs pdo.png", 8, 4, units="in", res=300)
par(mar=c(1.5,2.5,1,0.5),  tcl=tc.l, mfrow=c(1,2), oma=c(0,0,0,0.2))

zlim <- range(regr.early, regr.late)

z <- rep(NA, ncol(mld))
z[!land] <- regr.early 
z <- t(matrix(z,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image.plot(x,y,z, col=new.col, zlim=c(zlim[1],-zlim[1]), ylab="", xlab="", yaxt="n", xaxt="n",legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("Ln(mld)-PDO 1960-1988", cex=0.8)

z <- rep(NA, ncol(mld))
z[!land] <- regr.late
z <- t(matrix(z,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image.plot(x,y,z, col=new.col, zlim=c(zlim[1],-zlim[1]), ylab="", xlab="", yaxt="n", xaxt="n",legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("Ln(mld)-PDO 1989-2011", cex=0.8)

dev.off()


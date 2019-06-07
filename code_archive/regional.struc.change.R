# code for testing for change in relationship between PDO/NPGO and regional-scale climate and biology variables

library(dplyr)
library(MuMIn)
library(tidyr)
library(ggplot2)

# load GOA climate data

dat <- read.csv("/Users/MikeLitzow 1/Documents/R/FATE2/data/GOA.climate.csv", row.names=1)

unique(dat$key)


unique(dat$key)

keep <- c("FMA.FW", "NDJ.grad", "MJJ.UW", "FMA.WS", "Papa", "FMA.SSH")

sdat <- filter(dat, key %in% keep)

x <- filter(dat, key=="FMA.PDO")

codes <- unique(sdat$key)

y <- filter(dat, key==codes[1])

identical(x$year, y$year)

nrow(x)

out <-matrix(nrow=63-26, ncol=length(keep))
dimnames(out) <- list(1963:1999, keep)


splits <- 1963:1999

for(j in 1:length(keep)){ # loop through each column
  
y <- filter(sdat, key==keep[j])
for(i in 1:length(splits)){ # loop through each split point
  #i <- 1
  x$era <- 1
  x$era[x$year>splits[i]] <- 2
  out[i,j] <- AICc(lm(y$value ~ x$value*x$era))
}
} 

# change to AICc
delta.out <- out
for(j in 1:ncol(out)){
  delta.out[,j] <- out[,j] - min(out[,j])
}


plot <- gather(as.data.frame(delta.out))
plot$year <- splits

ggplot(plot, aes(year, value)) + geom_line() + facet_wrap(~key, scales="free") + ylab("delta-AICc")

# producing different smoothings for winter PDO/NPGO for use in salmon analysis

library(tidyverse)
library(zoo)

# load pdo and npgo
download.file("http://jisao.washington.edu/pdo/PDO.latest", "~pdo")
names <- read.table("~pdo", skip=30, nrows=1, as.is = T)
pdo <- read.table("~pdo", skip=31, nrows=119, fill=T, col.names = names)
pdo$YEAR <- 1900:(1899+nrow(pdo)) # drop asterisks!
pdo <- pdo %>%
  gather(month, value, -YEAR) %>%
  arrange(YEAR)

# load npgo
download.file("http://www.oces.us/npgo/npgo.php", "~npgo")
npgo <- read.table("~npgo", skip=10, nrows=828, fill=T, col.names = c("Year", "month", "value"))

win.pdo <- pdo[pdo$month %in% c("NOV", "DEC", "JAN", "FEB", "MAR"),]
win.pdo$win.yr <- ifelse(win.pdo$month %in% c("NOV", "DEC"), win.pdo$YEAR+1, win.pdo$YEAR)
win.pdo <- tapply(win.pdo$value, win.pdo$win.yr, mean)

win.npgo <- npgo[npgo$month %in% c(11,12,1:3),]
win.npgo$win.yr <- ifelse(win.npgo$month %in% 11:12, win.npgo$Year+1, win.npgo$Year)

win.npgo <- tapply(win.npgo$value, win.npgo$win.yr, mean)

# limit both to 1951:2018
win.pdo <- win.pdo[names(win.pdo) %in% 1951:2018]
win.npgo <- win.npgo[names(win.npgo) %in% 1951:2018]


# put into a data frame with various smoothings
# 2a = winter before and winter of ocean entry year
# 2b = winter of and winter after ocean entry year

pdo.npgo <- data.frame(year=1951:2018, pdo1=win.pdo, pdo2a=rollmean(win.pdo, 2, align="right", fill=NA), pdo2b=rollmean(win.pdo, 2, align="left", fill=NA),
                       pdo3=rollmean(win.pdo, 3, fill=NA), npgo=win.npgo, npgo2a=rollmean(win.npgo, 2, align="right", fill=NA),
                       npgo2b=rollmean(win.npgo, 2, align="left", fill=NA), npgo3=rollmean(win.npgo, 3, fill=NA))

# clean up
row.names(pdo.npgo) <- 1:nrow(pdo.npgo)

# save 
write.csv(pdo.npgo, "winter pdo npgo various smoothings.csv")

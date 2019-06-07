library(tidyverse)
library(zoo)
library(nlme)
library(MARSS)
library(MuMIn)

# script to test for changing relationships between PDO/NPGO and regional biology

# I start by fitting best one-trend DFA model to each data set in each era, and 
# then fit regressions for shared trends onto PDO/NPGO

# download PDO / NPGO and process
download.file("http://jisao.washington.edu/pdo/PDO.latest", "~pdo")
names <- read.table("~pdo", skip=30, nrows=1, as.is = T)
pdo <- read.table("~pdo", skip=32, nrows=119, fill=T, col.names = names)
pdo$YEAR <- 1900:(1899+nrow(pdo)) # drop asterisks!
pdo <- pdo %>%
  gather(month, value, -YEAR) %>%
  arrange(YEAR)

download.file("http://www.oces.us/npgo/npgo.php", "~npgo")

npgo <- read.table("~npgo", skip=10, nrows=828, fill=T, col.names = c("Year", "month", "value"))

# calculate NDJFM means for each index
pdo$win.yr <- ifelse(pdo$month %in% c("NOV", "DEC"), pdo$YEAR+1, pdo$YEAR)
win.pdo <- tapply(pdo$value, pdo$win.yr, mean)

npgo$win.yr <- ifelse(npgo$month %in% 11:12, npgo$Year+1, npgo$Year)
win.npgo <- tapply(npgo$value, npgo$win.yr, mean)

# and smoothed (2yr) values of each
npgo2 <- rollapply(win.npgo, 2, mean, align="right", fill=NA)
pdo2 <- rollapply(win.pdo, 2, mean, align="right", fill=NA)
names(pdo2) <- 1900:2019


# set up forms of R matrices and control values for DFA
levels.R = c("diagonal and equal",
             "diagonal and unequal",
             "equalvarcov",
             "unconstrained")

cntl.list = list(minit=200, maxit=10000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)


# now load farallon seabirds data

# dat <- read.csv("CCE-nonsalmon-biology&climate.csv")
# 
# our.dat <- grep("SBRD", colnames(dat))
# year.dat <- grep("year", colnames(dat))
# 
# dat <- dat[,c(year.dat, our.dat)]
# head(dat)
# 
# # limit to 1971:2013
# dat <- dat[dat$year %in% 1971:2013,]
# 
# write.csv(dat, "farallon.sbrd.biol.csv")

dat <- read.csv("farallon.sbrd.biol.csv", row.names = 1)
# transform
run.dat <- t(dat[,2:8]) 
colnames(run.dat) <- dat$year


# # select best DFA model

# model.data = data.frame()

# for(R in levels.R) {
#   for(m in 1:1) {  # allowing only 1 trend!
#     dfa.model = list(A="zero", R=R, m=m)
#     
#     
#     
#     kemz = MARSS(run.dat[,colnames(run.dat) <= 1988], model=dfa.model, control=cntl.list,
#                  form="dfa", z.score=TRUE)
#     model.data = rbind(model.data,
#                        data.frame(R=R,
#                                   m=m,
#                                   logLik=kemz$logLik,
#                                   K=kemz$num.params,
#                                   AICc=kemz$AICc,
#                                   stringsAsFactors=FALSE))
#     assign(paste("kemz", m, R, sep="."), kemz)
#   } # end m loop
# } # end R loop
# 
#  model.data

# diagonal and equal is best model
# refit and save
dfa.model = list(A="zero", R="diagonal and equal", m=1)
fara.sbrd.era1 <- MARSS(run.dat[,colnames(run.dat) <= 1988], model=dfa.model, control=cntl.list,
                      form="dfa", z.score=TRUE)

# now era 2
# model.data = data.frame()
# 
# for(R in levels.R) {
#   for(m in 1:1) {  # allowing only 1 trend!
#     dfa.model = list(A="zero", R=R, m=m)
#     
#     
#     
#     kemz = MARSS(run.dat[,colnames(run.dat) > 1988], model=dfa.model, control=cntl.list,
#                  form="dfa", z.score=TRUE)
#     model.data = rbind(model.data,
#                        data.frame(R=R,
#                                   m=m,
#                                   logLik=kemz$logLik,
#                                   K=kemz$num.params,
#                                   AICc=kemz$AICc,
#                                   stringsAsFactors=FALSE))
#     assign(paste("kemz", m, R, sep="."), kemz)
#   } # end m loop
# } # end R loop
# 
# model.data

# diagonal and unequal is best model
# refit and save
dfa.model = list(A="zero", R="diagonal and unequal", m=1)
fara.sbrd.era2 <- MARSS(run.dat[,colnames(run.dat) > 1988], model=dfa.model, control=cntl.list,
                      form="dfa", z.score=TRUE)
#####################

### now calcofi ichthyo
# dat <- read.csv("updated 1804 1601 1704 1604 1501 1407 1311 ichthyoplankton by line and station.csv", row.names = 1)
# 
# 
# sorter <- colMeans(dat[,10:ncol(dat)])
# sorter <- sorter[order(-sorter)]
# 
# 
# keepers <- names(sorter)[1:12]
# 
# keep.dat <- colnames(dat) %in% keepers
# 
# dat <- cbind(dat[,1:9], dat[,keep.dat])
# 
# # log transform
# dat[,10:ncol(dat)] <- log(dat[,10:ncol(dat)]+1)
# 
# # mean seasonal catch
# try.dat <- dat[,9:ncol(dat)] %>%
#   gather(species, catch, -season)
# 
# tapply(try.dat$catch, try.dat$season, mean) # so spring (unsurprisingly) is the highest abundance season
# 
# # restrict to spring
# dat <- dat %>%
#   filter(season=="spring")
# 
# # restrict to best-sampled S_S ("stations"?)
# keep <- c(28, 30, 35, 37, 40, 45, 50, 51, 55, 60, 70, 80, 90)
# dat <- dat %>%
#   filter(S_S %in% keep)
# 
# dat <- dat[,c(8,10:ncol(dat))] %>%
#   gather(key, catch, -year) %>%
#   group_by(year, key) %>%
#   summarise(mean(catch))
# 
# names(dat)[3] <- "value"
# 
# dat <- dat %>%
#   spread(key, value)
# 
# # restrict to 1952:2013!
# # (this is the first year for which we have a 2-yr mean for npgo)
# dat <- dat %>%
#   filter(year %in% 1952:2013)
# 
# write.csv(dat, "calcofi.biol.csv")

dat <- read.csv("calcofi.biol.csv", row.names=1)
# and transpose
run.dat <- as.matrix(t(dat[,2:13]))
colnames(run.dat) <- dat$year

# # run each era!
# model.data = data.frame()
# 
# # select model
# for(R in levels.R) {
#   for(m in 1:1) {  # allowing only 1 trend!
#     dfa.model = list(A="zero", R=R, m=m)
#     
#     
#     
#     kemz = MARSS(run.dat[,colnames(run.dat) <= 1988], model=dfa.model, control=cntl.list,
#                  form="dfa", z.score=TRUE)
#     model.data = rbind(model.data,
#                        data.frame(R=R,
#                                   m=m,
#                                   logLik=kemz$logLik,
#                                   K=kemz$num.params,
#                                   AICc=kemz$AICc,
#                                   stringsAsFactors=FALSE))
#     assign(paste("kemz", m, R, sep="."), kemz)
#   } # end m loop
# } # end R loop
# 
# model.data

# equalvarcov is best model
# refit and save
dfa.model = list(A="zero", R="equalvarcov", m=1)
calcof.ichthyo.era1 <- MARSS(run.dat[,colnames(run.dat) <= 1988], model=dfa.model, control=cntl.list,
                        form="dfa", z.score=TRUE)



# now era 2
# model.data = data.frame()
# 
# for(R in levels.R) {
#   for(m in 1:1) {  # allowing only 1 trend!
#     dfa.model = list(A="zero", R=R, m=m)
#     
#     
#     
#     kemz = MARSS(run.dat[,colnames(run.dat) > 1988], model=dfa.model, control=cntl.list,
#                  form="dfa", z.score=TRUE)
#     model.data = rbind(model.data,
#                        data.frame(R=R,
#                                   m=m,
#                                   logLik=kemz$logLik,
#                                   K=kemz$num.params,
#                                   AICc=kemz$AICc,
#                                   stringsAsFactors=FALSE))
#     assign(paste("kemz", m, R, sep="."), kemz)
#   } # end m loop
# } # end R loop
# 
# model.data

# equalvarcov is best model
# refit and save
dfa.model = list(A="zero", R="equalvarcov", m=1)
calcof.ichth.era2 <- MARSS(run.dat[,colnames(run.dat) > 1988], model=dfa.model, control=cntl.list,
                        form="dfa", z.score=TRUE)



#############
# GOA biology
# dat <- read.csv("/Users/MikeLitzow 1/Documents/R/time and climes/salmon and non-salmon biology mar 28.csv", row.names = 1)
# 
# dat <- dat[,c(1, 3:10,15)]
# 
# write.csv(dat, "goa.biol.csv")

dat <- read.csv("goa.biol.csv", row.names = 1)

run.dat <- t(dat[rownames(dat) <= 2013,])

# # run each era!
# model.data = data.frame()
# 
# for(R in levels.R) {
#   for(m in 1:1) {  # allowing only 1 trend!
#     dfa.model = list(A="zero", R=R, m=m)
#     
#     
#     
#     kemz = MARSS(run.dat[,colnames(run.dat) <= 1988], model=dfa.model, control=cntl.list,
#                  form="dfa", z.score=TRUE)
#     model.data = rbind(model.data,
#                        data.frame(R=R,
#                                   m=m,
#                                   logLik=kemz$logLik,
#                                   K=kemz$num.params,
#                                   AICc=kemz$AICc,
#                                   stringsAsFactors=FALSE))
#     assign(paste("kemz", m, R, sep="."), kemz)
#   } # end m loop
# } # end R loop
# 
# model.data

# diagonal and unequal is best model
# refit and save
dfa.model = list(A="zero", R="diagonal and unequal", m=1)
goa.bio.era1 <- MARSS(run.dat[,colnames(run.dat) <= 1988], model=dfa.model, control=cntl.list,
                         form="dfa", z.score=TRUE)

# now era 2
# model.data = data.frame()
# 
# for(R in levels.R) {
#   for(m in 1:1) {  # allowing only 1 trend!
#     dfa.model = list(A="zero", R=R, m=m)
#     
#     
#     
#     kemz = MARSS(run.dat[,colnames(run.dat) > 1988], model=dfa.model, control=cntl.list,
#                  form="dfa", z.score=TRUE)
#     model.data = rbind(model.data,
#                        data.frame(R=R,
#                                   m=m,
#                                   logLik=kemz$logLik,
#                                   K=kemz$num.params,
#                                   AICc=kemz$AICc,
#                                   stringsAsFactors=FALSE))
#     assign(paste("kemz", m, R, sep="."), kemz)
#   } # end m loop
# } # end R loop
# 
# model.data

# diagonal and equal is best model
# refit and save
dfa.model = list(A="zero", R="diagonal and equal", m=1)
goa.bio.era2 <- MARSS(run.dat[,colnames(run.dat) > 1988], model=dfa.model, control=cntl.list,
                         form="dfa", z.score=TRUE)


#####################################################
# combine regression coeffs for a plot
#####################################################
# make a df to hold regression coeffs
regr.out <- data.frame(class=rep(NA, 12), model=NA, era=NA, index= NA, lower=NA, est=NA, upper=NA)

regr.out$class[1:4] <- "biol"
regr.out$model[1:4] <- "California seabirds"
regr.out$era[1:4] <- rep(c(1,2), each=2)
regr.out$index[1:4] <- rep(c("PDO", "NPGO"), 2)

# start with Farallon sbrds
# reload data
dat <- read.csv("farallon.sbrd.biol.csv", row.names = 1)

# transpose
run.dat <- t(dat[,2:8]) 
colnames(run.dat) <- dat$year

# identify the relevant DFA model
this.mod <- fara.sbrd.era1
keep <- colnames(run.dat) <= 1988

dat <- data.frame(year=as.numeric(as.character(colnames(run.dat)[keep])), trend=scale(as.vector(this.mod$states)))
dat$pdo <- win.pdo[match(dat$year, names(win.pdo))]
dat$pdo2 <- pdo2[match(dat$year, names(pdo2))]

dat$npgo <- win.npgo[match(dat$year, names(win.npgo))]
dat$npgo2 <- npgo2[match(dat$year, names(npgo2))]

# PDO regression
mod1 <- gls(trend ~ pdo, data=dat, correlation = corAR1())
mod2 <- gls(trend ~ pdo2, data=dat, correlation = corAR1())

# store the regression coefficient and CI for the best model
regr.out[1,5] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), intervals(mod1)$coef[2,1],
                                  intervals(mod2)$coef[2,1]))

regr.out[1,6] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), intervals(mod1)$coef[2,2],
                                  intervals(mod2)$coef[2,2]))

regr.out[1,7] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), intervals(mod1)$coef[2,3],
                                  intervals(mod2)$coef[2,3]))
# NPGO regression
mod1 <- gls(trend ~ npgo, data=dat, correlation = corAR1())
mod2 <- gls(trend ~ npgo2, data=dat, correlation = corAR1())

# store the regression coefficient and CI for the best model
regr.out[2,5] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), intervals(mod1)$coef[2,1],
                                  intervals(mod2)$coef[2,1]))

regr.out[2,6] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), intervals(mod1)$coef[2,2],
                                  intervals(mod2)$coef[2,2]))

regr.out[2,7] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), intervals(mod1)$coef[2,3],
                                  intervals(mod2)$coef[2,3]))

# second era
this.mod <- fara.sbrd.era2
keep <- colnames(run.dat) > 1988

dat <- data.frame(year=as.numeric(as.character(colnames(run.dat)[keep])), trend=scale(as.vector(this.mod$states)))
dat$pdo <- win.pdo[match(dat$year, names(win.pdo))]
dat$pdo2 <- pdo2[match(dat$year, names(pdo2))]

dat$npgo <- win.npgo[match(dat$year, names(win.npgo))]
dat$npgo2 <- npgo2[match(dat$year, names(npgo2))]

# PDO regression
mod1 <- gls(trend ~ pdo, data=dat, correlation = corAR1())
mod2 <- gls(trend ~ pdo2, data=dat, correlation = corAR1())

# store the regression coefficient and CI for the best model
regr.out[3,5] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), intervals(mod1)$coef[2,1],
                                  intervals(mod2)$coef[2,1]))

regr.out[3,6] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), intervals(mod1)$coef[2,2],
                                  intervals(mod2)$coef[2,2]))

regr.out[3,7] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), intervals(mod1)$coef[2,3],
                                  intervals(mod2)$coef[2,3]))
# NPGO regression
mod1 <- gls(trend ~ npgo, data=dat, correlation = corAR1())
mod2 <- gls(trend ~ npgo2, data=dat, correlation = corAR1())

# store the regression coefficient and CI for the best model
regr.out[4,5] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), intervals(mod1)$coef[2,1],
                                  intervals(mod2)$coef[2,1]))

regr.out[4,6] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), intervals(mod1)$coef[2,2],
                                  intervals(mod2)$coef[2,2]))

regr.out[4,7] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), intervals(mod1)$coef[2,3],
                                  intervals(mod2)$coef[2,3]))

##############################
### now calcofi ichthyo
dat <- read.csv("calcofi.biol.csv", row.names = 1)


# and transpose
run.dat <- as.matrix(t(dat[,2:13]))
colnames(run.dat) <- dat$year

# and save regressions
this.mod <- calcof.ichthyo.era1
keep <- colnames(run.dat) <= 1988

regr.out$class[5:8] <- "biol"
regr.out$model[5:8] <- "CalCOFI ichthyo"
regr.out$era[5:8] <- rep(c(1,2), each=2)
regr.out$index[5:8] <- rep(c("PDO", "NPGO"), 2)

dat <- data.frame(year=as.numeric(as.character(colnames(run.dat)[keep])), trend=scale(as.vector(this.mod$states)))
dat$pdo <- win.pdo[match(dat$year, names(win.pdo))]
dat$pdo2 <- pdo2[match(dat$year, names(pdo2))]

dat$npgo <- win.npgo[match(dat$year, names(win.npgo))]
dat$npgo2 <- npgo2[match(dat$year, names(npgo2))]

# PDO regression
# mod1 <- gls(trend ~ pdo, data=dat, correlation = corAR1())
# mod2 <- gls(trend ~ pdo2, data=dat, correlation = corAR1())

# switch to lm as phi is estimated ~1!
mod1 <- lm(trend ~ pdo, data=dat)
mod2 <- lm(trend ~ pdo2, data=dat)
# store the regression coefficient and CI for the best model
# regr.out[7,5] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), intervals(mod1)$coef[2,1],
#                                   intervals(mod2)$coef[2,1]))
# 
# regr.out[7,6] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), intervals(mod1)$coef[2,2],
#                                   intervals(mod2)$coef[2,2]))
# 
# regr.out[7,7] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), intervals(mod1)$coef[2,3],
#                                   intervals(mod2)$coef[2,3]))

# again, change to lm
# lower CI
regr.out[5,5] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), summary(mod1)$coef[2,1]-1.96*summary(mod1)$coef[2,2],
                                  summary(mod2)$coef[2,1]-1.96*summary(mod2)$coef[2,2]))
# estimate
regr.out[5,6] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), summary(mod1)$coef[2,1],
                                  summary(mod2)$coef[2,1]))

# upper CI
regr.out[5,7] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), summary(mod1)$coef[2,1]+1.96*summary(mod1)$coef[2,2],
                                  summary(mod2)$coef[2,1]+1.96*summary(mod2)$coef[2,2]))


# NPGO regression
# mod1 <- gls(trend ~ npgo, data=dat, correlation = corAR1())
# mod2 <- gls(trend ~ npgo2, data=dat, correlation = corAR1())
# change to lm
mod1 <- lm(trend ~ npgo, data=dat)
mod2 <- lm(trend ~ npgo2, data=dat)
# store the regression coefficient and CI for the best model
regr.out[6,5] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), summary(mod1)$coef[2,1]-1.96*summary(mod1)$coef[2,2],
                                  summary(mod2)$coef[2,1]-1.96*summary(mod2)$coef[2,2]))

regr.out[6,6] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), summary(mod1)$coef[2,1],
                                  summary(mod2)$coef[2,1]))

regr.out[6,7] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), summary(mod1)$coef[2,1]+1.96*summary(mod1)$coef[2,2],
                                  summary(mod2)$coef[2,1]+1.96*summary(mod2)$coef[2,2]))

# now the second era
this.mod <- calcof.ichth.era2
keep <- colnames(run.dat) > 1988

dat <- data.frame(year=as.numeric(as.character(colnames(run.dat)[keep])), trend=scale(as.vector(this.mod$states)))
dat$pdo <- win.pdo[match(dat$year, names(win.pdo))]
dat$pdo2 <- pdo2[match(dat$year, names(pdo2))]

dat$npgo <- win.npgo[match(dat$year, names(win.npgo))]
dat$npgo2 <- npgo2[match(dat$year, names(npgo2))]

# PDO regression
# mod1 <- gls(trend ~ pdo, data=dat, correlation = corAR1())
# mod2 <- gls(trend ~ pdo2, data=dat, correlation = corAR1())

# switch to lm as phi is estimated ~1!
mod1 <- lm(trend ~ pdo, data=dat)
mod2 <- lm(trend ~ pdo2, data=dat)
# store the regression coefficient and CI for the best model
# regr.out[7,5] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), intervals(mod1)$coef[2,1],
#                                   intervals(mod2)$coef[2,1]))
# 
# regr.out[7,6] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), intervals(mod1)$coef[2,2],
#                                   intervals(mod2)$coef[2,2]))
# 
# regr.out[7,7] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), intervals(mod1)$coef[2,3],
#                                   intervals(mod2)$coef[2,3]))

# again, change to lm
# lower CI
regr.out[7,5] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), summary(mod1)$coef[2,1]-1.96*summary(mod1)$coef[2,2],
                                  summary(mod2)$coef[2,1]-1.96*summary(mod2)$coef[2,2]))
# estimate
regr.out[7,6] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), summary(mod1)$coef[2,1],
                                  summary(mod2)$coef[2,1]))

# upper CI
regr.out[7,7] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), summary(mod1)$coef[2,1]+1.96*summary(mod1)$coef[2,2],
                                  summary(mod2)$coef[2,1]+1.96*summary(mod2)$coef[2,2]))


# NPGO regression
# mod1 <- gls(trend ~ npgo, data=dat, correlation = corAR1())
# mod2 <- gls(trend ~ npgo2, data=dat, correlation = corAR1())
# change to lm
mod1 <- lm(trend ~ npgo, data=dat)
mod2 <- lm(trend ~ npgo2, data=dat)
# store the regression coefficient and CI for the best model
regr.out[8,5] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), summary(mod1)$coef[2,1]-1.96*summary(mod1)$coef[2,2],
                                  summary(mod2)$coef[2,1]-1.96*summary(mod2)$coef[2,2]))

regr.out[8,6] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), summary(mod1)$coef[2,1],
                                  summary(mod2)$coef[2,1]))

regr.out[8,7] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), summary(mod1)$coef[2,1]+1.96*summary(mod1)$coef[2,2],
                                  summary(mod2)$coef[2,1]+1.96*summary(mod2)$coef[2,2]))

##################
# GOA biology!

regr.out$class[9:12] <- "biol"
regr.out$model[9:12] <- "GOA fish/crustaceans"
regr.out$era[9:12] <- rep(c(1,2), each=2)
regr.out$index[9:12] <- rep(c("PDO", "NPGO"), 2)

dat <- read.csv("goa.biol.csv", row.names = 1)

run.dat <- t(dat[rownames(dat) <= 2013,])

# identify the relevant DFA model
this.mod <- goa.bio.era1
keep <- colnames(run.dat) <= 1988

dat <- data.frame(year=as.numeric(as.character(colnames(run.dat)[keep])), trend=scale(as.vector(this.mod$states)))
dat$pdo <- win.pdo[match(dat$year, names(win.pdo))]
dat$pdo2 <- pdo2[match(dat$year, names(pdo2))]

dat$npgo <- win.npgo[match(dat$year, names(win.npgo))]
dat$npgo2 <- npgo2[match(dat$year, names(npgo2))]


# switch to lm as phi is estimated ~1!
mod1 <- lm(trend ~ pdo, data=dat)
mod2 <- lm(trend ~ pdo2, data=dat)
# store the regression coefficient and CI for the best model
# regr.out[7,5] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), intervals(mod1)$coef[2,1],
#                                   intervals(mod2)$coef[2,1]))
# 
# regr.out[7,6] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), intervals(mod1)$coef[2,2],
#                                   intervals(mod2)$coef[2,2]))
# 
# regr.out[7,7] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), intervals(mod1)$coef[2,3],
#                                   intervals(mod2)$coef[2,3]))

# again, change to lm
# lower CI
regr.out[9,5] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), summary(mod1)$coef[2,1]-1.96*summary(mod1)$coef[2,2],
                                  summary(mod2)$coef[2,1]-1.96*summary(mod2)$coef[2,2]))
# estimate
regr.out[9,6] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), summary(mod1)$coef[2,1],
                                  summary(mod2)$coef[2,1]))

# upper CI
regr.out[9,7] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), summary(mod1)$coef[2,1]+1.96*summary(mod1)$coef[2,2],
                                  summary(mod2)$coef[2,1]+1.96*summary(mod2)$coef[2,2]))


# NPGO regression
# mod1 <- gls(trend ~ npgo, data=dat, correlation = corAR1())
# mod2 <- gls(trend ~ npgo2, data=dat, correlation = corAR1())
# change to lm
mod1 <- lm(trend ~ npgo, data=dat)
mod2 <- lm(trend ~ npgo2, data=dat)
# store the regression coefficient and CI for the best model
regr.out[10,5] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), summary(mod1)$coef[2,1]-1.96*summary(mod1)$coef[2,2],
                                  summary(mod2)$coef[2,1]-1.96*summary(mod2)$coef[2,2]))

regr.out[10,6] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), summary(mod1)$coef[2,1],
                                  summary(mod2)$coef[2,1]))

regr.out[10,7] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), summary(mod1)$coef[2,1]+1.96*summary(mod1)$coef[2,2],
                                  summary(mod2)$coef[2,1]+1.96*summary(mod2)$coef[2,2]))

# second era
this.mod <- goa.bio.era2
keep <- colnames(run.dat) > 1988

dat <- data.frame(year=as.numeric(as.character(colnames(run.dat)[keep])), trend=scale(as.vector(this.mod$states)))
dat$pdo <- win.pdo[match(dat$year, names(win.pdo))]
dat$pdo2 <- pdo2[match(dat$year, names(pdo2))]

dat$npgo <- win.npgo[match(dat$year, names(win.npgo))]
dat$npgo2 <- npgo2[match(dat$year, names(npgo2))]

# switch to lm as phi is estimated ~1!
mod1 <- lm(trend ~ pdo, data=dat)
mod2 <- lm(trend ~ pdo2, data=dat)
# store the regression coefficient and CI for the best model
# regr.out[7,5] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), intervals(mod1)$coef[2,1],
#                                   intervals(mod2)$coef[2,1]))
# 
# regr.out[7,6] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), intervals(mod1)$coef[2,2],
#                                   intervals(mod2)$coef[2,2]))
# 
# regr.out[7,7] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), intervals(mod1)$coef[2,3],
#                                   intervals(mod2)$coef[2,3]))

# again, change to lm
# lower CI
regr.out[11,5] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), summary(mod1)$coef[2,1]-1.96*summary(mod1)$coef[2,2],
                                  summary(mod2)$coef[2,1]-1.96*summary(mod2)$coef[2,2]))
# estimate
regr.out[11,6] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), summary(mod1)$coef[2,1],
                                  summary(mod2)$coef[2,1]))

# upper CI
regr.out[11,7] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), summary(mod1)$coef[2,1]+1.96*summary(mod1)$coef[2,2],
                                  summary(mod2)$coef[2,1]+1.96*summary(mod2)$coef[2,2]))


# NPGO regression
# mod1 <- gls(trend ~ npgo, data=dat, correlation = corAR1())
# mod2 <- gls(trend ~ npgo2, data=dat, correlation = corAR1())
# change to lm
mod1 <- lm(trend ~ npgo, data=dat)
mod2 <- lm(trend ~ npgo2, data=dat)
# store the regression coefficient and CI for the best model
regr.out[12,5] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), summary(mod1)$coef[2,1]-1.96*summary(mod1)$coef[2,2],
                                   summary(mod2)$coef[2,1]-1.96*summary(mod2)$coef[2,2]))

regr.out[12,6] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), summary(mod1)$coef[2,1],
                                   summary(mod2)$coef[2,1]))

regr.out[12,7] <- as.vector(ifelse(AICc(mod1) < AICc(mod2), summary(mod1)$coef[2,1]+1.96*summary(mod1)$coef[2,2],
                                   summary(mod2)$coef[2,1]+1.96*summary(mod2)$coef[2,2]))



# plot
dodge <- position_dodge(width=0.9)

# load the color-blind palette
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
regr.out$era <- as.factor(regr.out$era)
regr.out$index.order <- ifelse(regr.out$index=="PDO", 1, 2)
regr.out$index <- reorder(regr.out$index, regr.out$index.order)

png("biol DFA PDO NPGO coefficients.png", 4,4, units="in", res=300)
ggplot(filter(regr.out, class=="biol"), aes(model, est, fill=era)) + geom_hline(yintercept = 0, color="dark grey") + 
  theme_linedraw() +
  geom_bar(position=dodge, stat="identity") + 
  geom_errorbar(aes(ymax=upper, ymin=lower), position=dodge, width=0.5) + xlab("") + 
  ylab("Coefficient") + 
  facet_wrap(~index, nrow=1, scales="free") + 
  theme(legend.position = 'top', axis.text.x = element_text(angle=45, hjust=1)) +
  scale_fill_manual(values=cb[c(2,6)]) 
dev.off()

######################################
# below plots the loadings...

m1 <- fara.sbrd.era1
m2 <- fara.sbrd.era2

p.load <- data.frame(earlyZ=MARSSparamCIs(m1)$par$Z, lateZ=MARSSparamCIs(m2)$par$Z, 
                     earlyUCI=MARSSparamCIs(m1)$par.upCI$Z, lateUCI=MARSSparamCIs(m2)$par.upCI$Z,
                     earlyLCI=MARSSparamCIs(m1)$par.lowCI$Z, lateLCI=MARSSparamCIs(m2)$par.lowCI$Z) 

fsbbd <- ggplot(p.load, aes(earlyZ, lateZ)) +
  theme_linedraw() + 
  geom_abline(slope=1) +
  geom_errorbar(aes(ymin=lateLCI, ymax=lateUCI), color="dark grey", width=0.05) + 
  geom_errorbarh(aes(xmin=earlyLCI, xmax=earlyUCI), color="dark grey") + 
  ylim(range(p.load$lateLCI, p.load$lateUCI, p.load$earlyLCI, p.load$earlyUCI)) +
  xlim(range(p.load$lateLCI, p.load$lateUCI, p.load$earlyLCI, p.load$earlyUCI)) +
  geom_point() +
  ylab("Loading after 1988/89") + 
  xlab("Loading before 1988/89") +
  ggtitle("Farallon seabirds")

m1 <- fara.sbrd.era1
m2 <- fara.sbrd.era2

p.load <- data.frame(earlyZ=MARSSparamCIs(m1)$par$Z, lateZ=MARSSparamCIs(m2)$par$Z, 
                     earlyUCI=MARSSparamCIs(m1)$par.upCI$Z, lateUCI=MARSSparamCIs(m2)$par.upCI$Z,
                     earlyLCI=MARSSparamCIs(m1)$par.lowCI$Z, lateLCI=MARSSparamCIs(m2)$par.lowCI$Z) 

fsbbd <- ggplot(p.load, aes(earlyZ, lateZ)) +
  theme_linedraw() + 
  geom_abline(slope=1) +
  geom_errorbar(aes(ymin=lateLCI, ymax=lateUCI), color="dark grey", width=0.05) + 
  geom_errorbarh(aes(xmin=earlyLCI, xmax=earlyUCI), color="dark grey") + 
  ylim(range(p.load$lateLCI, p.load$lateUCI, p.load$earlyLCI, p.load$earlyUCI)) +
  xlim(range(p.load$lateLCI, p.load$lateUCI, p.load$earlyLCI, p.load$earlyUCI)) +
  geom_point() +
  ylab("Loading after 1988/89") + 
  xlab("Loading before 1988/89") +
  ggtitle("Farallon seabirds")

m1 <- calcof.ichthyo.era1
m2 <- calcof.ichth.era2

p.load <- data.frame(earlyZ=MARSSparamCIs(m1)$par$Z, lateZ=MARSSparamCIs(m2)$par$Z, 
                     earlyUCI=MARSSparamCIs(m1)$par.upCI$Z, lateUCI=MARSSparamCIs(m2)$par.upCI$Z,
                     earlyLCI=MARSSparamCIs(m1)$par.lowCI$Z, lateLCI=MARSSparamCIs(m2)$par.lowCI$Z) 

ich <- ggplot(p.load, aes(earlyZ, lateZ)) +
  theme_linedraw() + 
  geom_abline(slope=1) +
  geom_errorbar(aes(ymin=lateLCI, ymax=lateUCI), color="dark grey", width=0.01) + 
  geom_errorbarh(aes(xmin=earlyLCI, xmax=earlyUCI), color="dark grey") + 
  ylim(range(p.load$lateLCI, p.load$lateUCI, p.load$earlyLCI, p.load$earlyUCI)) +
  xlim(range(p.load$lateLCI, p.load$lateUCI, p.load$earlyLCI, p.load$earlyUCI)) +
  geom_point() +
  ylab("Loading after 1988/89") + 
  xlab("Loading before 1988/89") +
  ggtitle("CalCOFI ichthyoplankton")

m1 <- goa.bio.era1
m2 <- goa.bio.era2

p.load <- data.frame(earlyZ=MARSSparamCIs(m1)$par$Z, lateZ=MARSSparamCIs(m2)$par$Z, 
                     earlyUCI=MARSSparamCIs(m1)$par.upCI$Z, lateUCI=MARSSparamCIs(m2)$par.upCI$Z,
                     earlyLCI=MARSSparamCIs(m1)$par.lowCI$Z, lateLCI=MARSSparamCIs(m2)$par.lowCI$Z) 

g.biol <- ggplot(p.load, aes(earlyZ, lateZ)) +
  theme_linedraw() + 
  geom_abline(slope=1) +
  geom_errorbar(aes(ymin=lateLCI, ymax=lateUCI), color="dark grey", width=0.01) + 
  geom_errorbarh(aes(xmin=earlyLCI, xmax=earlyUCI), color="dark grey") + 
  ylim(range(p.load$lateLCI, p.load$lateUCI, p.load$earlyLCI, p.load$earlyUCI)) +
  xlim(range(p.load$lateLCI, p.load$lateUCI, p.load$earlyLCI, p.load$earlyUCI)) +
  geom_point() +
  ylab("Loading after 1988/89") + 
  xlab("Loading before 1988/89") +
  ggtitle("GOA fish & crustaceans")

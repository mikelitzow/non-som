library(tidyverse)
library(zoo)
library(nlme)
library(MARSS)
library(MuMIn)
library(ggpubr)

# script testing for era differences (pre/post-1988/89) 
# in relationships between regional environmental time series and PDO/NPGO indices 

# there are four regions that we're looking at: CalCOFI/southern CCE, Farallon/central CCE, GOA, and EBS

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


# # load cce env data
# (blanking out my data processing steps here, just load the .csv file below)
# dat <- read.csv("cce.climate.csv")
# 
# # restrict to 1951:2013!
# # (this is the first year for which we have a winter mean for npgo)
# dat <- dat %>%
#   filter(year %in% 1951:2013)
# 
# # and drop the regional climate variables that are winter-spring means
# drop <- grep("win_spr", colnames(dat))
# dat <- dat[,-drop]
# 
# 
# # save in this form
# write.csv(dat, "cce.env.dat.csv")

dat <- read.csv("cce.env.dat.csv", row.names = 1)
# add era term
dat$era <- as.factor(ifelse(dat$year <= 1988, 1, 2))

# and pdo/npgo
dat$pdo <- win.pdo[match(dat$year, names(win.pdo))]
dat$npgo <- win.npgo[match(dat$year, names(win.npgo))]


# response variables
yy <- colnames(dat)[2:7]
system = "Farallon"

# and scale!
ff <- function(x) as.vector(scale(x))
dat[,2:7] <- apply(dat[,2:7], 2, ff)


# output object
mod.out <- data.frame()
for(i in 1:length(yy)){
  # i <- 1
  keep <- colnames(dat)==yy[i]
  temp <- cbind(dat[,keep], dat[,c(8:10)])
  colnames(temp)[1] <- yy[i]
  
  # PDO
  form <- as.formula(paste(colnames(temp)[1]," ~ pdo*era", sep="" ))
  mod <- gls(form, data=temp, correlation = corAR1())

  mod.out <-  rbind(mod.out, 
                    data.frame(system=system, 
                    variable=yy[i],
                    index="PDO",
                    era="early",
                    est=intervals(mod)$coef[2,2], 
                    lower=intervals(mod)$coef[2,1],
                    upper=intervals(mod)$coef[2,3],
                    P=NA),
                    data.frame(system=system, variable=yy[i],
                    index="PDO", 
                    era="late",
                    est=intervals(mod)$coef[2,2] + intervals(mod)$coef[4,2],
                    lower=intervals(mod)$coef[2,2] + intervals(mod)$coef[4,1],
                    upper=intervals(mod)$coef[2,2] + intervals(mod)$coef[4,3],
                    P=summary(mod)$tTable[4,4]
                    )) 
  
  # NPGO
  form <- as.formula(paste(colnames(temp)[1]," ~ npgo*era", sep="" ))
  mod <- gls(form, data=temp, correlation = corAR1())
  
  mod.out <-  rbind(mod.out, 
                    data.frame(system=system, 
                               variable=yy[i],
                               index="NPGO",
                               era="early",
                               est=intervals(mod)$coef[2,2], 
                               lower=intervals(mod)$coef[2,1],
                               upper=intervals(mod)$coef[2,3],
                               P=NA),
                    data.frame(system=system, variable=yy[i],
                               index="NPGO", 
                               era="late",
                               est=intervals(mod)$coef[2,2] + intervals(mod)$coef[4,2],
                               lower=intervals(mod)$coef[2,2] + intervals(mod)$coef[4,1],
                               upper=intervals(mod)$coef[2,2] + intervals(mod)$coef[4,3],
                               P=summary(mod)$tTable[4,4]
                    )) 
}

# plot
dodge <- position_dodge(width=0.9)

# load the color-blind palette
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

mod.out$era <- as.factor(mod.out$era)

ggplot(mod.out, aes(variable, est, fill=era)) + geom_hline(yintercept = 0, color="dark grey") + 
  theme_linedraw() +
  geom_bar(position=dodge, stat="identity") + 
  geom_errorbar(aes(ymax=upper, ymin=lower), position=dodge, width=0.5) + xlab("") + 
  ylab("Coefficient") + 
  facet_wrap(~index, nrow=2, scales="free") + 
  theme(legend.position = 'top', legend.title = element_blank()) +
  scale_fill_manual(values=cb[c(2,6)]) 


# # and CalCOFI environmental data!
# # again, blanking out my data processing steps
# cdat <- read.csv("CCE-nonsalmon-biology&climate.csv", row.names = 1)
# 
# # limit to Southern CCE
# our.dat <- grep("SCC", colnames(cdat))
# 
# env.dat <- cdat[,our.dat]
# head(env.dat)
# 
# # lots of TS w/o good pre-76/77 data - check this
# ff <- function(x) sum(!is.na(x))
# 
# # count of pre-76/77 data 
# keep <- apply(env.dat[rownames(env.dat) <= 1976,], 2, ff) >= 10
# env.dat <- env.dat[rownames(env.dat) <=2012,keep]
# 
# env.dat$year <- rownames(env.dat)
# 
# dat <- env.dat %>%
#   gather(key, value, -year)
# 
# # restrict to 1951:2013!
# dat <- dat %>%
#   filter(year %in% 1951:2013)  %>%
#   spread(key, value)
# 
# # save in this form
# write.csv(dat, "calcofi.env.dat.csv")

# note that Mary is looking at providing a better set of CalCOFI physical data
dat <- read.csv("calcofi.env.dat.csv", row.names = 1)

# add era term
dat$era <- as.factor(ifelse(dat$year <= 1988, 1, 2))

# and pdo/npgo
dat$pdo <- win.pdo[match(dat$year, names(win.pdo))]
dat$npgo <- win.npgo[match(dat$year, names(win.npgo))]


# response variables
yy <- colnames(dat)[2:9]
system = "CalCOFI"

# and scale!
ff <- function(x) as.vector(scale(x))
dat[,2:9] <- apply(dat[,2:9], 2, ff)


# output object

for(i in 1:length(yy)){
  # i <- 1
  keep <- colnames(dat)==yy[i]
  temp <- cbind(dat[,keep], dat[,c(10:12)])
  colnames(temp)[1] <- yy[i]
  
  # PDO
  form <- as.formula(paste(colnames(temp)[1]," ~ pdo*era", sep="" ))
  mod <- gls(form, data=temp, correlation = corAR1(), na.action="na.exclude")
  
  mod.out <-  rbind(mod.out, 
                    data.frame(system=system, 
                               variable=yy[i],
                               index="PDO",
                               era="early",
                               est=intervals(mod)$coef[2,2], 
                               lower=intervals(mod)$coef[2,1],
                               upper=intervals(mod)$coef[2,3],
                               P=NA),
                    data.frame(system=system, variable=yy[i],
                               index="PDO", 
                               era="late",
                               est=intervals(mod)$coef[2,2] + intervals(mod)$coef[4,2],
                               lower=intervals(mod)$coef[2,2] + intervals(mod)$coef[4,1],
                               upper=intervals(mod)$coef[2,2] + intervals(mod)$coef[4,3],
                               P=summary(mod)$tTable[4,4]
                    )) 
  
  # NPGO
  form <- as.formula(paste(colnames(temp)[1]," ~ npgo*era", sep="" ))
  mod <- gls(form, data=temp, correlation = corAR1(), na.action="na.exclude")
  
  mod.out <-  rbind(mod.out, 
                    data.frame(system=system, 
                               variable=yy[i],
                               index="NPGO",
                               era="early",
                               est=intervals(mod)$coef[2,2], 
                               lower=intervals(mod)$coef[2,1],
                               upper=intervals(mod)$coef[2,3],
                               P=NA),
                    data.frame(system=system, variable=yy[i],
                               index="NPGO", 
                               era="late",
                               est=intervals(mod)$coef[2,2] + intervals(mod)$coef[4,2],
                               lower=intervals(mod)$coef[2,2] + intervals(mod)$coef[4,1],
                               upper=intervals(mod)$coef[2,2] + intervals(mod)$coef[4,3],
                               P=summary(mod)$tTable[4,4]
                    )) 
}
mod.out$era <- as.factor(mod.out$era)

ggplot(filter(mod.out, system=="CalCOFI"), aes(variable, est, fill=era)) + geom_hline(yintercept = 0, color="dark grey") + 
  theme_linedraw() +
  geom_bar(position=dodge, stat="identity") + 
  geom_errorbar(aes(ymax=upper, ymin=lower), position=dodge, width=0.5) + xlab("") + 
  ylab("Coefficient") + 
  facet_wrap(~index, nrow=2, scales="free") + 
  theme(legend.position = 'top', legend.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values=cb[c(2,6)]) 


############
# and GOA env
# load data...
# dat <- read.csv("/Users/MikeLitzow 1/Documents/R/time and climes/clim vars for ordination.csv", row.names=1)
# 
# # put years into row names
# rownames(dat) <- dat[,1]
# dat <- dat[,-1]
# 
# # dropping AL and NPI as I think we will leave out these dynamics to keep the paper simple
# # the relevant atmosphere-ocean dynamics have been addequately addressed in PRSB ms. for our purposes!
# dat <- dat[,-c(1,2)]
# 
# head(dat)
# 
# # remove some extra  time series
# dat <- dat[,-c(5,6,8,9,10,11,14,17)]
# 
# 
# # change names to plot-friendly labels
# colnames(dat) <- c("PDO", "NPGO", "SLP gradient", "Freshwater", "Wind stress", "Downwelling", "SST", "Advection", "SSH")
# 
# # and remove PDO/NPGO
# dat <- dat[,-c(1,2)]
# 
# write.csv(dat, "goa.env.dat.csv")

dat <- read.csv("goa.env.dat.csv", row.names = 1)

# add era term
dat$era <- as.factor(ifelse(rownames(dat) <= 1988, 1, 2))

# and pdo/npgo
dat$pdo <- win.pdo[match(rownames(dat), names(win.pdo))]
dat$npgo <- win.npgo[match(rownames(dat), names(win.npgo))]

# clean uo column names

colnames(dat)[1] <- "SLP.gradient"
colnames(dat)[3] <- "Wind.stress"

# response variables
yy <- colnames(dat)[1:7]
system = "GOA"

# and scale!
ff <- function(x) as.vector(scale(x))
dat[,1:7] <- apply(dat[,1:7], 2, ff)

# output object

for(i in 1:length(yy)){
  # i <- 1
  keep <- colnames(dat)==yy[i]
  temp <- cbind(dat[,keep], dat[,c(8:10)])
  colnames(temp)[1] <- yy[i]
  
  # PDO
  form <- as.formula(paste(colnames(temp)[1]," ~ pdo*era", sep="" ))
  mod <- gls(form, data=temp, correlation = corAR1(), na.action="na.exclude")
  
  mod.out <-  rbind(mod.out, 
                    data.frame(system=system, 
                               variable=yy[i],
                               index="PDO",
                               era="early",
                               est=intervals(mod)$coef[2,2], 
                               lower=intervals(mod)$coef[2,1],
                               upper=intervals(mod)$coef[2,3],
                               P=NA),
                    data.frame(system=system, variable=yy[i],
                               index="PDO", 
                               era="late",
                               est=intervals(mod)$coef[2,2] + intervals(mod)$coef[4,2],
                               lower=intervals(mod)$coef[2,2] + intervals(mod)$coef[4,1],
                               upper=intervals(mod)$coef[2,2] + intervals(mod)$coef[4,3],
                               P=summary(mod)$tTable[4,4]
                    )) 
  
  # NPGO
  form <- as.formula(paste(colnames(temp)[1]," ~ npgo*era", sep="" ))
  mod <- gls(form, data=temp, correlation = corAR1(), na.action="na.exclude")
  
  mod.out <-  rbind(mod.out, 
                    data.frame(system=system, 
                               variable=yy[i],
                               index="NPGO",
                               era="early",
                               est=intervals(mod)$coef[2,2], 
                               lower=intervals(mod)$coef[2,1],
                               upper=intervals(mod)$coef[2,3],
                               P=NA),
                    data.frame(system=system, variable=yy[i],
                               index="NPGO", 
                               era="late",
                               est=intervals(mod)$coef[2,2] + intervals(mod)$coef[4,2],
                               lower=intervals(mod)$coef[2,2] + intervals(mod)$coef[4,1],
                               upper=intervals(mod)$coef[2,2] + intervals(mod)$coef[4,3],
                               P=summary(mod)$tTable[4,4]
                    )) 
}

mod.out$era <- as.factor(mod.out$era)

ggplot(filter(mod.out, system=="GOA"), aes(variable, est, fill=era)) + geom_hline(yintercept = 0, color="dark grey") + 
  theme_linedraw() +
  geom_bar(position=dodge, stat="identity") + 
  geom_errorbar(aes(ymax=upper, ymin=lower), position=dodge, width=0.5) + xlab("") + 
  ylab("Coefficient") + 
  facet_wrap(~index, nrow=2, scales="free") + 
  theme(legend.position = 'top', legend.title = element_blank()) +
  scale_fill_manual(values=cb[c(2,6)]) 


########################################################
# EBS environment

# now load EBS data
# setwd("/Users/MikeLitzow 1/Documents/R/climate scripts")
# 
# ewidata <- tibble::as.tibble(read.csv("ewidata.csv", row.names=1))
# setwd("/Users/MikeLitzow 1/Documents/R/FATE2 non-som")
# 
# keep <- grep("AKCLIM_EBS", ewidata$code)
# 
# ewidata <- ewidata[keep,]
# dat = ewidata
# 
# # clean up codes
# temp1 <- strsplit(as.character(dat$code),"AKCLIM_EBS_")
# temp2 <- matrix(unlist(temp1), ncol=2, byrow=TRUE)
# 
# dat$code <- temp2[,2]
# 
# # remove bottom temp as it is too short
# dat <- dat %>%
#   filter(code != "BOTTOM.TEMP") %>%
#   select(year, code, value)
# 
# names(dat)[2] <- "key"
# 
# # restrict to 1951:2013
# dat <- dat %>%
#   filter(year %in% 1951:2013) %>%
#   spread(key, value)
# 
# write.csv(dat, "ebs.env.dat.csv")

dat <- read.csv("ebs.env.dat.csv", row.names = 1)

# add era term
dat$era <- as.factor(ifelse(dat$year <= 1988, 1, 2))

# and pdo/npgo
dat$pdo <- win.pdo[match(dat$year, names(win.pdo))]
dat$npgo <- win.npgo[match(dat$year, names(win.npgo))]

# response variables
yy <- colnames(dat)[2:10]
system = "EBS"

# and scale!
ff <- function(x) as.vector(scale(x))
dat[,2:10] <- apply(dat[,2:10], 2, ff)

# output object

for(i in 1:length(yy)){
  # i <- 1
  keep <- colnames(dat)==yy[i]
  temp <- cbind(dat[,keep], dat[,c(11:13)])
  colnames(temp)[1] <- yy[i]
  
  # PDO
  form <- as.formula(paste(colnames(temp)[1]," ~ pdo*era", sep="" ))
  mod <- gls(form, data=temp, correlation = corAR1(), na.action="na.exclude")
  
  mod.out <-  rbind(mod.out, 
                    data.frame(system=system, 
                               variable=yy[i],
                               index="PDO",
                               era="early",
                               est=intervals(mod)$coef[2,2], 
                               lower=intervals(mod)$coef[2,1],
                               upper=intervals(mod)$coef[2,3],
                               P=NA),
                    data.frame(system=system, variable=yy[i],
                               index="PDO", 
                               era="late",
                               est=intervals(mod)$coef[2,2] + intervals(mod)$coef[4,2],
                               lower=intervals(mod)$coef[2,2] + intervals(mod)$coef[4,1],
                               upper=intervals(mod)$coef[2,2] + intervals(mod)$coef[4,3],
                               P=summary(mod)$tTable[4,4]
                    )) 
  
  # NPGO
  form <- as.formula(paste(colnames(temp)[1]," ~ npgo*era", sep="" ))
  mod <- gls(form, data=temp, correlation = corAR1(), na.action="na.exclude")
  
  mod.out <-  rbind(mod.out, 
                    data.frame(system=system, 
                               variable=yy[i],
                               index="NPGO",
                               era="early",
                               est=intervals(mod)$coef[2,2], 
                               lower=intervals(mod)$coef[2,1],
                               upper=intervals(mod)$coef[2,3],
                               P=NA),
                    data.frame(system=system, variable=yy[i],
                               index="NPGO", 
                               era="late",
                               est=intervals(mod)$coef[2,2] + intervals(mod)$coef[4,2],
                               lower=intervals(mod)$coef[2,2] + intervals(mod)$coef[4,1],
                               upper=intervals(mod)$coef[2,2] + intervals(mod)$coef[4,3],
                               P=summary(mod)$tTable[4,4]
                    )) 
}

mod.out$era <- as.factor(mod.out$era)

ggplot(filter(mod.out, system=="EBS"), aes(variable, est, fill=era)) + geom_hline(yintercept = 0, color="dark grey") + 
  theme_linedraw() +
  geom_bar(position=dodge, stat="identity") + 
  geom_errorbar(aes(ymax=upper, ymin=lower), position=dodge, width=0.5) + xlab("") + 
  ylab("Coefficient") + 
  facet_wrap(~index, nrow=2, scales="free") + 
  theme(legend.position = 'top', legend.title = element_blank()) +
  scale_fill_manual(values=cb[c(2,6)]) 
dev.off()

#############

# now, trying to make sense of these results across the board
# start by just plotting results for time series with 95% CI
# for regression coefficient that does not include zero

mod.out$incl.zero <- NA
  
for(i in seq(2,nrow(mod.out), by=2)){
 # i <- 2
  mod.out$incl.zero[i] <- ifelse(between(0, mod.out[i,6], mod.out[i,7]) ==FALSE |  between(0, mod.out[(i-1),6], mod.out[(i-1),7])==FALSE, FALSE, TRUE)
  mod.out$incl.zero[(i-1)] <- mod.out$incl.zero[i]
}

# reorder regions
mod.out$reg.ord <- ifelse(mod.out$system=="EBS", 1, ifelse(mod.out$system=="GOA",2,3))
mod.out$system <- reorder(mod.out$system, mod.out$reg.ord)

# now I'm plotting these results together with p-values for era differences (in red text)

pdo.plot <- ggplot(filter(mod.out, incl.zero==FALSE, index=="PDO"), aes(variable, est, fill=era)) + geom_hline(yintercept = 0, color="dark grey") + 
  theme_linedraw() +
  geom_bar(position=dodge, stat="identity") + 
  geom_errorbar(aes(ymax=upper, ymin=lower), position=dodge, width=0.5) + xlab("") + 
  ylab("Coefficient") + 
  facet_wrap(~system, ncol=1, scales="free") + 
  theme(legend.position = 'top', legend.title = element_blank()) +
  scale_fill_manual(values=cb[c(2,6)]) +
  geom_text(aes(y=1.5,label = as.character(round(P,3))), color="red") 
            

npgo.plot <- ggplot(filter(mod.out, incl.zero==FALSE, index=="NPGO"), aes(variable, est, fill=era)) + geom_hline(yintercept = 0, color="dark grey") + 
  theme_linedraw() +
  geom_bar(position=dodge, stat="identity") + 
  geom_errorbar(aes(ymax=upper, ymin=lower), position=dodge, width=0.5) + xlab("") + 
  ylab("Coefficient") + 
  facet_wrap(~system, ncol=1, scales="free") + 
  theme(legend.position = 'top', legend.title = element_blank()) +
  scale_fill_manual(values=cb[c(2,6)]) +
  geom_text(aes(y=1,label = as.character(round(P,3))), color="red") 

# and the combined PDO-NPGO plot
png("combined env regression on pdo npgo.png", 8,8, units="in", res=300)
ggarrange(pdo.plot, npgo.plot, labels = c("a) PDO", "b) NPGO"), ncol=2, widths=c(1,0.5))
dev.off()

# so era differences are weak outside GOA
# generally noisy relationships w/ PDO & NPGO
# but it also appears that the era differences are similar across region/index combinations
# i.e., weaker or stronger in a given era
# look at distributions of coefficients between eras...


pdf.pdo <- ggplot(filter(mod.out, index=="PDO"), aes(abs(est), fill=era)) + 
  theme_linedraw() +
  facet_wrap(~system, ncol=1, scales="free_y") +
  geom_density(alpha=0.8) +
  xlim(0,1.3) + 
  scale_fill_manual(values=cb[c(2,6)]) + 
  theme(legend.position = 'top', legend.title = element_blank()) +
  xlab("Regression coefficient (absolute value)") 


pdf.npgo <- ggplot(filter(mod.out, index=="NPGO"), aes(abs(est), fill=era)) + 
  theme_linedraw() +
  facet_wrap(~system, ncol=1, scales="free_y") +
  geom_density(alpha=0.8) +
  xlim(0,1.3) + 
  scale_fill_manual(values=cb[c(2,6)]) + 
  theme(legend.position = 'top', legend.title = element_blank()) +
  xlab("Regression coefficient (absolute value)")  

# this is the plot in the outline
png("env regression density pdo npgo.png", 8,8, units="in", res=300)
ggarrange(pdf.pdo, pdf.npgo, labels = c("a) PDO", "b) NPGO"), ncol=2)
dev.off()

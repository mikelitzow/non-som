library(tidyverse)
library(zoo)
library(nlme)
library(MARSS)
library(MuMIn)
library(ggpubr)
library(reshape2)
library(glmmTMB)
library(rstanarm)
library(lme4)
library(rstan)

# load pdo/npgo 
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

# load three "non-salmon" data sets: GOA crustaceans/fish, Farallon seabirds, CalCOFI ichthyo
dat <- read.csv("farallon.sbrd.biol.csv", row.names = 1)
# add era term
dat$era <- as.factor(ifelse(dat$year <= 1988, 1, 2))

# and pdo/npgo
dat$pdo <- pdo2[match(dat$year, names(pdo2))]
dat$npgo <- npgo2[match(dat$year, names(npgo2))]

# reshape with year, era, and pdo and npgo as the grouping variables
melted <- melt(dat, id.vars = c("year","pdo","era","npgo"))
melted$variable_era = paste0(melted$era,melted$variable)
# standardize all the time series by variable -- so slopes are on same scale
m1 = dplyr::group_by(melted, variable) %>%
  mutate(scale_x = scale(value))
m1$system <- "Farallon seabirds"

####
# CalCOFI
dat <- read.csv("calcofi.biol.csv", row.names=1)

dat$era <- as.factor(ifelse(dat$year <= 1988, 1, 2))

# and pdo/npgo
dat$pdo <- pdo2[match(dat$year, names(pdo2))]
dat$npgo <- npgo2[match(dat$year, names(npgo2))]

# reshape with year, era, and pdo and npgo as the grouping variables
melted <- melt(dat, id.vars = c("year","pdo","era","npgo"))
melted$variable_era = paste0(melted$era,melted$variable)
# standardize all the time series by variable -- so slopes are on same scale
m2 = dplyr::group_by(melted, variable) %>%
  mutate(scale_x = scale(value))
m2$system <- "CalCOFI ichthyoplankton"

#########

dat <- read.csv("goa.biol.csv")
colnames(dat)[1] <- "year"
dat$era <- as.factor(ifelse(dat$year <= 1988, 1, 2))

# and pdo/npgo
dat$pdo <- win.pdo[match(dat$year, names(win.pdo))]
dat$npgo <- win.npgo[match(dat$year, names(win.npgo))]

# reshape with year, era, and pdo and npgo as the grouping variables
melted <- melt(dat, id.vars = c("year","pdo","era","npgo"))
melted$variable_era = paste0(melted$era,melted$variable)
# standardize all the time series by variable -- so slopes are on same scale
m3 = dplyr::group_by(melted, variable) %>%
  mutate(scale_x = scale(value))
m3$system <- "GOA fish & crustaceans"

#######
# combine 

melted <- rbind(m1, m2, m3)
melted$year <- as.numeric(melted$year)
melted$variable <- as.factor(melted$variable)
melted$variable_era <- as.factor(melted$variable_era)

model.data = data.frame()

levels.syst <- as.factor(unique(melted$system))

for(s in levels.syst) {
  
  # s <- levels.syst[1]
  
  temp <- melted %>%
    filter(system==s)
  temp <- na.omit(temp)
  # create data for stan
  stan_data = list(era = as.numeric(temp$era),
                   y = temp$pdo,
                   variable = as.numeric(temp$variable),
                   n = nrow(temp),
                   n_levels = max(as.numeric(temp$variable)),
                   x = temp$scale_x)
  
  mod = stan(file="mod.stan", data=stan_data, chains=3, warmup=4000, iter=6000,thin=2,
             pars = c("beta","mu_beta","ratio","mu_ratio","sigma_beta","sigma_ratio"),
             control=list(adapt_delta=0.99, max_treedepth=20))
  
  pars = rstan::extract(mod,permuted=TRUE)
  
  model.data = rbind(model.data,
                     data.frame(system=s, ratio=100*exp(pars$mu_ratio)))
  
}

# order the systems north-south
model.data$order <- ifelse(model.data$system=="GOA fish & crustaceans", 1, 
                           ifelse(model.data$system=="Farallon seabirds", 2, 3))
model.data$system <- reorder(model.data$system, model.data$order)

cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

pdo.data <- model.data

pdo.plot <- ggplot(pdo.data, aes(ratio/100)) + # just removing % for now
  theme_linedraw() +
  geom_density(fill=cb[3]) + xlab("Avg ratio: Era 1 slope / Era 2 slope") +
  facet_wrap(~system, ncol=1) +
  xlim(c(0,2)) +
  geom_vline(xintercept = 1)+ 
  ggtitle("a) PDO")

#################
## and the same thing for npgo
model.data <- data.frame()
for(s in levels.syst) {
  
  # s <- levels.syst[1]
  
  temp <- melted %>%
    filter(system==s)
  temp <- na.omit(temp)
  # create data for stan
  stan_data = list(era = as.numeric(temp$era),
                   y = temp$npgo,
                   variable = as.numeric(temp$variable),
                   n = nrow(temp),
                   n_levels = max(as.numeric(temp$variable)),
                   x = temp$scale_x)
  
  mod = stan(file="mod.stan", data=stan_data, chains=3, warmup=4000, iter=6000,thin=2,
             pars = c("beta","mu_beta","ratio","mu_ratio","sigma_beta","sigma_ratio"),
             control=list(adapt_delta=0.99, max_treedepth=20))
  
  pars = rstan::extract(mod,permuted=TRUE)
  
  model.data = rbind(model.data,
                     data.frame(system=s, ratio=100*exp(pars$mu_ratio)))
  
}

# order the systems north-south
model.data$order <- ifelse(model.data$system=="GOA fish & crustaceans", 1, 
                           ifelse(model.data$system=="Farallon seabirds", 2, 3))
model.data$system <- reorder(model.data$system, model.data$order)

npgo.data <- model.data

npgo.plot <- ggplot(npgo.data, aes(ratio/100)) +
  theme_linedraw() +
  geom_density(fill=cb[3]) + xlab("Avg ratio: Era 1 slope / Era 2 slope") +
  facet_wrap(~system, ncol=1) +
  geom_vline(xintercept=1) + 
  xlim(c(0,4.5)) +
  ggtitle("b) NPGO")

png("biol regression change pdo-npgo slope.png", 7, 7, units="in", res=300)
ggarrange(pdo.plot, npgo.plot, ncol=2)
dev.off()
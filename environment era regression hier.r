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

# there are five regions that we're looking at: CalCOFI/southern CCE, Farallon/central CCE, northern CCE, GOA, and EBS
# load 'em!
dat <- read.csv("ebs.env.dat.csv", row.names = 1)
# add era term
dat$era <- as.factor(ifelse(dat$year <= 1988, 1, 2))

# and pdo/npgo
dat$pdo <- win.pdo[match(dat$year, names(win.pdo))]
dat$npgo <- win.npgo[match(dat$year, names(win.npgo))]

# reshape with year, era, and pdo and npgo as the grouping variables
melted <- melt(dat, id.vars = c("year","pdo","era","npgo"))
melted$variable_era = paste0(melted$era,melted$variable)
# standardize all the time series by variable -- so slopes are on same scale
m1 = dplyr::group_by(melted, variable) %>%
  mutate(scale_x = scale(value))
m1$system <- "EBS"

####

dat <- read.csv("goa.env.dat.csv", row.names = 1)
dat$year <- rownames(dat)
dat$era <- as.factor(ifelse(dat$year <= 1988, 1, 2))
dat$year <- as.numeric(dat$year)

# and pdo/npgo
dat$pdo <- win.pdo[match(dat$year, names(win.pdo))]
dat$npgo <- win.npgo[match(dat$year, names(win.npgo))]

# reshape with year, era, and pdo and npgo as the grouping variables
melted <- melt(dat, id.vars = c("year","pdo","era","npgo"))
melted$variable_era = paste0(melted$era,melted$variable)
# standardize all the time series by variable -- so slopes are on same scale
m2 = dplyr::group_by(melted, variable) %>%
  mutate(scale_x = scale(value))
m2$system <- "GOA"

#########

dat <- read.csv("cce.env.dat.csv", row.names = 1)

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
m3$system <- "Central CCE"

###################

dat <- read.csv("calcofi.env.dat.csv", row.names = 1)

dat$era <- as.factor(ifelse(dat$year <= 1988, 1, 2))

# and pdo/npgo
dat$pdo <- win.pdo[match(dat$year, names(win.pdo))]
dat$npgo <- win.npgo[match(dat$year, names(win.npgo))]

# reshape with year, era, and pdo and npgo as the grouping variables
melted <- melt(dat, id.vars = c("year","pdo","era","npgo"))
melted$variable_era = paste0(melted$era,melted$variable)
# standardize all the time series by variable -- so slopes are on same scale
m4 = dplyr::group_by(melted, variable) %>%
  mutate(scale_x = scale(value))
m4$system <- "Southern CCE"

###########################

dat <- read.csv("ncc.env.dat.csv", row.names = 1)

dat$era <- as.factor(ifelse(dat$year <= 1988, 1, 2))

# and pdo/npgo
dat$pdo <- win.pdo[match(dat$year, names(win.pdo))]
dat$npgo <- win.npgo[match(dat$year, names(win.npgo))]

# reshape with year, era, and pdo and npgo as the grouping variables
melted <- melt(dat, id.vars = c("year","pdo","era","npgo"))
melted$variable_era = paste0(melted$era,melted$variable)
# standardize all the time series by variable -- so slopes are on same scale
m5 = dplyr::group_by(melted, variable) %>%
  mutate(scale_x = scale(value))
m5$system <- "Northern CCE"

##############

melted <- rbind(m1, m2, m3, m4, m5)
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
                       data.frame(system=s, ratio=pars$mu_ratio))

  }

# order the systems north-south
model.data$order <- ifelse(model.data$system=="EBS", 1, 
                           ifelse(model.data$system=="GOA", 2, 
                                  ifelse(model.data$system=="Northern CCE", 3,
                                         ifelse(model.data$system=="Central CCE", 4, 5))))
model.data$system <- reorder(model.data$system, model.data$order)

cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

pdo.data <- model.data

pdo.plot <- ggplot(pdo.data, aes(100*ratio)) +
  theme_linedraw() +
  geom_density(fill=cb[3]) + xlab("Avg change in slope (%) from Era 1 to Era 2") +
  facet_wrap(~system, ncol=1) +
  geom_vline(xintercept=0) + 
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
                     data.frame(system=s, ratio=pars$mu_ratio))
  
}

# order the systems north-south
model.data$order <- ifelse(model.data$system=="EBS", 1, 
                           ifelse(model.data$system=="GOA", 2, 
                                  ifelse(model.data$system=="Northern CCE", 3,
                                         ifelse(model.data$system=="Central CCE", 4, 5))))
model.data$system <- reorder(model.data$system, model.data$order)

npgo.data <- model.data

npgo.plot <- ggplot(npgo.data, aes(100*ratio)) +
  theme_linedraw() +
  geom_density(fill=cb[3]) + xlab("Avg change in slope (%) from Era 1 to Era 2") +
  facet_wrap(~system, ncol=1) +
  geom_vline(xintercept=0) + 
  ggtitle("b) NPGO")

png("env regression change pdo-npgo slope.png", 7, 9, units="in", res=300)
ggarrange(pdo.plot, npgo.plot, ncol=2)
dev.off()

########################
# everything below is older approach
########################

# trying a new approach:
# include the fixed effect as an overall measure of the degree of association across
# all variables.
# to do this we'll need to align the relationship between each individual variable/time series
# and each index as having a positive slope during the first era.

vars <- dat %>%
  select(-year, -era, -pdo, -npgo)
vars <- colnames(vars)

for(i in 1:length(vars)){
  # i <- 1
  temp <- dat %>%
    select(vars[i], pdo, era) %>%
    filter(era==1)
  form <- as.formula(paste("pdo ~", vars[i], sep=""))
  mod <- lm(form, data=temp, na.action="na.exclude")

  # if slope is negative, multiply time series by -1
  if(summary(mod)$coef[2,1] < 0) {dat[,colnames(dat)==vars[i]] <- -dat[,colnames(dat)==vars[i]]}
}






# fit model with random slopes by era:variable, same intercept
# mod = glmmTMB::glmmTMB(pdo ~ 1 + (-1+z_value)|variable_era, data=melted)

# coefs = data.frame(coef = ranef(mod)$cond$variable_era$z_value)
# coefs$era = substr(rownames(ranef(mod)$cond$variable_era),1,1)
# coefs$variable = substr(rownames(ranef(mod)$cond$variable_era),2,nchar(coefs$coef))
#
# ggplot(coefs, aes(variable, coef, fill=era)) +
#   theme_linedraw() +
#   geom_bar(position="dodge", stat="identity")

# era effect is changing by variable
# mod = glmmTMB::glmmTMB(pdo ~ (-1 + (era:z_value)|variable), data=melted,
                       # control=glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000)))

#mod = rstanarm::stan_lmer(pdo ~ (1+ z_value + era:z_value|variable), data=melted)

# changing to include fixed effects...
mod = rstanarm::stan_lmer(pdo ~ 1 + z_value + era:z_value + (1 + z_value + era:z_value|variable), data=melted)

summary(mod)
fit <- as.data.frame(mod)
head(fit)
plotfit <- data.frame(slope=c(fit[,2], fit[,2] + fit[,3]), era=as.factor(rep(c(1,2), each=nrow(fit))))

ggplot(plotfit, aes(slope, fill=era)) +
  theme_linedraw() +
  geom_density(alpha=0.8) +
  theme(legend.position = c(0.2, 0.8))

ggplot(plotfit, aes(slope, fill=era)) +
  theme_linedraw() +
  geom_density(alpha=0.8) +
  scale_fill_manual(values=cb[c(2,6)]) +
  theme(legend.position = 'top', legend.title = element_blank()) +
  xlab("Regression coefficient (absolute value)")

# make plot-friendly label names
names <- c("Spring SST", "Winter SST", "Summer NW wind", "Summer SE wind", "Winter NW wind", "Winter SE wind",
           "Ice 57.75N 164.5W", "Ice 58.75N 164.5W", "Ice 64.25N 163.5W")

# make a df to hold results for plotting
plot <- data.frame(coef=c(coef(mod)$variable[,1], coef(mod)$variable[,1]+coef(mod)$variable[,2]), names=names,
                   era=rep(c("Before 1988/89", "After 1988/89"), each=nrow(coef(mod)$variable)))

plot$names <- reorder(plot$names, rep(plot$coef[plot$era=="Before 1988/89"],2))
plot$index <- "PDO"
plot$system <- "EBS"

ggplot(plot, aes(names, coef, fill=era)) +
  theme_linedraw() +
  geom_bar(position="dodge", stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title.x = element_blank(),
        legend.position = c(0.2, 0.8), legend.title = element_blank())

# now npgo
mod = rstanarm::stan_lmer(npgo ~ (1+ z_value + era:z_value|variable), data=melted)
temp <- data.frame(coef=c(coef(mod)$variable[,1], coef(mod)$variable[,1]+coef(mod)$variable[,2]), names=names,
                   era=rep(c("Before 1988/89", "After 1988/89"), each=nrow(coef(mod)$variable)))

temp$names <- reorder(plot$names, rep(plot$coef[plot$era=="Before 1988/89"],2))
temp$index <- "NPGO"
temp$system <- "EBS"

plot <- rbind(plot, temp)

# and goa
system <- "GOA"
dat <- read.csv("goa.env.dat.csv", row.names = 1)
names <- colnames(dat)[1:7]

dat$year <- rownames(dat)
# add era term
dat$era <- as.factor(ifelse(dat$year <= 1988, 1, 2))

# and pdo/npgo
dat$pdo <- win.pdo[match(dat$year, names(win.pdo))]
dat$npgo <- win.npgo[match(dat$year, names(win.npgo))]

# reshape with year, era, and pdo and npgo as the grouping variables
melted <- melt(dat, id.vars = c("year","pdo","era","npgo"))
melted$variable_era = paste0(melted$era,melted$variable)


# # standardize all the time series?
melted = dplyr::group_by(melted, variable) %>%
  mutate(z_value = scale(value))

#pdo model
mod = rstanarm::stan_lmer(pdo ~ (1+ z_value + era:z_value|variable), data=melted)
temp <- data.frame(coef=c(coef(mod)$variable[,1], coef(mod)$variable[,1]+coef(mod)$variable[,2]), names=names,
                   era=rep(c("Before 1988/89", "After 1988/89"), each=nrow(coef(mod)$variable)))

temp$index <- "PDO"
temp$system <- system
plot <- rbind(plot, temp)

# #npgo model
mod = rstanarm::stan_lmer(npgo ~ (1+ z_value + era:z_value|variable), data=melted)
temp <- data.frame(coef=c(coef(mod)$variable[,1], coef(mod)$variable[,1]+coef(mod)$variable[,2]), names=names,
                   era=rep(c("Before 1988/89", "After 1988/89"), each=nrow(coef(mod)$variable)))

temp$index <- "NPGO"
temp$system <- system
plot <- rbind(plot, temp)

###############
# Farallon region
system <- "Farallon"
dat <- read.csv("cce.env.dat.csv", row.names = 1)
names <- c("Spring upwelling", "Winter upwelling", "Spring SST", "Winter SST", "Spring SSH", "Winter SSH")
# add era term
dat$era <- as.factor(ifelse(dat$year <= 1988, 1, 2))

# add era term
dat$era <- as.factor(ifelse(dat$year <= 1988, 1, 2))

# and pdo/npgo
dat$pdo <- win.pdo[match(dat$year, names(win.pdo))]
dat$npgo <- win.npgo[match(dat$year, names(win.npgo))]

# reshape with year, era, and pdo and npgo as the grouping variables
melted <- melt(dat, id.vars = c("year","pdo","era","npgo"))
melted$variable_era = paste0(melted$era,melted$variable)


# # standardize all the time series?
melted = dplyr::group_by(melted, variable) %>%
  mutate(z_value = scale(value))

#pdo model
mod = rstanarm::stan_lmer(pdo ~ (1+ z_value + era:z_value|variable), data=melted)
temp <- data.frame(coef=c(coef(mod)$variable[,1], coef(mod)$variable[,1]+coef(mod)$variable[,2]), names=names,
                   era=rep(c("Before 1988/89", "After 1988/89"), each=nrow(coef(mod)$variable)))

temp$index <- "PDO"
temp$system <- system
plot <- rbind(plot, temp)

# #npgo model
mod = rstanarm::stan_lmer(npgo ~ (1+ z_value + era:z_value|variable), data=melted)
temp <- data.frame(coef=c(coef(mod)$variable[,1], coef(mod)$variable[,1]+coef(mod)$variable[,2]), names=names,
                   era=rep(c("Before 1988/89", "After 1988/89"), each=nrow(coef(mod)$variable)))

temp$index <- "NPGO"
temp$system <- system
plot <- rbind(plot, temp)

# now calcofi
system <- "CalCOFI"
dat <- read.csv("calcofi.env.dat.csv", row.names = 1)
names <- colnames(dat)[2:9]

# add era term
dat$era <- as.factor(ifelse(dat$year <= 1988, 1, 2))

# and pdo/npgo
dat$pdo <- win.pdo[match(dat$year, names(win.pdo))]
dat$npgo <- win.npgo[match(dat$year, names(win.npgo))]

# reshape with year, era, and pdo and npgo as the grouping variables
melted <- melt(dat, id.vars = c("year","pdo","era","npgo"))
melted$variable_era = paste0(melted$era,melted$variable)


# # standardize all the time series?
melted = dplyr::group_by(melted, variable) %>%
  mutate(z_value = scale(value))

#pdo model
mod = rstanarm::stan_lmer(pdo ~ (1+ z_value + era:z_value|variable), data=melted)
temp <- data.frame(coef=c(coef(mod)$variable[,1], coef(mod)$variable[,1]+coef(mod)$variable[,2]), names=names,
                   era=rep(c("Before 1988/89", "After 1988/89"), each=nrow(coef(mod)$variable)))

temp$index <- "PDO"
temp$system <- system
plot <- rbind(plot, temp)

# #npgo model
mod = rstanarm::stan_lmer(npgo ~ (1+ z_value + era:z_value|variable), data=melted)
temp <- data.frame(coef=c(coef(mod)$variable[,1], coef(mod)$variable[,1]+coef(mod)$variable[,2]), names=names,
                   era=rep(c("Before 1988/89", "After 1988/89"), each=nrow(coef(mod)$variable)))

temp$index <- "NPGO"
temp$system <- system
plot <- rbind(plot, temp)


#################
# and combined plot
dodge <- position_dodge(width=0.9)
# load the color-blind palette
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

plot$order <- ifelse(plot$era=="Before 1988/89", 1, 2)
plot$names <- reorder(plot$names, plot$order)


pdo.plot <- ggplot(filter(plot, index=="PDO"), aes(names, coef, fill=era)) + geom_hline(yintercept = 0, color="dark grey") +
  theme_linedraw() +
  geom_bar(position=dodge, stat="identity") +
  xlab("") +
  ylab("Coefficient") +
  facet_wrap(~system, ncol=1, scales="free") +
  theme(legend.position = 'top', legend.title = element_blank()) +
  scale_fill_manual(values=cb[c(2,6)])


npgo.plot <- ggplot(filter(mod.out, incl.zero==FALSE, index=="NPGO"), aes(variable, est, fill=era)) + geom_hline(yintercept = 0, color="dark grey") +
  theme_linedraw() +
  geom_bar(position=dodge, stat="identity") +
  geom_errorbar(aes(ymax=upper, ymin=lower), position=dodge, width=0.5) + xlab("") +
  ylab("Coefficient") +
  facet_wrap(~system, ncol=1, scales="free") +
  theme(legend.position = 'top', legend.title = element_blank()) +
  scale_fill_manual(values=cb[c(2,6)]) +
  geom_text(aes(y=1,label = as.character(round(P,3))), color="red")

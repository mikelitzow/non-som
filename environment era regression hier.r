dat <- read.csv("ebs.env.dat.csv", row.names = 1)

# add era term
dat$era <- as.factor(ifelse(dat$year <= 1988, 1, 2))

# and pdo/npgo
dat$pdo <- win.pdo[match(dat$year, names(win.pdo))]
dat$npgo <- win.npgo[match(dat$year, names(win.npgo))]

# EW added this code for the hierarchical model
library(reshape2)
library(glmmTMB)
library(rstanarm)
library(lme4)

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
  

# reshape with year, era, and pdo and npgo as the grouping variables
melted <- melt(dat, id.vars = c("year","pdo","era","npgo"))
melted$variable_era = paste0(melted$era,melted$variable)
# # standardize all the time series? 
melted = dplyr::group_by(melted, variable) %>% 
  mutate(z_value = scale(value))

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
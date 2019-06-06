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

dat <- read.csv("calcofi.phys.gam.csv", row.names = 1)

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
                     data.frame(system=s, ratio=100*exp(pars$mu_ratio)))
  
}

# order the systems north-south
model.data$order <- ifelse(model.data$system=="EBS", 1, 
                           ifelse(model.data$system=="GOA", 2, 
                                  ifelse(model.data$system=="Northern CCE", 3,
                                         ifelse(model.data$system=="Central CCE", 4, 5))))
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
model.data$order <- ifelse(model.data$system=="EBS", 1, 
                           ifelse(model.data$system=="GOA", 2, 
                                  ifelse(model.data$system=="Northern CCE", 3,
                                         ifelse(model.data$system=="Central CCE", 4, 5))))
model.data$system <- reorder(model.data$system, model.data$order)

npgo.data <- model.data

npgo.plot <- ggplot(npgo.data, aes(ratio/100)) +
  theme_linedraw() +
  geom_density(fill=cb[3]) + xlab("Avg ratio: Era 1 slope / Era 2 slope") +
  facet_wrap(~system, ncol=1) +
  geom_vline(xintercept=1) + 
  xlim(c(0,4.5)) +
  ggtitle("b) NPGO")

png("env regression change pdo-npgo slope.png", 7, 9, units="in", res=300)
ggarrange(pdo.plot, npgo.plot, ncol=2)
dev.off()

# Caterpillar Plot ===============================
# Helper Functions
q.50 <- function(x) { return(quantile(x, probs=c(0.25,0.75))) }
q.95 <- function(x) { return(quantile(x, probs=c(0.025,0.975))) }

head(npgo.data)
head(pdo.data)

# Combine dataframes
npgo.data$var <- "NPGO"
pdo.data$var <- "PDO"
all.data <- rbind(pdo.data, npgo.data)

cat.plt <- ggplot(all.data, aes(x=system, y=ratio/100, fill=system)) +
  theme_linedraw() +
  # scale_fill_colorblind() +
  scale_fill_tableau() +
  # scale_fill_brewer(c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
  # geom_eye() +
  
  geom_violin(alpha = 0.75, lwd=0.1, scale='width') +
  stat_summary(fun.y="q.95", colour="black", geom="line", lwd=0.75) +
  stat_summary(fun.y="q.50", colour="black", geom="line", lwd=1.5) +
  stat_summary(fun.y="median", colour="black", size=2, geom="point", pch=21) +
  facet_wrap(~var, ncol=1) +
  ylab("Avg ratio: Era 1 slope / Era 2 slope") +
  theme(axis.text.y = element_blank()) +
  coord_flip(ylim=c(0,7))

cat.plt
ggsave("env regression change pdo-npgo slope_cater.png", plot=cat.plt, 
       height=7, width=7, units="in", dpi=300)

# Separate by 

cat.plt.pdo <- all.data %>% filter(var=='PDO') %>% ggplot(aes(x=system, y=ratio/100, fill=system)) +
  theme_linedraw() +
  scale_fill_tableau() +
  # geom_eye() +
  geom_violin(alpha = 0.75, lwd=0.1, scale='width') +
  stat_summary(fun.y="q.95", colour="black", geom="line", lwd=0.75) +
  stat_summary(fun.y="q.50", colour="black", geom="line", lwd=1.5) +
  stat_summary(fun.y="median", colour="black", size=2, geom="point", pch=21) +
  facet_wrap(~var, ncol=1) +
  ylab("Avg ratio: Era 1 slope / Era 2 slope") +
  theme(axis.text.y = element_blank(), legend.position='top') +
  coord_flip(ylim=c(0,3))

cat.plt.npgo <- all.data %>% filter(var=='NPGO') %>% ggplot(aes(x=system, y=ratio/100, fill=system)) +
  theme_linedraw() +
  scale_fill_tableau() +
  # geom_eye() +
  geom_violin(alpha = 0.75, lwd=0.1, scale='width') +
  stat_summary(fun.y="q.95", colour="black", geom="line", lwd=0.75) +
  stat_summary(fun.y="q.50", colour="black", geom="line", lwd=1.5) +
  stat_summary(fun.y="median", colour="black", size=2, geom="point", pch=21) +
  facet_wrap(~var, ncol=1) +
  ylab("Avg ratio: Era 1 slope / Era 2 slope") +
  theme(axis.text.y = element_blank(), legend.position="none") +
  coord_flip(ylim=c(0,6))

# Plot Combined with Sepearte 
cat.plt.2 <- plot_grid(cat.plt.pdo, cat.plt.npgo, ncol=1, rel_heights = c(1.1,1))
ggsave("env regression change pdo-npgo slope_cater2.png", plot=cat.plt.2, 
       height=7, width=7, units="in", dpi=300)
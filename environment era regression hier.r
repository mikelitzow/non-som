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

# reshape with year, era, and pdo and npgo as the grouping variables
melted <- melt(dat, id.vars = c("year","pdo","era","npgo"))
melted$variable_era = paste0(melted$era,melted$variable)
# # standardize all the time series? 
melted = dplyr::group_by(melted, variable) %>% 
  mutate(z_value = scale(value))

# fit model with random slopes by era:variable, same intercept
mod = glmmTMB::glmmTMB(pdo ~ 1 + (-1+z_value)|variable_era, data=melted)


coefs = data.frame(coef = ranef(mod)$cond$variable_era$z_value)
coefs$era = substr(rownames(ranef(mod)$cond$variable_era),1,1)
coefs$variable = substr(rownames(ranef(mod)$cond$variable_era),2,nchar(coefs$coef))

ggplot(coefs, aes(coef, fill=era)) + 
  theme_linedraw() +
  geom_density(alpha=0.8)

# era effect is changing by variable
mod = glmmTMB::glmmTMB(pdo ~ (-1 + (era:z_value)|variable), data=melted,
                       control=glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000)))

mod = rstanarm::stan_lmer(pdo ~ (1+ z_value + era:z_value|variable), data=melted)

coef(mod)
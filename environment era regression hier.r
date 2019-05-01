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

# make a df to hold results for plotting
plot <- as.data.frame(coef(mod)$variable[,1:2]) %>%
  gather()

plot$era <- ifelse(plot$key=="z_value", "Era.1", "Era.2")

# make plot-friendly label names
names <- c("Spring SST", "Winter SST", "Summer NW wind", "Summer SE wind", "Winter NW wind", "Winter SE wind",
           "Ice 57.75N 164.5W", "Ice 58.75N 164.5W", "Ice 64.25N 163.5W")


plot$names <- names
plot$names <- reorder(plot$names, rep(plot$value[plot$era=="Era.1"],2))


ggplot(plot, aes(value, fill=era)) + 
  theme_linedraw() +
  geom_density(alpha=0.8)

ggplot(plot, aes(names, value, fill=era)) + 
  theme_linedraw() +
  geom_bar(position="dodge", stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title.x = element_blank(),
        legend.position = c(0.2, 0.8), legend.title = element_blank())


# compare with lme results
mod <- lmer(pdo ~ -1 + era:z_value + (-1 + era:z_value | variable),  data=melted)
ranef(mod)

# combine and plot
plot$method <- "stan_lmer"

plot2 <- as.data.frame(coef(mod)$variable) %>%
  gather()

plot2$era <- ifelse(plot2$key=="era1:z_value", "Era.1", "Era.2")
plot2$names <- names
plot2$method <- "lmer"

plot <- rbind(plot, plot2)

png("method comparison", 4, 4, units="in", res=300)
ggplot(plot, aes(names, value, fill=era)) + 
  theme_linedraw() +
  geom_bar(position="dodge", stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title.x = element_blank(),
        legend.position = c(0.1, 0.9), legend.title = element_blank()) +
  facet_wrap(~method, ncol=1)
dev.off()

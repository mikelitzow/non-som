
library(MuMIn)
library(dplyr)
library(ggplot2)
library(nlme)
library(lmtest)
library(tidyverse)
library(ggpubr)
options(max.print = 99999)

#run.dat <- read.csv("/Users/MikeLitzow 1/Documents/R/coastwide-salmon/coastwide salmon data with different sst groupings.csv", row.names = 1)

#write.csv(run.dat, "salmon run dat.csv")

run.dat <- read.csv("salmon run dat.csv", row.names = 1)

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

# smooth over various scales
# winter before and winter of ocean entry, winter of and winter after ocean entry,
# and winter before, winter of, and winter after ocean entry
npgo2 <- rollmean(win.npgo, 2, align="right", fill=NA)
npgo2a <- rollmean(win.npgo, 2, align="left", fill=NA)
npgo3 <- rollmean(win.npgo, 3, fill=NA)

pdo2 <- rollmean(win.pdo, 2, align="right", fill=NA)
pdo2a <- rollmean(win.pdo, 2, align="left", fill=NA)
pdo3 <- rollmean(win.pdo, 3, fill=NA)

run.dat$era <- as.factor(run.dat$era)

run.dat$pdo <- win.pdo[match(run.dat$entry.yr, names(win.pdo))]
run.dat$pdo2 <- pdo2[match(run.dat$entry.yr, names(pdo2))]
run.dat$pdo2a <- pdo2a[match(run.dat$entry.yr, names(pdo2a))]
run.dat$pdo3 <- pdo3[match(run.dat$entry.yr, names(pdo3))]

run.dat$npgo <- win.npgo[match(run.dat$entry.yr, names(win.npgo))]
run.dat$npgo2 <- npgo2[match(run.dat$entry.yr, names(npgo2))]
run.dat$npgo2a <- npgo2a[match(run.dat$entry.yr, names(npgo2a))]
run.dat$npgo3 <- npgo3[match(run.dat$entry.yr, names(npgo3))]

# start with Pink

# compare three versions of pdo

# random intercept for stocks
p1 <- lme(log(recruits/spawners) ~ 1 + stock:spawners +  pdo:region + pdo:region:era, random = ~1 | stock,
method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Pink"))

p2 <- lme(log(recruits/spawners) ~ 1 + stock:spawners +  pdo2:region + pdo2:region:era, random = ~1 | stock,
method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Pink"))

p2a <- lme(log(recruits/spawners) ~ 1 + stock:spawners +  pdo2a:region + pdo2a:region:era, random = ~1 | stock,
          method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Pink"))

p3 <- lme(log(recruits/spawners) ~ 1 + stock:spawners +  pdo3:region + pdo3:region:era, random = ~1 | stock,
method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Pink"))

AIC <- AICc(p1, p2, p2a, p3)
AIC <- AIC %>%
mutate(model=rownames(AIC)) %>%
mutate(dAICc=AICc-min(AICc)) %>%
arrange(AICc)
AIC

# so pdo winter before/after entry is the best
# compare with the stationary model

p.red <- lme(log(recruits/spawners) ~ 1 + stock:spawners +  pdo2a:region, random = ~1 | stock,
method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Pink"))

anova(p2a, p.red)

# save best model to plot
pink.pdo <- p2a


# now pink NPGO
p1 <- lme(log(recruits/spawners) ~ 1 + stock:spawners +  npgo:region + npgo:region:era, random = ~1 | stock,
          method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Pink"))

p2 <- lme(log(recruits/spawners) ~ 1 + stock:spawners +  npgo2:region + npgo2:region:era, random = ~1 | stock,
          method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Pink"))

p2a <- lme(log(recruits/spawners) ~ 1 + stock:spawners +  npgo2a:region + npgo2a:region:era, random = ~1 | stock,
          method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Pink"))

p3 <- lme(log(recruits/spawners) ~ 1 + stock:spawners +  npgo3:region + npgo3:region:era, random = ~1 | stock,
          method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Pink"))

AIC <- AICc(p1, p2, p2a, p3)

AIC <- AIC %>%
  mutate(model=rownames(AIC)) %>%
  mutate(dAICc=AICc-min(AICc)) %>%
  arrange(AICc)
AIC

# so npgo2a is the best
# compare with the stationary model

p.red <- lme(log(recruits/spawners) ~ 1 + stock:spawners +  npgo2a:region, random = ~1 | stock,
             method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Pink"))

anova(p2a, p.red)

pink.npgo <- p2a


######
# now sockeye

# compare three versions of pdo

# random intercept for stocks
p1 <- lme(log(recruits/spawners) ~ 1 + stock:spawners +  pdo:region + pdo:region:era, random = ~1 | stock,
          method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Sockeye"))

p2 <- lme(log(recruits/spawners) ~ 1 + stock:spawners +  pdo2:region + pdo2:region:era, random = ~1 | stock,
          method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Sockeye"))

p2a <- lme(log(recruits/spawners) ~ 1 + stock:spawners +  pdo2a:region + pdo2a:region:era, random = ~1 | stock,
          method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Sockeye"))

p3 <- lme(log(recruits/spawners) ~ 1 + stock:spawners +  pdo3:region + pdo3:region:era, random = ~1 | stock,
          method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Sockeye"))

AIC <- AICc(p1, p2, p2a, p3)
AIC <- AIC %>%
  mutate(model=rownames(AIC)) %>%
  mutate(dAICc=AICc-min(AICc)) %>%
  arrange(AICc)
AIC

# so pdo2 is the best
# compare with the stationary model

p.red <- lme(log(recruits/spawners) ~ 1 + stock:spawners +  pdo2:region, random = ~1 | stock,
             method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Sockeye"))

anova(p2, p.red)

sock.pdo <- p2

# now sockeyee NPGO
p1 <- lme(log(recruits/spawners) ~ 1 + stock:spawners +  npgo:region + npgo:region:era, random = ~1 | stock,
          method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Sockeye"))

p2 <- lme(log(recruits/spawners) ~ 1 + stock:spawners +  npgo2:region + npgo2:region:era, random = ~1 | stock,
          method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Sockeye"))

p2a <- lme(log(recruits/spawners) ~ 1 + stock:spawners +  npgo2a:region + npgo2a:region:era, random = ~1 | stock,
          method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Sockeye"))

p3 <- lme(log(recruits/spawners) ~ 1 + stock:spawners +  npgo3:region + npgo3:region:era, random = ~1 | stock,
          method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Sockeye"))

AIC <- AICc(p1, p2, p2a, p3)

AIC <- AIC %>%
  mutate(model=rownames(AIC)) %>%
  mutate(dAICc=AICc-min(AICc)) %>%
  arrange(AICc)
AIC

# so npgo1 is the best
sock.npgo <- p1

# compare with the stationary model

p.red <- lme(log(recruits/spawners) ~ 1 + stock:spawners +  npgo:region, random = ~1 | stock,
             method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Sockeye"))

anova(p1, p.red)

############
# and chum!

# compare three versions of pdo

# random intercept for stocks
p1 <- lme(log(recruits/spawners) ~ 1 + stock:spawners +  pdo:region + pdo:region:era, random = ~1 | stock,
          method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Chum"))

p2 <- lme(log(recruits/spawners) ~ 1 + stock:spawners +  pdo2:region + pdo2:region:era, random = ~1 | stock,
          method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Chum"))

p2a <- lme(log(recruits/spawners) ~ 1 + stock:spawners +  pdo2a:region + pdo2a:region:era, random = ~1 | stock,
          method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Chum"))

p3 <- lme(log(recruits/spawners) ~ 1 + stock:spawners +  pdo3:region + pdo3:region:era, random = ~1 | stock,
          method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Chum"))

AIC <- AICc(p1, p2, p2a, p3)
AIC <- AIC %>%
  mutate(model=rownames(AIC)) %>%
  mutate(dAICc=AICc-min(AICc)) %>%
  arrange(AICc)
AIC

# so pdo2a is the best
chum.pdo <- p2a
# compare with the stationary model

p.red <- lme(log(recruits/spawners) ~ 1 + stock:spawners +  pdo2a:region, random = ~1 | stock,
             method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Chum"))

anova(p2a, p.red)


# now Chum NPGO
p1 <- lme(log(recruits/spawners) ~ 1 + stock:spawners +  npgo:region + npgo:region:era, random = ~1 | stock,
          method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Chum"))

p2 <- lme(log(recruits/spawners) ~ 1 + stock:spawners +  npgo2:region + npgo2:region:era, random = ~1 | stock,
          method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Chum"))

p2a <- lme(log(recruits/spawners) ~ 1 + stock:spawners +  npgo2a:region + npgo2a:region:era, random = ~1 | stock,
          method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Chum"))

p3 <- lme(log(recruits/spawners) ~ 1 + stock:spawners +  npgo3:region + npgo3:region:era, random = ~1 | stock,
          method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Chum"))

AIC <- AICc(p1, p2, p2a, p3)
AIC <- AIC %>%
  mutate(model=rownames(AIC)) %>%
  mutate(dAICc=AICc-min(AICc)) %>%
  arrange(AICc)
AIC

# so npgo3 is the best
chum.npgo <- p3
# compare with the stationary model

p.red <- lme(log(recruits/spawners) ~ 1 + stock:spawners +  npgo3:region, random = ~1 | stock,
             method = "ML", correlation=corAR1(form= ~ 1 | stock), data = filter(run.dat, species=="Chum"))

anova(p3, p.red)

cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

lm <- 0.1
png("aicc for salmon SR.png", 2,3, units="in", res=300)
ggplot(plot.out, aes(species, AICc, fill=index)) +
  geom_bar(position="dodge", stat="identity", color="dark grey") +
  theme_linedraw() +
  ylim(c(-max(plot.out$AICc), max(plot.out$AICc))) +
  theme(axis.title.x = element_blank(), legend.title = element_blank(),
        legend.position = c(0.65, 0.9), legend.text = element_text(size=8), legend.key.size = unit(0.4, "cm"),
        legend.margin = margin(lm,lm,lm,lm,"cm")) +
  geom_hline(yintercept = 0, color="dark grey") +
  ylab("delta-AICc") +
  scale_fill_manual(values=cb[c(3,7)])
dev.off()

## plot the regression coefficients
plot <- data.frame(species=rep(c("Pink", "Sockeye", "Chum"), each=12), region=rep(c("EBS", "GOA", "South"), each=4), index=rep(c("PDO", "NPGO"), each=2),
                   era=c("Before 1988/89", "After 1988/89"), LCI=NA, effect=NA, UCI=NA)

plot[1,5:7] <- intervals(pink.pdo, which="fixed")$fixed[37,] # pink PDO EBS early
plot[2,5:7] <- intervals(pink.pdo, which="fixed")$fixed[37,2]+
  intervals(pink.pdo, which="fixed")$fixed[40,] # pink PDO EBS late

plot[3,5:7] <- intervals(pink.npgo, which="fixed")$fixed[37,] # pink NPGO EBS early
plot[4,5:7] <- intervals(pink.npgo, which="fixed")$fixed[37,2]+
  intervals(pink.npgo, which="fixed")$fixed[40,] # pink NPGO EBS late

plot[5,5:7] <- intervals(pink.pdo, which="fixed")$fixed[38,] # pink PDO GOA early
plot[6,5:7] <- intervals(pink.pdo, which="fixed")$fixed[38,2]+
  intervals(pink.pdo, which="fixed")$fixed[41,] # pink PDO GOA late

plot[7,5:7] <- intervals(pink.npgo, which="fixed")$fixed[38,] # pink NPGO GOA early
plot[8,5:7] <- intervals(pink.npgo, which="fixed")$fixed[38,2]+
  intervals(pink.npgo, which="fixed")$fixed[41,] # pink NPGO GOA late

plot[9,5:7] <- intervals(pink.pdo, which="fixed")$fixed[39,] # pink PDO South early
plot[10,5:7] <- intervals(pink.pdo, which="fixed")$fixed[39,2]+
  intervals(pink.pdo, which="fixed")$fixed[42,] # pink PDO South late

plot[11,5:7] <- intervals(pink.npgo, which="fixed")$fixed[39,] # pink PDO South early
plot[12,5:7] <- intervals(pink.npgo, which="fixed")$fixed[39,2]+
  intervals(pink.npgo, which="fixed")$fixed[42,] # pink PDO South late

# now sockeye
plot[13,5:7] <- intervals(sock.pdo, which="fixed")$fixed[32,] # sockeye PDO EBS early
plot[14,5:7] <- intervals(sock.pdo, which="fixed")$fixed[32,2] +
  intervals(sock.pdo, which="fixed")$fixed[35,] # sockeye PDO EBS late

plot[15,5:7] <- intervals(sock.npgo, which="fixed")$fixed[32,] # sockeye NPGO EBS early
plot[16,5:7] <- intervals(sock.npgo, which="fixed")$fixed[32,2] +
  intervals(sock.npgo, which="fixed")$fixed[35,] # sockeye NPGO EBS late

plot[17,5:7] <- intervals(sock.pdo, which="fixed")$fixed[33,] # sockeye PDO GOA early
plot[18,5:7] <- intervals(sock.pdo, which="fixed")$fixed[33,2] +
  intervals(sock.pdo, which="fixed")$fixed[36,] # sockeye PDO GOA late

plot[19,5:7] <- intervals(sock.npgo, which="fixed")$fixed[33,] # sockeye NPGO GOA early
plot[20,5:7] <- intervals(sock.npgo, which="fixed")$fixed[33,2] +
  intervals(sock.npgo, which="fixed")$fixed[36,] # sockeye NPGO GOA late

plot[21,5:7] <- intervals(sock.pdo, which="fixed")$fixed[34,] # sockeye PDO South early
plot[22,5:7] <- intervals(sock.pdo, which="fixed")$fixed[34,2] +
  intervals(sock.pdo, which="fixed")$fixed[37,] # sockeye PDO South late

plot[23,5:7] <- intervals(sock.npgo, which="fixed")$fixed[34,] # sockeye NPGO South early
plot[24,5:7] <- intervals(sock.npgo, which="fixed")$fixed[34,2] +
  intervals(sock.npgo, which="fixed")$fixed[37,] # sockeye NPGO South late

# and chum!
plot[25,5:7] <- intervals(chum.pdo, which="fixed")$fixed[23,] # chum PDO EBS early
plot[26,5:7] <- intervals(chum.pdo, which="fixed")$fixed[23,2] +
  intervals(chum.pdo, which="fixed")$fixed[26,] # chum PDO EBS late

plot[27,5:7] <- intervals(chum.npgo, which="fixed")$fixed[23,] # chum NPGO EBS early
plot[28,5:7] <- intervals(chum.npgo, which="fixed")$fixed[23,2] +
  intervals(chum.npgo, which="fixed")$fixed[26,] # chum NPGO EBS late

plot[29,5:7] <- intervals(chum.pdo, which="fixed")$fixed[24,] # chum PDO GOA early
plot[30,5:7] <- intervals(chum.pdo, which="fixed")$fixed[24,2] +
  intervals(chum.pdo, which="fixed")$fixed[27,] # chum PDO GOA late

plot[31,5:7] <- intervals(chum.npgo, which="fixed")$fixed[24,] # chum NPGO GOA early
plot[32,5:7] <- intervals(chum.npgo, which="fixed")$fixed[24,2] +
  intervals(chum.npgo, which="fixed")$fixed[27,] # chum NPGO GOA late

plot[33,5:7] <- intervals(chum.pdo, which="fixed")$fixed[25,] # chum PDO South early
plot[34,5:7] <- intervals(chum.pdo, which="fixed")$fixed[25,2] +
  intervals(chum.pdo, which="fixed")$fixed[28,] # chum PDO South late

plot[35,5:7] <- intervals(chum.npgo, which="fixed")$fixed[25,] # chum NPGO South early
plot[36,5:7] <- intervals(chum.npgo, which="fixed")$fixed[25,2] +
  intervals(chum.npgo, which="fixed")$fixed[28,] # chum NPGO South late

plot$era <- reorder(plot$era, -order(plot$era))

plot$species <- reorder(plot$species, rep(c(1,2,3), each=12))
plot$plot.era <- as.factor(ifelse(plot$era=="Before 1988/89", 1, 2))

dodge <- position_dodge(width=0.9)

# load the color-blind palette
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


png("salmon SR pdo npgo.png", 6,8, units="in", res=300)
ggplot(plot, aes(region, effect, fill=plot.era)) + geom_hline(yintercept = 0, color="dark grey") +
  theme_linedraw() +
  geom_bar(position=dodge, stat="identity") +
  geom_errorbar(aes(ymax=UCI, ymin=LCI), position=dodge, width=0.5) + xlab("") +
  ylab("Coefficient") +
  facet_wrap(index~species, scales="free_y") +
  theme(legend.position = 'top', legend.title = element_blank()) +
  scale_fill_manual(values=cb[c(2,6)], breaks=c("Before 1988/89", "After 1988/89"))
dev.off()

head(plot)

# order era factors
plot$era <- reorder(plot$era, as.numeric(plot$plot.era))


pdo.sr <- ggplot(filter(plot, index=="PDO"), aes(species, effect, fill=era)) + geom_hline(yintercept = 0, color="dark grey") +
  theme_linedraw() +
  geom_bar(position=dodge, stat="identity") +
  geom_errorbar(aes(ymax=UCI, ymin=LCI), position=dodge, width=0.5) + xlab("") +
  ylab("Coefficient") +
  facet_wrap(~region, scales="free_y", ncol=1) +
  theme(legend.position = 'top', legend.title = element_blank()) +
  scale_fill_manual(values=cb[c(2,6)]) +
  ggtitle("a) PDO")

npgo.sr <- ggplot(filter(plot, index=="NPGO"), aes(species, effect, fill=era)) + geom_hline(yintercept = 0, color="dark grey") +
  theme_linedraw() +
  geom_bar(position=dodge, stat="identity") +
  geom_errorbar(aes(ymax=UCI, ymin=LCI), position=dodge, width=0.5) + xlab("") +
  ylab("Coefficient") +
  facet_wrap(~region, scales="free_y", ncol=1) +
  theme(legend.position = 'top', legend.title = element_blank()) +
  scale_fill_manual(values=cb[c(2,6)]) +
  ggtitle("b) NPGO")


png("alternate salmon SR pdo.png", 6,6, units="in", res=300)
ggarrange(pdo.sr, npgo.sr, ncol=2)
dev.off()


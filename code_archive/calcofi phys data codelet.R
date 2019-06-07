
dat <- read.csv("CalCOFI.phys.data.csv", row.names = 1)

all <- dat %>%
  spread(variable, pred)

all <- all %>%
  gather(key, value, -season, -year) %>%
  filter(season=="winter") %>%
  spread(key, value) %>%
  select(-season)

write.csv(all, "calcofi.phys.gam.csv")

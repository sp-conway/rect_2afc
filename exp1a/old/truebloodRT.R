# calculating median trueblood RT to get an idea of presentation times for exps 1a & 1b
rm(list=ls())
library(tidyverse)
library(fs)

files <- dir_ls("/Users/seanconway/Manuscripts/Repulsion_Effect/trueblood2013/Data/E1a")
dat <- map_dfr(files, read_delim, skip=6)

# filtered out very slow and very fast RTs, though it doesn't really matter 
dat_attraction <- filter(dat, str_detect(Effect, "Att")) %>%
  mutate(RTz=(RT-mean(RT))/sd(RT)) %>%
  filter(RTz<3 & RTz>-3)
 
ggplot(dat_attraction,aes(RT))+
  geom_histogram(bins=40)+
  scale_x_continuous(breaks=seq(0,20000,2000))
dat_attraction %>%
  summarise(median(RT),
            mean(RT),
            sd(RT),
            min(RT),
            max(RT))
t.test(dat_attraction$RT)


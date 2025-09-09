rm(list=ls())
library(here)
library(tidyverse)
load(here("analyses/stan/m14/m14_summary.RData"))
options(pillar.sigfig = 3)
t <- as_tibble(fit_summary$summary) %>%
  mutate(variable=rownames(fit_summary$summary)) %>%
  relocate(variable,.before=everything()) %>%
  filter(str_detect(variable,"\\[",negate=T)) %>%
  slice(1:13) %>% 
  select(-se_mean)

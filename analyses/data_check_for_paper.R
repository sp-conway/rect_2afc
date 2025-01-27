# data checks
rm(list=ls())

# packages
library(tidyverse)
library(here)

d <- here("exp1b","data","aggregated","rect_exp1b_aggregated.csv") %>%
  read_csv() %>%
  mutate(display=recode(display, spektor="triangle", trueblood="horizontal"),
         across(c(sub_n,tdd,tdo,probe,display,choice),as.factor)) %>%
  mutate(block=case_when(
    trial_number >= 1 & trial_number<=90  ~ 1,
    trial_number >= 91 & trial_number <= 180 ~ 2,
    trial_number >=181 & trial_number <= 270 ~ 3,
    trial_number >= 271 & trial_number <= 360 ~ 4,
    trial_number >= 361 & trial_number <= 450 ~ 5
  ))

# check number of subjects
# should be 85 (removed 1 subject for failing catch trials)
length(unique(d$sub_n))

# rts should be >= 100 ms and < 10000 ms
range(d$rt)

dd <- d %>%
  group_by(sub_n,block,trial) %>%
  summarise(N=n()) %>%
  ungroup()

  
ddd <- dd %>%
  group_by(sub_n,block) %>%
  mutate(NN=sum(N)) %>%
  ungroup()

max(ddd$NN)

dddd <- d %>%
  filter(trial=="critical") %>%
  group_by(sub_n,block,trial,tdd) %>%
  summarise(N=n()) %>%
  ungroup()

ddddd <- d %>%
  filter(trial=="critical") %>%
  group_by(sub_n,block,trial,tdo,tdd,display,probe) %>%
  summarise(N=n()) %>%
  ungroup()
max(ddddd$N)

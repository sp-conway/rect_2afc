rm(list=ls())
library(here)
library(tidyverse)
dat <- read_csv("/Users/seanconway/Downloads/rect_exp1_SeanTest - 2023-02-28T150626.626.csv") %>%
  mutate(across(c(display,tdo,tdd,probe),as.factor),
         choice=case_when(
           response==1~r1,
           response==2~r2,
           response==3~r3
         ),
         choice=as.factor(choice)) 

dat %>% 
  filter(trial=="critical" & !is.na(choice)) %>%
  count(probe,choice,.drop=F) %>%
  filter(choice %in% c("c","d","t"))

# data checks
rm(list=ls())

# packages
library(tidyverse)
library(here)

# read in datasets
exp1a <- here("exp1a","data","aggregated","rect_exp1a_FIXEDSPEKTOR_aggregated.csv") %>%
  read_csv() %>%
  mutate(exp="1a",
         display=recode(display, spektor="triangle", trueblood="horizontal"),
         across(c(sub_n,exp,tdd,tdo,probe,display,choice),as.factor))
exp1b <- here("exp1b","data","aggregated","rect_exp1b_aggregated.csv") %>%
  read_csv() %>%
  mutate(exp="1b",
         display=recode(display, spektor="triangle", trueblood="horizontal"),
         across(c(sub_n,exp,tdd,tdo,probe,display,choice),as.factor))

# number of subjects
n_subs_exp1a <- length(unique(exp1a$sub_n))
n_subs_exp1b <- length(unique(exp1b$sub_n))

# separate data by catch/critical 
catch <- bind_rows(
  filter(exp1a,trial=="catch"),
  filter(exp1b,trial=="catch")
) %>%
  select(-c(tdd,tdo,tw,th,cw,ch,dw,dh,ta,ca,da,probe)) %>%
  mutate(correct=case_when(
    choice=="l"~1,
    choice %in% c("sw","sh")~0
  ))

critical <- bind_rows(
  filter(exp1a,trial=="critical"),
  filter(exp1b,trial=="critical")
)
## checks ==================================================================
critical %>%
  filter(tdd==0) %>%
  filter(probe=="td") %>%
  count(exp,display,choice)

sub_0_props <- critical %>%
  filter(tdd==0) %>%
  filter(probe=="td") %>%
  count(sub_n,exp,display,choice) %>% 
  group_by(sub_n, display,exp) %>% 
  mutate(prop=n/sum(n)) %>% 
  ungroup()

sub_0_props %>% 
  filter(choice=="t") %>%
  ggplot(aes(prop))+
  geom_histogram(bins = 10,fill="lightblue",col="black")+
  geom_vline(xintercept=.5,linetype="dashed")+
  facet_wrap(vars(exp,display))+
  labs(x='p(target) on tdd=0.00 trials')+
  ggthemes::theme_few()

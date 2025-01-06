rm(list=ls())
library(tidyverse)
library(fs)

file_dir <- here("exp1a","test_data")
files <- dir_ls(file_dir)

read_data <- function(f){
  read_csv(f,col_types=cols()) %>%
    filter(screen_id=="make choice") %>%
    select(-platform_version) %>%
    mutate(sub_n=as.character(sub_n),
           rchoice=case_when(
      response==1~r1,
      response==2~r2,
      response==3~r3
    )) %>%
    mutate(across(c(sub_n,tdd,display,rchoice,probe),as.factor))
}

data <- map_dfr(files, read_data)
data %>%
  filter(trial=="critical" & !is.null(tdd)) %>%
  mutate(tdd=fct_drop(tdd,"null")) %>%
  count(tdd,display,probe,rchoice,.drop=F) %>%
  mutate(probe1=str_sub(probe,1,1),
         probe2=str_sub(probe,2,2)) %>%
  rowwise() %>%
  filter(rchoice==probe1 | rchoice==probe2) %>%
  ungroup() %>%
  group_by(tdd,display,probe) %>%
  mutate(prop=n/sum(n),
         prop=case_when(
           is.nan(prop)~0,
           TRUE~prop
         )) %>%
  ungroup() %>%
  ggplot(aes(tdd,prop,fill=rchoice))+
  geom_col(position="dodge")+
  facet_grid(probe~display,scales="free_x",labeller = label_both)+
  scale_y_continuous(limits=c(0,1),breaks=c(0,.5,1))+
  scale_fill_discrete(name="choice")+
  labs(x="target decoy distance",y="choice proportion")+
  ggthemes::theme_few()

ggsave(here("exp1a","testplot.pdf"),
       width=5,height=5)

data %>%
  filter(trial=="catch") %>%
  count(rchoice)

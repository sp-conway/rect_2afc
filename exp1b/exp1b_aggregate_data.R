rm(list=ls())
library(tidyverse)
library(fs)
library(glue)
library(here)

# specify directories ahead of time 
data_dir <- here("exp1b","data")

# read in data ====================================================
data_files <- dir_ls(data_dir,type="file",regexp="csv")

read_data <- function(f){
  print(f)
  f %>%
    data.table::fread() %>%
    as_tibble() %>%
    select(-c(success,trial_type,time_elapsed,internal_node_id,timeout,failed_audio,failed_images,failed_video,stimulus,view_history)) %>%
    filter(screen_id=="make choice") %>%
    arrange(trial_index) %>%
    mutate(across(c(display,tdd,tdo,probe,response),as.factor),
           choice=case_when(
             response==1~r1,
             response==2~r2,
             response==3~r3
           ),
           trial_number=1:n(),
           across(c(tdd,tw,th,cw,ch,dw,dh,ta,ca,da,
                    csls,csss,clls,clss,sa,la),~na_if(.x,"null")),
           across(c(csls,csss,clls,clss,sa,la,h1,w1,h2,w2,h3,w3,tw,th,cw,ch,dw,dh),as.numeric),
           across(c(h1,w1,h2,w2,h3,w3,tw,th,cw,ch,dw,dh),round),
           ta=tw*th,
           ca=cw*ch,
           da=dw*dh,
           probe=na_if(probe,""),
           rt=as.numeric(rt)) %>%
    select(-c(trial_index,screen_id,r1yloc,r2yloc,r3yloc)) %>%
    relocate(trial_number,.after=sub_n) %>%
    relocate(response,.before=choice) %>%
    relocate(rt,.before=response) %>%
    filter(rt>=100 & rt<=10000) # filter out very fast and very slow rts
}

# read in data
data_all <- map_dfr(data_files,read_data)

n_subs_all <- length(unique(data_all$sub_n))

catch <- filter(data_all, trial=="catch")

sub_ns_filtered <- catch %>%
  mutate(correct=case_when(
    choice=="l"~1,
    TRUE~0
  )) %>%
  count(sub_n,correct) %>%
  group_by(sub_n) %>%
  mutate(prop=n/sum(n)) %>%
  filter(correct==1) %>%
  filter(prop>=.8) %>%
  pull(sub_n)

data_all_filtered <- data_all %>%
  filter(sub_n %in% sub_ns_filtered) 

n_filtered <- n_subs_all-length(sub_ns_filtered)

glue("filtered out {n_filtered} subjects based on catch performance\nalso filtered out trials with rt <100ms and >10000ms") %>%
  write_lines(file=here("analyses","exp1b_filtered_subs.txt"))

# write aggregated data file 
glue("{data_dir}/aggregated/rect_exp1b_aggregated.csv") %>%
  write_csv(data_all_filtered,file=.)



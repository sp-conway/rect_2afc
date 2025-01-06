## NOTE ## 
## AN EXPERIMENT CODE ERROR ON THE SPEKTOR TRIALS CAUSED H & W TO BE SWITCHED AROUND
## THIS CHANGED ORIENTATION FOR ALL RECTANGLES ##
## THIS WILL ONLY CHANGE CONCLUSIONS IF TRIALS ARE SPLIT BY TD ORIENTATION ##
## THIS HAS BEEN FIXED IN THE AGGREGATED FILES (AND THUS IN THIS ANALYSIS) ##
## BUT NOT IN THE INDIVIDUAL DATA FILES (OR THE EXPERIMENT CODE) ##

# setup ==================================================================
rm(list=ls())
library(tidyverse)
library(fs)
library(here)
library(glue)

# specify directories ahead of time 
data_dir <- here("exp1a","data")

data_files <- dir_ls(data_dir, type = "file",regexp = "csv")


# unwieldy function to read in data files
read_data <- function(f){
  f %>%
    data.table::fread() %>%
    as_tibble() %>%
    select(-c(success,trial_type,time_elapsed,internal_node_id,timeout,failed_audio,failed_images,failed_video,stimulus,view_history)) %>%
    filter(screen_id=="make choice") %>%
    arrange(trial_index) %>%
    mutate(across(c(sub_n,display,tdd,tdo,response),as.factor),
           choice=case_when(
             response==1~r1,
             response==2~r2,
             response==3~r3
           ),
           trial_number=1:n(),
           across(c(tdd,tw,th,cw,ch,dw,dh,ta,ca,da,
                    csls,csss,clls,clss,sa,la),~na_if(.x,"null")),
           across(c(csls,csss,clls,clss,sa,la),as.numeric)) %>%
    mutate(w1_temp=case_when(
      display=="spektor"~h1,
      TRUE~w1
    ),
    w2_temp=case_when(
      display=="spektor"~h2,
      TRUE~w2
    ),
    w3_temp=case_when(
      display=="spektor"~h3,
      TRUE~w3
    ),
    h1_temp=case_when(
      display=="spektor"~w1,
      TRUE~h1
    ),
    h2_temp=case_when(
      display=="spektor"~w2,
      TRUE~h2
    ),
    h3_temp=case_when(
      display=="spektor"~w3,
      TRUE~h3
    )) %>%
    mutate(tdo=case_when(
      display=="spektor" & trial=="critical" & tdo=="h"~"w",
      display=="spektor" & trial=="critical" & tdo=="w"~"h",
      display=="trueblood" & trial=="critical" & tdo=="h"~"h",
      display=="trueblood" & trial=="critical" & tdo=="w"~"w",
      trial=="catch" ~ NA_character_
    ),
    probe=na_if(probe,"")) %>%
    mutate(h1=h1_temp,
           w1=w1_temp,
           h2=h2_temp,
           w2=w2_temp,
           h3=h3_temp,
           w3=w3_temp,
           across(c(h1,w1,h2,w2,h3,w3),round)) %>% # round all dim values since probably not frac of pixel) 
    mutate(dw=case_when(
      r1=="d" & display=="spektor" ~ w1,
      r2=="d" & display=="spektor" ~ w2,
      r3=="d" & display=="spektor" ~ w3,
      display=="trueblood" & trial=="catch" ~ NA_real_,
      display=="trueblood" & trial=="critical" & r1=="d" ~ w1,
      display=="trueblood" & trial=="critical" & r2=="d" ~ w2,
      display=="trueblood" & trial=="critical" & r3=="d" ~ w3,
    ),
    dh=case_when(
      r1=="d" & display=="spektor" ~ h1,
      r2=="d" & display=="spektor" ~ h2,
      r3=="d" & display=="spektor" ~ h3,
      display=="trueblood" & trial=="catch" ~ NA_real_,
      display=="trueblood" & trial=="critical" & r1=="d" ~ h1,
      display=="trueblood" & trial=="critical" & r2=="d" ~ h2,
      display=="trueblood" & trial=="critical" & r3=="d" ~ h3,
    ),
    tw=case_when(
      r1=="t" & display=="spektor" ~ w1,
      r2=="t" & display=="spektor" ~ w2,
      r3=="t" & display=="spektor" ~ w3,
      display=="trueblood" & trial=="catch" ~ NA_real_,
      display=="trueblood" & trial=="critical" & r1=="t" ~ w1,
      display=="trueblood" & trial=="critical" & r2=="t" ~ w2,
      display=="trueblood" & trial=="critical" & r3=="t" ~ w3
    ),
    th=case_when(
      r1=="t" & display=="spektor" ~ h1,
      r2=="t" & display=="spektor" ~ h2,
      r3=="t" & display=="spektor" ~ h3,
      display=="trueblood" & trial=="catch" ~ NA_real_,
      display=="trueblood" & trial=="critical" & r1=="t" ~ h1,
      display=="trueblood" & trial=="critical" & r2=="t" ~ h2,
      display=="trueblood" & trial=="critical" & r3=="t" ~ h3
    ),
    cw=case_when(
      r1=="c" & display=="spektor" ~ w1,
      r2=="c" & display=="spektor" ~ w2,
      r3=="c" & display=="spektor" ~ w3,
      display=="trueblood" & trial=="catch" ~ NA_real_,
      display=="trueblood" & trial=="critical" & r1=="c" ~ w1,
      display=="trueblood" & trial=="critical" & r2=="c" ~ w2,
      display=="trueblood" & trial=="critical" & r3=="c" ~ w3
    ),
    ch=case_when(
      r1=="c" & display=="spektor" ~ h1,
      r2=="c" & display=="spektor" ~ h2,
      r3=="c" & display=="spektor" ~ h3,
      display=="trueblood" & trial=="catch" ~ NA_real_,
      display=="trueblood" & trial=="critical" & r1=="c" ~ h1,
      display=="trueblood" & trial=="critical" & r2=="c" ~ h2,
      display=="trueblood" & trial=="critical" & r3=="c" ~ h3
    )
    ) %>%
    mutate(ta=tw*th, # recomputing decoy height and width with rounded values
           ca=cw*ch,
           da=dw*dh) %>%
    select(-c(trial_index,screen_id,w1_temp,w2_temp,w3_temp,h1_temp,h2_temp,h3_temp)) %>%
    relocate(trial_number,.after=sub_n) %>%
    relocate(response,.before=choice) %>%
    relocate(rt,.before=response) %>%
    mutate(rt=as.numeric(rt)) %>%
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
  write_lines(file=here("exp1a","data","exp1a_filtered_subs.txt"))

# write aggregated data file 
here(data_dir,"aggregated","rect_exp1a_FIXEDSPEKTOR_aggregated.csv") %>%
  write_csv(data_all_filtered, file=.)

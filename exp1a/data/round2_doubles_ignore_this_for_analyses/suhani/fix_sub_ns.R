rm(list=ls())
setwd("/Users/seanconway/Manuscripts/Repulsion_Effect/Experiments/exp1a/data/round2_doubles_ignore_this_for_analyses/suhani/corrected subs")
library(tidyverse)
files <- list.files()
read_and_fix <- function(f){
  d <- read_csv(f)
  sub_n_new <- str_replace_all(f, c("rect_exp1_"="",".csv"="")) %>%
    as.numeric()
  d$sub_n <- sub_n_new
  write_csv(d,f)
}
map(files, read_and_fix)

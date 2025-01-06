rm(list=ls())
# creating trials for experiments 1a and 1b
library(tidyverse)
library(jsonlite)
library(here)

# per spektor, large side will be 240, small side will be 164
high <- 240 
low <- 164

# target decoy distance
tdd <- c(0,.02,.05,.09,.14) 

# target decoy orientation
tdo <- c("w","h") 

# display (triangle is spektor, line is trueblood)
display <- c("spektor","trueblood")

# which two stimuli are we asking about?
probe <- c("td","tc","dc")

# critical trials
crit_trials <- crossing(
  tdd,
  tdo,
  display,
  probe
) %>%
  mutate(
    # figure out target and competitor dimensions
    tw=case_when(
      tdo=="h"~low,
      tdo=="w"~high
    ),
    th=case_when(
      tdo=="h"~high,
      tdo=="w"~low
    ),
    cw=case_when(
      tdo=="h"~high,
      tdo=="w"~low
    ),
    ch=case_when(
      tdo=="h"~low,
      tdo=="w"~high
    ),
    ta=th*tw,
    ca=ch*cw,
    threl=th/(th+tw), # target rel h/w
    twrel=tw/(th+tw), # target rel w/h
    dw=tw - (tw*twrel*tdd), 
    dh=th - (th*threl*tdd),
    da=dh*dw,
    dtar=da/ta,
    trial="critical",
    csls=NULL,
    csss=NULL,
    clls=NULL,
    clss=NULL,
    sa=NULL,
    la=NULL
  )

# catch trials 
catch_small_large_side <- 180
catch_small_small_side <- 120
catch_large_large_side <- 260
catch_large_small_side <- 200
n_catch <- 30
catch_trials <- tibble(
  trial="catch",
  csls=rep(catch_small_large_side,n_catch),
  csss=rep(catch_large_small_side,n_catch),
  clls=rep(catch_large_large_side,n_catch),
  clss=rep(catch_large_small_side,n_catch),
  sa=csls*csss,
  la=clls*clss,
  display=c(rep("spektor",n_catch/2),
            rep("trueblood",n_catch/2)),
  tdd=NULL,
  tdo=NULL,
  probe=NULL,
  tw=NULL,
  th=NULL,
  cw=NULL,
  ch=NULL,
  ta=NULL,
  ca=NULL,
  threl=NULL,
  twrel=NULL,
  dw=NULL,
  dh=NULL,
  dtar=NULL
)

all_trials <- select(crit_trials,
       trial,display,tdd,tdo,tw,th,cw,ch,dw,dh,ta,ca,da,probe) %>%
  bind_rows(
    catch_trials
  )
all_trials %>%
  write_json(path=here("exp1a","experiment_code","rect_exp1_trials.json"))

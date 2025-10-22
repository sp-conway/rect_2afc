# setup ==============================================================
rm(list=ls())
library(tidyverse)
library(glue)
library(ggsci)
library(tidybayes)
library(fs)
library(rstan)
library(posterior)
library(bayesplot)
library(here)
library(patchwork)
options(digits=5)

read_dat <- function(f){
  dd <- read_csv(f) %>%
    mutate(display=recode(display, spektor="triangle", trueblood="horizontal"),
           across(c(sub_n,tdo,display,choice),as.factor),
           display=fct_relevel(display,c("triangle","horizontal")))
  return(dd)
}

d <- here("exp1b","data","aggregated","rect_exp1b_aggregated.csv") %>%
  read_dat()


# number of ppts
get_subs <- function(d){
  p <- length(unique(d$sub_n))
  print(p)
  return(p)
}

n_subs <- get_subs(d)

# separate data by catch/critical 
catch <- d %>%
  filter(trial=="catch") %>%
  select(-c(tdd,tdo,tw,th,cw,ch,dw,dh,ta,ca,da,probe)) %>%
  mutate(correct=case_when(
    choice=="l"~1,
    choice %in% c("sw","sh")~0
  ))

critical <- d %>%
  filter(trial=="critical")

# checking props on all critical trials ========================================================================================
critical_mean_props_all <- critical %>%
  mutate(tdd=as.factor(tdd),
         probe=as.factor(probe)) %>%
  group_by(sub_n,tdd,probe,choice) %>%
  summarise(N=n()) %>%
  group_by(sub_n,tdd,probe) %>%
  mutate(prop=N/sum(N)) %>%
  ungroup() %>%
  group_by(tdd,probe,choice) %>%
  summarise(m=mean(prop),
            s=sd(prop),
            N=n(),
            se=s/sqrt(N),
            ci_lower=m-qt(.975,n()-1)*se,
            ci_upper=m+qt(.975,n()-1)*se) %>%
  ungroup() 
critical_mean_props_all

# critical props on all trial types, all tdd value
critical_mean_props_all %>%
  ggplot(aes(tdd,m,col=choice,group=choice))+
  geom_point()+
  geom_errorbar(aes(ymin=ci_lower,ymax=ci_upper),width=.2)+
  geom_line()+
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.2))+
  facet_grid(probe~.)+
  ggthemes::theme_few()
ggsave(filename=here("analyses","plots","crit_props_all.jpeg"),width=4,height=6)

# critical trial analyses - td and cd ========================================================================================
# NOTE THAT IN THE PAPER WE LABEL COMPETITOR-DECOY TRIALS CD BUT HERE I USE DC
# CHANGED LATER ON IN FIGURE

# baseline display level should be triangle for model
# baseline probe level is cd (dc)
# baseline distance is 2%
critical_1 <- critical %>%
  filter(tdd!=0 & probe!="tc") %>% # filter out tc trials
  mutate(tdd=as.factor(tdd),
         probe=factor(probe,levels=c("dc","td"))) %>%
  mutate(display=fct_relevel(display,c("triangle","horizontal")),
         discrim=case_when(# whether or not ppt can discriminate from decoy
           choice=="t" & probe=="td"~1,
           choice=="c" & probe=="dc"~1, 
           choice=="d"~0
         )) %>%
  arrange(sub_n,trial_number)

# Compute mean choice proportions, with CIs
# Not using this figure in paper
critical_mean_discrim <- critical_1 %>%
  group_by(sub_n,tdd,probe,display) %>%
  summarise(p_discrim=mean(discrim)) %>%
  group_by(tdd, probe, display) %>%
  summarise(mean_discrim=mean(p_discrim),
            se=sd(p_discrim)/sqrt(n()),
            ci_lower=mean_discrim-qt(.975,n()-1)*se,
            ci_upper=mean_discrim+qt(.975,n()-1)*se) %>%
  ungroup() %>%
  mutate(tdd=as.numeric(as.character(tdd))*100)

critical_mean_discrim %>%
  ggplot(aes(tdd, mean_discrim,col=probe))+
  geom_point(alpha=.8)+
  geom_line(alpha=.8)+
  geom_errorbar(aes(ymin=ci_lower,ymax=ci_upper),width=.75,alpha=.8)+
  facet_grid(.~display)+
  scale_x_continuous(breaks=c(2,5,9,14))+
  scale_y_continuous(limits=c(.5,1))+
  ggsci::scale_color_startrek(name="choice")+
  labs(x="target-decoy distance",y="proability discriminate")+
  ggthemes::theme_few()+
  theme(text=element_text(size=18),
        plot.title=element_text(hjust=0.5),
        legend.position = "inside",
        legend.position.inside = c(.1,.85))
ggsave(filename=here("analyses","plots","critical_discrim_CIs.jpeg"),width=6,height=5)

# Critical TD trials - Statistical modeling ========================================================================
# renumbering subject numbers to go sequentially
# this key specifies which original subject # goes with which new subject #
subs_key <- tibble(
  sub_n = sort(unique(critical_1$sub_n)),
  sub_n_new = seq(1,n_subs,1)
)
critical_for_model <- tibble(
  sub_n=critical_1$sub_n,
  trial_number=critical_1$trial_number,
  tdo=if_else(critical_1$tdo=="w",1,0),
  display=if_else(critical_1$display=="horizontal",1,0),
  probe=if_else(critical_1$probe=="td",1,0),
  tdd_5=if_else(critical_1$tdd==0.05,1,0),
  tdd_9=if_else(critical_1$tdd==0.09,1,0),
  tdd_14=if_else(critical_1$tdd==0.14,1,0),
  discrim=critical_1$discrim
) %>%
  left_join(subs_key) %>% ## assigning "new" subject numbers. Just going sequentially from 1-n subs for easier indexing in stan
  relocate(sub_n_new,.after=sub_n)
critical_unique <- critical_for_model %>% # ALL UNIQUE TRIALS FOR PREDICTION IN MODEL
  distinct(sub_n_new,tdo,display,probe,tdd_5,tdd_9,tdd_14)

# number of unique trials
n_unique <- nrow(critical_unique)

# data for stan
stan_data <- list(
  n_subs=n_subs, # N
  n_trials=nrow(critical_for_model),
  sub_n_new=critical_for_model$sub_n_new, 
  tdo=critical_for_model$tdo,
  display=critical_for_model$display,
  probe=critical_for_model$probe,
  tdd_5=critical_for_model$tdd_5,
  tdd_9=critical_for_model$tdd_9,
  tdd_14=critical_for_model$tdd_14,
  discrim=critical_for_model$discrim,
  n_trials_unique=n_unique,
  sub_n_new_unique=critical_unique$sub_n_new,
  tdo_unique=critical_unique$tdo,
  display_unique=critical_unique$display,
  probe_unique=critical_unique$probe,
  tdd_5_unique=critical_unique$tdd_5,
  tdd_9_unique=critical_unique$tdd_9,
  tdd_14_unique=critical_unique$tdd_14
)

# controls
debug_model <- F # whether or not we're testing the model to make sure stan code works
prefix <- "m14" # which model iteration
model_dir <- path(here("analyses","stan",prefix)) # stan files / save directory
stan_model_code <- path(model_dir,glue("{prefix}.stan")) # model code
fit_file <- path(model_dir,glue("{prefix}_fit.RData")) # name of our resulting fit object


if(debug_model){
  n_iter <- 100
  n_core <- 1
  n_chain <- 2
  n_warm <- 10
  to_save <- F
}else{
  # number of iterations for each model per core (for MCMC)
  n_iter <- 2500
  n_chain <- 5
  n_core <- 10
}

if(!file_exists(fit_file) | debug_model){
  model_compiled <- stan_model(stan_model_code)
  fit <- sampling(model_compiled,
                  data=stan_data,
                  chains=n_chain,
                  cores=n_core,
                  iter=n_iter,
                  control=list(adapt_delta=.95)) 
  if(!debug_model) save(fit,file=fit_file)
  diagnostics <- get_sampler_params(fit)
  if(!debug_model) save(diagnostics, file=path(model_dir,glue("{prefix}_diag.RData")))
  fit_summary <- summary(fit,probs=c(.025,.975))
  if(!debug_model) save(fit_summary, file=path(model_dir,glue("{prefix}_summary.RData")))
}else{
  load(fit_file)
}

# check modeling ===============================================================================
check_energy(fit)


color_scheme_set("red")

try({
  p <- mcmc_trace(fit,c("lp__"))
  p
  if(!debug_model) ggsave(p, filename=path(model_dir,glue("{prefix}_lp___trace.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_trace(fit,c("b_0"))
  p
  if(!debug_model) ggsave(p, filename=path(model_dir,glue("{prefix}_b_0_trace.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_trace(fit,c("b_tdd_5","b_tdd_9","b_tdd_14"))
  if(!debug_model) ggsave(p, filename=path(model_dir,glue("{prefix}_b_tdd_trace.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_trace(fit,c("b_td"))
  if(!debug_model) ggsave(p, filename=path(model_dir,glue("{prefix}_b_td_trace.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_trace(fit,c("b_w"))
  if(!debug_model) ggsave(p, filename=path(model_dir,glue("{prefix}_b_w_trace.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_trace(fit,c("b_h"))
  if(!debug_model) ggsave(p,filename=path(model_dir,glue("{prefix}_b_h_trace.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_trace(fit,c("b_tdd_5_X_td","b_tdd_9_X_td","b_tdd_14_X_td"))
  if(!debug_model) ggsave(p,filename=path(model_dir,glue("{prefix}_b_tdd_X_td_trace.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_trace(fit,c("sigma_b_0_s","sigma_b_s"))
  if(!debug_model) ggsave(p,filename=path(model_dir,glue("{prefix}_sigma_trace.jpeg")),width=7,height=6)
  rm(p)
})

try({
  p <- mcmc_hist(fit,c("lp__"))
  if(!debug_model) ggsave(p,filename=path(model_dir,glue("{prefix}_lp___hist.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_hist(fit,c("b_0"))
  if(!debug_model) ggsave(p,filename=path(model_dir,glue("{prefix}_b_0_hist.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_hist(fit,c("b_tdd_5","b_tdd_9","b_tdd_14"))
  if(!debug_model) ggsave(p,filename=path(model_dir,glue("{prefix}_b_tdd_hist.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_hist(fit,c("b_td"))
  if(!debug_model) ggsave(p,filename=path(model_dir,glue("{prefix}_b_td_hist.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_hist(fit,c("b_w"))
  if(!debug_model) ggsave(p,filename=path(model_dir,glue("{prefix}_b_w_hist.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_hist(fit,c("b_h"))
  if(!debug_model) ggsave(p,filename=path(model_dir,glue("{prefix}_b_h_hist.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_hist(fit,c("b_tdd_5_X_td","b_tdd_9_X_td","b_tdd_14_X_td"))
  if(!debug_model) ggsave(p,filename=path(model_dir,glue("{prefix}_b_tdd_X_td_hist.jpeg")),width=5,height=4)
  rm(p)
})


try({
  p <- mcmc_hist(fit,c("sigma_b_0_s","sigma_b_s"))
  if(!debug_model) ggsave(p,filename=path(model_dir,glue("{prefix}_sigma_hist.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_dens(fit,c("lp__"))
  if(!debug_model) ggsave(p,filename=path(model_dir,glue("{prefix}_lp___dens.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_dens(fit,c("b_0"))
  if(!debug_model) ggsave(p,filename=path(model_dir,glue("{prefix}_b_0_dens.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_dens(fit,c("b_tdd_5","b_tdd_9","b_tdd_14"))
  if(!debug_model) ggsave(p,filename=path(model_dir,glue("{prefix}_b_tdd_dens.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_dens(fit,c("b_td"))
  if(!debug_model) ggsave(p,filename=path(model_dir,glue("{prefix}_b_td_dens.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_dens(fit,c("b_w"))
  if(!debug_model) ggsave(p,filename=path(model_dir,glue("{prefix}_b_w_dens.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_dens(fit,c("b_h"))
  if(!debug_model) ggsave(p,filename=path(model_dir,glue("{prefix}_b_h_dens.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_dens(fit,c("b_tdd_5_X_td","b_tdd_9_X_td","b_tdd_14_X_td"))
  if(!debug_model) ggsave(p,filename=path(model_dir,glue("{prefix}_b_tdd_X_td_dens.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_dens(fit,c("sigma_b_0_s","sigma_b_s"))
  if(!debug_model) ggsave(p,filename=path(model_dir,glue("{prefix}_sigma_dens.jpeg")),width=5,height=4)
  rm(p)
})


try({
  p <- mcmc_dens_chains(fit,c("lp__"))
  if(!debug_model) ggsave(p,filename=path(model_dir,glue("{prefix}_lp___dens_chains.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_dens_chains(fit,c("b_0"))
  if(!debug_model) ggsave(p,filename=path(model_dir,glue("{prefix}_b_0_dens_chains.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_dens_chains(fit,c("b_tdd_5","b_tdd_9","b_tdd_14"))
  if(!debug_model) ggsave(p,filename=path(model_dir,glue("{prefix}_b_tdd_dens_chains.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_dens_chains(fit,c("b_td"))
  if(!debug_model) ggsave(p,filename=path(model_dir,glue("{prefix}_b_td_dens_chains.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_dens_chains(fit,c("b_w"))
  if(!debug_model) ggsave(p,filename=path(model_dir,glue("{prefix}_b_w_dens_chains.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_dens_chains(fit,c("b_h"))
  if(!debug_model) ggsave(p,filename=path(model_dir,glue("{prefix}_b_h_dens_chains.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_dens_chains(fit,c("b_tdd_5_X_td","b_tdd_9_X_td","b_tdd_14_X_td"))
  if(!debug_model) ggsave(p,filename=path(model_dir,glue("{prefix}_b_tdd_X_td_dens_chains.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_dens_chains(fit,c("sigma_b_0_s","sigma_b_s"))
  if(!debug_model) ggsave(p,filename=path(model_dir,glue("{prefix}_sigma_dens_chains.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_pairs(fit,c("b_tdd_5","b_tdd_9","b_tdd_14"))
  if(!debug_model) ggsave(p,filename=path(model_dir,glue("{prefix}_b_tdd_pairs.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_pairs(fit,pars = c("b_0","b_w","b_h","b_td","b_tdd_5","b_tdd_9","b_tdd_14","b_tdd_5_X_td","b_tdd_9_X_td","b_tdd_14_X_td"),
    off_diag_args = list(size = .5, alpha = 0.5))
  ggsave(p,filename=path(model_dir,glue("{prefix}_b_pairs_fixed.jpeg")),width=10,height=10)
  rm(p)
})

try({
  p <- mcmc_pairs(fit,pars = c("b_0","b_w","b_h","b_td"),
                  off_diag_args = list(size = .5, alpha = 0.5))
  p
  ggsave(p,filename=path(model_dir,glue("{prefix}_b_pairs_noTDD_fixed.jpeg")),width=10,height=10)
  rm(p)
})

try({
  p <- mcmc_pairs(fit,pars = c("b_tdd_5","b_tdd_9","b_tdd_14","b_tdd_5_X_td","b_tdd_9_X_td","b_tdd_14_X_td"),
                  off_diag_args = list(size = .5, alpha = 0.5))
  p
  ggsave(p,filename=path(model_dir,glue("{prefix}_b_tdd_pairs_fixed.jpeg")),width=10,height=10)
  rm(p)
})

try({
  p <- mcmc_pairs(fit,regex_pars="sigma",
                  off_diag_args = list(size = .5, alpha = 0.5))
  if(!debug_model) ggsave(p,filename=path(model_dir,glue("{prefix}_sigma_pairs.jpeg")),width=10,height=10)
  rm(p)
})

try({
  p <- mcmc_intervals(fit, regex_pars = "^b_0_s",prob=.5,prob_outer = .95)+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())+
    labs(y="participant",x="estimate")
  if(!debug_model) ggsave(p,filename=path(model_dir,glue("{prefix}_b_0_s.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_intervals(fit, regex_pars = "^b_h_s",prob=.5,prob_outer = .95)+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())+
    labs(y="participant",x="estimate")
  if(!debug_model) ggsave(p,filename=path(model_dir,glue("{prefix}_b_h_s.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_intervals(fit, regex_pars = "^b_w_s",prob=.5,prob_outer = .95)+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())+
    labs(y="participant",x="estimate")
  if(!debug_model) ggsave(p,filename=path(model_dir,glue("{prefix}_b_w_s.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_intervals(fit, regex_pars = "^b_td_s",prob=.5,prob_outer = .95)+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())+
    labs(y="participant",x="estimate")
  if(!debug_model) ggsave(p,filename=path(model_dir,glue("{prefix}_b_td_s.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_intervals(fit, regex_pars = "^b_tdd_5_s",prob=.5,prob_outer = .95)+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())+
    labs(y="participant",x="estimate")
  if(!debug_model) ggsave(p,filename=path(model_dir,glue("{prefix}_b_tdd_5_s.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_intervals(fit, regex_pars = "^b_tdd_9_s",prob=.5,prob_outer = .95)+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())+
    labs(y="participant",x="estimate")
  if(!debug_model) ggsave(p,filename=path(model_dir,glue("{prefix}_b_tdd_9_s.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_intervals(fit, regex_pars = "^b_tdd_14_s",prob=.5,prob_outer = .95)+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())+
    labs(y="participant",x="estimate")
  if(!debug_model) ggsave(p,filename=path(model_dir,glue("{prefix}_b_tdd_14_s.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_intervals(fit, regex_pars = "^b_tdd_5_X_td_s",prob=.5,prob_outer = .95)+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())+
    labs(y="participant",x="estimate")
  if(!debug_model) ggsave(p,filename=path(model_dir,glue("{prefix}_b_tdd_5_X_td_s.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_intervals(fit, regex_pars = "^b_tdd_9_X_td_s",prob=.5,prob_outer = .95)+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())+
    labs(y="participant",x="estimate")
  if(!debug_model) ggsave(p,filename=path(model_dir,glue("{prefix}_b_tdd_9_X_td_s.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_intervals(fit, regex_pars = "^b_tdd_14_X_td_s",prob=.5,prob_outer = .95)+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())+
    labs(y="participant",x="estimate")
  if(!debug_model) ggsave(p,filename=path(model_dir,glue("{prefix}_b_tdd_14_X_td_s.jpeg")),width=5,height=4)
  rm(p)
})

try({
  p <- mcmc_intervals(fit,pars=c("b_0","b_w","b_h","b_td","b_tdd_5","b_tdd_9","b_tdd_14",
                              "b_tdd_5_X_td","b_tdd_9_X_td","b_tdd_14_X_td"),prob=.5,prob_outer = .95)
  rm(p)
})


# get model predictions ================================================================================================================================================
p <- extract(fit, pars="p")$p
if(!debug_model) save(p,file=path(model_dir,glue("{prefix}_p.RData")))

# convert p to data frame and clean
pd <- as.data.frame(t(p))
colnames(pd) <- paste0("iter_",1:ncol(pd))
pdd <- bind_cols(critical_unique,pd) %>%
  pivot_longer(cols=contains("iter"),values_to = "p",names_to = "iter") %>%
  mutate(tdd=case_when(
    tdd_5==1~5,
    tdd_9==1~9,
    tdd_14==1~14,
    T~2
  ),
  tdo=case_when(
    tdo==1~"w",
    tdo==0~"h"
  ),
  display=case_when(
    display==1~"horizontal",
    display==0~"triangle"
  ),
  probe=case_when(
    probe==1~"td",
    probe==0~"dc"
  )) %>%
  group_by(iter,display,probe,tdd) %>%
  summarise(mp=mean(p)) %>% # mean predictions across subjects for each iteration
  group_by(display,probe,tdd) %>%
  summarise(p=mean(mp), # mean predictions, HDI across iterations
            hdi_lower=hdi(mp)[,1],
            hdi_upper=hdi(mp)[,2]) %>%
  ungroup() %>%
  mutate(source="model")

# join to data and clean
pddd <- critical_1 %>%
  mutate(tdd=as.numeric(as.character(tdd))*100) %>%
  group_by(sub_n,display,probe,tdd) %>%
  summarise(p=mean(discrim)) %>%
  group_by(display,probe,tdd) %>%
  summarise(p=mean(p)) %>%
  ungroup() %>%
  mutate(source="data") %>%
  bind_rows(pdd) %>%
  mutate(probe=case_when(
    probe=="dc"~"cd", # IMPORTANT RELABELING DC TO CD
    T~probe,
  ),
  display=factor(display,levels=c("triangle","horizontal"))) # RELEVEL DISPLAY SO TRIANGLE IS ON TOP
write_csv(pddd, file=path(model_dir,glue("{prefix}_2afc_preds_v_data.csv")))

# PLOTTING MODEL PREDICTIONS VS. DATA IMPORTANT
p <- pddd %>%
  ggplot(aes(tdd,p,shape=source,col=probe))+
  geom_point(alpha=.65,size=2.5)+
  geom_line(data = filter(pddd, source == "data"), aes(tdd, p, group = probe), alpha = 0.8) +  
  geom_errorbar(aes(ymin=hdi_lower,ymax=hdi_upper),width=.25,alpha=.5,col="black")+
  scale_x_continuous(breaks=c(2,5,9,14),limits=c(1.5,14.5),labels=c("2%","5%","9%","14%"))+
  scale_y_continuous(limits=c(.5,1))+
  ggsci::scale_color_startrek(name="comparison")+
  scale_shape_manual(values=c(1,4))+
  # scale_color_manual(values=c("grey","black"))+
  labs(x="tdd",y="p(correct)")+
  facet_grid(display~.)+
  ggthemes::theme_few()+
  theme(text=element_text(size=18))
p
ggsave(p,filename=path(model_dir,glue("{prefix}_model_preds_v_data.jpeg")),
       width=5,height=4,dpi=500)

# check random effect constraints ================================================================================================================================================
b_0_s <- extract(fit,pars="b_0_s")$b_0_s
hist(rowSums(b_0_s))
mean(rowSums(b_0_s))

b_h_s <- extract(fit,pars="b_h_s")$b_h_s
hist(rowSums(b_h_s))

b_td_s <- extract(fit,pars="b_td_s")$b_td_s
hist(rowSums(b_td_s))

b_w_s <- extract(fit,pars="b_w_s")$b_w_s
hist(rowSums(b_w_s))

b_tdd_5_s <- extract(fit,pars="b_tdd_5_s")$b_tdd_5_s
hist(rowSums(b_tdd_5_s))

b_tdd_9_s <- extract(fit,pars="b_tdd_9_s")$b_tdd_9_s
hist(rowSums(b_tdd_9_s))

b_tdd_14_s <- extract(fit,pars="b_tdd_14_s")$b_tdd_14_s
hist(rowSums(b_tdd_14_s))

b_tdd_5_X_td_s <- extract(fit,pars="b_tdd_5_X_td_s")$b_tdd_5_X_td_s
hist(rowSums(b_tdd_5_X_td_s))

b_tdd_9_X_td_s <- extract(fit,pars="b_tdd_9_X_td_s")$b_tdd_9_X_td_s
hist(rowSums(b_tdd_9_X_td_s))

b_tdd_14_X_td_s <- extract(fit,pars="b_tdd_14_X_td_s")$b_tdd_14_X_td_s
hist(rowSums(b_tdd_14_X_td_s))
# tc trials =================================================================================================================================
critical_tc <- critical %>%
  filter(probe=="tc" & tdd!=0) %>%
  mutate(choose_t=choice=="t")
morey_correction <- function(X){
  J <- length(unique(X$display))*length(unique(X$tdd))
  corr_factor <- J/(J-1)
  X %>%
    group_by(sub_n) %>%
    mutate(pp=mean(prop)) %>%
    ungroup() %>%
    mutate(z=prop-pp+mean(prop)) %>%
    group_by(tdd,display) %>%
    summarise(sd=sqrt(var(z)*corr_factor),
              t=qt(.975,n()-1)*(sd/sqrt(n())),
              m=mean(prop),
              ci_lower=m-t,
              ci_upper=m+t)
}
critical_tc <- critical %>%
  filter(probe=="tc" & tdd!=0) %>%
  mutate(choose_t=choice=="t") %>%
  group_by(sub_n,tdd,display) %>%
  summarise(prop=mean(choose_t)) %>%
  ungroup()
critical_tc_means <- critical_tc %>%
  morey_correction()
  # 
  # group_by(sub_n, tdd, display) %>%
  # mutate(prop=mean(choose_t)) %>%
  # ungroup() %>%
  # group_by(tdd, display) %>%
  # summarise(mean_prop=mean(prop),
  #           se_prop=sd(prop)/sqrt(n()),
  #           ci_lower=mean_prop-qt(.975,n()-1)*se_prop,
  #           ci_upper=mean_prop+qt(.975,n()-1)*se_prop) %>%
  # ungroup() %>%
  

critical_tc_means %>%
  mutate(tdd=as.numeric(as.character(tdd))*100) %>%
  ggplot(aes(tdd, m))+
  geom_point(alpha=.8)+
  geom_line(alpha=.8)+
  geom_errorbar(aes(ymin=ci_lower,ymax=ci_upper),width=.75,alpha=.8)+
  geom_hline(yintercept = .5, linetype='dashed',alpha=.25)+
  facet_grid(display~.)+
  scale_x_continuous(breaks=c(2,5,9,14),limits=c(1.5,14.5),labels=c("2%","5%","9%","14%"))+
  scale_y_continuous(limits=c(.4,.6))+
  ggsci::scale_color_startrek(name="choice")+
  labs(y="p(t)")+
  ggthemes::theme_few()
ggsave(filename=here("analyses","plots","2afc_tc_choices.jpeg"),width=4,height=5)

# catch analyses ========================================================================================
# nothing too substantive here, just checking how people did on catch trials
# poor performing ppts (<80% correct) already removed prior to this analysis
catch_correct_indiv <- catch  %>% 
  group_by(sub_n,display) %>%
  summarise(correct_prop=sum(correct)/n()) %>%
  ungroup()

catch_correct_indiv_hist <- catch_correct_indiv %>%
  ggplot(aes(correct_prop))+
  geom_histogram(fill="lightblue",col="black")+
  facet_grid(.~display,scales="free_y")+
  ggthemes::theme_few()+
  theme(text=element_text(size=18))
catch_correct_indiv_hist

catch %>%
  group_by(display) %>%
  summarise(correct_prop=sum(correct)/n()) %>%
  ungroup()

catch %>%
  group_by(sub_n,display) %>%
  summarise(correct=mean(correct)) %>%
  ungroup() %>%
  summarise(m_correct=mean(correct),
            s_correct=sd(correct),.by = display)



# if(!file_exists(fit_file)){
#   m <- stan_glmer(discrim~tdd*probe+display+tdo+(1+tdd*probe+display+tdo|sub_n),
#                            family=binomial,
#                            data=critical_1,
#                            refresh=250,
#                            cores=n_core,
#                            iter=n_iter)
#   save(m,file=fit_file)
#   m_summary <- summary(m, probs = c(.025,.975))
#   save(m_summary,file=here("analyses","stats","m_summary.RData"))
# }else{
#   load(fit_file)
# }
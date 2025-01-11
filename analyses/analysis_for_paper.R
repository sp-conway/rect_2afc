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

read_dat <- function(f){
  dd <- read_csv(f) %>%
    mutate(display=recode(display, spektor="triangle", trueblood="horizontal"),
           across(c(sub_n,tdo,display,choice),as.factor))
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

# critical TD trial analyses ========================================================================================
# baseline display level should be triangle for model
critical_1 <- critical %>%
  filter(tdd!=0 & probe!="tc") %>%
  mutate(tdd=as.factor(tdd),
         probe=factor(probe,levels=c("dc","td"))) %>%
  mutate(display=fct_relevel(display,c("triangle","horizontal")),
         discrim=case_when(
           choice=="t" & probe=="td"~1,
           choice=="c" & probe=="dc"~1,
           choice=="d"~0
         )) %>%
  arrange(sub_n,trial_number)

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
critical_unique <- critical_for_model %>%
  distinct(sub_n_new,tdo,display,probe,tdd_5,tdd_9,tdd_14)
n_unique <- nrow(critical_unique)

stan_data <- list(
  n_subs=n_subs,
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
debug_model <- T
prefix <- "m1"
model_dir <- path(here("analyses","stan",prefix))
stan_model_code <- path(model_dir,glue("{prefix}.stan"))
fit_file <- path(model_dir,glue("{prefix}_fit.RData"))


if(debug_model){
  n_iter <- 100
  n_core <- 1
  n_chain <- 2
  n_warm <- 10
  to_save <- F
}else{
  # number of iterations for each model per core (for MCMC)
  n_iter <- 4000
  n_core <- n_chain <- 4
  n_warm <- 500
  to_save <- T
}

if(!file_exists(fit_file) | debug_model){
  model_compiled <- stan_model(stan_model_code)
  fit <- sampling(model_compiled,
                  data=stan_data,
                  chains=n_chain,
                  warmup=n_warm,
                  iter=n_iter)
  if(!debug_model) save(fit,file=fit_file)
  fit_summary <- summary(fit,probs=c(.025,.975))
  if(!debug_model) save(summary, file=path(model_dir,glue("{prefix}_summary.RData")))
  diagnostics <- get_sampler_params(fit)
  if(!debug_model) save(diagnostics, file=path(model_dir,glue("{prefix}_diag.RData")))
}else{
  load(fit_file)
}

# check modeling ===============================================================================
m_fit <- as.matrix(fit)
color_scheme_set("red")

try({
  mcmc_trace(m_fit, c("lp__"))
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_lp___trace.jpeg")),width=5,height=4)
})

try({
  mcmc_trace(m_fit, c("b_0"))
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_0_trace.jpeg")),width=5,height=4)
})

try({
  mcmc_trace(m_fit, c("b_tdd_5","b_tdd_9","b_tdd_14"))
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_tdd_trace.jpeg")),width=5,height=4)
})

try({
  mcmc_trace(m_fit, c("b_td"))
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_td_trace.jpeg")),width=5,height=4)
})

try({
  mcmc_trace(m_fit, c("b_w"))
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_w_trace.jpeg")),width=5,height=4)
})

try({
  mcmc_trace(m_fit, c("b_h"))
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_h_trace.jpeg")),width=5,height=4)
})

try({
  mcmc_trace(m_fit, c("b_tdd_5_X_td","b_tdd_9_X_td","b_tdd_14_X_td"))
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_tdd_X_td_trace.jpeg")),width=5,height=4)
})

try({
  mcmc_trace(m_fit, c("sigma_b_0_s","sigma_b_w_s","sigma_b_h_s","sigma_b_td_s","sigma_b_tdd_5_s","sigma_b_tdd_9_s","sigma_b_tdd_14_s","sigma_b_tdd_5_X_td_s","sigma_b_tdd_9_X_td_s","sigma_b_tdd_14_X_td_s"))
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_tdd_X_td_trace.jpeg")),width=5,height=4)
})

try({
  mcmc_hist(m_fit, c("lp__"))
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_lp___hist.jpeg")),width=5,height=4)
})

try({
  mcmc_hist(m_fit, c("b_0"))
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_0_hist.jpeg")),width=5,height=4)
})

try({
  mcmc_hist(m_fit, c("b_tdd_5","b_tdd_9","b_tdd_14"))
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_tdd_hist.jpeg")),width=5,height=4)
})

try({
  mcmc_hist(m_fit, c("b_td"))
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_td_hist.jpeg")),width=5,height=4)
})

try({
  mcmc_hist(m_fit, c("b_w"))
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_w_hist.jpeg")),width=5,height=4)
})

try({
  mcmc_hist(m_fit, c("b_h"))
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_h_hist.jpeg")),width=5,height=4)
})

try({
  mcmc_hist(m_fit, c("b_tdd_5_X_td","b_tdd_9_X_td","b_tdd_14_X_td"))
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_tdd_X_td_hist.jpeg")),width=5,height=4)
})


try({
  mcmc_hist(m_fit, c("sigma_b_0_s","sigma_b_w_s","sigma_b_h_s","sigma_b_td_s","sigma_b_tdd_5_s","sigma_b_tdd_9_s","sigma_b_tdd_14_s","sigma_b_tdd_5_X_td_s","sigma_b_tdd_9_X_td_s","sigma_b_tdd_14_X_td_s"))
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_tdd_X_td_hist.jpeg")),width=5,height=4)
})

try({
  mcmc_dens(m_fit, c("lp__"))
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_lp___dens.jpeg")),width=5,height=4)
})

try({
  mcmc_dens(m_fit, c("b_0"))
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_0_dens.jpeg")),width=5,height=4)
})

try({
  mcmc_dens(m_fit, c("b_tdd_5","b_tdd_9","b_tdd_14"))
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_tdd_dens.jpeg")),width=5,height=4)
})

try({
  mcmc_dens(m_fit, c("b_td"))
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_td_dens.jpeg")),width=5,height=4)
})

try({
  mcmc_dens(m_fit, c("b_w"))
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_w_dens.jpeg")),width=5,height=4)
})

try({
  mcmc_dens(m_fit, c("b_h"))
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_h_dens.jpeg")),width=5,height=4)
})

try({
  mcmc_dens(m_fit, c("b_tdd_5_X_td","b_tdd_9_X_td","b_tdd_14_X_td"))
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_tdd_X_td_dens.jpeg")),width=5,height=4)
})

try({
  mcmc_dens(m_fit, c("sigma_b_0_s","sigma_b_w_s","sigma_b_h_s","sigma_b_td_s","sigma_b_tdd_5_s","sigma_b_tdd_9_s","sigma_b_tdd_14_s","sigma_b_tdd_5_X_td_s","sigma_b_tdd_9_X_td_s","sigma_b_tdd_14_X_td_s"))
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_tdd_X_td_dens.jpeg")),width=5,height=4)
})


try({
  mcmc_dens_chains(m_fit, c("lp__"))
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_lp___dens_chains.jpeg")),width=5,height=4)
})

try({
  mcmc_dens_chains(m_fit, c("b_0"))
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_0_dens_chains.jpeg")),width=5,height=4)
})

try({
  mcmc_dens_chains(m_fit, c("b_tdd_5","b_tdd_9","b_tdd_14"))
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_tdd_dens_chains.jpeg")),width=5,height=4)
})

try({
  mcmc_dens_chains(m_fit, c("b_td"))
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_td_dens_chains.jpeg")),width=5,height=4)
})

try({
  mcmc_dens_chains(m_fit, c("b_w"))
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_w_dens_chains.jpeg")),width=5,height=4)
})

try({
  mcmc_dens_chains(m_fit, c("b_h"))
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_h_dens_chains.jpeg")),width=5,height=4)
})

try({
  mcmc_dens_chains(m_fit, c("b_tdd_5_X_td","b_tdd_9_X_td","b_tdd_14_X_td"))
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_tdd_X_td_dens_chains.jpeg")),width=5,height=4)
})

try({
  mcmc_dens_chains(m_fit, c("sigma_b_0_s","sigma_b_w_s","sigma_b_h_s","sigma_b_td_s","sigma_b_tdd_5_s","sigma_b_tdd_9_s","sigma_b_tdd_14_s","sigma_b_tdd_5_X_td_s","sigma_b_tdd_9_X_td_s","sigma_b_tdd_14_X_td_s"))
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_tdd_X_td_dens_chains.jpeg")),width=5,height=4)
})

try({
  mcmc_pairs(m_fit, c("b_tdd_5","b_tdd_9","b_tdd_14"))
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_tdd_pairs.jpeg")),width=5,height=4)
})

try({
  mcmc_pairs(m_fit, c("sigma_b_0_s","sigma_b_w_s","sigma_b_h_s","sigma_b_td_s","sigma_b_tdd_5_s","sigma_b_tdd_9_s","sigma_b_tdd_14_s","sigma_b_tdd_5_X_td_s","sigma_b_tdd_9_X_td_s","sigma_b_tdd_14_X_td_s"))
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_tdd_X_td_pairs.jpeg")),width=5,height=4)
})

try({
  mcmc_pairs(m_fit, c("b_tdd_5_X_td","b_tdd_9_X_td","b_tdd_14_X_td"))
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_tdd_X_td_pair.jpeg")),width=5,height=4)
})

try({
  mcmc_intervals(m_fit, regex_pars = "^b_0_s",prob=.5,prob_outer = .95)+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())+
    labs(y="participant",x="estimate")
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_0_s.jpeg")),width=5,height=4)
})

try({
  mcmc_intervals(m_fit, regex_pars = "^b_h_s",prob=.5,prob_outer = .95)+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())+
    labs(y="participant",x="estimate")
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_h_s.jpeg")),width=5,height=4)
})

try({
  mcmc_intervals(m_fit, regex_pars = "^b_w_s",prob=.5,prob_outer = .95)+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())+
    labs(y="participant",x="estimate")
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_w_s.jpeg")),width=5,height=4)
})

try({
  mcmc_intervals(m_fit, regex_pars = "^b_td_s",prob=.5,prob_outer = .95)+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())+
    labs(y="participant",x="estimate")
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_td_s.jpeg")),width=5,height=4)
})

try({
  mcmc_intervals(m_fit, regex_pars = "^b_tdd_5_s",prob=.5,prob_outer = .95)+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())+
    labs(y="participant",x="estimate")
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_tdd_5_s.jpeg")),width=5,height=4)
})

try({
  mcmc_intervals(m_fit, regex_pars = "^b_tdd_9_s",prob=.5,prob_outer = .95)+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())+
    labs(y="participant",x="estimate")
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_tdd_9_s.jpeg")),width=5,height=4)
})

try({
  mcmc_intervals(m_fit, regex_pars = "^b_tdd_14_s",prob=.5,prob_outer = .95)+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())+
    labs(y="participant",x="estimate")
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_tdd_14_s.jpeg")),width=5,height=4)
})

try({
  mcmc_intervals(m_fit, regex_pars = "^b_tdd_5_X_td_s",prob=.5,prob_outer = .95)+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())+
    labs(y="participant",x="estimate")
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_tdd_5_X_td_s.jpeg")),width=5,height=4)
})

try({
  mcmc_intervals(m_fit, regex_pars = "^b_tdd_9_X_td_s",prob=.5,prob_outer = .95)+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())+
    labs(y="participant",x="estimate")
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_tdd_9_X_td_s.jpeg")),width=5,height=4)
})

try({
  mcmc_intervals(m_fit, regex_pars = "^b_tdd_14_X_td_s",prob=.5,prob_outer = .95)+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())+
    labs(y="participant",x="estimate")
  if(!debug_model) ggsave(filename=path(model_dir,glue("{prefix}_b_tdd_14_X_td_s.jpeg")),width=5,height=4)
})

try({
  mcmc_intervals(m_fit,pars=c("b_0","b_w","b_h","b_td","b_tdd_5","b_tdd_9","b_tdd_14",
                              "b_tdd_5_X_td","b_tdd_9_X_td","b_tdd_14_X_td"),prob=.5,prob_outer = .95)
})


# get model predictions ================================================================================================================================================
p <- extract(fit, pars="p")$p
if(!debug) save(file=path(model_dir,glue("{prefix}_p.RData")))
# tc trials =================================================================================================================================
critical_tc_means <- critical %>%
  filter(probe=="tc") %>%
  mutate(choose_t=choice=="t") %>%
  group_by(sub_n, tdd, display) %>%
  mutate(prop=mean(choose_t)) %>%
  ungroup() %>%
  group_by(tdd, display) %>%
  summarise(mean_prop=mean(prop),
            se_prop=sd(prop)/sqrt(n()),
            ci_lower=mean_prop-qt(.975,n()-1)*se_prop,
            ci_upper=mean_prop+qt(.975,n()-1)*se_prop) %>%
  ungroup() %>%
  mutate(tdd=as.numeric(as.character(tdd))*100)

critical_tc_means %>%
  ggplot(aes(tdd, mean_prop))+
  geom_point(alpha=.8)+
  geom_line(alpha=.8)+
  geom_errorbar(aes(ymin=ci_lower,ymax=ci_upper),width=.75,alpha=.8)+
  geom_hline(yintercept = .5, linetype='dashed',alpha=.25)+
  facet_grid(display~exp)+
  scale_x_continuous(breaks=c(0,2,5,9,14))+
  scale_y_continuous(limits=c(.33,.66))+
  ggsci::scale_color_startrek(name="choice")+
  labs(x="target-decoy distance",y="prop choose target")+
  ggthemes::theme_few()+
  theme(text=element_text(size=18),
        plot.title=element_text(hjust=0.5),
        legend.position = "inside",
        legend.position.inside = c(.1,.85))
ggsave(filename=here("analyses","plots","tc_choices.jpeg"),width=5,height=6)

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
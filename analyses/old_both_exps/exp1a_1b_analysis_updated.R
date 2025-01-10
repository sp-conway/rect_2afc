# exp1a_1b_analysis.R
# Sean Conway

# setup ==============================================================
rm(list=ls())

# packages
library(tidyverse)
library(glue)
library(ggsci)
library(tidybayes)
library(fs)
library(rstanarm)
library(posterior)
library(bayesplot)
library(here)
library(patchwork)

read_dat <- function(f,exp){
  dd <- read_csv(f) %>%
    mutate(exp=exp,
           display=recode(display, spektor="triangle", trueblood="horizontal"),
           across(c(sub_n,exp,tdo,display,choice),as.factor))
  return(dd)
}
# read in datasets
exp1a <- here("exp1a","data","aggregated","rect_exp1a_FIXEDSPEKTOR_aggregated.csv") %>%
  read_dat(exp="1a")
exp1b <- here("exp1b","data","aggregated","rect_exp1b_aggregated.csv") %>%
  read_dat(exp="1b")

# number of ppts
get_subs <- function(d){
  p <- length(unique(d$sub_n))
  print(p)
  return(p)
}

# number of subjects
n_subs_exp1a <- get_subs(exp1a)
n_subs_exp1b <- get_subs(exp1b)

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

# catch analyses ====================================================
# nothing too substantive here, just checking how people did on catch trials
# poor performing ppts (<80% correct) already removed prior to this analysis
catch_correct_indiv <- catch  %>% 
  group_by(exp,sub_n,display) %>%
  summarise(correct_prop=sum(correct)/n()) %>%
  ungroup()

catch_correct_indiv_hist <- catch_correct_indiv %>%
  ggplot(aes(correct_prop))+
  geom_histogram(fill="lightblue",col="black")+
  facet_grid(exp~display,scales="free_y")+
  ggthemes::theme_few()+
  theme(text=element_text(size=18))
catch_correct_indiv_hist

catch %>%
  group_by(exp,display) %>%
  summarise(correct_prop=sum(correct)/n()) %>%
  ungroup()

# critical analysis ==================================================================
# td discrimn, target decoy distance by orientation
critical_agg_choice_display_tdd_tdo <- critical %>%
  count(exp,display,probe,tdd,tdo,choice) %>%
  group_by(exp,display,probe,tdo,tdd) %>%
  mutate(prop=n/sum(n)) %>%
  ungroup() %>%
  mutate(probe1=str_sub(probe,1,1),
         probe2=str_sub(probe,2,2)) %>%
  rowwise() %>%
  filter(choice==probe1 | choice==probe2) %>%
  ungroup()

# discriminability, COLLAPSED ACROSS ORIENTATION
# this is the crucial analysis and the only one we report in paper
critical_agg_choice_display_tdd <- critical %>%
  count(exp,display,probe,tdd,choice) %>%
  group_by(exp,display,probe,tdd) %>%
  mutate(prop=n/sum(n)) %>%
  ungroup() %>%
  mutate(probe1=str_sub(probe,1,1),
         probe2=str_sub(probe,2,2)) %>%
  rowwise() %>%
  filter(choice==probe1 | choice==probe2) %>%
  ungroup()

# critical trial analysis =============================================================================================================
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
         ))

critical_1a <- filter(critical_1, exp=="1a")
critical_1b <- filter(critical_1, exp=="1b")

critical_1b_mean_discrim <- critical_1b %>%
  group_by(sub_n,tdd,probe,display) %>%
  summarise(p_discrim=mean(discrim)) %>%
  group_by(tdd, probe, display) %>%
  summarise(mean_discrim=mean(p_discrim),
            se=sd(p_discrim)/sqrt(n()),
            ci_lower=mean_discrim-qt(.975,n()-1)*se,
            ci_upper=mean_discrim+qt(.975,n()-1)*se) %>%
  ungroup() %>%
  mutate(exp="1b")
critical_1a_mean_discrim <- critical_1a %>%
  group_by(sub_n,tdd,probe,display) %>%
  summarise(p_discrim=mean(discrim)) %>%
  group_by(tdd, probe, display) %>%
  summarise(mean_discrim=mean(p_discrim),
            se=sd(p_discrim)/sqrt(n()),
            ci_lower=mean_discrim-qt(.975,n()-1)*se,
            ci_upper=mean_discrim+qt(.975,n()-1)*se) %>%
  ungroup() %>%
  mutate(exp="1a")

critical_mean_discrim <- bind_rows(critical_1a_mean_discrim,critical_1b_mean_discrim) %>%
  mutate(tdd=as.numeric(as.character(tdd))*100)

critical_mean_discrim %>%
  ggplot(aes(tdd, mean_discrim,col=probe))+
  geom_point(alpha=.8)+
  geom_line(alpha=.8)+
  geom_errorbar(aes(ymin=ci_lower,ymax=ci_upper),width=.75,alpha=.8)+
  facet_grid(display~exp)+
  scale_x_continuous(breaks=c(2,5,9,14))+
  scale_y_continuous(limits=c(.5,1))+
  ggsci::scale_color_startrek(name="choice")+
  labs(x="target-decoy distance",y="proability discriminate")+
  ggthemes::theme_few()+
  theme(text=element_text(size=18),
        plot.title=element_text(hjust=0.5),
        legend.position = "inside",
        legend.position.inside = c(.1,.85))
ggsave(filename=here("analyses","plots","critical_discrim_1a1b.jpeg"),width=6,height=5)


# Statistical modeling ========================================================================

# number of iterations for each model per core (for MCMC)
n_iter <- 2000
n_core <- 4

fit_file <- here("analyses","stats","m1b_intxn.RData")
if(!file_exists(fit_file)){
  m_1b_intxn <- stan_glmer(discrim~tdd*probe+display+tdo+(1|sub_n),
                           family=binomial,
                           data=critical_1b,
                           refresh=250,
                           cores=n_core,
                           iter=n_iter)
  save(m_1b_intxn,file=fit_file)
}else{
  load(fit_file)
}

# check modeling ===============================================================================
m_1b_intxn_fit <- as.matrix(m_1b_intxn)
color_scheme_set("red")
prefix <- "m1b_intxn"

mcmc_trace(m_1b_intxn, c("tdd0.05","tdd0.09","tdd0.14"))
ggsave(filename=here("analyses","stats",glue("{prefix}_tdd_trace.jpeg")),width=5,height=4)

mcmc_trace(m_1b_intxn, c("probetd"))
ggsave(filename=here("analyses","stats",glue("{prefix}_probetd_trace.jpeg")),width=5,height=4)

mcmc_trace(m_1b_intxn, c("(Intercept)"))
ggsave(filename=here("analyses","stats",glue("{prefix}_fixed_int_trace.jpeg")),width=5,height=4)

mcmc_trace(m_1b_intxn, c("displayhorizontal"))
ggsave(filename=here("analyses","stats",glue("{prefix}_displayhorizontal_trace.jpeg")),width=5,height=4)

mcmc_trace(m_1b_intxn, c("tdd0.05:probetd","tdd0.09:probetd","tdd0.14:probetd"))
ggsave(filename=here("analyses","stats",glue("{prefix}_tddXprobe_trace.jpeg")),width=5,height=4)

mcmc_hist(m_1b_intxn, c("tdd0.05","tdd0.09","tdd0.14"))
ggsave(filename=here("analyses","stats",glue("{prefix}_tdd_hist.jpeg")),width=5,height=4)

mcmc_hist(m_1b_intxn, c("probetd"))
ggsave(filename=here("analyses","stats",glue("{prefix}_probetd_hist.jpeg")),width=5,height=4)

mcmc_hist(m_1b_intxn, c("(Intercept)"))
ggsave(filename=here("analyses","stats",glue("{prefix}_fixed_int_hist.jpeg")),width=5,height=4)

mcmc_hist(m_1b_intxn, c("displayhorizontal"))
ggsave(filename=here("analyses","stats",glue("{prefix}_displayhorizontal_hist.jpeg")),width=5,height=4)

mcmc_hist(m_1b_intxn, c("tdd0.05:probetd","tdd0.09:probetd","tdd0.14:probetd"))
ggsave(filename=here("analyses","stats",glue("{prefix}_tddXprobe_hist.jpeg")),width=5,height=4)

x=m_1b_intxn %>% 
  as_tibble() %>%
  pivot_longer(contains("sub_n")) %>%
  mutate(sub_n=str_extract(name,"(?<=sub_n:)[:digit:]{1,100}")) %>%
  select(-name) %>% 
# critical tc trials ====================================================================================
critical_tc_means <- critical %>%
  filter(probe=="tc") %>%
  mutate(choose_t=choice=="t") %>%
  group_by(exp,sub_n, tdd, display) %>%
  mutate(prop=mean(choose_t)) %>%
  ungroup() %>%
  group_by(exp,tdd, display) %>%
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

# # compute means for later ========================================
# here("analyses","output")
# comp_means <- function(data,fpath){
#   coll <- data %>%
#     group_by(exp,sub_n,display,tdd) %>%
#     summarise(p=mean(tchoice)) %>%
#     group_by(exp,display,tdd) %>%
#     summarise(m_pt=mean(p)) %>%
#     ungroup()
#   write_csv(coll,path(fpath,"data_means_discrim_TD_collapsed.csv"))
#   sep <- data %>%
#     group_by(exp,sub_n,display,tdd,tdo) %>%
#     summarise(p=mean(tchoice)) %>%
#     group_by(exp,display,tdd,tdo) %>%
#     summarise(m_pt=mean(p)) %>%
#     ungroup()
#   write_csv(sep,path(fpath,"data_means_discrim_TD.csv"))
# }
# comp_means(bind_rows(critical_td_1a,critical_td_1b),
#            here("analyses","output"))

# m_1a_intxn <- stan_glmer(discrim~tdd*probe+display+tdo+(1|sub_n),
#                          family=binomial,
#                          data=critical_1a,
#                          cores=n_core,
#                          iter=n_iter)
# MORE OLD
# summary(m_1a_intxn,probs = c(.025, .975))


# m_1a_intxn_fit <- as.matrix(m_1a_intxn)
# save(m_1a_intxn, m_1b_intxn,
#      m_1a, m_1b, file=here("tmp.RData"))
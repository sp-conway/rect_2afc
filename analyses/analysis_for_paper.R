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

# catch analyses ====================================================
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

# Statistical modeling ========================================================================

# number of iterations for each model per core (for MCMC)
n_iter <- 4000
n_core <- 4

fit_file <- here("analyses","stats","m1.RData")
if(!file_exists(fit_file)){
  m <- stan_glmer(discrim~tdd*probe+display+tdo+(1+tdd*probe+display+tdo|sub_n),
                           family=binomial,
                           data=critical_1,
                           refresh=250,
                           cores=n_core,
                           iter=n_iter)
  save(m,file=fit_file)
  m_summary <- summary(m, probs = c(.025,.975))
  save(m_summary,file=here("analyses","stats","m_summary.RData"))
}else{
  load(fit_file)
}

# check modeling ===============================================================================
m_fit <- as.matrix(m)
color_scheme_set("red")
prefix <- "m"

mcmc_trace(m, c("tdd0.05","tdd0.09","tdd0.14"))
ggsave(filename=here("analyses","stats",glue("{prefix}_tdd_trace.jpeg")),width=5,height=4)

mcmc_trace(m, c("probetd"))
ggsave(filename=here("analyses","stats",glue("{prefix}_probetd_trace.jpeg")),width=5,height=4)

mcmc_trace(m, c("(Intercept)"))
ggsave(filename=here("analyses","stats",glue("{prefix}_fixed_int_trace.jpeg")),width=5,height=4)

mcmc_trace(m, c("displayhorizontal"))
ggsave(filename=here("analyses","stats",glue("{prefix}_displayhorizontal_trace.jpeg")),width=5,height=4)

mcmc_trace(m, c("tdd0.05:probetd","tdd0.09:probetd","tdd0.14:probetd"))
ggsave(filename=here("analyses","stats",glue("{prefix}_tddXprobe_trace.jpeg")),width=5,height=4)

mcmc_hist(m, c("tdd0.05","tdd0.09","tdd0.14"))
ggsave(filename=here("analyses","stats",glue("{prefix}_tdd_hist.jpeg")),width=5,height=4)

mcmc_hist(m, c("probetd"))
ggsave(filename=here("analyses","stats",glue("{prefix}_probetd_hist.jpeg")),width=5,height=4)

mcmc_hist(m, c("(Intercept)"))
ggsave(filename=here("analyses","stats",glue("{prefix}_fixed_int_hist.jpeg")),width=5,height=4)

mcmc_hist(m, c("displayhorizontal"))
ggsave(filename=here("analyses","stats",glue("{prefix}_displayhorizontal_hist.jpeg")),width=5,height=4)

mcmc_hist(m, c("tdd0.05:probetd","tdd0.09:probetd","tdd0.14:probetd"))
ggsave(filename=here("analyses","stats",glue("{prefix}_tddXprobe_hist.jpeg")),width=5,height=4)

# ppc =====================================================================================================================
mppc <- posterior_predict(m)

ppc_dens_overlay(critical_1$discrim,mppc)



# critical tc trials =================================================================================================================================
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

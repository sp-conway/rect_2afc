rm(list=ls())
library(tidyverse)
library(here)
library(fs)
library(glue)

n <- 14

fits <- vector("list",n)
for(model in 1:n){
  f <- here("analyses","stan",glue("m{model}"),glue("m{model}_summary.RData"))
  load(f)
  tmp <- as.data.frame(fit_summary$summary)
  tmp$param <- rownames(tmp)
  tmp <- filter(tmp, str_detect(param,"_s|p\\[",negate=T))
  tmp$m <- model
  fits[[model]] <- tmp
  rm(fit_summary)
  rm(tmp)
}
fits_df <- list_rbind(fits)
pl <- function(f,regex){
  f %>%
    filter(str_detect(param,regex)) %>%
    select(param,mean,`2.5%`,`97.5%`,m) %>%
    mutate(m=as.factor(m))%>%
    ggplot(aes(m,mean))+
    geom_point()+
    # geom_path(linetype="dashed")+
    geom_errorbar(aes(ymin=`2.5%`,ymax=`97.5%`))+
    labs(x="model iteration")+
    facet_grid(param~.,scales="free_y")+
    ggthemes::theme_few()
}
pl(fits_df,"tdd_[:digit:]{1,2}_X_td")
pl(fits_df,"b_0$|b_h$|b_w$|b_td$")
pl(fits_df,"b_tdd_")
pl(fits_df, "b_w$")
pl(fits_df,"b_h$")
pl(fits_df,"b_td$")


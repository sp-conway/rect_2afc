# get all covariance matrices from model =======================================================s <- c(1,1,1)
rm(list=ls())
library(tidyverse)
library(mvtnorm)
library(here)

# control params
N <- 50000
which_choice <- c("td","dc")
tdd <- 1:4

# correlations
cors <- matrix(c(1, .55, .68,
                 .55, 1, .55,
                 .68, .55, 1),
               nrow=3,byrow=T)

# means
dec <- c(.8,.7,.6,.5)
tar_comp <- 1

# standar deviations
s <- c(.8,.8,.8)

# compute variance covariance matrices
cv <- cors * (s %*% t(s) )

# function to draw and simulate choice
sim_choice <- function(N,mu,cv,which_choice){
  samp <- rmvnorm(N,mu,cv)
  if(which_choice=="td"){
    p <- sum(samp[,1]>samp[,3])/N
  }else if(which_choice=="dc"){
    p <- sum(samp[,2]>samp[,3])/N
  }
  return(p)
}

p <- array(NA, dim=c(2,4))
for(i in 1:length(which_choice)){
  for(j in 1:length(tdd)){
    p[i,j] <- sim_choice(N, 
                         c(tar_comp, tar_comp, dec[j]),
                         cv,
                         which_choice[i])
  }
}

tibble(
  p=as.vector(p),
  tdd=c(2,2,5,5,9,9,14,14),
  probe=c("td","cd","td","cd","td","cd","td","cd")
) %>%
  ggplot(aes(tdd, p,col=probe,shape=probe))+
  geom_point()+
  geom_line()+
  scale_x_continuous(breaks=c(2,5,9,14))+
  scale_y_continuous(limits=c(.5,.8))+
  scale_color_manual(values=c("gray","black"),name="pair")+
  scale_shape_manual(values=c(17,15),name="pair")+
  labs(x="tdd",y="p(discriminate)")+
  ggthemes::theme_few()+
  theme(text=element_text(size=18),
        plot.title=element_text(hjust=0.5),
        legend.position="inside",
        legend.position.inside = c(.8,.2))
ggsave(filename=here("analyses","plots","model_sim_2afc.jpeg"),width=4,height=4)

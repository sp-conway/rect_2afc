data {
  int<lower=0> n_subs; // number of subjects
  int<lower=0> n_trials; // number of TOTAL trials across all subjects
  int<lower=0> n_trials_unique;
  array[n_trials] int sub_n_new;  // vector of all subject ids for indexing of random effects
  vector[n_trials] tdo; // target-decoy orientation. 1=w, 0=h
  vector[n_trials] display; // display. 1=horizontal, 0=triangle
  vector[n_trials] probe; // probe. 1=target-decoy, 0=competitor-decoy
  vector[n_trials] tdd_5; // dummy variable for target-decoy distance=5
  vector[n_trials] tdd_9; // dummy variable for target-decoy distance=9
  vector[n_trials] tdd_14; // dummy variable for target-decoy distance=14
  array[n_trials] int discrim; // whether participant discriminates, i.e., does not choose decoy
  
  // unique trials for prediction
  array[n_trials_unique] int sub_n_new_unique;  // vector of all subject ids for indexing of random effects
  vector[n_trials_unique] tdo_unique; // target-decoy orientation. 1=w, 0=h
  vector[n_trials_unique] display_unique; // display. 1=horizontal, 0=triangle
  vector[n_trials_unique] probe_unique; // probe. 1=target-decoy, 0=competitor-decoy
  vector[n_trials_unique] tdd_5_unique; // dummy variable for target-decoy distance=5
  vector[n_trials_unique] tdd_9_unique; // dummy variable for target-decoy distance=9
  vector[n_trials_unique] tdd_14_unique; // dummy variable for target-decoy distance=5
}

parameters {
  // fixed effects
  // fixed intercept
  real b_0;
  
  // fixed effect of wide
  real b_w;
  
  // fixed effect of horizontal
  real b_h;
  
  // fixed effect of td comparison
  real b_td;
  
  // fixed effect(s) of TDD
  real b_tdd_5;
  real b_tdd_9;
  real b_tdd_14;
  
  // fixed interaction(s) between tdd and probe
  real b_tdd_5_X_td;
  real b_tdd_9_X_td;
  real b_tdd_14_X_td;
  
  // random subject-level intercepts
  vector[n_subs] b_0_s;

  // Std. Dev. of random-subject level intercepts
  real<lower=0> sigma_b_0_s;
  
  // std dev of effect of wide
  real<lower=0> sigma_bw_s;
  
  // std dev of effect of horizontal
  real<lower=0> sigma_bh_s;
  
  // std dev of effect of td compariosn
  real<lower=0> sigma_btd_s;
  
  // std dev of effect of tdd=5
  real<lower=0> sigma_tdd5_s;
  
  // std dev of effect of tdd=9
  real<lower=0> sigma_tdd9_s;
  
  // std dev of effect of tdd=14
  real<lower=0> sigma_tdd14_s;
  
  // std dev of effect of tdd5 x td
  real<lower=0> sigma_tdd5xtd_s;
  
  // std dev of effect of tdd9 x td
  real<lower=0> sigma_tdd9xtd_s;
  
  // std dev of effect of tdd14 x td
  real<lower=0> sigma_tdd14xtd_s;
  
  // random effect of wide
  vector[n_subs] b_w_s;
  
  // random effect of horizontal
  vector[n_subs] b_h_s;
  
  // random effect of td comparison
  vector[n_subs] b_td_s;

  // random effect of distance=5
  vector[n_subs] b_tdd_5_s;
  
  // random effect of distance=9
  vector[n_subs] b_tdd_9_s;
  
  // random effect of distance=14
  vector[n_subs] b_tdd_14_s;

  // random interaction(s) between tdd and probe
  vector[n_subs] b_tdd_5_X_td_s;
  vector[n_subs] b_tdd_9_X_td_s;
  vector[n_subs] b_tdd_14_X_td_s;
  
}

model {
  // fixed effects
  b_0 ~ normal(0,5);
  b_w ~ normal(0,5);
  b_h ~ normal(0,5);
  b_td ~ normal(0,5);
  b_tdd_5 ~ normal(0,5);
  b_tdd_9 ~ normal(0,5);
  b_tdd_14 ~ normal(0,5);
  b_tdd_5_X_td ~ normal(0,2.5);
  b_tdd_9_X_td ~ normal(0,2.5);
  b_tdd_14_X_td~ normal(0,2.5);
  
  // random effects - standard deviations
  sigma_b_0_s ~ lognormal(0,.5);
  sigma_bw_s ~ lognormal(0,.5);
  sigma_btd_s ~ lognormal(0,.5);
  sigma_bh_s ~ lognormal(0,.5);
  sigma_btdd5_s ~ lognormal(0,.5);
  sigma_btdd9_s ~ lognormal(0,.5);
  sigma_btdd14_s ~ lognormal(0,.5);
  sigma_btdd5xtd_s ~ lognormal(0,.5);
  sigma_btdd9xtd_s ~ lognormal(0,.5);
  sigma_btdd14_s ~ lognormal(0,.5);  
  
  // random effects - coefficients
  b_0_s ~ normal(0,sigma_b_0_s);
  b_w_s ~ normal(0,sigma_bw_s);
  b_td_s ~ normal(0,sigma_btd_s);
  b_h_s ~ normal(0,sigma_bh_s);
  b_tdd_5_s ~ normal(0,sigma_btdd5_s);
  b_tdd_9_s ~ normal(0,sigma_btdd9_s);
  b_tdd_14_s ~ normal(0,sigma_btdd14_s);
  b_tdd_5_X_td_s ~ normal(0,sigma_btdd5xtd_s);
  b_tdd_9_X_td_s ~ normal(0,sigma_btdd9xtd_s);
  b_tdd_14_X_td_s ~ normal(0,sigma_btdd14xtd_s);

  for(i in 1:n_trials) {
    discrim[i] ~ bernoulli_logit( 
      (b_0 + b_0_s[sub_n_new[i]]) +  // intercept
      (b_w + b_w_s[sub_n_new[i]]) * tdo[i] +  // effect of td wide orientation 
      (b_h + b_h_s[sub_n_new[i]]) * display[i] + // effect of horizontal
      (b_td + b_td_s[sub_n_new[i]]) * probe[i] + // effect of td comparison
      (b_tdd_5 + b_tdd_5_s[sub_n_new[i]]) * tdd_5[i] + // effect of tdd=5
      (b_tdd_9 + b_tdd_9_s[sub_n_new[i]]) * tdd_9[i] + // effect of tdd=9
      (b_tdd_14 + b_tdd_14_s[sub_n_new[i]]) * tdd_14[i] + // effect of tdd=14
      (b_tdd_5_X_td + b_tdd_5_X_td_s[sub_n_new[i]]) * tdd_5[i] * probe[i] +  // effect of tdd=5 & probe=td
      (b_tdd_9_X_td + b_tdd_9_X_td_s[sub_n_new[i]]) * tdd_9[i] * probe[i] + // effect of tdd=9 & probe=td
      (b_tdd_14_X_td + b_tdd_14_X_td_s[sub_n_new[i]]) * tdd_14[i] * probe[i]  // effect of tdd=14 & probe=td
  );
}

}


generated quantities{
  vector[n_trials_unique] p;
  for(i in 1:n_trials_unique){
    p[i] = inv_logit((b_0+b_0_s[sub_n_new_unique[i]]) +  // intercept
                     (b_w+b_w_s[sub_n_new_unique[i]])*tdo_unique[i] +  // effect of td wide orientation 
                     (b_h+b_h_s[sub_n_new_unique[i]])*display_unique[i] + // effect of horizontal
                     (b_td+b_td_s[sub_n_new_unique[i]])*probe_unique[i] + // effect of td comparison
                     (b_tdd_5+b_tdd_5_s[sub_n_new_unique[i]])*tdd_5_unique[i] + // effect of tdd=5
                     (b_tdd_9+b_tdd_9_s[sub_n_new_unique[i]])*tdd_9_unique[i] + // effect of tdd=9
                     (b_tdd_14+b_tdd_14_s[sub_n_new_unique[i]])*tdd_14_unique[i] + // effect of tdd=14
                     (b_tdd_5_X_td + b_tdd_5_X_td_s[sub_n_new_unique[i]])*tdd_5_unique[i]*probe_unique[i] +  // effect of tdd=5 & probe=td
                     (b_tdd_9_X_td + b_tdd_9_X_td_s[sub_n_new_unique[i]])*tdd_9_unique[i]*probe_unique[i] + // effect of tdd=9 & probe=td
                     (b_tdd_14_X_td + b_tdd_14_X_td_s[sub_n_new_unique[i]])*tdd_14_unique[i]*probe_unique[i] );
  }
}

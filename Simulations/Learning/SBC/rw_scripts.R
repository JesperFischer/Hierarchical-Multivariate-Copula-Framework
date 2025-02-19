
learning_update = function(u,e,learningrate){
  return(e + learningrate * (u-e))
}


entropy = function(p){
  return(-(p*log(p)+(1-p)*log(1-p)))
}



generate_trialwise_data_rt_copula = function(u,learningrate,e0,rt_int,rt_beta,rt_ndt,
                                             rt_sd,rho,participant){
  
  
  expectation = array(NA,length(u)+1)
  
  expectation[1] = e0
  for(i in 1: length(u)){
    expectation[i+1] = learning_update(u[i], expectation[i], learningrate)
  }
  expectation = expectation[1:length(u)]
  
  mu_rts = rt_int + rt_beta * entropy(expectation)
  
  
  p_resp = expectation
  
  
  cov_matrix <- matrix(c(1, rho, rho, 1), nrow = 2)
  
  y_mix <- pnorm(matrix(rnorm(length(u)), ncol = 2) %*% chol(cov_matrix))
  
  #rts
  rts <- qlnorm(y_mix[,1], mean = mu_rts, sd = rt_sd) + rt_ndt
  
  #binary resps
  resp <- qbinom(y_mix[,2], size = 1, prob = expectation)
  
  
  return(data.frame(rts,resp,expectation,u,participant, trial = 1:length(expectation), expect = expectation, mu_rts = mu_rts))
}


#ekstra math functions

# A generator function should return a named list containing elements "variables" and "generated"

generator_single <- function(N = 60, S = 20){  # N is the number of data points we are generating
  
  # mu_learningrate = rnorm(1,-1,0.2)
  mu_learningrate = rnorm(1,-1,1)
  mu_e0 =     rnorm(1,0,0.2)
  mu_rt_int =    rnorm(1,-1,0.25)
  mu_rt_beta =   rnorm(1,1.5,0.5)
  mu_rt_sd =     rnorm(1,-1,0.25)
  
  
  tau_learningrate = abs(rnorm(1,0.5,0.1))
  tau_e0 =     abs(rnorm(1,0.5,0.1))
  tau_rt_int =    abs(rnorm(1,0.6,0.2))
  tau_rt_beta =   abs(rnorm(1,0.6,0.2))
  tau_rt_sd =     abs(rnorm(1,0.5,0.1))
  
  learningrate = brms::inv_logit_scaled(rnorm(S,mu_learningrate,tau_learningrate)) 
  e0 = brms::inv_logit_scaled(rnorm(S,mu_e0,tau_e0))
  rt_int = rnorm(S,mu_rt_int,tau_rt_int)  
  rt_beta = rnorm(S,mu_rt_beta,tau_rt_beta)  
  rt_sd = exp(rnorm(S,mu_rt_sd,tau_rt_sd)) 
  
  
  rho = ggdist::rlkjcorr_marginal(S,2,12)
  rt_ndt = rnorm(S,0.3,0.05)  
  
  
  u = c(1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1,
        0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0,
        0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1)
  
  
  threshold = rnorm(S,mu_learningrate,tau_learningrate)  
  slope = exp(rnorm(S,mu_e0,tau_e0))
  rt_int = rnorm(S,mu_rt_int,tau_rt_int)  
  rt_beta = rnorm(S,mu_rt_beta,tau_rt_beta)  
  rt_sd = exp(rnorm(S,mu_rt_sd,tau_rt_sd)) 
  
  
  trial_df = data.frame(learningrate = learningrate, e0 = e0,
                        rt_int = rt_int, rt_beta = rt_beta, rt_sd = rt_sd,
                        rt_ndt = rt_ndt,rho = rho, participant = 1:S) %>% rowwise() %>% 
    summarize((generate_trialwise_data_rt_copula(u, learningrate = learningrate,e0 = e0,participant = participant,
                                                 rt_int = rt_int,rt_beta = rt_beta,rt_sd = rt_sd,rt_ndt = rt_ndt, rho = rho)))
  
  # trial_df %>% ggplot(aes(x = x, y = expectation))+facet_wrap(~participant)+geom_line()
  # trial_df %>% ggplot(aes(x = x, y = rts))+facet_wrap(~participant)+geom_line()
  
  ranges <- data.frame(
    start = seq(1, by = N, length.out = S),
    end = seq(N, by = N, length.out = S)
  )
  
  if(max(trial_df$rts) > 8){
    "Error"
  }
  
  list(
    variables = list(
      mu_learningrate = mu_learningrate,
      mu_e0 = mu_e0,
      mu_rt_int = mu_rt_int,
      mu_rt_beta = mu_rt_beta,
      mu_rt_sd = mu_rt_sd,
      tau_learningrate = tau_learningrate,
      tau_e0 = tau_e0,
      tau_rt_int = tau_rt_int,
      tau_rt_beta = tau_rt_beta,
      tau_rt_sd = tau_rt_sd,
      learningrate = learningrate,
      e0 = e0,
      rt_int = rt_int,
      rt_beta = rt_beta,
      rt_sd = rt_sd,
      rt_ndt = rt_ndt,
      rho = rho
    ),
    generated = list(
      N = nrow(trial_df),
      expect = trial_df$expect,
      mu_rts = trial_df$mu_rts,
      x = trial_df$u,
      binom_y = trial_df$resp,
      RT = (trial_df$rts),
      trial = trial_df$trial,
      minRT = unique(trial_df %>% group_by(participant) %>% summarize(minRT = min(rts)) %>% .$minRT),
      minRT_t = rep(unique(trial_df %>% group_by(participant) %>% summarize(minRT = min(rts)) %>% .$minRT),each = N),
      S_id = trial_df$participant,
      id_ind = ifelse(!duplicated(trial_df$participant), 1, 0),
      starts = ranges$start,
      ends = ranges$end,
      S = length(unique(trial_df$participant))
    )
  )
  
}

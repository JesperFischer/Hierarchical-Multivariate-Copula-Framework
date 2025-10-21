

logistic_psychometric = function(x,threshold,slope,lapse){
  return(lapse+(1-2*lapse)*(1 / (1+exp(-slope * (x - threshold)))))
}

entropy = function(p){
  return(-(p*log(p)+(1-p)*log(1-p)))
}

#shannon_entropy
generate_trial_shannon_entropy_rt = function(x,threshold,slope,lapse,participant,
                                             rt_int,rt_beta, rt_sd, rt_ndt, rho){
  
  expectation = logistic_psychometric(x,threshold, slope, lapse)
  
  mu_rts = rt_int + rt_beta * entropy(expectation)
  
  
  cov_matrix <- matrix(c(1, rho, rho, 1), nrow = 2)
  
  y_mix <- pnorm(matrix(rnorm(length(x)), ncol = 2) %*% chol(cov_matrix))
  
  #rts
  rts <- qlnorm(y_mix[,1], mean = mu_rts, sd = rt_sd) + rt_ndt
  
  #binary resps
  resp <- qbinom(y_mix[,2], size = 1, prob = expectation)
  
  
  return(data.frame(rts,resp,expectation,x,participant, trial = 1:length(expectation),mu_rts = mu_rts))
}



#ekstra math functions

# A generator function should return a named list containing elements "variables" and "generated"

generator_single <- function(N = 60, S = 20){  # N is the number of data points we are generating
  repeat {
      
    mu_threshold = rnorm(1,0,5)
    mu_slope =     rnorm(1,-1.5,0.5)
    mu_lapse =     rnorm(1,-4,1)
    mu_rt_int =    rnorm(1,-1.2,0.25)
    mu_rt_beta =   rnorm(1,1.5,0.5)
    mu_rt_sd =     rnorm(1,-1.2,0.25)
    
    
    tau_threshold = abs(rnorm(1,5,2))
    tau_slope =     abs(rnorm(1,0.6,0.2))
    tau_lapse =     abs(rnorm(1,0.6,0.2))
    tau_rt_int =    abs(rnorm(1,0.6,0.2))
    tau_rt_beta =   abs(rnorm(1,0.6,0.2))
    tau_rt_sd =     abs(rnorm(1,0.4,0.1))
    
    mu_vec <- c(mu_rt_int, mu_rt_sd)
    
    # SDs
    tau_vec <- c(tau_rt_int, tau_rt_sd)
    
    # Correlation matrix with negative correlation, e.g., -0.5
    rho <- -0.5
    Sigma <- matrix(c(
      tau_vec[1]^2,              rho * tau_vec[1] * tau_vec[2],
      rho * tau_vec[1] * tau_vec[2],  tau_vec[2]^2
    ), nrow = 2)
    
    # Joint draws
    joint <- rmvnorm(S, mu_vec, Sigma)
    
  
    rt_int <- joint[,1]
    rt_sd  <- exp(joint[,2]) 
    
    
    threshold = rnorm(S,mu_threshold,tau_threshold)  
    slope = exp(rnorm(S,mu_slope,tau_slope))
    lapse = brms::inv_logit_scaled(rnorm(S,mu_lapse,tau_lapse)) / 2 
    # rt_int = rnorm(S,mu_rt_int,tau_rt_int)  
    rt_beta = rnorm(S,mu_rt_beta,tau_rt_beta)  
    # rt_sd = exp(rnorm(S,mu_rt_sd,tau_rt_sd)) 
    
    
    rho = ggdist::rlkjcorr_marginal(S,2,12)
    rt_ndt = rnorm(S,0.3,0.05)  
    
    x = seq(-40,40, length.out = N)
    
    threshold = rnorm(S,mu_threshold,tau_threshold)  
    slope = exp(rnorm(S,mu_slope,tau_slope))
    lapse = brms::inv_logit_scaled(rnorm(S,mu_lapse,tau_lapse)) / 2 
    rt_int = rnorm(S,mu_rt_int,tau_rt_int)  
    rt_beta = rnorm(S,mu_rt_beta,tau_rt_beta)  
    rt_sd = exp(rnorm(S,mu_rt_sd,tau_rt_sd)) 
    
    
    trial_df = data.frame(threshold = threshold, slope = slope, lapse = lapse,
                          rt_int = rt_int, rt_beta = rt_beta, rt_sd = rt_sd,
                          rt_ndt = rt_ndt,rho = rho, participant = 1:S) %>% rowwise() %>% 
      summarize((generate_trial_shannon_entropy_rt(x, threshold = threshold,slope = slope,lapse = lapse,participant = participant,
                                                   rt_int = rt_int,rt_beta = rt_beta,rt_sd = rt_sd,rt_ndt = rt_ndt, rho = rho)))
    
    # trial_df %>% ggplot(aes(x = x, y = expectation))+facet_wrap(~participant)+geom_line()
    # trial_df %>% ggplot(aes(x = x, y = rts))+facet_wrap(~participant)+geom_line()
    
    ranges <- data.frame(
      start = seq(1, by = N, length.out = S),
      end = seq(N, by = N, length.out = S)
    )
    if (max(trial_df$rts) <= 8) break
  }
    
  list(
    variables = list(
      mu_threshold = mu_threshold,
      mu_slope = mu_slope,
      mu_lapse = mu_lapse,
      mu_rt_int = mu_rt_int,
      mu_rt_beta = mu_rt_beta,
      mu_rt_sd = mu_rt_sd,
      tau_threshold = tau_threshold,
      tau_slope = tau_slope,
      tau_lapse = tau_lapse,
      tau_rt_int = tau_rt_int,
      tau_rt_beta = tau_rt_beta,
      tau_rt_sd = tau_rt_sd,
      threshold = threshold,
      slope = slope,
      lapse = lapse,
      rt_int = rt_int,
      rt_beta = rt_beta,
      rt_sd = rt_sd,
      rt_ndt = rt_ndt,
      rho = rho
    ),
    generated = list(
      N = nrow(trial_df),
      expect = trial_df$expectation,
      mu_rts = trial_df$mu_rts,
      x = trial_df$x,
      binom_y = trial_df$resp,
      RT = (trial_df$rts),
      minRT = unique(trial_df %>% group_by(participant) %>% summarize(minRT = min(rts)) %>% .$minRT),
      S_id = trial_df$participant,
      starts = ranges$start,
      ends = ranges$end,
      S = length(unique(trial_df$participant))
    )
  )
  
}

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


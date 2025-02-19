
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




generate_trial_ddm = function(x,threshold,slope,lapse,participant,
                              alpha,beta, delta, rt_ndt){
  
  expectation = logistic_psychometric(x,threshold, slope, lapse)
  
  deltat = (expectation - (1-expectation)) * delta
  
  data = data.frame()
  for(i in 1:length(x)){
    resp = RWiener::rwiener(1,alpha,rt_ndt,beta,deltat[i])
    data = rbind(data,resp)
  }
  
  
  return(data %>% mutate(expectation = expectation,x = x,participant = participant, trial = 1:length(expectation)))
}




generate_psycho = function(x,threshold,slope,lapse,participant){
  
  expectation = logistic_psychometric(x,threshold, slope, lapse)
  
  #binary resps
  resp <- rbinom(length(expectation),size = 1, prob = expectation)
  
  
  return(data.frame(resp,expectation,x,participant, trial = 1:length(expectation)))
}


#ekstra math functions

# A generator function should return a named list containing elements "variables" and "generated"

generator_single_rt <- function(N = 60, S = 20){  # N is the number of data points we are generating
  
  mu_threshold = rnorm(1,0,5)
  mu_slope =     rnorm(1,-1.5,0.5)
  mu_lapse =     rnorm(1,-3,1)
  mu_rt_int =    rnorm(1,-1,0.25)
  mu_rt_beta =   rnorm(1,1.5,0.5)
  mu_rt_sd =     rnorm(1,-1,0.25)
  
  
  tau_threshold = abs(rnorm(1,5,2))
  tau_slope =     abs(rnorm(1,0.6,0.2))
  tau_lapse =     abs(rnorm(1,0.6,0.2))
  tau_rt_int =    abs(rnorm(1,0.6,0.2))
  tau_rt_beta =   abs(rnorm(1,0.6,0.2))
  tau_rt_sd =     abs(rnorm(1,0.5,0.1))
  
  threshold = rnorm(S,mu_threshold,tau_threshold)  
  slope = exp(rnorm(S,mu_slope,tau_slope))
  lapse = brms::inv_logit_scaled(rnorm(S,mu_lapse,tau_lapse)) / 2 
  rt_int = rnorm(S,mu_rt_int,tau_rt_int)  
  rt_beta = rnorm(S,mu_rt_beta,tau_rt_beta)  
  rt_sd = exp(rnorm(S,mu_rt_sd,tau_rt_sd)) 
  
  
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
  
  if(max(trial_df$rts) > 8){
    "Error"
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

generator_single_ddm <- function(N = 60, S = 20){  # N is the number of data points we are generating
  
  mu_threshold = rnorm(1,0,5)
  mu_slope =     rnorm(1,-1.5,0.5)
  mu_lapse =     rnorm(1,-3,1)
  mu_alpha =    rnorm(1,0.75,0.25)
  mu_beta =   rnorm(1,0,0.1)
  mu_delta =     rnorm(1,2,0.25)
  
  
  tau_threshold = abs(rnorm(1,5,2))
  tau_slope =     abs(rnorm(1,0.6,0.2))
  tau_lapse =     abs(rnorm(1,0.6,0.2))
  tau_alpha =    abs(rnorm(1,0.4,0.1))
  tau_beta =   abs(rnorm(1,0.4,0.1))
  tau_delta =     abs(rnorm(1,0.5,0.1))
  
  threshold = rnorm(S,mu_threshold,tau_threshold)  
  slope = exp(rnorm(S,mu_slope,tau_slope))
  lapse = brms::inv_logit_scaled(rnorm(S,mu_lapse,tau_lapse)) / 2 
  alpha = exp(rnorm(S,mu_alpha,tau_alpha))
  beta = brms::inv_logit_scaled(rnorm(S,mu_beta,tau_beta))
  delta = rnorm(S,mu_delta,tau_delta)
  
  rt_ndt = rnorm(S,0.3,0.05)  
  
  x = seq(-40,40, length.out = N)
  
  
  trial_df = data.frame(threshold = threshold, slope = slope, lapse = lapse,
                        alpha = alpha, beta = beta, delta = delta,
                        rt_ndt = rt_ndt, participant = 1:S) %>% rowwise() %>% 
    summarize((generate_trial_ddm(x, threshold = threshold,slope = slope,lapse = lapse,participant = participant,
                                  alpha = alpha,beta = beta,delta = delta,rt_ndt = rt_ndt))) %>% rename(rts = q) %>% 
    mutate(resp = ifelse(resp == "upper",1,0))
  
  # trial_df %>% ggplot(aes(x = x, y = expectation, col = participant))+facet_wrap(~participant)+geom_point()
  # trial_df %>% ggplot(aes(x = x, y = rts, col = participant))+facet_wrap(~participant)+geom_point()
  # 
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
      mu_threshold = mu_threshold,
      mu_slope = mu_slope,
      mu_lapse = mu_lapse,
      mu_alpha = mu_alpha,
      mu_beta = mu_beta,
      mu_delta = mu_delta,
      tau_threshold = tau_threshold,
      tau_slope = tau_slope,
      tau_lapse = tau_lapse,
      tau_alpha = tau_alpha,
      tau_beta = tau_beta,
      tau_delta = tau_delta,
      threshold = threshold,
      slope = slope,
      lapse = lapse,
      alpha = alpha,
      beta = beta,
      delta = delta,
      rt_ndt = rt_ndt
    ),
    generated = list(
      N = nrow(trial_df),
      expect = trial_df$expectation,
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
# 
# priors_rt <- do.call(rbind, lapply(1:30, function(i) data.frame(generator_single_rt(60, 10)[[2]])))
# 
# priors_rt = priors_rt %>% mutate(sim_id = rep(1:30,each = 600))
# 
# priors_ddm <- do.call(rbind, lapply(1:30, function(i) data.frame(generator_single_ddm(60, 10)[[2]])))
# 
# priors_ddm = priors_ddm %>% mutate(sim_id = rep(1:30,each = 600))
# 
# priors_ddm %>% group_by(x) %>% summarize(RT = mean(RT), expect = mean(expect)) %>% ggplot(aes(x = x, y = RT))+geom_point()
# 
# priors_ddm %>% filter(RT < 20) %>%  group_by(sim_id,x) %>% summarize(RT = mean(RT), expect = mean(expect)) %>% ggplot(aes(x = x, y = RT,group = sim_id))+geom_line()
# 
# priors_rt %>% group_by(x) %>% summarize(RT = mean(RT), expect = mean(expect)) %>% ggplot(aes(x = x, y = RT))+geom_point()+
#   scale_y_continuous(limits = c(1,2))
# 
# priors_rt %>% filter(RT < 20) %>% group_by(sim_id,x) %>% summarize(RT = mean(RT), expect = mean(expect)) %>% ggplot(aes(x = x, y = RT,group = sim_id))+geom_line()


fitter_cop_rt = function(samples){
  
  # fit and simualate rts model:
  
  
  sim_data = generator_single_rt(60,10)
  
  sim_data_x = sim_data[[2]]
  
  rt_model <- cmdstanr::cmdstan_model(here::here("Simulations","psychometric_scripts","Simulation_based_calibration","copulas","hier","subj_ndt","model_recovery_estimation","subj_ndt.stan"))
  nort_model <- cmdstanr::cmdstan_model(here::here("Simulations","psychometric_scripts","Simulation_based_calibration","copulas","hier","subj_ndt","model_recovery_estimation","no_rt.stan"))
  ddm_model <- cmdstanr::cmdstan_model(here::here("Simulations","psychometric_scripts","Simulation_based_calibration","copulas","hier","subj_ndt","model_recovery_estimation","hier_ddm.stan"))
  
  fit_rt = rt_model$sample(data = sim_data_x,
                           iter_warmup = samples,
                           iter_sampling = samples,
                           chains = 4,
                           refresh = 100,
                           max_treedepth = 10,
                           parallel_chains = 4,
                           adapt_delta = 0.95)
  
  rt_loo_bin = fit_rt$loo("log_lik_bin")
  rt_loo_full = fit_rt$loo("log_lik")
  
  
  divs_rt = fit_rt$diagnostic_summary()
  
  
  group_means_names = names(sim_data[[1]])[grepl("mu_",names(sim_data[[1]]))]
  
  groupmeans_rt = data.frame(fit_rt$summary(group_means_names)) %>% mutate(max_div = max(divs_rt$num_divergent),
                                                                           max_tree = max(divs_rt$num_max_treedepth),
                                                                           estimation_time = fit_rt$time()$total) %>% 
    mutate(simulated = as.numeric(data.frame(sim_data[[1]][grepl("mu_",names(sim_data[[1]]))]))) %>% 
    mutate(trials = 60, subjects = 10)
  
  
  group_var_names = names(sim_data[[1]])[grepl("tau_",names(sim_data[[1]]))]
  
  groupvar_rt = data.frame(fit_rt$summary(group_var_names)) %>% mutate(max_div = max(divs_rt$num_divergent),
                                                                       max_tree = max(divs_rt$num_max_treedepth),
                                                                       estimation_time = fit_rt$time()$total) %>% 
    mutate(simulated = as.numeric(data.frame(sim_data[[1]][grepl("tau_",names(sim_data[[1]]))])))%>% 
    mutate(trials = 60, subjects = 10)
  
  
  
  first = names(sim_data[[1]])[!grepl("mu_",names(sim_data[[1]]))]
  subj_names = first[!grepl("tau_",(first))]
  
  subj_parameters = data.frame(fit_rt$summary(subj_names)) %>% mutate(max_div = max(divs_rt$num_divergent),
                                                                      max_tree = max(divs_rt$num_max_treedepth),
                                                                      estimation_time = fit_rt$time()$total) %>% 
    mutate(
      subj_id = str_extract(variable, "\\[\\d+\\]") %>% str_remove_all("[\\[\\]]"),
      variable = str_remove(variable, "\\[\\d+\\]"),
    )%>% 
    mutate(trials = 60, subjects = 10)
  
  
  simulated_subj_param = data.frame(sim_data[[1]][subj_names]) %>% mutate(subj_id = as.character(1:10)) %>% pivot_longer(-subj_id, values_to = "simulated",names_to = "variable")
  
  subj_parameters_rt = inner_join(subj_parameters,simulated_subj_param)
  
  rt_list = list(groupmeans_rt,groupvar_rt,subj_parameters_rt)
  
  
  ########### fiitting no response time
  
  fit_nort = nort_model$sample(data = sim_data_x,
                               iter_warmup = samples,
                               iter_sampling = samples,
                               chains = 4,
                               refresh = 100,
                               max_treedepth = 10,
                               parallel_chains = 4,
                               adapt_delta = 0.95)
  
  
  nort_loo = fit_nort$loo()
  
  
  divs_nort = fit_nort$diagnostic_summary()
  
  
  
  group_means_names = c("mu_threshold","mu_slope","mu_lapse")
  
  groupmeans_nort = data.frame(fit_nort$summary(group_means_names)) %>% mutate(max_div = max(divs_nort$num_divergent),
                                                                               max_tree = max(divs_nort$num_max_treedepth),
                                                                               estimation_time = fit_nort$time()$total) %>% 
    mutate(simulated = as.numeric(data.frame(sim_data[[1]][group_means_names]))) %>% 
    mutate(simulated_model = "RT")%>% 
    mutate(trials = 60, subjects = 10)
  
  
  
  group_var_names = c("tau_threshold","tau_slope","tau_lapse")
  
  groupvar_nort = data.frame(fit_nort$summary(group_var_names)) %>% mutate(max_div = max(divs_nort$num_divergent),
                                                                           max_tree = max(divs_nort$num_max_treedepth),
                                                                           estimation_time = fit_nort$time()$total) %>% 
    mutate(simulated = as.numeric(data.frame(sim_data[[1]][group_var_names]))) %>% 
    mutate(simulated_model = "RT")%>% 
    mutate(trials = 60, subjects = 10)
  
  
  
  subj_names = c("threshold","slope","lapse")
  
  subj_parameters = data.frame(fit_nort$summary(subj_names)) %>% mutate(max_div = max(divs_nort$num_divergent),
                                                                        max_tree = max(divs_nort$num_max_treedepth),
                                                                        estimation_time = fit_nort$time()$total) %>% 
    mutate(
      subj_id = str_extract(variable, "\\[\\d+\\]") %>% str_remove_all("[\\[\\]]"),
      variable = str_remove(variable, "\\[\\d+\\]"),
    ) %>% mutate(simulated_model = "RT")%>% 
    mutate(trials = 60, subjects = 10)
  
  
  
  simulated_subj_param = data.frame(sim_data[[1]][subj_names]) %>% mutate(subj_id = as.character(1:10)) %>% pivot_longer(-subj_id, values_to = "simulated",names_to = "variable")
  
  subj_parameters_nort = inner_join(subj_parameters,simulated_subj_param)
  
  
  nort_list = list(groupmeans_nort,groupvar_nort,subj_parameters_nort)
  
  
  ####### fit ddm
  
  
  fit_ddm = ddm_model$sample(data = sim_data_x,
                             iter_warmup = samples,
                             iter_sampling = samples,
                             chains = 4,
                             refresh = 100,
                             max_treedepth = 10,
                             parallel_chains = 4,
                             adapt_delta = 0.95)
  
  
  ddm_loo = fit_ddm$loo()
  
  divs_ddm = fit_ddm$diagnostic_summary()
  
  
  full_loo = data.frame(loo::loo_compare(list(rt = rt_loo_full, ddm = ddm_loo))) %>% 
    rownames_to_column("models") %>% arrange(models) %>% 
    mutate(div = c(max(divs_ddm$num_divergent), max(divs_rt$num_divergent)),
           tree = c(max(divs_ddm$num_max_treedepth), max(divs_rt$num_max_treedepth)),
           pareto_over_07 = c(sum(ddm_loo$diagnostics$pareto_k > 0.7), sum(rt_loo_full$diagnostics$pareto_k > 0.7)))%>% 
    mutate(trials = 60, subjects = 10,simulated_model = "RT")
  
  
  
  bin_loo = data.frame(loo::loo_compare(list(rt = rt_loo_bin, nort = nort_loo))) %>% 
    rownames_to_column("models") %>% arrange(models) %>% 
    mutate(div = c(max(divs_nort$num_divergent), max(divs_rt$num_divergent)),
           tree = c(max(divs_nort$num_max_treedepth), max(divs_rt$num_max_treedepth)),
           pareto_over_07 = c(sum(nort_loo$diagnostics$pareto_k > 0.7), sum(rt_loo_bin$diagnostics$pareto_k > 0.7)))%>% 
    mutate(trials = 60, subjects = 10,simulated_model = "RT")
  
  
  
  return(list(full_loo, bin_loo,rt_list,nort_list))
  
}


fitter_ddm = function(samples){
  
  # fit and simualate rts model:
  
  sim_data = generator_single_ddm(60,10)
  
  sim_data_x = sim_data[[2]]
  
  rt_model <- cmdstanr::cmdstan_model(here::here("Simulations","psychometric_scripts","Simulation_based_calibration","copulas","hier","subj_ndt","model_recovery_estimation","subj_ndt.stan"))
  nort_model <- cmdstanr::cmdstan_model(here::here("Simulations","psychometric_scripts","Simulation_based_calibration","copulas","hier","subj_ndt","model_recovery_estimation","no_rt.stan"))
  ddm_model <- cmdstanr::cmdstan_model(here::here("Simulations","psychometric_scripts","Simulation_based_calibration","copulas","hier","subj_ndt","model_recovery_estimation","hier_ddm.stan"))
  
  fit_rt = rt_model$sample(data = sim_data_x,
                           iter_warmup = samples,
                           iter_sampling = samples,
                           chains = 4,
                           refresh = 100,
                           max_treedepth = 10,
                           parallel_chains = 4,
                           adapt_delta = 0.95)
  
  rt_loo_bin = fit_rt$loo("log_lik_bin")
  rt_loo_full = fit_rt$loo("log_lik")
  
  divs_rt = fit_rt$diagnostic_summary()
  
  
  
  ########### fiitting no response time
  
  fit_nort = nort_model$sample(data = sim_data_x,
                               iter_warmup = samples,
                               iter_sampling = samples,
                               chains = 4,
                               refresh = 100,
                               max_treedepth = 10,
                               parallel_chains = 4,
                               adapt_delta = 0.95)
  
  
  nort_loo = fit_nort$loo()
  
  
  divs_nort = fit_nort$diagnostic_summary()
  
  
  
  
  ####### fit ddm
  
  fit_ddm = ddm_model$sample(data = sim_data_x,
                             iter_warmup = samples,
                             iter_sampling = samples,
                             chains = 4,
                             refresh = 100,
                             max_treedepth = 10,
                             parallel_chains = 4,
                             adapt_delta = 0.95)
  
  
  ddm_loo = fit_ddm$loo()
  divs_ddm = fit_ddm$diagnostic_summary()
  
  
  group_means_names = names(sim_data[[1]])[grepl("mu_",names(sim_data[[1]]))]
  
  groupmeans_ddm = data.frame(fit_ddm$summary(group_means_names)) %>% mutate(max_div = max(divs_rt$num_divergent),
                                                                             max_tree = max(divs_rt$num_max_treedepth),
                                                                             estimation_time = fit_ddm$time()$total) %>% 
    mutate(simulated = as.numeric(data.frame(sim_data[[1]][grepl("mu_",names(sim_data[[1]]))])))%>% 
    mutate(trials = 60, subjects = 10,simulated_model = "DDM")
  
  group_var_names = names(sim_data[[1]])[grepl("tau_",names(sim_data[[1]]))]
  
  groupvar_ddm = data.frame(fit_ddm$summary(group_var_names)) %>% mutate(max_div = max(divs_rt$num_divergent),
                                                                         max_tree = max(divs_rt$num_max_treedepth),
                                                                         estimation_time = fit_ddm$time()$total) %>% 
    mutate(simulated = as.numeric(data.frame(sim_data[[1]][grepl("tau_",names(sim_data[[1]]))])))%>% 
    mutate(trials = 60, subjects = 10,simulated_model = "DDM")
  
  
  first = names(sim_data[[1]])[!grepl("mu_",names(sim_data[[1]]))]
  subj_names = first[!grepl("tau_",(first))]
  
  subj_parameters = data.frame(fit_ddm$summary(subj_names)) %>% mutate(max_div = max(divs_rt$num_divergent),
                                                                       max_tree = max(divs_rt$num_max_treedepth),
                                                                       estimation_time = fit_ddm$time()$total) %>% 
    mutate(
      subj_id = str_extract(variable, "\\[\\d+\\]") %>% str_remove_all("[\\[\\]]"),
      variable = str_remove(variable, "\\[\\d+\\]"),
    )%>% 
    mutate(trials = 60, subjects = 10,simulated_model = "DDM")
  
  
  simulated_subj_param = data.frame(sim_data[[1]][subj_names]) %>% mutate(subj_id = as.character(1:10)) %>% pivot_longer(-subj_id, values_to = "simulated",names_to = "variable")
  
  subj_parameters_ddm = inner_join(subj_parameters,simulated_subj_param)
  
  ddm_list = list(groupmeans_ddm,groupvar_ddm,subj_parameters_ddm)
  
  
  
  
  
  full_loo = data.frame(loo::loo_compare(list(rt = rt_loo_full, ddm = ddm_loo))) %>% 
    rownames_to_column("models") %>% arrange(models) %>% 
    mutate(div = c(max(divs_ddm$num_divergent), max(divs_rt$num_divergent)),
           tree = c(max(divs_ddm$num_max_treedepth), max(divs_rt$num_max_treedepth)),
           pareto_over_07 = c(sum(ddm_loo$diagnostics$pareto_k > 0.7), sum(rt_loo_full$diagnostics$pareto_k > 0.7)))%>% 
    mutate(trials = 60, subjects = 10,simulated_model = "DDM")
  
  
  
  bin_loo = data.frame(loo::loo_compare(list(rt = rt_loo_bin, nort = nort_loo))) %>% 
    rownames_to_column("models") %>% arrange(models) %>% 
    mutate(div = c(max(divs_nort$num_divergent), max(divs_rt$num_divergent)),
           tree = c(max(divs_nort$num_max_treedepth), max(divs_rt$num_max_treedepth)),
           pareto_over_07 = c(sum(nort_loo$diagnostics$pareto_k > 0.7), sum(rt_loo_bin$diagnostics$pareto_k > 0.7)))%>% 
    mutate(trials = 60, subjects = 10,simulated_model = "DDM")
  
  
  return(list(full_loo, bin_loo,ddm_list))
  
}

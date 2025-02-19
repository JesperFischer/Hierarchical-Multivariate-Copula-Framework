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


generate_learn = function(u,learningrate,e0,participant){
  
  
  expectation = array(NA,length(u)+1)
  
  expectation[1] = e0
  for(i in 1: length(u)){
    expectation[i+1] = learning_update(u[i], expectation[i], learningrate)
  }
  expectation = expectation[1:length(u)]
  
  
  

  #binary resps
  resp <- rbinom(length(expectation), size = 1, prob = expectation)
  
  
  return(data.frame(resp,expectation,u,participant, trial = 1:length(expectation), expect = expectation))
}


#ekstra math functions

# A generator function should return a named list containing elements "variables" and "generated"

generate_trial_ddm = function(u,learningrate,e0,participant,
                              alpha,beta, delta, rt_ndt){
  
  expectation = array(NA,length(u)+1)
  
  expectation[1] = e0
  for(i in 1: length(u)){
    expectation[i+1] = learning_update(u[i], expectation[i], learningrate)
  }
  expectation = expectation[1:length(u)]
  deltat = (expectation - (1-expectation)) * delta
  
  data = data.frame()
  for(i in 1:length(u)){
    resp = RWiener::rwiener(1,alpha,rt_ndt,beta,deltat[i])
    data = rbind(data,resp)
  }
  
  
  return(data %>% mutate(expect = expectation,u = u,participant = participant, trial = 1:length(expectation)))
}



# A generator function should return a named list containing elements "variables" and "generated"

generator_single_rt <- function(N = 60, S = 20){  # N is the number of data points we are generating
  
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


generator_single_ddm <- function(N = 60, S = 20){  # N is the number of data points we are generating
  
  mu_learningrate = rnorm(1,-1,1)
  mu_e0 =     rnorm(1,0,0.2)
  mu_alpha =    rnorm(1,0.75,0.25)
  mu_beta =   rnorm(1,0,0.1)
  mu_delta =     rnorm(1,2,0.25)
  
  
  tau_learningrate = abs(rnorm(1,0.5,0.1))
  tau_e0 =     abs(rnorm(1,0.5,0.1))
  tau_alpha =    abs(rnorm(1,0.4,0.1))
  tau_beta =   abs(rnorm(1,0.4,0.1))
  tau_delta =     abs(rnorm(1,0.5,0.1))
  
  learningrate = brms::inv_logit_scaled(rnorm(S,mu_learningrate,tau_learningrate))
  e0 = brms::inv_logit_scaled(rnorm(S,mu_e0,tau_e0))
  alpha = exp(rnorm(S,mu_alpha,tau_alpha))
  beta = brms::inv_logit_scaled(rnorm(S,mu_beta,tau_beta))
  delta = rnorm(S,mu_delta,tau_delta)
  
  rt_ndt = rnorm(S,0.3,0.05)  
  
  
  u = c(1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1,
        0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0,
        0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1)
  
  
  trial_df = data.frame(learningrate = learningrate, e0 = e0,
                        alpha = alpha, beta = beta, delta = delta,
                        rt_ndt = rt_ndt, participant = 1:S) %>% rowwise() %>% 
    summarize((generate_trial_ddm(u = u, learningrate = learningrate,e0 = e0,participant = participant,
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
      mu_learningrate = mu_learningrate,
      mu_e0 = mu_e0,
      mu_alpha = mu_alpha,
      mu_beta = mu_beta,
      mu_delta = mu_delta,
      tau_learningrate = tau_learningrate,
      tau_e0 = tau_e0,
      tau_alpha = tau_alpha,
      tau_beta = tau_beta,
      tau_delta = tau_delta,
      learningrate = learningrate,
      e0 = e0,
      alpha = alpha,
      beta = beta,
      delta = delta,
      rt_ndt = rt_ndt
    ),
    generated = list(
      N = nrow(trial_df),
      expect = trial_df$expect,
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


fitter_cop_rt = function(samples){
  
  # fit and simualate rts model:
  
  sim_data = generator_single_rt(60,10)
  
  sim_data_x = sim_data[[2]]
  
  
  
  rt_model <- cmdstanr::cmdstan_model(here::here("Simulations","learning_scripts","Copulas","Simulation_based_calibration","hier","ddm_compar","rw_rt.stan"))
  nort_model <- cmdstanr::cmdstan_model(here::here("Simulations","learning_scripts","Copulas","Simulation_based_calibration","hier","ddm_compar","rw_nort.stan"))
  ddm_model <- cmdstanr::cmdstan_model(here::here("Simulations","learning_scripts","Copulas","Simulation_based_calibration","hier","ddm_compar","rw_ddm.stan"))
  
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
  
  
  
  group_means_names = c("mu_learningrate","mu_e0")
  
  groupmeans_nort = data.frame(fit_nort$summary(group_means_names)) %>% mutate(max_div = max(divs_nort$num_divergent),
                                                                               max_tree = max(divs_nort$num_max_treedepth),
                                                                               estimation_time = fit_nort$time()$total) %>% 
    mutate(simulated = as.numeric(data.frame(sim_data[[1]][group_means_names]))) %>% 
    mutate(simulated_model = "RT")%>% 
    mutate(trials = 60, subjects = 10)
  
  
  
  group_var_names = c("tau_learningrate","tau_e0")
  
  groupvar_nort = data.frame(fit_nort$summary(group_var_names)) %>% mutate(max_div = max(divs_nort$num_divergent),
                                                                           max_tree = max(divs_nort$num_max_treedepth),
                                                                           estimation_time = fit_nort$time()$total) %>% 
    mutate(simulated = as.numeric(data.frame(sim_data[[1]][group_var_names]))) %>% 
    mutate(simulated_model = "RT")%>% 
    mutate(trials = 60, subjects = 10)
  
  
  
  subj_names = c("learningrate","e0")
  
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
  
  data.frame(sim_data_x) %>% ggplot(aes(x = expectation))
  
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
  
  rt_model <- cmdstanr::cmdstan_model(here::here("Simulations","learning_scripts","Copulas","Simulation_based_calibration","hier","ddm_compar","rw_rt.stan"))
  nort_model <- cmdstanr::cmdstan_model(here::here("Simulations","learning_scripts","Copulas","Simulation_based_calibration","hier","ddm_compar","rw_nort.stan"))
  ddm_model <- cmdstanr::cmdstan_model(here::here("Simulations","learning_scripts","Copulas","Simulation_based_calibration","hier","ddm_compar","rw_ddm.stan"))
  
  
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
packages = c("brms","tidyverse","bayesplot","pracma","here",
             "patchwork","posterior","HDInterval","loo", "furrr", "SBC","future")

do.call(pacman::p_load, as.list(packages))


df <- read_csv(here::here("Real data","Learning","data","Data_as_csv.csv"))%>% 
  group_by(Subject) %>% 
  mutate(trial = 1:n(),
         rts  = exp(rts) / 1000) %>% drop_na()

df %>% drop_na() %>% 
  group_by(Subject) %>% mutate(trial = 1:n()) %>% group_by(trial) %>% 
  summarize(u = mean(u), resp = mean(bin_resp), rt = mean(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>% 
  ggplot(aes(x = trial))+geom_pointrange(aes(y = rt, ymin = rt-se_rts, ymax = rt+se_rts))



df %>% drop_na() %>% 
  group_by(Subject) %>% mutate(trial = 1:n()) %>% group_by(trial) %>% 
  summarize(u = mean(u), resp = mean(bin_resp), rt = mean(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>% 
  ggplot(aes(x = trial))+geom_pointrange(aes(y = resp, ymin = resp-se_resp, ymax = resp+se_resp))+geom_line(aes(y = u),col = "red")



df %>% drop_na() %>% 
  group_by(Subject) %>% mutate(trial = 1:n()) %>% group_by(trial) %>% 
  summarize(u = mean(u), resp = mean(bin_resp), rt = mean(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>% 
  ggplot(aes(x = resp))+geom_pointrange(aes(y = rt, ymin = rt-se_rts, ymax = rt+se_rts))


df_filtered <- df %>%
  group_by(Subject) %>%
  filter(rts > min(rts))

df = df_filtered

N = df %>% group_by(Subject) %>% summarize(N = n()) %>% .$N

S = df %>% group_by(Subject) %>% summarize(N = n()) %>% ungroup() %>% summarize(N = n()) %>% .$N



ranges <- data.frame(
  start = c(1, head(cumsum(N) + 1, -1)),
  end = cumsum(N)
)



minRT_repeated <- df %>%
  group_by(Subject) %>%
  summarize(minRT = min(rts), .groups = "drop") %>%
  pull(minRT) %>%
  rep(times = N)  # Repeat for each subject's trials



dflist = list(
  N = nrow(df),
  x = df$u,
  binom_y = df$bin_resp,
  RT = (df$rts),
  trial = df$trial,
  minRT = df %>% group_by(Subject) %>% summarize(minRT = min(rts)) %>% .$minRT,
  minRT_t = minRT_repeated,
  S_id = df$Subject,
  id_ind = ifelse(!duplicated(df$Subject), 1, 0),
  starts = ranges$start,
  ends = ranges$end,
  S = length(unique(df$Subject))
)

# rt_model_test <- cmdstanr::cmdstan_model(here::here("Simulations","learning_scripts","Copulas","Simulation_based_calibration","hier","ddm_compar","rw_rt_test.stan"), force_recompile = T)
# rw_rt_test  = rt_model_test$sample(data = dflist,
#                                    iter_warmup = 250,
#                                    iter_sampling = 250,
#                                    chains = 4,
#                                    refresh = 50,
#                                    max_treedepth = 10,
#                                    parallel_chains = 4,
#                                    adapt_delta = 0.95)
# 
# rw_rt_test$save_object(here::here("real data","reinforcement learning","rw_rt_test.rds"))


rt_model <- cmdstanr::cmdstan_model(here::here("real data","reinforcement learning","rw_rt.stan"), force_recompile = T)
nort_model <- cmdstanr::cmdstan_model(here::here("real data","reinforcement learning","rw_nort.stan"), force_recompile = T)
ddm_model <- cmdstanr::cmdstan_model(here::here("real data","reinforcement learning","rw_ddm.stan"), force_recompile = T)



rw_rt  = rt_model$sample(data = dflist,
                         iter_warmup = 1000,
                         iter_sampling = 1000,
                         chains = 4,
                         refresh = 50,
                         max_treedepth = 10,
                         parallel_chains = 4,
                         adapt_delta = 0.95)

rw_rt_nolow$save_object(here::here("real data","reinforcement learning","rw_rt_nolow.rds"))



rw_nort_nolow = nort_model$sample(data = dflist,
                             iter_warmup = 1000,
                             iter_sampling = 1000,
                             chains = 4,
                             refresh = 100,
                             max_treedepth = 10,
                             parallel_chains = 4,
                             adapt_delta = 0.95)

rw_nort_nolow$save_object(here::here("real data","reinforcement learning","rw_nort_nolow.rds"))



rw_ddm_nolow = ddm_model$sample(data = dflist,
                           iter_warmup = 1000,
                           iter_sampling = 1000,
                           chains = 4,
                           refresh = 50,
                           max_treedepth = 10,
                           parallel_chains = 4,
                           adapt_delta = 0.95)


rw_ddm_nolow$save_object(here::here("real data","reinforcement learning","rw_ddm_nolow.rds"))


rt_loo_bin = rw_rt_nolow$loo("log_lik_bin")
rt_loo_full = rw_rt_nolow$loo("log_lik")

# rt_loo_bin = rw_rt $loo("log_lik_bin", moment_match = T)
# rt_loo_full = rw_rt $loo("log_lik", moment_match = T)

divs_rt = rw_rt_nolow$diagnostic_summary()
nort_loo = rw_nort_nolow$loo()

# nort_loo = fit_nort$loo(moment_match = T)
divs_nort = rw_nort_nolow$diagnostic_summary()



ddm_loo = rw_ddm_nolow$loo()
divs_ddm = rw_ddm_nolow$diagnostic_summary()


full_loo = data.frame(loo::loo_compare(list(rt = rt_loo_full, ddm = ddm_loo))) %>% 
  rownames_to_column("models") %>% arrange(models) %>% 
  mutate(pareto_over_07 = c(sum(ddm_loo$diagnostics$pareto_k > 0.7), sum(rt_loo_full$diagnostics$pareto_k > 0.7)))



bin_loo = data.frame(loo::loo_compare(list(rt = rt_loo_bin, nort = nort_loo))) %>% 
  rownames_to_column("models") %>% arrange(models) %>% 
  mutate(pareto_over_07 = c(sum(nort_loo$diagnostics$pareto_k > 0.7), sum(rt_loo_bin$diagnostics$pareto_k > 0.7)))

full_loo

bin_loo

draw_id = sample(1:4000,500)

# subject level 
big_df_ddm = data.frame()

for(s in unique(df$Subject)){
  print(s)
  us = df %>% filter(Subject == s) %>% .$u
  
  minRT = min(df %>% filter(Subject == s) %>% .$rts)
  
  parameters = paste0(c("learningrate","e0","alpha","beta","delta","rt_ndt"), "[",s,"]")
  
  dfq = as_draws_df(rw_ddm_nolow$draws(parameters))%>% select(-contains(".")) %>% mutate(draw = 1:n()) %>% 
    filter(draw %in% draw_id) %>% 
    rename_with(~c("learningrate","e0","alpha","beta","delta","rt_ndt","draw")) %>% 
    group_by(draw) %>% 
    summarize(list(generate_trial_ddm(u = us, learningrate = learningrate,e0 = e0,participant = s,
                                      alpha = alpha,beta = beta,delta = delta,rt_ndt = rt_ndt))) %>% unnest() %>% 
    rename(rts = q) %>% 
    mutate(resp = ifelse(resp == "upper",1,0))
  
  big_df_ddm = rbind(big_df_ddm,dfq)
  
}

qq_summar_ddm = big_df_ddm %>%  ungroup() %>%  mutate(rts = (rts),bin_resp = resp) %>% drop_na() %>%
  group_by(trial,draw) %>% 
  summarize(resp = mean(bin_resp), rt = median(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n()))


qq_summar_sum_ddm = qq_summar_ddm %>%   group_by(trial) %>% 
  summarize(resps = mean(resp), rts = median(rt), se_resp =  (mean(resp) * (1- mean(resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n()))


df %>%   ungroup() %>%  mutate(rts = (rts),bin_resp = bin_resp) %>% drop_na() %>%
  group_by(trial) %>% 
  summarize(resp = mean(bin_resp), rt = median(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>%
  ggplot(aes(x = trial))+
  geom_line(data = qq_summar_ddm, aes(x = trial, y = rt, group = draw), alpha = 0.01, col = "orange")+
  geom_pointrange(aes(y = rt, ymin = rt-2*se_rts, ymax = rt+2*se_rts))+
  geom_line(data = qq_summar_sum_ddm, aes(x = trial, y = rts), col = "red")+
  scale_y_continuous(limits = c(0,1))+
  theme_classic()


df %>%   ungroup() %>%  mutate(rts = (rts),bin_resp = bin_resp) %>% drop_na() %>%
  group_by(trial) %>% 
  summarize(resp = mean(bin_resp), rt = median(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>%
  ggplot(aes(x = trial))+
  geom_line(data = qq_summar_ddm, aes(x = trial, y = resp, group = draw), alpha = 0.01, col = "orange")+
  geom_line(data = qq_summar_sum_ddm, aes(x = trial, y = resps), col = "red")+
  geom_pointrange(aes(y = resp, ymin = resp-2*se_resp, ymax = resp+2*se_resp))+
  theme_classic()


# subject level 
big_df_rt = data.frame()

for(s in unique(df$Subject)){
  print(s)
  us = df %>% filter(Subject == s) %>% .$u
  
  minRT = min(df %>% filter(Subject == s) %>% .$rts)
  
  parameters = paste0(c("learningrate","e0","rt_int","rt_beta","rt_ndt","rt_sd","rho"), "[",s,"]")
  
  dfq = as_draws_df(rw_rt_nolow$draws(parameters))%>% select(-contains(".")) %>% mutate(draw = 1:n()) %>% 
    filter(draw %in% draw_id) %>% 
    rename_with(~c("learningrate","e0","rt_int","rt_beta","rt_ndt","rt_sd","rho","draw")) %>% 
    group_by(draw) %>% 
    summarize(list(generate_trialwise_data_rt_copula(u = us, learningrate = learningrate,e0 = e0,participant = s,
                                                     rt_int = rt_int,rt_beta = rt_beta,rt_sd = rt_sd,rt_ndt = rt_ndt,rho = rho))) %>% unnest()
  
  big_df_rt = rbind(big_df_rt,dfq)
}

# big_df_rt %>% ggplot(aes(x = trial, y = expect))+geom_line(aes(group = draw))+facet_wrap(~participant)+geom_point(aes(y = u))


qq_summar_rt = big_df_rt %>%  ungroup() %>%  mutate(rts = (rts),bin_resp = resp) %>% drop_na() %>%
  group_by(trial,draw) %>% 
  summarize(resp = mean(bin_resp), rt = median(rts))

qq_summar_sum_rt = qq_summar_rt %>% 
  group_by(trial) %>% 
  summarize(resps = mean(resp), rts = median(rt))


df %>%   ungroup() %>%  mutate(rts = (rts),bin_resp = bin_resp) %>% drop_na() %>%
  group_by(trial) %>% 
  summarize(resp = mean(bin_resp), rt = median(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>%
  ggplot(aes(x = trial))+
  geom_line(data = qq_summar_rt, aes(x = trial, y = rt, group = draw), alpha = 0.01, col = "orange")+
  geom_pointrange(aes(y = rt, ymin = rt-2*se_rts, ymax = rt+2*se_rts))+
  geom_line(data = qq_summar_sum_rt, aes(x = trial, y = rts), col = "red")+
  scale_y_continuous(limits = c(0,1))+
  theme_classic()


df %>%   ungroup() %>%  mutate(rts = (rts),bin_resp = bin_resp) %>% drop_na() %>%
  group_by(trial) %>% 
  summarize(resp = mean(bin_resp), rt = median(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>%
  ggplot(aes(x = trial))+
  geom_line(data = qq_summar_rt, aes(x = trial, y = resp, group = draw), alpha = 0.01, col = "orange")+
  geom_pointrange(aes(y = resp, ymin = resp-2*se_resp, ymax = resp+2*se_resp))+
  geom_line(data = qq_summar_sum_rt, aes(x = trial, y = resps), col = "red")+
  theme_classic()


df %>% drop_na() %>% 
  ggplot(aes(x = trial))+
  geom_point(aes(y = bin_resp))+
  geom_point(aes(y = u-0.05),col = "red")+
  geom_line(data = big_df_rt %>% rename(Subject = participant), aes(x = trial, y = expect, group = draw), alpha = 0.05)+
  facet_wrap(~Subject)



# subject level  no rts
big_df_nort = data.frame()

for(s in unique(df$Subject)){
  print(s)
  us = df %>% filter(Subject == s) %>% .$u
  
  parameters = paste0(c("learningrate","e0"), "[",s,"]")
  
  dfq = as_draws_df(rw_nort_nolow$draws(parameters))%>% select(-contains(".")) %>% mutate(draw = 1:n()) %>% 
    filter(draw %in% draw_id) %>% 
    rename_with(~c("learningrate","e0","draw")) %>% 
    group_by(draw) %>% 
    summarize(list(generate_learn(u = us, learningrate = learningrate,e0 = e0,participant = s))) %>% unnest()
  
  big_df_nort = rbind(big_df_nort,dfq)
}


qq_summar_nort = big_df_nort %>% mutate(bin_resp = resp) %>% 
  ungroup() %>% drop_na() %>% 
  group_by(trial,draw) %>% 
  summarize(resp = mean(bin_resp))


qq_summar_sum_nort = qq_summar_nort %>%   
  group_by(trial) %>% 
  summarize(resps = mean(resp))



df %>%   ungroup() %>% drop_na() %>%
  group_by(trial) %>% 
  summarize(resp = mean(bin_resp), rt = median(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>%
  ggplot(aes(x = trial))+
  geom_line(data = qq_summar_nort, aes(x = trial, y = resp, group = draw), alpha = 0.01, col = "#00C853")+
  geom_line(data = qq_summar_sum_nort, aes(x = trial, y = resps), col = "green")+
  geom_pointrange(aes(y = resp, ymin = resp-se_resp, ymax = resp+se_resp))+
  theme_classic()




## combine

alpha = 0.1

rw_rw = df %>%   ungroup() %>% drop_na()  %>%
  group_by(trial)  %>% 
  summarize(resp = mean(bin_resp), rt = median(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>%
  ggplot(aes(x = trial))+
  geom_line(data = qq_summar_rt %>% filter(draw %in% draw_id), aes(x = trial, y = resp, group = draw), alpha = alpha, col = "#6CEEF8", show.legend = FALSE)+
  # geom_line(data = qq_summar_nort, aes(x = Difficulty_bin, y = resp, group = draw), alpha = 0.01, col = "black")+
  # geom_line(data = qq_summar_sum_nort, aes(x = trial, y = resps, col = "No RT"))+
  geom_line(data = qq_summar_ddm%>% filter(draw %in% draw_id), aes(x = trial, y = resp, group = draw), alpha = alpha, col = "orange", show.legend = FALSE)+
  geom_pointrange(aes(y = resp, ymin = resp-2*se_resp, ymax = resp+2*se_resp), show.legend = FALSE)+
  geom_line(data = qq_summar_sum_ddm , aes(x = trial, y = resps, col = "DDM"), show.legend = FALSE, linewidth = 1.1)+
  # geom_line(data = qq_summar_sum_nort, aes(x = trial, y = resps, col = "No RT"), show.legend = FALSE, linewidth = 1.1)+
  geom_line(data = qq_summar_sum_rt, aes(x = trial, y = resps, col = "RT"), show.legend = FALSE, linewidth = 1.1)+
  theme_classic()+
  scale_color_manual(name = "Fitted model",values = c("red","blue"))+
  scale_x_continuous("trials", breaks = scales::pretty_breaks(n = 3))+
  scale_y_continuous("P(Response == 1)", breaks = scales::pretty_breaks(n = 3))+ 
  ggtitle("Probabilistic learning")+
  theme(plot.title = element_text(hjust = 0.5))

rw_rw
## rts

rts_rw = df %>%   ungroup() %>% drop_na() %>%
  group_by(trial) %>% 
  summarize(resp = mean(bin_resp), rt = median(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>%
  ggplot(aes(x = trial))+
  geom_line(data = qq_summar_rt%>% filter(draw %in% draw_id), aes(x = trial, y = rt, group = draw), alpha = alpha, col = "#6CEEF8")+
  geom_line(data = qq_summar_ddm%>% filter(draw %in% draw_id), aes(x = trial, y = rt, group = draw), alpha = alpha, col = "orange")+
  geom_pointrange(aes(y = rt, ymin = rt-2*se_rts, ymax = rt+2*se_rts))+
  geom_line(data = qq_summar_sum_rt, aes(x = trial, y = rts), col = "blue", linewidth = 1.1)+
  geom_line(data = qq_summar_sum_ddm, aes(x = trial, y = rts), col = "red", linewidth = 1.1)+
  scale_color_manual(name = "Fitted model",values = c("red","blue"))+
  theme_classic()+
  scale_x_continuous("Trial", breaks = scales::pretty_breaks(n = 3))+
  scale_y_continuous("Response time (S)", limits = c(0.2,0.7), breaks = scales::pretty_breaks(n = 3))


reinforcement = rw_rw/rts_rw+
  plot_layout(axis_titles = "collect")&
  theme(legend.position = "top")



reinforcement








# as rolling:

library(slider)

# Define rolling window size (adjust as needed)
window_size <- 10  # Example: 10-trial moving average

df_rolled <- df %>%
  drop_na() %>%
  group_by(Subject) %>%
  mutate(
    bin_resp_smooth = slide_dbl(bin_resp, mean, .before = window_size, .complete = TRUE),
    u_smooth = slide_dbl(u, mean, .before = window_size, .complete = TRUE)
  )

# Plot with rolling averages
df_rolled %>%
  ggplot(aes(x = trial)) +
  geom_line(data = big_df_rt %>% rename(Subject = participant), 
            aes(x = trial, y = expect, group = draw), alpha = 0.05) +
  geom_line(aes(y = bin_resp_smooth), col = "blue") +  # Smoothed bin_resp
  geom_line(aes(y = u_smooth), col = "red") +  # Smoothed u
  facet_wrap(~Subject)


# rts

df %>% drop_na() %>% 
  ggplot(aes(x = trial))+
  geom_point(aes(y = rts))+
  geom_line(data = big_df_rt %>% rename(Subject = participant), aes(x = trial, y = rts  , group = draw), alpha = 0.05)+
  coord_cartesian(ylim = c(0,2.5))+
  facet_wrap(~Subject)



## convergence!

# trace and prior posterior updates: (DDM)

parameters = paste0("mu_",c("learningrate","e0","alpha","beta","delta"))

traces_rw_ddm = as_draws_df(rw_ddm_nolow$draws(parameters)) %>% pivot_longer(all_of(parameters), names_to = "variables") %>% 
  mutate(.chain = as.factor(.chain)) %>% 
  ggplot(aes(x = .iteration, y = value, col = .chain)) + facet_wrap(~variables, nrow = 1, scales = "free")+geom_line()+
  theme_minimal()+
  scale_color_manual(values = c(bayesplot::color_scheme_get("red")[[1]],
                                bayesplot::color_scheme_get("red")[[2]],
                                bayesplot::color_scheme_get("red")[[3]],
                                bayesplot::color_scheme_get("red")[[4]]))+
  theme(legend.position = "top")+
  scale_x_continuous("Iterations", breaks = scales::pretty_breaks(n = 3))+
  scale_y_continuous("Value", breaks = scales::pretty_breaks(n = 3))



posteriors_rw_ddm = as_draws_df(rw_ddm_nolow$draws(parameters)) %>% pivot_longer(all_of(parameters), names_to = "variables") %>% 
  mutate(distribution = "posterior")

priors_rw_ddm = data.frame(variables = parameters) %>% mutate(value = list(rnorm(4000,-1,1),
                                                                           rnorm(4000,0,1),
                                                                           rnorm(4000,0.5,0.5),
                                                                           rnorm(4000,0,0.5),
                                                                           rnorm(4000,2,1))) %>% unnest()%>% 
  mutate(distribution = "prior")


pp_update_rw_ddm = rbind(priors_rw_ddm,posteriors_rw_ddm %>% select(variables,value,distribution)) %>% ggplot(aes(x = value, fill = distribution))+
  geom_histogram(col = "black", alpha = 0.5)+
  facet_wrap(~variables, nrow = 1, scales = "free")+
  theme_minimal()+
  theme(legend.position = "top")+
  scale_x_continuous("Value", breaks = scales::pretty_breaks(n = 3))+
  scale_y_continuous("Frequency", breaks = scales::pretty_breaks(n = 3))+
  theme(legend.position = "top")



convergence_rw_ddm = traces_rw_ddm / pp_update_rw_ddm
convergence_rw_ddm




parameters = paste0("mu_",c("learningrate","e0","rt_int","rt_beta","rt_sd"))

traces_rw_rt = as_draws_df(rw_rt_nolow$draws(parameters)) %>% pivot_longer(all_of(parameters), names_to = "variables") %>% 
  mutate(.chain = as.factor(.chain)) %>% 
  ggplot(aes(x = .iteration, y = value, col = .chain)) + facet_wrap(~variables, nrow = 1, scales = "free")+geom_line()+
  theme_minimal()+
  scale_color_manual(values = c(bayesplot::color_scheme_get("red")[[1]],
                                bayesplot::color_scheme_get("red")[[2]],
                                bayesplot::color_scheme_get("red")[[3]],
                                bayesplot::color_scheme_get("red")[[4]]))+
  theme(legend.position = "top")+
  scale_x_continuous("Iterations", breaks = scales::pretty_breaks(n = 3))+
  scale_y_continuous("Value", breaks = scales::pretty_breaks(n = 3))


posteriors_rw_rt = as_draws_df(rw_rt_nolow$draws(parameters)) %>% pivot_longer(all_of(parameters), names_to = "variables") %>% 
  mutate(distribution = "posterior")

priors_rw_rt = data.frame(variables = parameters) %>% mutate(value = list(rnorm(4000,-1,1),
                                                                           rnorm(4000,0,1),
                                                                           rnorm(4000,-1,1),
                                                                           rnorm(4000,1.5,1),
                                                                           rnorm(4000,-1,1))) %>% unnest()%>% 
  mutate(distribution = "prior")


pp_update_rw_rt = rbind(priors_rw_rt,posteriors_rw_rt %>% select(variables,value,distribution)) %>% ggplot(aes(x = value, fill = distribution))+
  geom_histogram(col = "black", alpha = 0.5)+
  facet_wrap(~variables, nrow = 1, scales = "free")+
  theme_minimal()+
  theme(legend.position = "top")+
  scale_x_continuous("Value", breaks = scales::pretty_breaks(n = 3))+
  scale_y_continuous("Frequency", breaks = scales::pretty_breaks(n = 3))+
  theme(legend.position = "top")



convergence_rw_rt = traces_rw_rt / pp_update_rw_rt


convergence_rw_rt
convergence_rw_ddm




save.image("~/Modeling-the-experiment/real data/reinforcement learning/plotting.RData")






















# # group level
# subject_level <- big_df %>%
#   group_by(participant, trial) %>%
#   summarize(
#     mean_rts = mean(rts),
#     mean_expect = mean(expect),
#     rts_hdi_low = HDInterval::hdi(rts, 0.2)[1],  # 95% credible interval
#     rts_hdi_high = HDInterval::hdi(rts, 0.2)[2],
#     .groups = "drop"
#   )
# 
# # Step 2: Compute Group-Level Aggregates
# group_estimates <- subject_level %>%
#   group_by(trial) %>%
#   summarize(
#     mean_rts = mean(mean_rts),
#     rts_hdi_low = mean(rts_hdi_low),
#     rts_hdi_high = mean(rts_hdi_high),
#     mean_expect = mean(mean_expect),
#     se_expect = sd(mean_expect) / sqrt(n()),
#     expect_hdi_low = HDInterval::hdi(mean_expect, 0.95)[1],
#     expect_hdi_high = HDInterval::hdi(mean_expect, 0.95)[2],
#     .groups = "drop"
#   )
# 
# df %>% mutate(rts = (rts)) %>% drop_na() %>% 
#   group_by(Subject) %>% mutate(trial = 1:n()) %>% group_by(trial) %>% 
#   summarize(u = mean(u), resp = mean(bin_resp), rt = median(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>% 
#   ggplot(aes(x = trial))+geom_pointrange(aes(y = rt, ymin = rt-se_rts, ymax = rt+se_rts))+
#   geom_ribbon(data = group_estimates, aes(x = trial, y = mean_rts, ymin = rts_hdi_low, ymax = rts_hdi_high), alpha = 0.2, fill = "black")





# group level shit

params = c("mu_learningrate",
"mu_e0",
"mu_alpha",
"mu_beta",
"mu_delta")

minRT = mean(df %>% group_by(Subject) %>% summarize(minRT = min(rts)) %>% .$minRT)

us = df %>% filter(Subject == 1) %>% .$u

draw_id = sample(1:2000,500)

ddm_group = as_draws_df(rw_ddm$draws(params)) %>% select(-contains(".")) %>% 
  rename_with(~c("learningrate","e0","alpha","beta","delta")) %>% mutate(draw = 1:n()) %>% 
  mutate(learningrate = brms::inv_logit_scaled(learningrate),
         e0 = brms::inv_logit_scaled(e0),
         alpha = exp(alpha),
         beta = brms::inv_logit_scaled(beta)
  )

tester = ddm_group %>% 
  filter(draw %in% draw_id) %>%
  mutate(rt_ndt = minRT+0.05) %>% group_by(draw) %>% 
  summarize(list(generate_trial_ddm(u = us, learningrate = learningrate,e0 = e0,participant = 1,
                                alpha = alpha,beta = beta,delta = delta,rt_ndt = rt_ndt))) %>% unnest() %>% 
  rename(rts = q) %>% 
  mutate(resp = ifelse(resp == "upper",1,0))

tester_sum = tester %>% group_by(trial) %>% summarize(mean = median(rts),
                                                      q50_low = HDInterval::hdi(rts,0.5)[1],
                                                      q50_high = HDInterval::hdi(rts,0.5)[2],
                                                      q80_low = HDInterval::hdi(rts,0.8)[1],
                                                      q80_high = HDInterval::hdi(rts,0.8)[2])

df %>% mutate(rts = (rts)) %>% drop_na() %>% 
  group_by(Subject) %>% mutate(trial = 1:n()) %>% group_by(trial) %>% 
  summarize(u = mean(u), resp = mean(bin_resp), rt = median(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>% 
  ggplot(aes(x = trial))+geom_pointrange(aes(y = rt, ymin = rt-se_rts, ymax = rt+se_rts))+
  geom_ribbon(data = tester_sum, aes(x = trial, y = mean, ymin = q80_low, ymax = q80_high), alpha = 0.2, fill = "black")+
  geom_ribbon(data = tester_sum, aes(x = trial, y = mean, ymin = q50_low, ymax = q50_high), alpha = 0.5, fill = "black")


tester_sum2 = tester %>% group_by(trial) %>% summarize(mean = median(expect),
                                                      q5 = HDInterval::hdi(expect)[1],
                                                      q95 = HDInterval::hdi(expect)[2])



df %>% drop_na() %>% 
  group_by(Subject) %>% mutate(trial = 1:n()) %>% group_by(trial) %>% 
  summarize(u = mean(u), resp = mean(bin_resp), rt = mean(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>% 
  ggplot(aes(x = trial))+
  geom_pointrange(aes(y = resp, ymin = resp-se_resp, ymax = resp+se_resp))+
  # geom_line(aes(y = u),col = "black", alpha = 0.2)+
  geom_line(data = tester, aes(x = trial, y = expect, group = draw), alpha = 0.5, col = "orange")


df %>% mutate(rts = (rts)) %>% drop_na() %>% 
  group_by(Subject) %>% mutate(trial = 1:n()) %>% group_by(trial) %>% 
  summarize(u = mean(u), resp = mean(bin_resp), rt = mean(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>% 
  ggplot(aes(x = trial))+geom_pointrange(aes(y = rt, ymin = rt-se_rts, ymax = rt+se_rts))+
  geom_ribbon(data = tester_sum, aes(x = trial, y = mean, ymin = q5, ymax = q95), alpha = 0.2)

tester_sum2 = tester %>% group_by(trial) %>% summarize(expect = median(expect),
                                                       rts = median(rts))



df %>% ggplot(aes(x = rts))+
  geom_histogram(aes(y = after_stat(density)), col = "black")+
  geom_line(data = tester, aes(x = rts, group = draw), alpha = 0.1, stat = "density")




### rts

# subject level:


# subject level 
big_df_rt = data.frame()

for(s in unique(df$Subject)){
  print(s)
  draw_id = sample(1:1000,50)
  us = df %>% filter(Subject == s) %>% .$u
  
  minRT = min(df %>% filter(Subject == s) %>% .$rts)
  
  parameters = paste0(c("learningrate","e0","rt_int","rt_beta","rt_ndt","rt_sd","rho"), "[",s,"]")
  
  dfq = as_draws_df(rw_rt$draws(parameters))%>% select(-contains(".")) %>% mutate(draw = 1:n()) %>% 
    filter(draw %in% draw_id) %>% 
    rename_with(~c("learningrate","e0","rt_int","rt_beta","rt_ndt","rt_sd","rho","draw")) %>% 
    group_by(draw) %>% 
    summarize(list(generate_trialwise_data_rt_copula(u = us, learningrate = learningrate,e0 = e0,participant = s,
                                      rt_int = rt_int,rt_beta = rt_beta,rt_sd = rt_sd,rt_ndt = rt_ndt,rho = rho))) %>% unnest()
  
  big_df_rt = rbind(big_df_rt,dfq)
  
}


df %>% drop_na() %>% 
  ggplot(aes(x = trial))+
  geom_point(aes(y = bin_resp))+
  geom_point(aes(y = u-0.05),col = "red")+
  geom_line(data = big_df %>% rename(Subject = participant), aes(x = trial, y = expect, group = draw), alpha = 0.05)+
  facet_wrap(~Subject)


# as rolling:

library(slider)

# Define rolling window size (adjust as needed)
window_size <- 10  # Example: 10-trial moving average

df_rolled <- df %>%
  drop_na() %>%
  group_by(Subject) %>%
  mutate(
    bin_resp_smooth = slide_dbl(bin_resp, mean, .before = window_size, .complete = TRUE),
    u_smooth = slide_dbl(u, mean, .before = window_size, .complete = TRUE)
  )

# Plot with rolling averages
df_rolled %>%
  ggplot(aes(x = trial)) +
  geom_line(data = big_df %>% rename(Subject = participant), 
            aes(x = trial, y = expect, group = draw), alpha = 0.05) +
  geom_line(aes(y = bin_resp_smooth), col = "blue") +  # Smoothed bin_resp
  geom_line(aes(y = u_smooth), col = "red") +  # Smoothed u
  facet_wrap(~Subject)




# rts

df %>% drop_na() %>% 
  ggplot(aes(x = trial))+
  geom_point(aes(y = rts))+
  geom_line(data = big_df %>% rename(Subject = participant), aes(x = trial, y = rts  , group = draw), alpha = 0.05)+
  coord_cartesian(ylim = c(0,2.5))+
  facet_wrap(~Subject)


# tryng

### group elvel



params = c("mu_learningrate",
           "mu_e0",
           "mu_rt_int",
           "mu_rt_beta",
           "mu_rt_sd")


minRT = mean(df %>% group_by(Subject) %>% summarize(minRT = min(rts)) %>% .$minRT)

us = df %>% filter(Subject == 1) %>% .$u

draw_id = sample(1:2000,2000)

data.frame(rw_rt$summary("rho")) %>% mutate(id =1:n())%>% ggplot(aes(x = id, y = mean, ymin = q5,ymax = q95))+geom_pointrange()

rt_group = as_draws_df(rw_rt$draws(params)) %>% select(-contains(".")) %>% 
  rename_with(~c("learningrate","e0","rt_int","rt_beta","rt_sd")) %>% mutate(draw = 1:n()) %>% 
  mutate(learningrate = brms::inv_logit_scaled(learningrate),
         e0 = brms::inv_logit_scaled(e0),
         rt_sd = exp(rt_sd)
  )

tester = rt_group %>% filter(draw %in% draw_id) %>% mutate(rt_ndt = minRT, rho = 0) %>% group_by(draw) %>% 
  summarize(list(generate_trialwise_data_rt_copula(u = us, learningrate = learningrate,e0 = e0,participant = 1,
                                    rt_int = rt_int,rt_beta = rt_beta,rt_ndt = rt_ndt,rt_sd = rt_sd,rho = 0))) %>% unnest()

tester_sum = tester %>% group_by(trial) %>% summarize(mean = median(rts),
                                                      q50_low = HDInterval::hdi(rts,0.5)[1],
                                                      q50_high = HDInterval::hdi(rts,0.5)[2],
                                                      q80_low = HDInterval::hdi(rts,0.8)[1],
                                                      q80_high = HDInterval::hdi(rts,0.8)[2])


df %>% mutate(rts = (rts)) %>% drop_na() %>% 
  group_by(Subject) %>% mutate(trial = 1:n()) %>% group_by(trial) %>% 
  summarize(u = mean(u), resp = mean(bin_resp), rt = mean(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>% 
  ggplot(aes(x = trial))+geom_pointrange(aes(y = rt, ymin = rt-se_rts, ymax = rt+se_rts))+
  geom_ribbon(data = tester_sum, aes(x = trial, y = mean, ymin = q80_low, ymax = q80_high), alpha = 0.2, fill = "black")+
  geom_ribbon(data = tester_sum, aes(x = trial, y = mean, ymin = q50_low, ymax = q50_high), alpha = 0.5, fill = "black")
  





library(tidyverse)
library(HDInterval)

# Generate HDIs at intervals of 1% from 50% to 95%
hdi_levels <- seq(0.50, 0.95, by = 0.05)

hdi_df <- map_dfr(hdi_levels, function(level) {
  tester %>% 
    group_by(trial) %>% 
    summarize(
      mean = median(rts),
      hdi_low = HDInterval::hdi(rts, level)[1],
      hdi_high = HDInterval::hdi(rts, level)[2]
    ) %>%
    mutate(interval = level) # Track the interval size
})

# Plot with continuous shading
df %>% 
  drop_na() %>% 
  group_by(Subject) %>% 
  mutate(trial = row_number()) %>% 
  group_by(trial) %>% 
  summarize(
    u = mean(u), 
    resp = mean(bin_resp), 
    rt = mean(rts), 
    se_resp = (mean(bin_resp) * (1 - mean(bin_resp))) / sqrt(n()), 
    se_rts = sd(rts) / sqrt(n())
  ) %>% 
  ggplot(aes(x = trial)) +
    # Mean RT with error bars
  geom_pointrange(aes(y = rt, ymin = rt - se_rts, ymax = rt + se_rts)) +
  # HDI intervals with increasing transparency
  geom_ribbon(data = hdi_df, aes(x = trial, y = mean, ymin = hdi_low, ymax = hdi_high, alpha = as.factor(interval)), fill = "lightblue2") +
  scale_alpha_discrete(range = c(0.01, 0.1), guide = "none") + # Adjust transparency
  theme_minimal()












tester_sum2 = tester %>% group_by(trial) %>% summarize(mean = median(expect),
                                                       q5 = HDInterval::hdi(expect)[1],
                                                       q95 = HDInterval::hdi(expect)[2])



df %>% drop_na() %>% 
  group_by(Subject) %>% mutate(trial = 1:n()) %>% group_by(trial) %>% 
  summarize(u = mean(u), resp = mean(bin_resp), rt = mean(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>% 
  ggplot(aes(x = trial))+
  geom_pointrange(aes(y = resp, ymin = resp-se_resp, ymax = resp+se_resp))+
  # geom_line(aes(y = u),col = "black", alpha = 0.2)+
  geom_line(data = tester, aes(x = trial, y = expect, group = draw), alpha = 0.5, col = "orange")


df %>% mutate(rts = (rts)) %>% drop_na() %>% 
  group_by(Subject) %>% mutate(trial = 1:n()) %>% group_by(trial) %>% 
  summarize(u = mean(u), resp = mean(bin_resp), rt = mean(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>% 
  ggplot(aes(x = trial))+geom_pointrange(aes(y = rt, ymin = rt-se_rts, ymax = rt+se_rts))+
  geom_ribbon(data = tester_sum, aes(x = trial, y = mean, ymin = q5, ymax = q95), alpha = 0.2)

tester_sum2 = tester %>% group_by(trial) %>% summarize(expect = median(expect),
                                                       rts = median(rts))



df %>% drop_na() %>% 
  group_by(Subject) %>% mutate(trial = 1:n()) %>% group_by(trial) %>% 
  summarize(u = mean(u), resp = mean(bin_resp), rt = mean(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>% 
  ggplot(aes(x = resp))+geom_pointrange(aes(y = rt, ymin = rt-se_rts, ymax = rt+se_rts))+
  geom_line(data = tester, aes(x = expect, y = rts, group = draw), alpha = 0.2)


df %>% ggplot(aes(x = rts))+
  geom_histogram(aes(y = after_stat(density)), col = "black")+
  geom_line(data = tester, aes(x = rts, group = draw), alpha = 0.1, stat = "density")


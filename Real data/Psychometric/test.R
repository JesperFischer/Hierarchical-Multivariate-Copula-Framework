packages = c("brms","tidyverse","bayesplot","pracma","here",
             "patchwork","posterior","HDInterval","loo", "furrr", "SBC","future")

do.call(pacman::p_load, as.list(packages))

df <- read_csv(here::here("Real data","Psychometric","data","data_Bang_2019_Exp1.csv"))


df %>%
  group_by(Condition) %>% filter(SN != max(SN) & SN != min(SN)) %>%
  mutate(ACC = ifelse(Stimulus == Response,1,0),
         Difficulty  = ifelse(Stimulus == 2, SN , -SN),
         Response = Response-1,
         Confidence = Confidence/max(Confidence,na.rm = T)) %>% 
  
  filter(Difficulty != min(Difficulty)) %>% 
  ungroup() %>% 
  pivot_longer(cols = c(Response,Confidence,RT_dec))%>% 
  mutate(Difficulty_bin = cut(Difficulty, breaks = 20, labels = FALSE)) %>% 
  group_by(Difficulty_bin ,Stimulus,name,Condition) %>% 
  summarize(mean = mean(value,na.rm = T),
            se = sd(value,na.rm = T)/sqrt(n())) %>% 
  ggplot(aes(x = Difficulty_bin , col = as.factor(name)))+
  geom_pointrange(aes(y = mean,ymin = mean-2*se, ymax = mean+2*se))+
  facet_wrap(name~Condition, scales = "free")+theme_minimal()+
  theme(legend.position = "top",
        text = element_text(size = 16),        # General text size
        axis.title = element_text(size = 18), # Axis titles
        axis.text = element_text(size = 14),  # Axis tick labels
        plot.title = element_text(size = 20)  # Plot title)
  )


df %>%
  group_by(Condition) %>% 
  mutate(ACC = ifelse(Stimulus == Response,1,0),
         Difficulty  = ifelse(Stimulus == 2, SN , -SN),
         Response = Response-1,
         Confidence = Confidence/max(Confidence,na.rm = T)) %>% 
  ungroup() %>% 
  pivot_longer(cols = c(Response,Confidence,RT_conf,RT_dec))%>% 
  mutate(Difficulty_bin = cut(Difficulty, breaks = 20, labels = FALSE)) %>% 
  group_by(Difficulty_bin ,name,Condition,ACC) %>%  
  summarize(mean = mean(value,na.rm = T),
            se = sd(value,na.rm = T)/sqrt(n())) %>% 
  ggplot(aes(x = Difficulty_bin , col = as.factor(ACC)))+
  geom_pointrange(aes(y = mean,ymin = mean-2*se, ymax = mean+2*se))+
  facet_wrap(Condition~name, scales = "free")+theme_minimal()+theme(legend.position = "top")



df %>% filter(Condition == 1) %>% filter(SN != max(SN) & SN != min(SN)) %>%
  mutate(ACC = ifelse(Stimulus == Response,1,0),
         Difficulty  = ifelse(Stimulus == 2, SN , -SN),
         Response = Response-1,
         Confidence = Confidence/max(Confidence,na.rm = T)) %>% 
  
  filter(Difficulty != min(Difficulty)) %>% 
  ungroup() %>% 
  pivot_longer(cols = c(Response,RT_dec))%>% 
  mutate(Difficulty_bin = cut(Difficulty, breaks = 20, labels = FALSE)) %>% 
  group_by(Difficulty_bin ,Stimulus,name) %>% 
  summarize(mean = mean(value,na.rm = T),
            se = sd(value,na.rm = T)/sqrt(n())) %>% 
  ggplot(aes(x = Difficulty_bin , col = as.factor(name)))+
  geom_pointrange(aes(y = mean,ymin = mean-2*se, ymax = mean+2*se))+
  facet_wrap(~name, scales = "free")+theme_minimal()+
  theme(legend.position = "top",
        text = element_text(size = 16),        # General text size
        axis.title = element_text(size = 18), # Axis titles
        axis.text = element_text(size = 14),  # Axis tick labels
        plot.title = element_text(size = 20)  # Plot title)
  )

df %>% filter(Condition == 1, Day == 1) %>% filter(SN != max(SN) & SN != min(SN)) %>%
  filter(RT_dec < 5) %>% 
  mutate(ACC = ifelse(Stimulus == Response,1,0),
         Difficulty  = ifelse(Stimulus == 2, SN , -SN),
         Response = Response-1,
         Confidence = Confidence/max(Confidence,na.rm = T)) %>% 
  # filter(Difficulty != min(Difficulty)) %>% 
  ungroup() %>% 
  pivot_longer(cols = c(Response,RT_dec))%>% 
  mutate(Difficulty_bin = cut(Difficulty, breaks = 20, labels = FALSE)) %>% 
  group_by(Difficulty_bin ,Stimulus,name,Subj_idx) %>% 
  summarize(mean = mean(value,na.rm = T),
            se = sd(value,na.rm = T)/sqrt(n())) %>% 
  ggplot(aes(x = Difficulty_bin , col = as.factor(name)))+
  geom_pointrange(aes(y = mean,ymin = mean-2*se, ymax = mean+2*se))+
  facet_wrap(Subj_idx~name, scales = "free")+theme_minimal()+
  theme(legend.position = "top",
        text = element_text(size = 16),        # General text size
        axis.title = element_text(size = 18), # Axis titles
        axis.text = element_text(size = 14),  # Axis tick labels
        plot.title = element_text(size = 20)  # Plot title)
  )



df = df %>% filter(Condition == 1, Day == 1) %>% filter(SN != max(SN) & SN != min(SN)) %>%
  mutate(ACC = ifelse(Stimulus == Response,1,0),
         Difficulty  = ifelse(Stimulus == 2, SN , -SN),
         Response = Response-1,
         Confidence = Confidence/max(Confidence,na.rm = T)) %>% rename(participant = Subj_idx) %>% 
  group_by(participant) %>% mutate(trial = 1:n())


df_filtered <- df %>%
  group_by(participant) %>%
  filter(RT_dec > min(RT_dec))

df = df_filtered

N = df %>% group_by(participant) %>% summarize(N = n()) %>% .$N

S = df %>% group_by(participant) %>% summarize(N = n()) %>% ungroup() %>% summarize(N = n()) %>% .$N

ranges <- data.frame(
  start = c(1, head(cumsum(N) + 1, -1)),
  end = cumsum(N)
)


dflist = list(
  N = nrow(df),
  x = df$Difficulty,
  binom_y = df$Response,
  RT = (df$RT_dec),
  minRT = unique(df %>% group_by(participant) %>% summarize(minRT = min(RT_dec)) %>% .$minRT),
  S_id = df$participant,
  starts = ranges$start,
  ends = ranges$end,
  S = length(unique(df$participant))
)



rt_model <- cmdstanr::cmdstan_model(here::here("Real data","Psychometric","stanmodels","subj_ndt.stan"), force_recompile = T)
nort_model <- cmdstanr::cmdstan_model(here::here("Real data","Psychometric","stanmodels","no_rt.stan"), force_recompile = T)
ddm_model <- cmdstanr::cmdstan_model(here::here("Real data","Psychometric","stanmodels","hier_ddm.stan"), force_recompile = T)



psy_rt  = rt_model$sample(data = dflist,
                         iter_warmup = 1000,
                         iter_sampling = 1000,
                         chains = 4,
                         refresh = 50,
                         max_treedepth = 10,
                         parallel_chains = 4,
                         adapt_delta = 0.95)

psy_rt$save_object(here::here("Real data","Psychometric","psy_rt1.rds"))
divs_rt = psy_rt$diagnostic_summary()



psy_nort = nort_model$sample(data = dflist,
                             iter_warmup = 1000,
                             iter_sampling = 1000,
                             chains = 4,
                             refresh = 100,
                             max_treedepth = 10,
                             parallel_chains = 4,
                             adapt_delta = 0.99)

psy_nort$save_object(here::here("Real data","Psychometric","psy_nort1.rds"))
divs_nort = psy_nort$diagnostic_summary()


psy_ddm = ddm_model$sample(data = dflist,
                          iter_warmup = 1000,
                          iter_sampling = 1000,
                          chains = 4,
                          refresh = 50,
                          max_treedepth = 10,
                          parallel_chains = 4,
                          adapt_delta = 0.95)
psy_ddm$save_object(here::here("Real data","Psychometric","psy_ddm1.rds"))


rt_loo_bin = psy_rt$loo("log_lik_bin")
rt_loo_full = psy_rt$loo("log_lik")

nort_loo = psy_nort$loo()

psy_loo = psy_ddm$loo()


full_loo = data.frame(loo::loo_compare(list(rt_loo_full,psy_loo)))

full_loo %>% filter(elpd_diff != 0) %>% mutate(ratio = elpd_diff/se_diff)

bin_loo = loo::loo_compare(list(rt_loo_bin,nort_loo))

## plotting:

draw_id = sample(1:4000,100)

# subject level 
big_df_ddm = data.frame()
subject_means = data.frame()

for(s in unique(df$participant)){
  print(s)
  xs = seq(-0.3,0.3,by = 0.01)
  
  minRT = min(df %>% filter(participant == s) %>% .$RT_dec)
  
  parameters = paste0(c("threshold","slope","lapse","alpha","beta","delta","rt_ndt"), "[",s,"]")
  
  dfq = as_draws_df(psy_ddm$draws(parameters))%>% select(-contains(".")) %>% mutate(draw = 1:n()) %>% 
    filter(draw %in% draw_id) %>% 
    rename_with(~c("threshold","slope","lapse","alpha","beta","delta","rt_ndt","draw")) %>% 
    group_by(draw) %>% 
    summarize(list(generate_trial_ddm(x = xs, threshold = threshold,slope = slope, lapse = lapse,participant = s,
                                      alpha = alpha,beta = beta,delta = delta,rt_ndt = rt_ndt))) %>% unnest() %>% 
    rename(RT_dec = q) %>% 
    mutate(resp = ifelse(resp == "upper",1,0))

  
  means = as_draws_df(psy_ddm$summary(parameters))%>% select(-contains(".")) %>% select(median,variable) %>% 
    pivot_wider(names_from = variable,values_from = median) %>% 
    rename_with(~c("threshold","slope","lapse","alpha","beta","delta","rt_ndt")) %>% 
    summarize(list(generate_trial_ddm(x = xs, threshold = threshold,slope = slope, lapse = lapse,participant = s, 
                                      alpha = alpha,beta = beta,delta = delta,rt_ndt = rt_ndt))) %>% unnest() %>% 
    rename(RT_dec = q) %>% 
    mutate(resp = ifelse(resp == "upper",1,0))
  
  subject_means = rbind(subject_means,means)
  
  big_df_ddm = rbind(big_df_ddm,dfq)
  
}

subject_means %>% group_by(participant,x) %>% summarize(resp = mean(resp), RT = median(RT_dec)) %>% 
  pivot_longer(cols = c("resp","RT")) %>% ggplot(aes(x = x, y = value))+
  facet_wrap(participant~name, scales = "free")+geom_point()

big_df_ddm = big_df_ddm %>% mutate(resp = ifelse(resp == 0, 1, 0))

subject_means = subject_means %>% mutate(name = "RT_dec")


df %>%   ungroup() %>%  mutate(rts = (RT_dec),bin_resp = Response) %>% drop_na() %>%
  mutate(Difficulty_bin = cut((Difficulty), breaks = 5, labels = FALSE),
         Difficulty = abs(Difficulty),
         Difficulty_bin = (Difficulty_bin - min(Difficulty_bin)) / (max(Difficulty_bin) - min(Difficulty_bin)) * (0.3 - (-0.3)) + (-0.3)
  ) %>%
  group_by(Difficulty_bin,participant ) %>% 
  summarize(resp = mean(bin_resp), rt = median(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>%
  ggplot()+
  geom_pointrange(aes(x = Difficulty_bin, y = resp, ymin = resp-se_resp, ymax = resp+se_resp))+
  # geom_line(data = subject_means_rt, aes(x = x, y = expectation))+
  geom_line(data = big_df_ddm %>% filter(draw %in% 1:100), aes(x = x, y = expectation, group = draw), alpha = 0.05)+
  facet_wrap(~participant)+
  theme_minimal()

df %>%   ungroup() %>%  mutate(rts = (RT_dec),bin_resp = Response) %>% drop_na() %>%
  mutate(Difficulty_bin = cut((Difficulty), breaks = 5, labels = FALSE),
         Difficulty = abs(Difficulty),
         Difficulty_bin = (Difficulty_bin - min(Difficulty_bin)) / (max(Difficulty_bin) - min(Difficulty_bin)) * (0.3 - (-0.3)) + (-0.3)
  ) %>%
  group_by(Difficulty_bin,participant ) %>% 
  summarize(resp = mean(bin_resp), rt = median(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>%
  ggplot()+
  geom_pointrange(aes(x = Difficulty_bin, y = rt, ymin = rt-2*se_rts, ymax = rt+2*se_rts))+
  # geom_line(data = subject_means_rt, aes(x = x, y = rts))+
  geom_line(data = big_df_ddm %>% filter(draw %in% 1:100), aes(x = x, y = RT_dec, group = draw), alpha = 0.05)+
  scale_y_continuous(limits = c(0,5))+
  facet_wrap(~participant)+
  theme_minimal()




qq_summar = big_df_ddm %>%  ungroup() %>%  mutate(rts = (RT_dec),bin_resp = resp, Difficulty = x) %>% drop_na() %>%
  mutate(Difficulty_bin = cut((Difficulty), breaks = 5, labels = FALSE),
         Difficulty = abs(Difficulty),
         Difficulty_bin = (Difficulty_bin - min(Difficulty_bin)) / (max(Difficulty_bin) - min(Difficulty_bin)) * (0.3 - (-0.3)) + (-0.3)
  ) %>%
  group_by(Difficulty_bin,draw) %>% 
  summarize(resp = mean(bin_resp), rt = median(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n()))



qq_summar_sum = qq_summar %>%   group_by(Difficulty_bin) %>% 
  summarize(resps = mean(resp), rts = median(rt), se_resp =  (mean(resp) * (1- mean(resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n()))


df %>%   ungroup() %>%  mutate(rts = (RT_dec),bin_resp = Response) %>% drop_na() %>%
  mutate(Difficulty_bin = cut((Difficulty), breaks = 5, labels = FALSE),
         Difficulty = abs(Difficulty),
         Difficulty_bin = (Difficulty_bin - min(Difficulty_bin)) / (max(Difficulty_bin) - min(Difficulty_bin)) * (0.3 - (-0.3)) + (-0.3)
  ) %>%
  group_by(Difficulty_bin) %>% 
  summarize(resp = mean(bin_resp), rt = median(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>%
  ggplot(aes(x = Difficulty_bin))+
  geom_line(data = qq_summar, aes(x = Difficulty_bin, y = rt, group = draw), alpha = 0.01, col = "orange")+
  geom_line(data = qq_summar_sum, aes(x = Difficulty_bin, y = rts), col = "red")+
  geom_pointrange(aes(y = rt, ymin = rt-2*se_rts, ymax = rt+2*se_rts))+
  theme_classic()


df %>%   ungroup() %>%  mutate(rts = (RT_dec),bin_resp = Response) %>% drop_na() %>%
  mutate(Difficulty_bin = cut((Difficulty), breaks = 5, labels = FALSE),
         Difficulty = abs(Difficulty),
         Difficulty_bin = (Difficulty_bin - min(Difficulty_bin)) / (max(Difficulty_bin) - min(Difficulty_bin)) * (0.3 - (-0.3)) + (-0.3)
  ) %>%
  group_by(Difficulty_bin) %>% 
  summarize(resp = mean(bin_resp), rt = median(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>%
  ggplot(aes(x = Difficulty_bin))+
  geom_line(data = qq_summar, aes(x = Difficulty_bin, y = resp, group = draw), alpha = 0.01, col = "orange")+
  geom_line(data = qq_summar_sum, aes(x = Difficulty_bin, y = resps), col = "red")+
  geom_pointrange(aes(y = resp, ymin = resp-se_resp, ymax = resp+se_resp))+
  theme_classic()



# rts
# subject level 
big_df_rt = data.frame()
subject_means_rt = data.frame()
for(s in unique(df$participant)){
  print(s)
  xs = seq(-0.3,0.3,by = 0.01)
  
  minRT = min(df %>% filter(participant == s) %>% .$RT_dec)
  
  parameters = paste0(c("threshold","slope","lapse","rt_int","rt_beta","rt_ndt","rt_sd","rho"), "[",s,"]")
  
  dfq = as_draws_df(psy_rt$draws(parameters))%>% select(-contains(".")) %>% mutate(draw = 1:n()) %>% 
    filter(draw %in% draw_id) %>% 
    rename_with(~c("threshold","slope","lapse","rt_int","rt_beta","rt_ndt","rt_sd","rho","draw")) %>% 
    group_by(draw) %>% 
    summarize(list(generate_trial_shannon_entropy_rt(x = xs, threshold = threshold,slope = slope, lapse = lapse,participant = s,
                                                     rt_int = rt_int,rt_beta = rt_beta,rt_ndt = rt_ndt,rt_sd = rt_sd,rho = rho))) %>% unnest()
  
  means = as_draws_df(psy_rt$summary(parameters))%>% select(-contains(".")) %>% select(median,variable) %>% 
    pivot_wider(names_from = variable,values_from = median) %>% 
    rename_with(~c("threshold","slope","lapse","rt_int","rt_beta","rt_ndt","rt_sd","rho")) %>% 
    summarize(list(generate_trial_shannon_entropy_rt(x = xs, threshold = threshold,slope = slope, lapse = lapse,participant = s,
                                                     rt_int = rt_int,rt_beta = rt_beta,rt_ndt = rt_ndt,rt_sd = rt_sd,rho = rho))) %>% unnest()
  
  subject_means_rt = rbind(subject_means_rt,means)
  
  big_df_rt = rbind(big_df_rt,dfq)
  
}

subject_means_rt = subject_means_rt %>% mutate(name = "RT_dec")

df %>%   ungroup() %>%  mutate(rts = (RT_dec),bin_resp = Response) %>% drop_na() %>%
  mutate(Difficulty_bin = cut((Difficulty), breaks = 5, labels = FALSE),
         Difficulty = abs(Difficulty),
         Difficulty_bin = (Difficulty_bin - min(Difficulty_bin)) / (max(Difficulty_bin) - min(Difficulty_bin)) * (0.3 - (-0.3)) + (-0.3)
  ) %>%
  group_by(Difficulty_bin,participant ) %>% 
  summarize(resp = mean(bin_resp), rt = median(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>%
  ggplot()+
  geom_pointrange(aes(x = Difficulty_bin, y = resp, ymin = resp-se_resp, ymax = resp+se_resp))+
  # geom_line(data = subject_means_rt, aes(x = x, y = expectation))+
  geom_line(data = big_df_rt %>% filter(draw %in% 1:100), aes(x = x, y = expectation, group = draw), alpha = 0.05)+
  facet_wrap(~participant)+
  theme_minimal()

df %>%   ungroup() %>%  mutate(rts = (RT_dec),bin_resp = Response) %>% drop_na() %>%
  mutate(Difficulty_bin = cut((Difficulty), breaks = 5, labels = FALSE),
         Difficulty = abs(Difficulty),
         Difficulty_bin = (Difficulty_bin - min(Difficulty_bin)) / (max(Difficulty_bin) - min(Difficulty_bin)) * (0.3 - (-0.3)) + (-0.3)
  ) %>%
  group_by(Difficulty_bin,participant ) %>% 
  summarize(resp = mean(bin_resp), rt = median(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>%
  ggplot()+
  geom_pointrange(aes(x = Difficulty_bin, y = rt, ymin = rt-2*se_rts, ymax = rt+2*se_rts))+
  # geom_line(data = subject_means_rt, aes(x = x, y = rts))+
  geom_line(data = big_df_rt %>% filter(draw %in% 1:100), aes(x = x, y = rts, group = draw), alpha = 0.05)+
  scale_y_continuous(limits = c(0,5))+
  facet_wrap(~participant)+
  theme_minimal()





qq_summar_rt = big_df_rt %>% rename(Response = expectation,RT_dec = rts) %>%   
  ungroup() %>%  
  mutate(rts = (RT_dec),bin_resp = resp, Difficulty = x) %>% drop_na() %>%
  mutate(Difficulty_bin = cut((Difficulty), breaks = 5, labels = FALSE),
         Difficulty = abs(Difficulty),
         Difficulty_bin = (Difficulty_bin - min(Difficulty_bin)) / (max(Difficulty_bin) - min(Difficulty_bin)) * (0.3 - (-0.3)) + (-0.3)
  ) %>%
  group_by(Difficulty_bin,draw) %>% 
  summarize(resp = mean(bin_resp), rt = median(rts))



qq_summar_sum_rt = qq_summar_rt %>%   group_by(Difficulty_bin) %>% 
  summarize(resps = mean(resp), rts = median(rt))


df %>%   ungroup() %>%  mutate(rts = (RT_dec),bin_resp = Response) %>% drop_na() %>%
  mutate(Difficulty_bin = cut((Difficulty), breaks = 5, labels = FALSE),
         Difficulty = abs(Difficulty),
         Difficulty_bin = (Difficulty_bin - min(Difficulty_bin)) / (max(Difficulty_bin) - min(Difficulty_bin)) * (0.3 - (-0.3)) + (-0.3)
  ) %>%
  group_by(Difficulty_bin) %>% 
  summarize(resp = mean(bin_resp), rt = median(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>%
  ggplot(aes(x = Difficulty_bin))+
  geom_line(data = qq_summar_rt, aes(x = Difficulty_bin, y = rt, group = draw), alpha = 0.01, col = "#6CEEF8")+
  geom_line(data = qq_summar_sum_rt, aes(x = Difficulty_bin, y = rts), col = "blue")+
  geom_pointrange(aes(y = rt, ymin = rt-se_rts, ymax = rt+se_rts))+
  theme_classic()


df %>%   ungroup() %>%  mutate(rts = (RT_dec),bin_resp = Response) %>% drop_na() %>%
  mutate(Difficulty_bin = cut((Difficulty), breaks = 5, labels = FALSE),
         Difficulty = abs(Difficulty),
         Difficulty_bin = (Difficulty_bin - min(Difficulty_bin)) / (max(Difficulty_bin) - min(Difficulty_bin)) * (0.3 - (-0.3)) + (-0.3)
  ) %>%
  group_by(Difficulty_bin) %>% 
  summarize(resp = mean(bin_resp), rt = median(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>%
  ggplot(aes(x = Difficulty_bin))+
  geom_line(data = qq_summar_rt, aes(x = Difficulty_bin, y = resp, group = draw), alpha = 0.01, col = "#6CEEF8")+
  geom_line(data = qq_summar_sum_rt, aes(x = Difficulty_bin, y = resps), col = "blue")+
  geom_pointrange(aes(y = resp, ymin = resp-se_resp, ymax = resp+se_resp))+
  theme_classic()



# no rts:

big_df_nort = data.frame()
subject_means_nort = data.frame()
for(s in unique(df$participant)){
  print(s)
  xs = seq(-0.3,0.3,by = 0.01)
  
  
  parameters = paste0(c("threshold","slope","lapse"), "[",s,"]")
  
  dfq = as_draws_df(psy_nort$draws(parameters))%>% select(-contains(".")) %>% mutate(draw = 1:n()) %>% 
    filter(draw %in% draw_id) %>% 
    rename_with(~c("threshold","slope","lapse","draw")) %>% 
    group_by(draw) %>% 
    summarize(list(generate_psycho(x = xs, threshold = threshold,slope = slope, lapse = lapse,participant = s))) %>% unnest()
  
  means = as_draws_df(psy_nort$summary(parameters))%>% select(-contains(".")) %>% select(median,variable) %>% 
    pivot_wider(names_from = variable,values_from = median) %>% 
    rename_with(~c("threshold","slope","lapse")) %>% 
    summarize(list(generate_psycho(x = xs, threshold = threshold,slope = slope, lapse = lapse,participant = s))) %>% unnest()
  
  big_df_nort = rbind(big_df_nort,dfq)
  
  subject_means_nort = rbind(subject_means_nort,means)
  
}

qq_summar_nort = big_df_nort %>% 
  ungroup() %>%  mutate(bin_resp = resp, Difficulty = x) %>% drop_na() %>%
  mutate(Difficulty_bin = cut((Difficulty), breaks = 5, labels = FALSE),
         Difficulty = abs(Difficulty),
         Difficulty_bin = (Difficulty_bin - min(Difficulty_bin)) / (max(Difficulty_bin) - min(Difficulty_bin)) * (0.3 - (-0.3)) + (-0.3)
  ) %>%
  group_by(Difficulty_bin,draw) %>% 
  summarize(resp = mean(bin_resp))


qq_summar_sum_nort = qq_summar_nort %>%   
  group_by(Difficulty_bin) %>% 
  summarize(resps = mean(resp))



df %>%   ungroup() %>%  mutate(rts = (RT_dec),bin_resp = Response) %>% drop_na() %>%
  mutate(Difficulty_bin = cut((Difficulty), breaks = 5, labels = FALSE),
         Difficulty = abs(Difficulty),
         Difficulty_bin = (Difficulty_bin - min(Difficulty_bin)) / (max(Difficulty_bin) - min(Difficulty_bin)) * (0.3 - (-0.3)) + (-0.3)
  ) %>%
  group_by(Difficulty_bin) %>% 
  summarize(resp = mean(bin_resp), rt = median(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>%
  ggplot(aes(x = Difficulty_bin))+
  geom_line(data = qq_summar_nort, aes(x = Difficulty_bin, y = resp, group = draw), alpha = 0.01, col = "#00C853")+
  geom_line(data = qq_summar_sum_nort, aes(x = Difficulty_bin, y = resps), col = "green")+
  geom_pointrange(aes(y = resp, ymin = resp-se_resp, ymax = resp+se_resp))+
  theme_classic()



## Both

## psychometric fits
alpha = 0.1

psycho_psycho = df %>%   ungroup() %>%  mutate(rts = (RT_dec),bin_resp = Response) %>% drop_na() %>%
  mutate(Difficulty_bin = cut((Difficulty), breaks = 10, labels = FALSE),
         Difficulty = abs(Difficulty),
         Difficulty_bin = (Difficulty_bin - min(Difficulty_bin)) / (max(Difficulty_bin) - min(Difficulty_bin)) * (0.3 - (-0.3)) + (-0.3)
  ) %>%
  group_by(Difficulty_bin) %>% 
  summarize(resp = mean(bin_resp), rt = median(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>%
  ggplot(aes(x = Difficulty_bin))+
  geom_line(data = qq_summar_rt %>% filter(draw %in% draw_id), aes(x = Difficulty_bin, y = resp, group = draw), alpha = alpha, col = "#6CEEF8")+

  # geom_line(data = qq_summar_nort, aes(x = Difficulty_bin, y = resp, group = draw), alpha = 0.01, col = "black")+

  geom_line(data = qq_summar%>% filter(draw %in% draw_id), aes(x = Difficulty_bin, y = resp, group = draw), alpha = alpha, col = "orange")+
  geom_line(data = qq_summar_sum, aes(x = Difficulty_bin, y = resps, col = "DDM"), linewidth = 1.1)+
  geom_line(data = qq_summar_sum_rt, aes(x = Difficulty_bin, y = resps, col = "RT"), linewidth = 1.1)+
  geom_line(data = qq_summar_sum_nort, aes(x = Difficulty_bin, y = resps, col = "No RT"), linewidth = 1.1)+
  geom_pointrange(aes(y = resp, ymin = resp-2*se_resp, ymax = resp+2*se_resp))+
  theme_minimal()+
  scale_color_manual(name = "Fitted model",values = c("red","black","blue"))+
  scale_x_continuous("Binned stimulus intensity", breaks = scales::pretty_breaks(n = 3))+
  scale_y_continuous("P(Response == 1)", breaks = scales::pretty_breaks(n = 3))+  
  ggtitle("Psychophysics")+
  theme(plot.title = element_text(hjust = 0.5))


psycho_psycho
## rts

rts_psycho = df %>%   ungroup() %>%  mutate(rts = (RT_dec),bin_resp = Response) %>% drop_na() %>%
  mutate(Difficulty_bin = cut((Difficulty), breaks = 10, labels = FALSE),
         Difficulty = abs(Difficulty),
         Difficulty_bin = (Difficulty_bin - min(Difficulty_bin)) / (max(Difficulty_bin) - min(Difficulty_bin)) * (0.3 - (-0.3)) + (-0.3)
  ) %>%
  group_by(Difficulty_bin) %>% 
  summarize(resp = mean(bin_resp), rt = median(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>%
  ggplot(aes(x = Difficulty_bin))+
  geom_line(data = qq_summar_rt%>% filter(draw %in% draw_id), aes(x = Difficulty_bin, y = rt, group = draw), alpha = alpha, col = "#6CEEF8")+
  geom_line(data = qq_summar%>% filter(draw %in% draw_id), aes(x = Difficulty_bin, y = rt, group = draw), alpha = alpha, col = "orange")+
  geom_line(data = qq_summar_sum_rt, aes(x = Difficulty_bin, y = rts), col = "blue", linewidth = 1.1)+
  geom_line(data = qq_summar_sum, aes(x = Difficulty_bin, y = rts), col = "red", linewidth = 1.1)+
  geom_pointrange(aes(y = rt, ymin = rt-2*se_rts, ymax = rt+2*se_rts))+
  scale_color_manual(name = "Fitted model",values = c("red","black","blue"))+
  theme_minimal()+
  scale_x_continuous("Binned stimulus intensity", breaks = scales::pretty_breaks(n = 3))+
  scale_y_continuous("Response time (S)", breaks = scales::pretty_breaks(n = 3))


psycho_psycho/rts_psycho+
  plot_layout(axis_titles = "collect")


# trace and prior posterior updates: (DDM)

parameters = paste0("mu_",c("threshold","slope","lapse","alpha","beta","delta"))

traces_psycho_ddm = as_draws_df(psy_ddm$draws(parameters)) %>% pivot_longer(all_of(parameters), names_to = "variables") %>% 
  mutate(.chain = as.factor(.chain)) %>% 
  ggplot(aes(x = .iteration, y = value, col = .chain)) + facet_wrap(~variables, nrow = 1, scales = "free")+geom_line()+
  theme_minimal()+
  scale_color_manual(values = c(bayesplot::color_scheme_get("red")[[1]],
                                bayesplot::color_scheme_get("red")[[2]],
                                bayesplot::color_scheme_get("red")[[3]],
                                bayesplot::color_scheme_get("red")[[4]]))+
  scale_x_continuous("Iterations", breaks = scales::pretty_breaks(n = 3))+
  scale_y_continuous("Value", breaks = scales::pretty_breaks(n = 3))+
  theme(legend.position = "top")



posteriors_psycho_ddm = as_draws_df(psy_ddm$draws(parameters)) %>% pivot_longer(all_of(parameters), names_to = "variables") %>% 
  mutate(distribution = "posterior")

priors_psycho_ddm = data.frame(variables = parameters) %>% mutate(value = list(rnorm(4000,0,5),
                                                                     rnorm(4000,1,1),
                                                                     rnorm(4000,-3,1),
                                                                     rnorm(4000,0.5,1),
                                                                     rnorm(4000,0,0.5),
                                                                     rnorm(4000,5,1))) %>% unnest()%>% 
  mutate(distribution = "prior")

   
pp_update_psycho_ddm = rbind(priors_psycho_ddm,posteriors_psycho_ddm %>% select(variables,value,distribution)) %>% ggplot(aes(x = value, fill = distribution))+
  geom_histogram(col = "black", alpha = 0.5)+
  facet_wrap(~variables, nrow = 1, scales = "free")+
  theme_minimal()+
  scale_x_continuous("Value", breaks = scales::pretty_breaks(n = 3))+
  scale_y_continuous("Frequency", breaks = scales::pretty_breaks(n = 3))+
  theme(legend.position = "top")



convergence_psycho_ddm = traces_psycho_ddm / pp_update_psycho_ddm
convergence_psycho_ddm


# response


parameters = paste0("mu_",c("threshold","slope","lapse","rt_int","rt_beta","rt_sd"))

traces_psycho_rt = as_draws_df(psy_rt$draws(parameters)) %>% pivot_longer(all_of(parameters), names_to = "variables") %>% 
  mutate(.chain = as.factor(.chain)) %>% 
  ggplot(aes(x = .iteration, y = value, col = .chain)) + facet_wrap(~variables, nrow = 1, scales = "free")+geom_line()+
  theme_minimal()+
  scale_color_manual(values = c(bayesplot::color_scheme_get("red")[[1]],
                                bayesplot::color_scheme_get("red")[[2]],
                                bayesplot::color_scheme_get("red")[[3]],
                                bayesplot::color_scheme_get("red")[[4]]))+
  scale_x_continuous("Iterations", breaks = scales::pretty_breaks(n = 3))+
  scale_y_continuous("Value", breaks = scales::pretty_breaks(n = 3))+
  theme(legend.position = "top")



posteriors_psycho_rt = as_draws_df(psy_rt$draws(parameters)) %>% pivot_longer(all_of(parameters), names_to = "variables") %>% 
  mutate(distribution = "posterior")

priors_psycho_rt = data.frame(variables = parameters) %>% mutate(value = list(rnorm(4000,0,5),
                                                                               rnorm(4000,1,1),
                                                                               rnorm(4000,-3,1),
                                                                               rnorm(4000,-1,1),
                                                                               rnorm(4000,1.5,1),
                                                                               rnorm(4000,-1,1))) %>% unnest()%>% 
  mutate(distribution = "prior")


pp_update_psycho_rt = rbind(priors_psycho_rt,posteriors_psycho_rt %>% select(variables,value,distribution)) %>% ggplot(aes(x = value, fill = distribution))+
  geom_histogram(col = "black", alpha = 0.5)+
  facet_wrap(~variables, nrow = 1, scales = "free")+
  theme_minimal()+
  scale_x_continuous("Value", breaks = scales::pretty_breaks(n = 3))+
  scale_y_continuous("Frequency", breaks = scales::pretty_breaks(n = 3))+
  theme(legend.position = "top")



convergence_psycho_rt = traces_psycho_rt / pp_update_psycho_rt
convergence_psycho_rt
convergence_psycho_ddm

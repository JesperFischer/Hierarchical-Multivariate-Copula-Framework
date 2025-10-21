
packages = c("brms","tidyverse","bayesplot","pracma","here", "patchwork","posterior","HDInterval","loo", "furrr","cmdstanr")

do.call(pacman::p_load, as.list(packages))

df = read_csv(here::here("Real data","Psychometric","data","VMP","raw_hrd.csv"))

df = df %>% filter(session == 1) %>% 
  mutate(DecisionRT = as.numeric(DecisionRT),
         Confidence = as.numeric(Confidence),
         ResponseCorrect = as.numeric(ifelse(ResponseCorrect == "1",1,
                                             ifelse(ResponseCorrect == "1.0",1,
                                                    ifelse(ResponseCorrect == "True",1,
                                                           ifelse(ResponseCorrect == "0",0,
                                                                  ifelse(ResponseCorrect == "0.0",0,
                                                                         ifelse(ResponseCorrect == "False",0,NA))))))),
         Decision = ifelse(Decision == "More",1,ifelse(Decision == "Less",0,NA)))

df$participant_id = as.numeric(as.factor(df$participant_id))

df = df  %>% filter(DecisionRT < 8 & DecisionRT > 0.1)%>% 
  select(ResponseCorrect,participant_id,DecisionRT,Decision, Confidence, ConfidenceRT,Modality, Alpha)  %>% drop_na()

minRT = df %>% group_by(participant_id) %>% summarize(minRT = min(DecisionRT))

# idq = (1:length(unique(df$participant_id)))[1:50]


df1 = inner_join(df,minRT) %>% 
  filter(Modality == "Extero")

df1$participant_id = as.numeric(as.factor(df1$participant_id))



df1 %>% mutate(DecisionRT = log(DecisionRT), Confidence = Confidence/100) %>% 
  pivot_longer(cols = c("Decision","DecisionRT","Confidence")) %>% 
  # filter(abs(Alpha)<25) %>% 
  group_by(Alpha, name) %>% 
  summarize(mean = mean(value), se = sd(value)/sqrt(n())) %>% 
  ggplot(aes(x = Alpha, y = mean, col = name))+
  geom_pointrange(aes(ymin = mean-2*se,ymax = mean+2*se))+
  theme_minimal()
# facet_wrap(~participant_id)

df1 %>% mutate(DecisionRT = log(DecisionRT), Confidence = Confidence/100) %>% 
  filter(participant_id < 5) %>% 
  pivot_longer(cols = c("Decision","DecisionRT","Confidence")) %>% 
  # filter(abs(Alpha)<25) %>% 
  group_by(Alpha, name,participant_id) %>% 
  summarize(mean = mean(value), se = sd(value)/sqrt(n())) %>% 
  ggplot(aes(x = Alpha, y = mean, col = name))+
  geom_point()+
  theme_minimal()+
  facet_wrap(~participant_id)+
  geom_smooth()

df1 %>% mutate(DecisionRT = log(DecisionRT), Confidence = Confidence/100) %>% 
  pivot_longer(cols = c("Decision","DecisionRT","Confidence")) %>% 
  filter(abs(Alpha)<25) %>%
  group_by(Alpha, name) %>% 
  summarize(mean = mean(value), se = sd(value)/sqrt(n())) %>% 
  ggplot(aes(x = Alpha, y = mean, col = name))+
  geom_pointrange(aes(ymin = mean-2*se,ymax = mean+2*se))+
  theme_minimal()+
  facet_grid(name~1,scales = "free")+
  theme(legend.position = "top",
        text = element_text(size = 16),        # General text size
        axis.title = element_text(size = 18), # Axis titles
        axis.text = element_text(size = 14),  # Axis tick labels
        plot.title = element_text(size = 20)  # Plot title)
  )

remover = function(df){
  avg_acc = df %>% group_by(participant_id) %>% summarize(mean = mean(ResponseCorrect)) %>% mutate(se = NA, measure = "correct") 
  avg_rt = df %>% group_by(participant_id) %>% summarize(mean = mean(DecisionRT), se = sd(DecisionRT)/sqrt(n())) %>% mutate(measure = "rt")
  avg_conf = df %>% group_by(participant_id) %>% summarize(mean = mean(Confidence), se = sd(Confidence)/sqrt(n())) %>% mutate(measure = "confidence")
  avg_confrt = df %>% mutate(ConfidenceRT = as.numeric(ConfidenceRT)) %>% group_by(participant_id) %>% summarize(mean = mean(ConfidenceRT, na.rm = T), se = sd(ConfidenceRT)/sqrt(n())) %>% mutate(measure = "confidencert")
  
  
  outliers = rbind(avg_acc,avg_rt,avg_conf,avg_confrt) %>% group_by(measure) %>% summarize(meann = mean(mean), se = sd(mean, na.rm = T)) %>%
    mutate(low_int = meann - 2*se,high_int = meann + 2*se) %>% mutate(se = NULL) 
  
  dd = inner_join(rbind(avg_acc,avg_rt,avg_conf,avg_confrt),outliers, by = "measure") %>% 
    mutate(outlier = ifelse(mean < low_int | mean > high_int,1,0))
  
  
  inner_join(rbind(avg_acc,avg_rt,avg_conf,avg_confrt),outliers, by = "measure") %>% 
    mutate(outlier = ifelse(mean < low_int | mean > high_int,1,0)) %>% 
    ggplot(aes(x = participant_id, y = mean, ymin = mean-se, ymax = mean+se, col = as.factor(outlier)))+
    facet_wrap(~measure, scales = "free")+theme_minimal()+geom_pointrange()
  
  return(badids = dd %>% filter(measure == "correct" & outlier == 1) %>% .$participant_id)
}


removes = remover(df)

df %>% filter(participant_id %in% removes) %>% 
  group_by(Alpha,participant_id) %>% summarize(pfaster = mean(Decision), se = sd(Decision)/sqrt(n())) %>% 
  ggplot(aes(x = Alpha, y = pfaster,ymin = pfaster-2*se, ymax = pfaster+2*se))+geom_pointrange()+facet_wrap(~participant_id)



df1 %>% filter(!participant_id %in% removes) %>% 
  mutate(DecisionRT = log(DecisionRT), Confidence = Confidence/100) %>% 
  pivot_longer(cols = c("Decision","DecisionRT","Confidence")) %>% 
  # filter(abs(Alpha)<25) %>% 
  group_by(Alpha, name) %>% 
  summarize(mean = mean(value), se = sd(value)/sqrt(n())) %>% 
  ggplot(aes(x = Alpha, y = mean, col = name))+
  geom_pointrange(aes(ymin = mean-2*se,ymax = mean+2*se))+
  theme_minimal()
# facet_wrap(~participant_id)


df1 %>% filter(!participant_id %in% removes )%>% 
  mutate(DecisionRT = log(DecisionRT), Confidence = Confidence/100) %>% 
  pivot_longer(cols = c("Decision","DecisionRT","Confidence")) %>% 
  filter(abs(Alpha)<25) %>%
  group_by(Alpha, name) %>% 
  summarize(mean = mean(value), se = sd(value)/sqrt(n())) %>% 
  ggplot(aes(x = Alpha, y = mean, col = name))+
  geom_pointrange(aes(ymin = mean-2*se,ymax = mean+2*se))+
  theme_minimal()+
  facet_grid(name~1,scales = "free")+
  theme(legend.position = "top",
        text = element_text(size = 16),        # General text size
        axis.title = element_text(size = 18), # Axis titles
        axis.text = element_text(size = 14),  # Axis tick labels
        plot.title = element_text(size = 20)  # Plot title)
  )


# qq = qq %>% filter(participant_id < 5)

qq = df1 %>% filter(participant_id %in% c(10))


qq %>% mutate(DecisionRT = log(DecisionRT), Confidence = Confidence/100) %>% 
  # filter(participant_id < 5) %>% 
  pivot_longer(cols = c("Decision","DecisionRT","Confidence")) %>% 
  # filter(abs(Alpha)<25) %>% 
  group_by(Alpha, name,participant_id) %>% 
  summarize(mean = mean(value), se = sd(value)/sqrt(n())) %>% 
  ggplot(aes(x = Alpha, y = mean, col = name))+
  geom_point()+
  theme_minimal()+
  facet_wrap(~participant_id, scales = "free")+
  geom_smooth()

# 
# copula_discrete_cor <- cmdstan_model(here::here("Single_Subject_play","Psychometric","stanmodels","ss.stan"))
# 
# cor <- copula_discrete_cor$sample(
#   data = list(N = nrow(qq),
#               minRT = unique(qq %>% group_by(participant_id) %>% summarize(minRT = min(DecisionRT)) %>% .$minRT),
#               X = qq$Alpha,
#               ACC = qq$ResponseCorrect,
#               RT = (qq$DecisionRT),
#               Conf = qq$Confidence/100,
#               binom_y = qq$Decision),
#   refresh = 100,
#   init = 0,
#   iter_sampling = 500,
#   iter_warmup = 500,
#   adapt_delta = 0.95,
#   parallel_chains = 4)


s = 3


sample = function(s){
 
   qq = df1 %>% filter(participant_id %in% s)
  
  
  no_cop <- cmdstan_model(here::here("Single_Subject_play","Psychometric","stanmodels","ss_nocop.stan"))
  
  nocor <- no_cop$sample(
    data = list(N = nrow(qq),
                minRT = unique(qq %>% group_by(participant_id) %>% summarize(minRT = min(DecisionRT)) %>% .$minRT),
                X = qq$Alpha,
                ACC = qq$ResponseCorrect,
                RT = (qq$DecisionRT),
                Conf = qq$Confidence/100,
                binom_y = qq$Decision),
    refresh = 100,
    init = 0,
    iter_sampling = 500,
    iter_warmup = 500,
    adapt_delta = 0.95,
    parallel_chains = 4)
  
  df_nocor = nocor$summary(c("gm","rt_ndt","c0","c11")) %>% mutate(variable = c("threshold","slope","lapse","rt_int","rt_beta",
                                              "rt_sd","rt_stim","conf_int","conf_acc","conf_entro","conf_acc_entro",
                                              "conf_prec","rt_ndt","c0","c11")) %>% 
    mutate(div = sum(nocor$diagnostic_summary()$num_divergent),
           tree = sum(nocor$diagnostic_summary()$num_max_treedepth),
           model = "pure",
           subject = s)
  
  no_cop_fluc <- cmdstan_model(here::here("Single_Subject_play","Psychometric","stanmodels","ss_nocop_fluc.stan"))
  
  nocor_fluc <- no_cop_fluc$sample(
    data = list(N = nrow(qq),
                minRT = unique(qq %>% group_by(participant_id) %>% summarize(minRT = min(DecisionRT)) %>% .$minRT),
                X = qq$Alpha,
                ACC = qq$ResponseCorrect,
                RT = (qq$DecisionRT),
                Conf = qq$Confidence/100,
                binom_y = qq$Decision),
    refresh = 100,
    init = 0,
    iter_sampling = 500,
    iter_warmup = 500,
    adapt_delta = 0.95,
    parallel_chains = 4)
  
  
  df_nocor_fluc = nocor_fluc$summary(c("gm","rt_ndt","c0","c11","threshold")) %>% mutate(variable = c("gm_alpha","slope","lapse","rt_int","rt_beta",
                                                                                "rt_sd","rt_stim","conf_int","conf_acc","conf_entro","conf_acc_entro",
                                                                                "conf_prec","prec_alpha","drift","rt_ndt","c0","c11","threshold")) %>% 
    mutate(div = sum(nocor_fluc$diagnostic_summary()$num_divergent),
           tree = sum(nocor_fluc$diagnostic_summary()$num_max_treedepth),
           model = "fluc",
           subject = s)
  
  copula_discrete_cor <- cmdstan_model(here::here("Single_Subject_play","Psychometric","stanmodels","ss.stan"))

  cor <- copula_discrete_cor$sample(
    data = list(N = nrow(qq),
                minRT = unique(qq %>% group_by(participant_id) %>% summarize(minRT = min(DecisionRT)) %>% .$minRT),
                X = qq$Alpha,
                ACC = qq$ResponseCorrect,
                RT = (qq$DecisionRT),
                Conf = qq$Confidence/100,
                binom_y = qq$Decision),
    refresh = 100,
    init = 0,
    iter_sampling = 500,
    iter_warmup = 500,
    adapt_delta = 0.95,
    parallel_chains = 4)

  
  df_cop = cor$summary(c("gm","rt_ndt","c0","c11","rho_p_rt","rho_rt_conf","rho_p_conf")) %>% mutate(variable = c("alpha","slope","lapse","rt_int","rt_beta",
                                                                                                      "rt_sd","rt_stim","conf_int","conf_acc","conf_entro","conf_acc_entro",
                                                                                                      "conf_prec","rt_ndt","c0","c11","rho_p_rt","rho_rt_conf","rho_p_conf")) %>% 
    mutate(div = sum(cor$diagnostic_summary()$num_divergent),
           tree = sum(cor$diagnostic_summary()$num_max_treedepth),
           model = "cop",
           subject = s)
  
  
  
  
  fluc = nocor_fluc$loo()
  pure = nocor$loo()
  cop = cor$loo()
  
  loos = data.frame(loo::loo_compare(list(fluc = fluc, pure = pure,cop = cop))) %>% rownames_to_column("model") %>% 
    mutate(pareto_k = c(sum(fluc$diagnostics$pareto_k > 0.7),sum(pure$diagnostics$pareto_k > 0.7), sum(cop$diagnostics$pareto_k > 0.7)),
           subject = s)
  
  loglik_fluc = nocor_fluc$summary("log_lik") %>%   mutate(div = sum(nocor_fluc$diagnostic_summary()$num_divergent),
                                                           tree = sum(nocor_fluc$diagnostic_summary()$num_max_treedepth),
                                                           subject = s)
  
  loglik_nocor = nocor$summary("log_lik") %>%  mutate(div = sum(nocor$diagnostic_summary()$num_divergent),
                                                  tree = sum(nocor$diagnostic_summary()$num_max_treedepth),
                                                  subject = s)
  
  loglik_cop = cor$summary("log_lik")  %>% mutate(div = sum(cor$diagnostic_summary()$num_divergent),
                                              tree = sum(cor$diagnostic_summary()$num_max_treedepth),
                                              subject = s)
  
  
  return(list(loos,df_nocor_fluc,df_nocor,df_cop,loglik_fluc,loglik_nocor,loglik_cop))
  
}

fits = sample(1)

library(future)
library(future.apply)

library(furrr)

plan(multisession, workers = 10)  # adjust workers

# subjects vector
subjects <- unique(df1$participant_id)[1:20]

# run in parallel with progress bar
results <- future_map(subjects, function(s) {
  sample(s)
}, .progress = TRUE)


save.image("~/Hierarchical-Multivariate-Copula-Framework/Single_Subject_play/Psychometric/all_copula.RData")

loos_all <- do.call(rbind, lapply(results, `[[`, 1))
df_fluc_all <- do.call(rbind, lapply(results, `[[`, 2))
df_nocor_all <- do.call(rbind, lapply(results, `[[`, 3))
df_cop_all <- do.call(rbind, lapply(results, `[[`, 4))


df_fluc_all  %>% select(subject,div,tree) %>% distinct() %>% ggplot(aes(x = subject, y = div))+geom_point()
df_nocor_all %>% select(subject,div,tree) %>% distinct() %>% ggplot(aes(x = subject, y = div))+geom_point()
df_cop_all %>% select(subject,div,tree) %>% distinct() %>% ggplot(aes(x = subject, y = div))+geom_point()

df_fluc_all %>% filter(div != 0)%>% select(subject,div,tree) %>% distinct() %>% ggplot(aes(x = subject, y = div))+geom_point()
df_nocor_all%>% filter(div != 0) %>% select(subject,div,tree) %>% distinct() %>% ggplot(aes(x = subject, y = div))+geom_point()
df_cop_all%>% filter(div != 0) %>% select(subject,div,tree) %>% distinct() %>% ggplot(aes(x = subject, y = div))+geom_point()


loos_all %>% filter(elpd_diff != 0) %>% ggplot(aes(x = subject, y = elpd_diff, col = model))+geom_point()
loos_all  %>% ggplot(aes(x = subject, y = pareto_k))+geom_point()+facet_wrap(~ model)


df_fluc_all %>% ggplot(aes(x = subject,y = mean))+geom_point()+facet_wrap(~variable, scales = "free")
df_fluc_all %>% filter(div == 0) %>% ggplot(aes(x = subject,y = mean))+geom_point()+facet_wrap(~variable, scales = "free")

df_nocor_all %>% ggplot(aes(x = subject,y = mean))+geom_point()+facet_wrap(~variable, scales = "free")
df_nocor_all %>% filter(div == 0) %>% ggplot(aes(x = subject,y = mean))+geom_point()+facet_wrap(~variable, scales = "free")

rbind(df_fluc_all,df_nocor_all) %>%
  mutate(div_bin = ifelse(div == 0,"div","no_div")) %>% 
  filter(variable == "slope") %>% 
  select(mean,model,subject,div_bin) %>% pivot_wider(names_from = "model",values_from = "mean") %>% 
  ggplot(aes(x = fluc, y= pure, col = div_bin))+geom_point()+geom_abline()+
  geom_smooth(method = "lm")

m1= rbind(df_fluc_all,df_nocor_all) %>% filter(div == 0 & variable == "slope")%>% 
  select(mean,model,subject) %>% pivot_wider(names_from = "model",values_from = "mean") %>% 
#  filter(fluc > -3.5)  %>% 
  lm(fluc ~1+pure, data = .)

summary(m1)    

probs = function(lapse,beta,x,alpha){lapse + (1 - 2 * lapse) * (tanh(beta*(x-alpha)) / 2 + 0.5)}
x = seq(-40,40,by = 0.1)

plot(x, probs(brms::inv_logit_scaled(nocor_fluc$summary("gm[3]")%>% .$mean),exp(nocor_fluc$summary("gm[2]")%>% .$mean),x,(nocor_fluc$summary("gm[1]") %>% .$mean)))

plot(x, probs(brms::inv_logit_scaled(nocor$summary("gm[3]")%>% .$mean),exp(nocor$summary("gm[2]")%>% .$mean),x,(nocor$summary("gm[1]") %>% .$mean)))


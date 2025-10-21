

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

extreme_ids = df1 %>% filter(Confidence == 100 | Confidence == 0)
extreme_ids = unique(extreme_ids$participant_id)

qq = df1 %>% filter(!participant_id %in% extreme_ids)
qq$participant_id = as.numeric(as.factor(qq$participant_id))

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

qq = df1 %>% filter(participant_id < 200)


# qq = df1 %>% filter(!participant_id %in% removes)




qq$participant_id = as.numeric(as.factor(qq$participant_id))

t_p_s = qq %>% group_by(participant_id) %>% summarize(n = n())


ends <- cumsum(t_p_s$n)

# Calculate the start points
starts <- c(1, head(ends, -1) + 1)

# copula_discrete_cor <- cmdstan_model(here::here("corpulas","real_copula_hier_mixed.stan"))
# copula_discrete_cor <- cmdstan_model(here::here("corpulas","real_copula_hier_mixed_conf.stan"))
copula_discrete_cor <- cmdstan_model(here::here("Real data","Psychometric","stanmodels","test","all.stan"))

cor <- copula_discrete_cor$sample(
  data = list(N = nrow(qq),
              S = length(unique(qq$participant_id)),
              starts = starts,
              minRT = unique(qq %>% group_by(participant_id) %>% summarize(minRT = min(DecisionRT)) %>% .$minRT),
              ends = ends,
              X = qq$Alpha,
              ACC = qq$ResponseCorrect,
              S_id = qq$participant_id,
              RT = (qq$DecisionRT),
              Conf = qq$Confidence/100,
              binom_y = qq$Decision),
  refresh = 100,
  init = 0,
  iter_sampling = 500,
  iter_warmup = 500,
  adapt_delta = 0.95,
  parallel_chains = 4)

cor$save_object(here::here("VMP_100.rds"))

cor$summary(c("rho_p_rt"))
cor$summary(c("rho_p_conf"))
cor$summary(c("rho_rt_conf"))


entropy = function(p){
  return(-p * log(p) - (1-p) * log(1-p))
}


as_draws_df(cor$draws("gm")) %>% select(-contains(".")) %>% rename_with(~c("alpha","beta","lapse",
                                                                           "rt_int","rt_slope","rt_sd","rt_stim",
                                                                           "conf_int","conf_acc","conf_slope","conf_acc_slope","conf_prec")) %>% 
  mutate(x = list(seq(-20,20,1))) %>% 
  unnest() %>% 
  mutate(acc = list(c(0,1))) %>% 
  unnest() %>% 
  mutate(Decision = brms::inv_logit_scaled(lapse)/2 + (1-2*brms::inv_logit_scaled(lapse)/2) * (tanh(exp(beta)*(x-alpha)) / 2 + 0.5)) %>% 
  mutate(DecisionRT = exp(rt_int + rt_slope * entropy(Decision) + rt_stim * x) + 0.3) %>% 
  mutate(Confidence = brms::inv_logit_scaled(conf_int + conf_slope * entropy(Decision) + conf_acc * acc + conf_acc_slope * entropy(Decision) * acc)) %>% 
  pivot_longer(cols = c("Decision","DecisionRT","Confidence")) %>%
  group_by(x, name,acc) %>% 
  summarize(mean = mean(value),
            q5 = quantile(value,0.05),
            q95 = quantile(value,0.95)) %>% 
  ggplot(aes(x = x, y = mean, col = interaction(name,acc)))+
  geom_ribbon(aes(ymin = q5,ymax = q95), alpha = 0.1)+
  theme_minimal()+
  geom_pointrange(data = realdd, aes(x = Alpha, y = mean,ymin = q5,ymax = q95))+
  facet_grid(name~1,scales = "free")+
  theme(legend.position = "top",
        text = element_text(size = 16),        # General text size
        axis.title = element_text(size = 18), # Axis titles
        axis.text = element_text(size = 14),  # Axis tick labels
        plot.title = element_text(size = 20)  # Plot title)
  )

  






realdd = qq %>% 
  mutate(DecisionRT = (DecisionRT), Confidence = Confidence/100) %>% 
  pivot_longer(cols = c("Decision","DecisionRT","Confidence")) %>% 
  filter(abs(Alpha)<20) %>%
  group_by(Alpha, name,ResponseCorrect) %>% 
  summarize(mean = mean(value), se = sd(value)/sqrt(n())) %>% 
    mutate(q5 = mean-2*se, q95 = mean+2*se) %>% 
  rename(acc = ResponseCorrect) %>% 
  mutate(mean = ifelse(acc == 0 & name %in% c("Decision","DecisionRT"),NA,mean))%>% 
  mutate(q5 = ifelse(acc == 0 & name %in% c("Decision","DecisionRT"),NA,q5))%>% 
  mutate(q95 = ifelse(acc == 0 & name %in% c("Decision","DecisionRT"),NA,q95))




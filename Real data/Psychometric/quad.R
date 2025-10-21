

packages = c("brms","tidyverse","bayesplot","pracma","here", "patchwork","posterior","HDInterval","loo", "furrr","cmdstanr")

do.call(pacman::p_load, as.list(packages))


dat <- bind_rows(read_csv(here::here("Real data","Psychometric","data","quad","filtered_concurrent.csv"), col_types = cols()) %>% mutate(exp = "con"), 
                 read_csv(here::here("Real data","Psychometric","data","quad","filtered_delayed.csv"), col_types = cols()) %>% mutate(exp = "delay"))



dat %>% filter(trial_type == "decision") %>% 
  mutate(answer = recode(answer, 'short' = 0, 'long' = 1)) %>% 
  mutate(confidence = recode(confidence, 'low' = 0, 'high' = 1)) %>% 
  pivot_longer(cols = c("confidence","rt","answer")) %>% 
  group_by(name,bias_source,bias_direction,target_length,accuracy,exp) %>% 
  summarize(mean = mean(value,na.rm = T), se = sd(value,na.rm = T)/sqrt(n())) %>% 
  ggplot(aes(x = target_length, y = mean,ymin = mean-2*se, ymax = mean+2*se, col = interaction(accuracy,bias_direction)))+
  geom_pointrange()+facet_grid(name~bias_source+exp, scales = "free")+
  theme(legend.position = "top")


qq = dat %>% filter(trial_type == "decision" & bias_direction == "short",bias_source == "mullerlyer") %>% 
  mutate(answer = recode(answer, 'short' = 0, 'long' = 1)) %>% 
  mutate(confidence = recode(confidence, 'low' = 0, 'high' = 1)) %>% 
  filter(exp == "delay" & trial_type == "decision") %>% 
  rename(DecisionRT = rt, Confidence = confidence, Alpha = target_length, ResponseCorrect = accuracy, Decision = answer) %>% 
  mutate(Alpha = Alpha-400,
         ResponseCorrect = ifelse(ResponseCorrect == T, 1,0))

qq$participant_id = as.numeric(as.factor(qq$participant))

qq = qq %>% filter(participant_id < 100)

t_p_s = qq %>% group_by(participant_id) %>% summarize(n = n())


ends <- cumsum(t_p_s$n)

# Calculate the start points
starts <- c(1, head(ends, -1) + 1)



copula_discrete_cor <- cmdstan_model(here::here("Real data","Psychometric","stanmodels","test","all_bin.stan"))

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
              conf_y = qq$Confidence,
              binom_y = qq$Decision),
  refresh = 100,
  init = 0,
  iter_sampling = 500,
  iter_warmup = 500,
  adapt_delta = 0.95,
  parallel_chains = 4)






as_draws_df(cor$draws("gm")) %>% select(-contains(".")) %>% rename_with(~c("alpha","beta","lapse",
                                                                           "rt_int","rt_slope","rt_sd","rt_stim",
                                                                           "conf_int","conf_acc","conf_slope","conf_acc_slope")) %>% 
  mutate(x = list(seq(-30,30,10))) %>% 
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





realdd %>% ggplot()+geom_pointrange(data = realdd, aes(x = Alpha, y = mean,ymin = q5,ymax = q95, col = as.factor(acc)))+facet_wrap(~name)


realdd = qq %>% 
  mutate(DecisionRT = (DecisionRT), Confidence = Confidence) %>% 
  pivot_longer(cols = c("Decision","DecisionRT","Confidence")) %>% 
  # filter(abs(Alpha)<20) %>%
  group_by(Alpha, name,ResponseCorrect) %>% 
  summarize(mean = mean(value), se = sd(value)/sqrt(n())) %>% 
  mutate(q5 = mean-2*se, q95 = mean+2*se) %>% 
  rename(acc = ResponseCorrect) %>% 
  mutate(mean = ifelse(acc == 0 & name %in% c("Decision","DecisionRT"),NA,mean))%>% 
  mutate(q5 = ifelse(acc == 0 & name %in% c("Decision","DecisionRT"),NA,q5))%>% 
  mutate(q95 = ifelse(acc == 0 & name %in% c("Decision","DecisionRT"),NA,q95))



realdd_sub = qq %>% 
  mutate(DecisionRT = (DecisionRT), Confidence = Confidence) %>% 
  pivot_longer(cols = c("Decision","DecisionRT","Confidence")) %>% 
  # filter(abs(Alpha)<20) %>%
  group_by(Alpha, name,ResponseCorrect, participant_id) %>% 
  summarize(mean = mean(value), se = sd(value)/sqrt(n())) %>% 
  mutate(q5 = mean-2*se, q95 = mean+2*se) %>% 
  rename(acc = ResponseCorrect) %>% 
  mutate(mean = ifelse(acc == 0 & name %in% c("Decision","DecisionRT"),NA,mean))%>% 
  mutate(q5 = ifelse(acc == 0 & name %in% c("Decision","DecisionRT"),NA,q5))%>% 
  mutate(q95 = ifelse(acc == 0 & name %in% c("Decision","DecisionRT"),NA,q95))



names  = c("alpha","beta","lapse",
          "rt_int","rt_slope","rt_stim",
          "conf_int","conf_ACC","conf_entropy","conf_entropy_ACC")

as_draws_df(cor$draws(names)) %>% select(-contains(".")) %>% mutate(draw = 1:n()) %>% 
  pivot_longer(-draw) %>%
  mutate(
    participant_id = str_extract(name, "\\[(\\d+)\\]") %>% str_remove_all("\\[|\\]"),
    parameter = str_remove(name, "\\[\\d+\\]"),
    name = NULL
  ) %>% pivot_wider(names_from = "parameter") %>% 
  mutate(x = list(seq(-30,30,10))) %>% 
  unnest() %>% 
  mutate(acc = list(c(0,1))) %>% 
  unnest() %>% 
  mutate(Decision = (lapse) + (1-2*(lapse)) * (tanh((beta)*(x-alpha)) / 2 + 0.5)) %>% 
  mutate(DecisionRT = exp(rt_int + rt_slope * entropy(Decision) + rt_stim * x) + 0.3) %>% 
  mutate(Confidence = brms::inv_logit_scaled(conf_int + conf_entropy  * entropy(Decision) + conf_ACC * acc + conf_entropy_ACC * entropy(Decision) * acc)) %>% 
  pivot_longer(cols = c("Decision","DecisionRT","Confidence")) %>%
  group_by(x, name,acc,participant_id) %>% 
  summarize(mean = mean(value),
            q5 = quantile(value,0.05),
            q95 = quantile(value,0.95)) %>% 
  ggplot(aes(x = x, y = mean, col = interaction(name,acc)))+
  geom_ribbon(aes(ymin = q5,ymax = q95), alpha = 0.1)+
  theme_minimal()+
  geom_pointrange(data = realdd_sub, aes(x = Alpha, y = mean,ymin = q5,ymax = q95))+
  facet_grid(name~participant_id,scales = "free")+
  theme(legend.position = "top",
        text = element_text(size = 16),        # General text size
        axis.title = element_text(size = 18), # Axis titles
        axis.text = element_text(size = 14),  # Axis tick labels
        plot.title = element_text(size = 20)  # Plot title)
  )


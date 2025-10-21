font = "sans"
font_size = 26
font_size_small = 16
axis_width = 1.5
tick_width = 1.5

text = ggplot2::theme(text = ggplot2::element_text(family = font, size = font_size))
theme = theme_classic()

#patch theme
patchtheme = ggplot2::theme(plot.tag = ggplot2::element_text(size = (font_size+4),       
                                                             family = "sans",     
                                                             face = "bold",            
                                                             hjust = 0.5,              
                                                             vjust = 0.5),
                            text=ggplot2::element_text(family=font))





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

## convergence!

# trace and prior posterior updates: (DDM)

parameters = paste0("mu_",c("learningrate","e0","alpha","beta","delta"))

traces_rw_ddm = as_draws_df(rw_ddm_nolow$draws(parameters)) %>% pivot_longer(all_of(parameters), names_to = "variables") %>% 
  mutate(.chain = as.factor(.chain),
         variables = factor(variables, levels = parameters)  # <-- set factor levels here
         ) %>% 
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


pp_update_rw_ddm = rbind(priors_rw_ddm,posteriors_rw_ddm %>% select(variables,value,distribution)) %>% 
  mutate(
    variables = factor(variables, levels = parameters)  # <-- set factor levels here
    
  ) %>% 
  ggplot(aes(x = value, fill = distribution))+
  geom_histogram(col = "black", alpha = 0.5, position = "identity")+
  facet_wrap(~variables, nrow = 1, scales = "free")+
  theme_minimal()+
  theme(legend.position = "top")+
  scale_x_continuous("Value", breaks = scales::pretty_breaks(n = 3))+
  scale_y_continuous("Frequency", breaks = scales::pretty_breaks(n = 3))+
  theme(legend.position = "top")



convergence_rw_ddm = traces_rw_ddm / pp_update_rw_ddm
convergence_rw_ddm

ggsave(here::here("Figures","supplementary3_convergence_rw_ddm.tiff"),convergence_rw_ddm,width = 7,height = 5,dpi = 300, units = "in")



parameters = paste0("mu_",c("learningrate","e0","rt_int","rt_beta","rt_sd"))

traces_rw_rt = as_draws_df(rw_rt_nolow$draws(parameters)) %>% pivot_longer(all_of(parameters), names_to = "variables") %>% 
  mutate(.chain = as.factor(.chain),
         variables = factor(variables, levels = parameters)  # <-- set factor levels here
         ) %>% 
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


pp_update_rw_rt = rbind(priors_rw_rt,posteriors_rw_rt %>% select(variables,value,distribution)) %>% 
  mutate(
    variables = factor(variables, levels = parameters)  # <-- set factor levels here
  ) %>% 
  ggplot(aes(x = value, fill = distribution))+
  geom_histogram(col = "black", alpha = 0.5, position = "identity")+
  facet_wrap(~variables, nrow = 1, scales = "free")+
  theme_minimal()+
  theme(legend.position = "top")+
  scale_x_continuous("Value", breaks = scales::pretty_breaks(n = 3))+
  scale_y_continuous("Frequency", breaks = scales::pretty_breaks(n = 3))+
  theme(legend.position = "top")

convergence_rw_rt = traces_rw_rt / pp_update_rw_rt
convergence_rw_rt
ggsave(here::here("Figures","supplementary4_convergence_rw_rt.tiff"),convergence_rw_rt,width = 7,height = 5,dpi = 300, units = "in")



pp_update_rw_rt / pp_update_rw_ddm


save.image("~/Modeling-the-experiment/real data/reinforcement learning/plotting.RData")



posteriors_rw_nort = as_draws_df(rw_nort_nolow$draws(parameters[1:2])) %>% 
  pivot_longer(all_of(parameters[1:2]), names_to = "variables") %>% 
  mutate(distribution = "posterior")



RW_parameters = rbind(
  rbind(priors_rw_rt,posteriors_rw_rt %>% select(variables,value,distribution)) %>% 
    mutate(
      variables = factor(variables, levels = parameters)  # <-- set factor levels here
    ) %>% filter(distribution == "posterior" & variables %in% c("mu_learningrate","mu_e0")) %>% mutate(model = "CBM")
  ,
  rbind(priors_rw_ddm,posteriors_rw_ddm %>% select(variables,value,distribution)) %>% 
    mutate(
      variables = factor(variables, levels = parameters)  # <-- set factor levels here
    )%>% filter(distribution == "posterior"& variables %in% c("mu_learningrate","mu_e0")) %>% mutate(model = "DDM")
  ,
  posteriors_rw_nort%>% select(variables,value,distribution) %>% 
  mutate(
    variables = factor(variables, levels = parameters)  # <-- set factor levels here
  ) %>%  mutate(model = "No_RT")
) %>% 
  ggplot(aes(x = value, fill = model))+
  geom_histogram(col = "black", alpha = 0.35, position = "identity", show.legend = F)+
  facet_wrap(~variables, nrow = 1, scales = "free")+
  theme_minimal()+
  scale_fill_manual(values = c("darkgreen","red","darkblue"))+
  theme(legend.position = "top")+
  scale_x_continuous("Parameter value", breaks = scales::pretty_breaks(n = 3))+
  scale_y_continuous("Frequency", breaks = scales::pretty_breaks(n = 3))+
  theme(legend.position = "top")+
  theme_classic()+
  theme+text+
  theme(strip.background = element_blank(),
        legend.position = "top",
        legend.box = "horizontal",
        legend.direction="horizontal",
        legend.justification='center',
        legend.key.height = unit(0.5, 'cm'),
        legend.key.width = unit(0.5, 'cm'),
        plot.title = element_text(hjust = 0.5, size = font_size), #change legend key height
        legend.text = element_text(size=font_size_small),
        axis.text=element_text(size=font_size_small),
        axis.title=element_text(size=font_size),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line=element_line(size=axis_width),
        axis.ticks=element_line(size=axis_width),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))






combined = (psychometric_parameters | RW_parameters)+
  plot_layout(tag_level = 'new', guides = "collect", axis_titles = "collect")&
  theme(legend.position = "bottom")

ggsave(here::here("Figures","parameter_values.tiff"),combined, dpi = 600,
       height = 15,width = 40, units = "cm")



parameters = paste0(c("learningrate","e0","alpha","beta","delta","rt_ndt"))

df1 = as_draws_df(rw_ddm_nolow$draws(parameters)) %>% 
  select(-contains(".")) %>% 
  mutate(draw = 1:n()) %>% 
  pivot_longer(-draw, names_to = "variables") %>% 
  mutate(
    id = as.integer(str_extract(variables, "(?<=\\[)\\d+(?=\\])")),  # extract number inside brackets
    variables = str_remove(variables, "\\[\\d+\\]")                  # remove [number]
  ) %>% pivot_wider(values_from = "value", names_from = "variables")


parameters = paste0(c("learningrate","e0","rt_int","rt_beta","rt_sd","rt_ndt","rho"))


df2 = as_draws_df(rw_rt_nolow$draws(parameters)) %>% 
  select(-contains(".")) %>% 
  mutate(draw = 1:n()) %>% 
  pivot_longer(-draw, names_to = "variables") %>% 
  mutate(
    id = as.integer(str_extract(variables, "(?<=\\[)\\d+(?=\\])")),  # extract number inside brackets
    variables = str_remove(variables, "\\[\\d+\\]")                  # remove [number]
  ) %>% pivot_wider(values_from = "value", names_from = "variables")

joined <- left_join(df1, df2, by = c("draw", "id"), suffix = c("_m1", "_m2")) %>% filter(draw %in% 1:1000)

params_m1 <- c("learningrate_m1", "e0_m1", "alpha_m1", "beta_m1", "delta_m1", "rt_ndt_m1")
params_m2 <- c("learningrate_m2", "e0_m2", "rt_int_m2", "rt_beta_m2", "rt_sd_m2", "rt_ndt_m2", "rho_m2")


colnames(joined) = c("draw" ,        "id"  ,         "learningrate_m1", "e0_m1","alpha_m1","beta_m1","delta_m1","rt_ndt_m1",
                     "learningrate_m2", "e0_m2","rt_int_m2","rt_beta_m2",      "rt_sd_m2",        "rt_ndt_m2","rho_m2")



# Generate all param pairs
param_grid <- expand_grid(param1 = params_m1, param2 = params_m2)

# Function to get correlation per draw for one pair
get_cor_by_draw <- function(p1, p2) {
  joined %>%
    group_by(draw) %>%
    summarize(cor = cor(.data[[p1]], .data[[p2]]), .groups = "drop") %>%
    mutate(var1 = p1, var2 = p2)
}

# Map over all pairs
all_correlations <- purrr::map2_dfr(param_grid$param1, param_grid$param2, get_cor_by_draw)


library(HDInterval)

summary_cor <- all_correlations %>%
  group_by(var1, var2) %>%
  summarize(
    median = median(cor),
    hdi_low = hdi(cor, credMass = 0.95)[1],
    hdi_high = hdi(cor, credMass = 0.95)[2],
    .groups = "drop"
  )


# Format text label with median + HDI
summary_cor <- summary_cor %>%
  mutate(label = sprintf("%.2f\n[%.2f, %.2f]", median, hdi_low, hdi_high))


ordered_m2 <- c(
  "learningrate_m2", "e0_m2", "rt_int_m2", "rt_beta_m2",
  "rt_ndt_m2", "rt_sd_m2", "rho_m2"
)

ordered_m1 <- c(
  "learningrate_m1", "e0_m1", "alpha_m1", "delta_m1",
  "rt_ndt_m1", "beta_m1"
)

param_labels_m1 <- c(
  "learningrate_m1" = "Learning rate",
  "e0_m1" = "E0",
  "alpha_m1" = "Boundary",
  "beta_m1" = "Bias",
  "delta_m1" = "Drift",
  "rt_ndt_m1" = "NDT"
)

param_labels_m2 <- c(
  "learningrate_m2" = "Learning rate",
  "e0_m2" = "E0",
  "rt_int_m2" = "RT Int",
  "rt_beta_m2" = "RT Beta",
  "rt_sd_m2" = "RT SD",
  "rt_ndt_m2" = "NDT",
  "rho_m2" = "Rho"
)


# Plot matrix-style heatmap with labels
corplot = summary_cor %>% mutate(
  var1_label = param_labels_m1[var1],
  var2_label = param_labels_m2[var2],
  var1_label = factor(var1_label, levels = param_labels_m1[ordered_m1]),
  var2_label = factor(var2_label, levels = param_labels_m2[ordered_m2])
) %>% 
  ggplot(aes(x = var2_label, y = var1_label, fill = median)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 2) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       name = "Median\nCorr", breaks = c(0.75,0,-0.50), labels = c("0.75","0","-0.50")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        panel.grid = element_blank()) +
  labs(x = "CBM parameters", y = "DDM parameters")
corplot

ggsave(here::here("Figures","subject level","subject_corplot_RW.tiff"),corplot,width = 7,height = 5,dpi = 300, units = "in")


deltas = rw_ddm_nolow$summary("delta") %>% .$mean

rbind(
  data.frame(rw_nort_nolow$summary("learningrate")) %>% mutate(model = "nort"),
  data.frame(rw_ddm_nolow$summary("learningrate")) %>% mutate(model = "ddm")) %>% select(median,q5,q95,model) %>% 
  pivot_wider(names_from = "model", values_from = c("q5","q95","median")) %>% 
  unnest() %>% ggplot(aes(col = deltas))+geom_pointrange(aes(x = median_ddm, y = median_nort,ymax = q95_nort,ymin = q5_nort))+
  geom_abline()


rbind(
  data.frame(rw_rt_nolow$summary("learningrate")) %>% mutate(model = "nort"),
  data.frame(rw_nort_nolow$summary("learningrate")) %>% mutate(model = "ddm")) %>% select(median,q5,q95,model) %>% pivot_wider(names_from = "model", values_from = c("q5","q95","median")) %>% 
  unnest() %>% ggplot()+geom_pointrange(aes(x = median_ddm, y = median_nort,ymax = q95_nort,ymin = q5_nort))+geom_abline()

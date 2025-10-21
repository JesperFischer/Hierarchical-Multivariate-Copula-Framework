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
ddm_model <- cmdstanr::cmdstan_model(here::here("Real data","Psychometric","stanmodels","hierddm_x.stan"), force_recompile = T)



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
psy_ddm$save_object(here::here("Real data","Psychometric","psy_ddm1_x.rds"))


rt_loo_bin = psy_rt$loo("log_lik_bin")
rt_loo_full = psy_rt$loo("log_lik")

nort_loo = psy_nort$loo()

psy_loo = psy_ddm$loo()


full_loo = data.frame(loo::loo_compare(list(rt_loo_full,psy_loo)))

full_loo %>% filter(elpd_diff != 0) %>% mutate(ratio = elpd_diff/se_diff)

bin_loo = loo::loo_compare(list(rt_loo_bin,nort_loo))


# trace and prior posterior updates: (DDM)

parameters = paste0("mu_",c("threshold","slope","lapse","alpha","beta","delta"))

traces_psycho_ddm = as_draws_df(psy_ddm$draws(parameters)) %>% pivot_longer(all_of(parameters), names_to = "variables") %>% 
  mutate(.chain = as.factor(.chain),
         variables = factor(variables, levels = parameters)  # <-- set factor levels here
         ) %>% 
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

   
pp_update_psycho_ddm = rbind(priors_psycho_ddm,posteriors_psycho_ddm %>% 
                               select(variables,value,distribution)) %>% 
  mutate(
    variables = factor(variables, levels = parameters)  # <-- set factor levels here
  ) %>% 
  ggplot(aes(x = value, fill = distribution))+
  geom_histogram(col = "black", alpha = 0.5, position = "identity")+
  facet_wrap(~variables, nrow = 1, scales = "free")+
  theme_minimal()+
  scale_x_continuous("Value", breaks = scales::pretty_breaks(n = 3))+
  scale_y_continuous("Frequency", breaks = scales::pretty_breaks(n = 3))+
  theme(legend.position = "top")



convergence_psycho_ddm = traces_psycho_ddm / pp_update_psycho_ddm
convergence_psycho_ddm


ggsave(here::here("Figures","supplementary1_convergence_psycho_ddm.tiff"),convergence_psycho_ddm,width = 7,height = 5,dpi = 300, units = "in")



# response


parameters = paste0("mu_",c("threshold","slope","lapse","rt_int","rt_beta","rt_sd"))

traces_psycho_rt = as_draws_df(psy_rt$draws(parameters)) %>% pivot_longer(all_of(parameters), names_to = "variables") %>% 
  mutate(.chain = as.factor(.chain),
         variables = factor(variables, levels = parameters)  # <-- set factor levels here
         ) %>% 
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


pp_update_psycho_rt = rbind(priors_psycho_rt,posteriors_psycho_rt %>% select(variables,value,distribution)) %>% 
  mutate(
    variables = factor(variables, levels = parameters)  # <-- set factor levels here
  ) %>% 
  ggplot(aes(x = value, fill = distribution))+
  geom_histogram(col = "black", alpha = 0.5, position = "identity")+
  facet_wrap(~variables, nrow = 1, scales = "free")+
  theme_minimal()+
  scale_x_continuous("Value", breaks = scales::pretty_breaks(n = 3))+
  scale_y_continuous("Frequency", breaks = scales::pretty_breaks(n = 3))+
  theme(legend.position = "top")



convergence_psycho_rt = traces_psycho_rt / pp_update_psycho_rt
convergence_psycho_rt

ggsave(here::here("Figures","supplementary2_convergence_psycho_rt.tiff"),convergence_psycho_rt,width = 7,height = 5,dpi = 300, units = "in")


posteriors_psy_nort = as_draws_df(psy_nort$draws(parameters[1:3])) %>% 
  pivot_longer(all_of(parameters[1:3]), names_to = "variables") %>% 
  mutate(distribution = "posterior")





psychometric_parameters = rbind(
  rbind(priors_psycho_rt,posteriors_psycho_rt %>% select(variables,value,distribution)) %>% 
    mutate(
      variables = factor(variables, levels = parameters)  # <-- set factor levels here
    ) %>% filter(distribution == "posterior" & variables %in% c("mu_threshold","mu_slope","mu_lapse")) %>% mutate(model = "CBM")
  ,
  pp_update_rw_ddm = rbind(priors_psycho_ddm,posteriors_psycho_ddm %>% select(variables,value,distribution)) %>% 
    mutate(
      variables = factor(variables, levels = parameters)  # <-- set factor levels here
    )%>% filter(distribution == "posterior"& variables %in% c("mu_threshold","mu_slope","mu_lapse")) %>% mutate(model = "DDM")
  ,
  posteriors_psy_nort%>% select(variables,value,distribution) %>% 
  mutate(model = "no_RT")) %>%
  filter(grepl("mu_",variables)) %>% 
  mutate(variables = str_replace(variables, "^mu_", "μ["),
         variables = paste0(variables, "]"),
         variables = str_replace(variables, "rt", "RT")) %>% 
  ggplot(aes(x = value, fill = model))+
  geom_histogram(col = "black", alpha = 0.35, position = "identity")+
  facet_wrap(~variables, nrow = 1, scales = "free")+
  theme_minimal()+
  theme(legend.position = "top")+
  scale_fill_manual(values = c("darkgreen","red","darkblue"))+
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








parameters = paste0(c("threshold","slope","lapse","rt_int","rt_beta","rt_sd","rt_ndt","rho"))

df1 = as_draws_df(psy_rt$draws(parameters)) %>% 
  select(-contains(".")) %>% 
  mutate(draw = 1:n()) %>% 
  pivot_longer(-draw, names_to = "variables") %>% 
  mutate(
    id = as.integer(str_extract(variables, "(?<=\\[)\\d+(?=\\])")),  # extract number inside brackets
    variables = str_remove(variables, "\\[\\d+\\]")                  # remove [number]
  ) %>% pivot_wider(values_from = "value", names_from = "variables")


parameters = paste0(c("threshold","slope","lapse","alpha","beta","delta","rt_ndt"))

df2 = as_draws_df(psy_ddm$draws(parameters)) %>% 
  select(-contains(".")) %>% 
  mutate(draw = 1:n()) %>% 
  pivot_longer(-draw, names_to = "variables") %>% 
  mutate(
    id = as.integer(str_extract(variables, "(?<=\\[)\\d+(?=\\])")),  # extract number inside brackets
    variables = str_remove(variables, "\\[\\d+\\]")                  # remove [number]
  ) %>% pivot_wider(values_from = "value", names_from = "variables")

joined <- left_join(df1, df2, by = c("draw", "id"), suffix = c("_m1", "_m2")) %>% filter(draw %in% 1:1000)

params_m1 <- c("threshold_m1", "slope_m1", "lapse_m1", "alpha_m1", "beta_m1", "delta_m1", "rt_ndt_m1")
params_m2 <- c("threshold_m2", "slope_m2", "lapse_m2", "rt_int_m2", "rt_beta_m2", "rt_sd_m2", "rt_ndt_m2", "rho_m2")


colnames(joined) = c("draw" ,        "id"  ,         "threshold_m1", "slope_m1"   ,  "lapse_m1" ,    "rt_int_m1",
                     "rt_beta_m1",      "rt_sd_m1",        "rt_ndt_m1","rho_m1",
                     "threshold_m2", "slope_m2"   ,  "lapse_m2",     "alpha_m2" ,       "beta_m2" ,        "delta_m2" ,       "rt_ndt_m2")


joined %>% group_by(draw) %>% 
  summarize(cor = cor.test(threshold_m1,threshold_m2)$estimate[[1]]) %>% 
  ggplot(aes(x = cor))+geom_histogram(col = "black")



params_m1 <- c("threshold_m1", "slope_m1", "lapse_m1", "rt_int_m1", "rt_beta_m1", "rt_sd_m1", "rt_ndt_m1","rho_m1")
params_m2 <- c("threshold_m2", "slope_m2", "lapse_m2", "alpha_m2", "beta_m2", "delta_m2", "rt_ndt_m2")

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



ordered_m1 <- c("threshold_m1", "slope_m1", "lapse_m1", "rt_int_m1", "rt_beta_m1", "rt_sd_m1", "rt_ndt_m1", "rho_m1")
ordered_m2 <- c("threshold_m2", "slope_m2", "lapse_m2", "alpha_m2", "beta_m2", "delta_m2", "rt_ndt_m2")

# Labels
param_labels_m1 <- c(
  "threshold_m1" = "Threshold",
  "slope_m1"     = "Slope",
  "lapse_m1"     = "Lapse",
  "rt_int_m1"    = "RT Intercept",
  "rt_beta_m1"   = "RT Slope",
  "rt_sd_m1"     = "RT SD",
  "rt_ndt_m1"    = "NDT",
  "rho_m1"       = "ρ"
)

param_labels_m2 <- c(
  "threshold_m2" = "Threshold",
  "slope_m2"     = "Slope",
  "lapse_m2"     = "Lapse",
  "alpha_m2"     = "Boundary",
  "beta_m2"      = "Bias",
  "delta_m2"     = "Drift rate",
  "rt_ndt_m2"    = "NDT"
)


# Plot matrix-style heatmap with labels
corplot = summary_cor %>% mutate(
  var1_label = param_labels_m1[var1],
  var2_label = param_labels_m2[var2],
  var1_label = factor(var1_label, levels = param_labels_m1[ordered_m1]),
  var2_label = factor(var2_label, levels = param_labels_m2[ordered_m2])
) %>% 
  ggplot(aes(x = var1_label, y = var2_label, fill = median)) +
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


ggsave(here::here("Figures","subject level","subject_corplot_PSY.tiff"),corplot,width = 7,height = 5,dpi = 300, units = "in")


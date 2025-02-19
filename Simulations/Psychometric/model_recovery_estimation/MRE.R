packages = c("brms","tidyverse","bayesplot","pracma","here",
             "patchwork","posterior","HDInterval","loo", "furrr", "SBC","future")

do.call(pacman::p_load, as.list(packages))

source(here::here("Simulations","Psychometric","model_recovery_estimation","scripts_mre.R"))


## with ddm as simulated  model:

cores = 20
n_sim = 1000
plan(multisession, workers = cores)

possfit_model = possibly(.f = fitter_ddm, otherwise = "Error")

results = future_map(rep(1000,n_sim), ~possfit_model(.x), .options = furrr_options(seed = T), .progress = T)

save.image(here::here("Simulations","Psychometric","model_recovery_estimation","results","1000_fitter_ddm_1000_samples1.RData"))



## with response time model as simulated  model:




cores = 20
n_sim = 1000
plan(multisession, workers = cores)

possfit_model = possibly(.f = fitter_cop_rt, otherwise = "Error")

results = future_map(rep(1000,n_sim), ~possfit_model(.x), .options = furrr_options(seed = T), .progress = T)

save.image(here::here("Simulations","Psychometric","model_recovery_estimation","results","1000_fitter_cop_rt_1000_samples1.RData"))


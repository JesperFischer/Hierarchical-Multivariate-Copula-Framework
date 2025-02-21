packages = c("brms","tidyverse","bayesplot","pracma","here",
             "patchwork","posterior","HDInterval","loo", "furrr", "SBC","future")

do.call(pacman::p_load, as.list(packages))


source(here::here("Simulations","Learning","model_recovery_estimation","sim_learn.R"))

## ddm

cores = 20
n_sim = 1000
plan(multisession, workers = cores)

qq = fitter_ddm(500)

possfit_model = possibly(.f = fitter_ddm, otherwise = "Error")

results <- future_map(rep(1000,n_sim), ~possfit_model(.x), .options = furrr_options(seed = TRUE), .progress = T)

save.image(here::here("Simulations","Learning","model_recovery_estimation","results","1000_rw_fitter_ddm_1000_samples1.RData"))

#### fitter rts

cores = 20
n_sim = 600
plan(multisession, workers = cores)

qq = fitter_cop_rt(1000)

possfit_model = possibly(.f = fitter_cop_rt, otherwise = "Error")

results <- future_map(rep(1000,n_sim), ~possfit_model(.x), .options = furrr_options(seed = TRUE), .progress = T)

save.image(here::here("Simulations","Learning","model_recovery_estimation","results","600_rw_fitter_rt_1000_samples2.RData"))



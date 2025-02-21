font = "sans"
font_size = 12
font_size_small = 8 
axis_width = 1 
tick_width = 1

text = ggplot2::theme(text = ggplot2::element_text(family = font, size = font_size))
theme = theme_classic()

#patch theme
patchtheme = ggplot2::theme(plot.tag = ggplot2::element_text(size = (font_size+4),       
                                                             family = "sans",     
                                                             face = "bold",            
                                                             hjust = 0.5,              
                                                             vjust = 0.5),
                            text=ggplot2::element_text(family=font))





plot_figure1 = function(rerun = F, seed = 1997){
  
  
  if(rerun){
    
    set.seed(seed)
    source(here::here("Simulations","Learning","SBC","rw_scripts.R"))
    
    cmdstan_model_rw <- cmdstanr::cmdstan_model(here::here("Simulations","Learning","SBC","stanmodels","plot_rw.stan"))
    
    sim_data_rw = generator_single(60,10)
    
    sim_data_data_rw = sim_data_rw[[2]]
    
    simulated_parameters_rw = sim_data_rw[[1]]
    
    
    fit_rw = cmdstan_model_rw$sample(data = sim_data_data_rw,
                               iter_warmup = 1000,
                               iter_sampling = 1000,
                               chains = 4,
                               seed = seed,
                               refresh = 100,
                               max_treedepth = 10,
                               parallel_chains = 4,
                               adapt_delta = 0.97)
    
    sim_data = data.frame(sim_data_data_rw)
    
    # plot it:
    sim_datax_rw = sim_data %>% 
      mutate(mu_rts = RT,
             expect = binom_y) %>% 
      group_by(trial) %>% 
      summarize(se_expect = (mean(expect) * (1-mean(expect))) /sqrt(n()),
                expect = mean(expect),
                se_murts = sd(mu_rts)/sqrt(n()),
                murts = median(mu_rts),
      ) %>% ungroup() %>% mutate(trialx = 1:n(),sim_id = NA)
    
    draw_idx_rw = sample(1:4000,100)
    
    source("~/Hierarchical-Multivariate-Copula-Framework/Plots/utility.R")
    # subject level 
    big_df_rw = data.frame()
    
    for(s in 1:10){
      print(s)
      us = sim_data %>% filter(S == S_id) %>% .$x
      
      minRT = min(sim_data %>% filter(S == S_id) %>% .$RT)
      
      parameters = paste0(c("learningrate","e0","rt_int","rt_beta","rt_ndt","rt_sd","rho"), "[",s,"]")
      
      dfq = as_draws_df(fit_rw$draws(parameters))%>% select(-contains(".")) %>% mutate(draw = 1:n()) %>% 
        filter(draw %in% draw_idx_rw) %>% 
        rename_with(~c("learningrate","e0","rt_int","rt_beta","rt_ndt","rt_sd","rho","draw")) %>% 
        group_by(draw) %>% 
        summarize(list(generate_trialwise_data_rt_copula(u = us, learningrate = learningrate,e0 = e0,participant = s,
                                                         rt_int = rt_int,rt_beta = rt_beta,rt_sd = rt_sd,rt_ndt = rt_ndt,rho = rho))) %>% unnest()

      
      big_df_rw = rbind(big_df_rw,dfq)
      
    }
    
  
    posterior_draws_rw = big_df_rw %>%  ungroup() %>%  mutate(rts = (rts),bin_resp = resp) %>% drop_na() %>%
      group_by(trial,draw) %>% 
      summarize(resp = mean(expectation), rt = median(rts)) %>% rename(trialx = trial)
    
    posterior_means_rw = posterior_draws_rw %>% 
      group_by(trialx) %>% 
      summarize(resps = mean(resp), rts = median(rt))
    
    
    
    # plots 50 priors
    priors_rw <- do.call(rbind, lapply(1:100, function(i) data.frame(generator_single(60, 10)[[2]])))
    
    priors_rw = priors_rw %>% mutate(sim_id = rep(1:100,each = 600))
    #subject level
    priors_rw %>% ggplot(aes(x = trial, y = expect, group = sim_id))+geom_line()+facet_wrap(~S_id,scales = "free", nrow = 2)
    priors_rw %>% ggplot(aes(x = trial, y = RT, group = sim_id))+geom_line()+facet_wrap(~S_id,scales = "free", nrow = 2)
    
    
    #group level:
    expect_rw = priors_rw %>% group_by(sim_id,trial) %>% summarize(expect = mean(expect)) %>% 
      mutate(paradigm = "Learning Paradigm") %>% 
      ggplot(aes(x = trial, y = expect, group = sim_id))+
      geom_line(aes(col = "Prior"),alpha = 0.5)+
      geom_line(data = posterior_draws_rw, aes(x = trialx, y = resp, group = draw,col = "Posterior"),,alpha = 0.3)+
      geom_pointrange(data = sim_datax_rw, aes(x = trialx, y = expect, ymin = expect-se_expect, ymax = expect+se_expect, col = "Data"),alpha = 0.5)+
      theme_classic()+
      scale_y_continuous(labels = c("0.0","0.5","1.0"), breaks = c(0,0.5,1))+
      scale_color_manual(name = "Model Predictions & Data",values = c("black","orange","lightblue3"))+
      coord_cartesian(ylim = c(0,1))+
      xlab(" ")+
      ylab("P(Response = 1)")+
      # facet_wrap(~paradigm)+
      ggtitle("Learning Paradigm")+
      theme+text+
      theme(strip.background = element_blank(),
            strip.text = element_text(size=font_size),
            legend.position = "center",
            legend.box = "horizontal",
            legend.direction="horizontal",
            legend.justification='center',
            legend.key.height = unit(0.5, 'cm'),
            legend.key.width = unit(0.5, 'cm'),
            legend.spacing.x = unit(0.35, "cm"),
            plot.title = element_text(hjust = 0.5, size = font_size), #change legend key height
            legend.title = element_blank(), #change legend title font size
            legend.text = element_text(size=font_size_small),
            axis.text=element_text(size=font_size_small),
            axis.title=element_text(size=font_size),
            axis.line=element_line(size=axis_width),
            axis.ticks=element_line(size=axis_width),
            axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))+
      ylab("")

    
    expect_rw
    
    rts_rw = priors_rw %>% rename(trialx = trial) %>% group_by(sim_id,trialx) %>% summarize(RT = median(exp(mu_rts))) %>% 
      ggplot(aes(x = trialx, y = RT, group = sim_id))+
      geom_line(aes(col = "Prior"),alpha = 0.5, show.legend = F)+
      geom_line(data = qq_summar_rt, aes(x = trialx, y = rt, group = draw,col = "Posterior"),alpha = 0.3, show.legend = F)+
      geom_pointrange(data = sim_datax_rw, aes(x = trialx, y = murts, ymin = murts-se_murts, ymax = murts+se_murts, col = "Data"),alpha = 0.5, show.legend = F)+
      theme_classic()+
      scale_color_manual(name = "Model Predictions & Data",values = c("black","orange","lightblue3"))+
      scale_y_continuous(limits = c(0.2,2)," ", breaks = scales::pretty_breaks(n = 3))+
      scale_x_continuous("Trial", breaks = scales::pretty_breaks(n = 3))+
      theme+text+
      theme(strip.background = element_blank(),
            legend.position = "center",
            strip.text = element_text(size=font_size),
            legend.box = "horizontal",
            legend.direction="horizontal",
            legend.justification='center',
            legend.key.height = unit(0.5, 'cm'),
            legend.key.width = unit(0.5, 'cm'),
            legend.spacing.x = unit(0.35, "cm"),
            plot.title = element_text(hjust = 0.5, size = font_size), #change legend key height
            legend.title = element_blank(), #change legend title font size
            legend.text = element_text(size=font_size_small),
            axis.text=element_text(size=font_size_small),
            axis.title=element_text(size=font_size),
            axis.line=element_line(size=axis_width),
            axis.ticks=element_line(size=axis_width),
            axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))
    
    rts_rw
    
    
    
    # psychometric paradigm:
    
    
    source(here::here("Simulations","Psychometric","SBC","scripts.R"))
    set.seed(seed)
    cmdstan_model <- cmdstanr::cmdstan_model(here::here("Simulations","Psychometric","SBC","stanmodels","plot.stan"))
    
    sim_data = generator_single(60,10)
    
    sim_data_data = sim_data[[2]]
    
    
    simulated_parameters = sim_data[[1]]
    
    
    fit = cmdstan_model$sample(data = sim_data_data,
                               iter_warmup = 1000,
                               iter_sampling = 1000,
                               chains = 4,
                               seed = seed,
                               refresh = 100,
                               max_treedepth = 10,
                               parallel_chains = 4,
                               adapt_delta = 0.97)
    
    # plot it:
    
    sim_datax = data.frame(sim_data_data) %>% 
      mutate(sim_id = NA,
             mu_rts = (RT)) %>% 
      group_by(x,sim_id) %>% 
      summarize(se_expect = sd(expect)/sqrt(n()),
                expect = mean(expect),
                se_murts = sd(mu_rts)/sqrt(n()),
                murts = median(mu_rts),
      ) %>% ungroup() %>% mutate(trialx = 1:n())
    
  
    draw_idx_psy = sample(1:4000,100)
    
    source("~/Hierarchical-Multivariate-Copula-Framework/Plots/utility.R")
    # subject level 
    big_df_psy = data.frame()
    
    for(s in 1:10){
      print(s)
      xs = seq(-40,40, length.out = 60)
      
      minRT = min(sim_data %>% filter(S == S_id) %>% .$RT)
      
      parameters = paste0(c("threshold","slope","lapse","rt_int","rt_beta","rt_ndt","rt_sd","rho"), "[",s,"]")
      
      dfq = as_draws_df(fit$draws(parameters))%>% select(-contains(".")) %>% mutate(draw = 1:n()) %>% 
        filter(draw %in% draw_idx_psy) %>% 
        rename_with(~c("threshold","slope","lapse","rt_int","rt_beta","rt_ndt","rt_sd","rho","draw")) %>% 
        group_by(draw) %>% 
        summarize(list(generate_trial_shannon_entropy_rt(x = xs, threshold = threshold,slope = slope, lapse = lapse,participant = s,
                                                         rt_int = rt_int,rt_beta = rt_beta,rt_ndt = rt_ndt,rt_sd = rt_sd,rho = rho))) %>% unnest()
      
      
      big_df_psy = rbind(big_df_psy,dfq)
      
    }
    
    
    posterior_draws_psy = big_df_psy %>%  ungroup() %>%  mutate(rts = (rts),bin_resp = resp) %>% drop_na() %>%
      group_by(x,draw) %>% 
      summarize(resp = mean(expectation), rt = median(rts))
    
    posterior_means_rw = posterior_draws_psy %>% 
      group_by(x) %>% 
      summarize(resps = mean(resp), rts = median(rt))
    
    
    
    
    
    # plots 50 priors
    priors <- do.call(rbind, lapply(1:100, function(i) data.frame(generator_single(60, 10)[[2]])))
    
    priors = priors %>% mutate(sim_id = rep(1:100,each = 600))
    #subject level
    priors %>% ggplot(aes(x = x, y = expect, group = sim_id))+geom_line()+facet_wrap(~S_id,scales = "free", nrow = 2)
    priors %>% ggplot(aes(x = x, y = RT, group = sim_id))+geom_line()+facet_wrap(~S_id,scales = "free", nrow = 2)
    
    
    #group level:
    expect_pyscho = priors %>% group_by(sim_id,x) %>% summarize(expect = mean(expect)) %>%
      mutate(paradigm = "Psychophysical Paradigm") %>% 
      ggplot(aes(x = x, y = expect, group = sim_id))+
      geom_line(aes(col = "Prior"),alpha = 0.5, show.legend = F)+
      geom_line(data = posterior_draws_psy, aes(x = x, y = resp, group = draw,col = "Posterior"),alpha = 0.3, show.legend = F)+
      geom_pointrange(data = sim_datax, aes(x = x, y = expect, ymin = expect-se_expect, ymax = expect+se_expect, col = "Data"),alpha = 0.5, show.legend = F)+
      theme_classic()+
      scale_y_continuous("P(Response = 1)", breaks = scales::pretty_breaks(n = 3))+
      scale_color_manual(name = "Model Predictions & Data",values = c("black","orange","lightblue3"))+
      coord_cartesian(ylim = c(0,1))+
      scale_x_continuous(" ", breaks = scales::pretty_breaks(n = 4))+
      ggtitle("Psychophysical Paradigm")+
      # facet_wrap(~paradigm)+
      theme+text+
      theme(strip.background = element_blank(),
            strip.text = element_text(size=font_size),
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
            axis.line=element_line(size=axis_width),
            axis.ticks=element_line(size=axis_width),
            axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))
    
    
    expect_pyscho
    
  
    rts_pyscho = priors %>% group_by(sim_id,x) %>% summarize(RT = median(exp(mu_rts))) %>% 
      ggplot(aes(x = x, y = RT, group = sim_id))+
      geom_line(aes(col = "Prior"),alpha = 0.5, show.legend = F)+
      geom_line(data = posterior_draws_psy, aes(x = x, y = rt, group = draw,col = "Posterior"),alpha = 0.3, show.legend = F)+
      geom_pointrange(data = sim_datax, aes(x = x, y = murts, ymin = murts-se_murts, ymax = murts+se_murts, col = "Data"),alpha = 0.5, show.legend = F)+
      theme_classic()+
      xlab("Stimulus Intensity")+
      scale_color_manual(name = "Model Predictions & Data", 
                         values = c("Prior" = "lightblue3", 
                                    "Posterior" = "orange", 
                                    "Data" = "black"))+
      scale_y_continuous(limits = c(0.2,2),"Response time (s)", breaks = scales::pretty_breaks(n = 3))+
      scale_x_continuous("Stimulus Intensity", breaks = scales::pretty_breaks(n = 4))+
      theme+text+
      theme(strip.background = element_blank(),
            strip.text = element_text(size=font_size),
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
            axis.line=element_line(size=axis_width),
            axis.ticks=element_line(size=axis_width),
            axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))
    
    
    rts_pyscho

    
    
    
    
  
    plot1_x = (expect_pyscho/rts_pyscho)
    
    plot1_y = (expect_rw/rts_rw)
    
    
  
    
    
    combined_plot = (plot1_x | plot1_y)+plot_layout(tag_level = 'new', guides = "collect")&
      theme(legend.position = "top",
            legend.direction="horizontal",
            legend.text = element_text(size = font_size),  # Increase legend text size
            legend.title = element_blank(),  # Increase title size if needed
            legend.key.size = unit(1, "cm"))& 
      guides(
        color = guide_legend(override.aes = list(
          linewidth = c(0.5, 3, 3)
        ))
      )
      
    combined_plot
    
    ggsave(here::here("Figures","plot_1.tiff"),combined_plot, dpi = 600,
           height = 6,width = 12)
    
    return(combined_plot)
  }else{
  
  load(here::here("Plots","workspace","figure_1.RData"))
  
  combined_plot
  
  # ggsave(here::here("Plots","plot_1.tiff"),combined_plot, dpi = 600,
  #        height = 5,width = 8)
  # 
  
  ggsave(here::here("Figures","plot_1.tiff"),combined_plot, dpi = 600,
         height = 6,width = 12)
  
  
  return(combined_plot)  
  }

}

plot_figure2 = function(){
  
  #####
  ### learning paradigm:
  
  load(here::here("Simulations","Learning","SBC","results","2000.RData"))
  
  results_learning = results
  
  sim_ids_to_keep <- results_learning$backend_diagnostics %>%
    dplyr::filter(n_divergent == 0 & n_max_treedepth == 0) %>%
    dplyr::pull(sim_id)
  
  
  sim_ids_to_exclude <- results_learning$backend_diagnostics %>%
    dplyr::filter(n_divergent != 0 | n_max_treedepth != 0) %>%
    dplyr::pull(sim_id)
  
  
  results_subset <- results_learning[sim_ids_to_keep]
  results_excluded <- results_learning[sim_ids_to_exclude]
  
  
  sim_ids_to_keep <- results_subset$default_diagnostics %>%
    dplyr::filter(min_ess_bulk  > 400)%>%
    dplyr::filter(min_ess_tail  > 400) %>%
    dplyr::pull(sim_id)
  
  results_subset <- results_subset[sim_ids_to_keep]
  
  print(paste0("Excluded based on divergent or treedepth = ",length(sim_ids_to_exclude)))
  
  
  grouplevel = results_excluded$stats %>% 
    filter(grepl("^mu_", results_excluded$stats$variable) | grepl("^tau_", results_excluded$stats$variable)) %>% 
    mutate(variable = gsub("\\[\\d+\\]", "", variable))%>% select(sim_id, simulated_value,variable) %>% distinct() %>% 
    mutate(sim_id = as.factor(sim_id), excluded = T) 
  
  grouplevel_inc = results_subset$stats %>% 
    filter(grepl("^mu_", results_subset$stats$variable) | grepl("^tau_", results_subset$stats$variable)) %>% 
    mutate(variable = gsub("\\[\\d+\\]", "", variable))%>% select(sim_id, simulated_value,variable) %>% distinct() %>% 
    mutate(sim_id = as.factor(sim_id), excluded = F)
  
  #rbind(grouplevel,grouplevel_inc)%>% ggplot(aes(x = simulated_value, fill = excluded))+geom_histogram(col = "black")+facet_wrap(~variable, scales = "free")
  
  
  subj_level = results_excluded$stats %>% 
    filter(!grepl("^mu_", results_excluded$stats$variable)& !grepl("^tau_", results_excluded$stats$variable)) %>% 
    mutate(variable = gsub("\\[\\d+\\]", "", variable))%>% select(sim_id, simulated_value,variable) %>% distinct() %>% 
    mutate(excluded = T)
  
  subj_level_inc = results_subset$stats %>% 
    filter(!grepl("^mu_", results_subset$stats$variable)& !grepl("^tau_", results_subset$stats$variable)) %>% 
    mutate(variable = gsub("\\[\\d+\\]", "", variable))%>% select(sim_id, simulated_value,variable) %>% distinct()  %>% 
    mutate(excluded = F)
  
  
  # rbind(subj_level,subj_level_inc) %>% ggplot(aes(x = simulated_value, fill = excluded))+
  #   geom_histogram(col = "black")+facet_wrap(~variable, scales = "free")
  
  
  
  
  
  # Extract subject-level parameters from results2$stats
  subject_level_params <- results_subset$stats %>% 
    filter(!grepl("^mu_", results_subset$stats$variable)& !grepl("^tau_", results_subset$stats$variable)) %>% 
    mutate(variable = gsub("\\[\\d+\\]", "", variable)) %>% 
    # filter(!sim_id %in% sims$sim_id) %>% 
    mutate()
  
  
  
  # Extract subject-level parameters from results2$stats
  group_level_params_learning <- results_subset$stats %>% 
    filter(grepl("^mu_", results_subset$stats$variable) | grepl("^tau_", results_subset$stats$variable)) %>% 
    # filter(!sim_id %in% sims$sim_id) %>% 
    mutate()
  
  # all group level parameters
  histograms_learning = plot_rank_hist(group_level_params_learning)+facet_wrap(~variable, nrow = 2)+theme_minimal()
  histograms_learning
  

  
  means_group_level_params_learning = group_level_params_learning %>% filter(grepl("^mu_", variable))
  
  datas = data.frame()
  parameters = unique(means_group_level_params_learning$variable)
  
  for (parameter in parameters){
    
    ranks = data.frame(means_group_level_params_learning %>% filter(variable == parameter)) %>% .$rank
    
    observed <- table(factor(ranks, levels = 1:399))  # Force levels to include all 1:399
    
    # Calculate expected frequencies (uniform distribution)
    n_ranks <- 399  # Total number of unique ranks (1:399)
    total_count <- length(ranks)      # Total number of observations
    expected <- rep(total_count / n_ranks, n_ranks)
    
    qq = chisq.test(x = observed, p = expected / total_count)
    print(paste0("parameter = ", parameter, " p = ",qq$p.value))
    
    dd = data.frame(variable = parameter, pvalue = round(qq$p.value,2), chisq = round(qq$statistic,2))
    
    datas = rbind(dd,datas)
  }

  datas$variable = c("μ[RT_sd]","μ[RT_slope]","μ[RT_int]","μ[e0]","μ[learning]")

  variables  = c("μ[learning]","μ[e0]","μ[RT_int]","μ[RT_slope]","μ[RT_sd]")
  
  #only group means
  means_group_level_params_learning = group_level_params_learning %>% 
    filter(grepl("mu",variable)) %>% 
    mutate(variable = ifelse(variable == "mu_rt_beta", "mu_rt_slope",variable),
           variable = ifelse(variable == "mu_learningrate","mu_learning",variable))
  
  plot_data = means_group_level_params_learning %>%
    mutate(variable = str_replace(variable, "^mu_", "μ["),
           variable = paste0(variable, "]"),
           variable = str_replace(variable, "rt", "RT"))%>% 
    mutate(variable = factor(variable, levels = variables))
  
  facet_labels <- as_labeller(setNames(paste0("italic('", plot_data$variable, "')"), plot_data$variable), default = label_parsed)
  
  #scatter plot:
  pp_recov_learning = plot_sim_estimated(plot_data, alpha = 0.05)+
    facet_wrap(~variable, nrow = 1, scales = "free", labeller = label_parsed)+theme_classic()+
    theme+text+
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          legend.position = "center",
          legend.box = "horizontal",
          legend.direction="horizontal",
          legend.justification='center',
          legend.key.height = unit(0.5, 'cm'),
          legend.key.width = unit(0.5, 'cm'),
          legend.spacing.x = unit(0.35, "cm"),
          plot.title = element_text(hjust = 0.5,size = font_size), #change legend key height
          legend.title = element_blank(), #change legend title font size
          legend.text = element_text(size=font_size_small),
          axis.text=element_text(size=font_size_small),
          axis.title=element_text(size=font_size),
          axis.line=element_line(size=axis_width),
          axis.ticks=element_line(size=axis_width),
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))+
    scale_x_continuous("Simulated value",breaks = scales::pretty_breaks(n = 3))+
    scale_y_continuous("",breaks = scales::pretty_breaks(n = 3))+
    facetted_pos_scales(
      y = list(
        variable == "μ[learning]" ~ scale_y_continuous(limits = c(-4-0.25,2), breaks = c(-4,-1,2), labels = c("-4","-1","2")),
        variable == "μ[e0]" ~ scale_y_continuous(limits = c(-0.5-0.25,0.5), breaks = c(-0.5,0,0.5), labels = c("-0.5","0","0.5")),
        variable == "μ[RT_int]" ~ scale_y_continuous(limits = c(-1.5-0.25,0), breaks = c(-1.5,-0.75,0), labels = c("-1.5","-0.75","0")),
        variable == "μ[RT_slope]" ~ scale_y_continuous(limits = c(0-0.25,3), breaks = c(0, 1.5,3), labels = c("0","1.5","3")),
        variable == "μ[RT_sd]" ~ scale_y_continuous(limits = c(-1.5-0.25,0), breaks = c(-1.6,-0.8,0), labels = c("-1.6","-0.8","0"))
      )
    ) +
    facetted_pos_scales(
      x = list(
        variable == "μ[learning]" ~ scale_x_continuous(limits = c(-4-0.25,2), breaks = c(-4,-1,2), labels = c("-4","-1","2")),
        variable == "μ[e0]" ~ scale_x_continuous(limits = c(-0.5-0.25,0.5), breaks = c(-0.5,0,0.5), labels = c("-0.5","0","0.5")),
        variable == "μ[RT_int]" ~ scale_x_continuous(limits = c(-1.5-0.25,0), breaks = c(-1.5,-0.75,0), labels = c("-1.5","-0.75","0")),
        variable == "μ[RT_slope]" ~ scale_x_continuous(limits = c(0-0.25,3), breaks = c(0, 1.5,3), labels = c("0","1.5","3")),
        variable == "μ[RT_sd]" ~ scale_x_continuous(limits = c(-1.5-0.25,0), breaks = c(-1.6,-0.8,0), labels = c("-1.6","-0.8","0"))
      )
    ) +
    ylab(" ")+
    xlab("Simulated value")
  
  
  pp_recov_learning
  
  
  #histograms
  histograms_learning_rw = plot_rank_hist(plot_data)+
    facet_wrap(~variable, nrow = 1, labeller = label_parsed)+theme_classic()+
    # geom_text(data = datas %>%  mutate(variable = factor(variable, levels = variables)),
    #           aes(x = 200, y = 24, label = paste0("\np = ", pvalue)), size = font_size_small/5)+
    ggtitle("Learning Paradigm")+
    scale_x_continuous(breaks = c(0,200,400), labels = c(0,200,400))+
    scale_y_continuous(" ",breaks = c(0,10,20), labels = c(0,10,20))+
    coord_cartesian(ylim = c(0,22.5))+
    theme+text+
    theme(strip.background = element_blank(),
          strip.text = element_text(font_size),
          legend.position = "center",
          legend.box = "horizontal",
          legend.direction="horizontal",
          legend.justification='center',
          legend.key.height = unit(0.5, 'cm'),
          legend.key.width = unit(0.5, 'cm'),
          legend.spacing.x = unit(0.35, "cm"),
          plot.title = element_text(hjust = 0.5,size = font_size), #change legend key height
          legend.title = element_blank(), #change legend title font size
          legend.text = element_text(size=font_size_small),
          axis.text=element_text(size=font_size_small),
          axis.title=element_text(size=font_size),
          axis.line=element_line(size=axis_width),
          axis.ticks=element_line(size=axis_width),
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))
    
  
  histograms_learning_rw
  
  comb_learning = (histograms_learning_rw) / pp_recov_learning
  
  comb_learning
  
  ###### 
  
  ## psychophysics:
  # loading
  
  load(here::here("Simulations","Psychometric","SBC","results","2000.RData"))
  
  results_psycho = results
  
  sim_ids_to_keep <- results_psycho$backend_diagnostics %>%
    dplyr::filter(n_divergent == 0 & n_max_treedepth == 0) %>%
    dplyr::pull(sim_id)
  
  
  sim_ids_to_exclude <- results_psycho$backend_diagnostics %>%
    dplyr::filter(n_divergent != 0 | n_max_treedepth != 0) %>%
    dplyr::pull(sim_id)
  
  print(paste0("Excluded based on divergent or treedepth = ",length(sim_ids_to_exclude)))
  
  # combing grouplevel from what was excluded and included to see if there are any emergent patterns
  
  results_subset <- results_psycho[sim_ids_to_keep]
  results_excluded <- results_psycho[sim_ids_to_exclude]
  
  
  grouplevel = results_excluded$stats %>% 
    filter(grepl("^mu_", results_excluded$stats$variable) | grepl("^tau_", results_excluded$stats$variable)) %>% 
    mutate(variable = gsub("\\[\\d+\\]", "", variable))%>% select(sim_id, simulated_value,variable) %>% distinct() %>% 
    mutate(sim_id = as.factor(sim_id), excluded = T) 
  
  grouplevel_inc = results_subset$stats %>% 
    filter(grepl("^mu_", results_subset$stats$variable) | grepl("^tau_", results_subset$stats$variable)) %>% 
    mutate(variable = gsub("\\[\\d+\\]", "", variable))%>% select(sim_id, simulated_value,variable) %>% distinct() %>% 
    mutate(sim_id = as.factor(sim_id), excluded = F)
  
  rbind(grouplevel,grouplevel_inc)%>% ggplot(aes(x = simulated_value, fill = excluded))+geom_histogram(col = "black")+facet_wrap(~variable, scales = "free")
  
  # combing subjeclevel from what was excluded and included to see if there are any emergent patterns
  
  subj_level = results_excluded$stats %>% 
    filter(!grepl("^mu_", results_excluded$stats$variable)& !grepl("^tau_", results_excluded$stats$variable)) %>% 
    mutate(variable = gsub("\\[\\d+\\]", "", variable))%>% select(sim_id, simulated_value,variable) %>% distinct() %>% 
    mutate(excluded = T)
  
  subj_level_inc = results_subset$stats %>% 
    filter(!grepl("^mu_", results_subset$stats$variable)& !grepl("^tau_", results_subset$stats$variable)) %>% 
    mutate(variable = gsub("\\[\\d+\\]", "", variable))%>% select(sim_id, simulated_value,variable) %>% distinct()  %>% 
    mutate(excluded = F)
  
  
  rbind(subj_level,subj_level_inc) %>% ggplot(aes(x = simulated_value, fill = excluded))+
    geom_histogram(col = "black")+facet_wrap(~variable, scales = "free")

  
  #extract the datasets
  
  sim_ids_to_keep <- results_subset$default_diagnostics %>%
    dplyr::filter(min_ess_bulk  > 400)%>%
    dplyr::filter(min_ess_tail  > 400) %>%
    dplyr::pull(sim_id)
  
  sim_ids_to_exclude <- results_subset$default_diagnostics %>%
    dplyr::filter(min_ess_bulk  < 400)%>%
    dplyr::filter(min_ess_tail  < 400) %>%
    dplyr::pull(sim_id)
  
  
  results_subset <- results_subset[sim_ids_to_keep]
  
  # Extract subject-level parameters from results$stats
  subject_level_params <- results_subset$stats %>% 
    filter(!grepl("^mu_", results_subset$stats$variable)& !grepl("^tau_", results_subset$stats$variable)) %>% 
    mutate(variable = gsub("\\[\\d+\\]", "", variable)) %>% 
    # filter(!sim_id %in% sims$sim_id) %>% 
    mutate()
  
  
  # Extract group-level parameters from results2$stats
  group_level_params_psycho <- results_subset$stats %>% 
    filter(grepl("^mu_", results_subset$stats$variable) | grepl("^tau_", results_subset$stats$variable)) %>% 
    # filter(!sim_id %in% sims$sim_id) %>% 
    mutate()
  
  # full histogram on all group level parameters
  plot_rank_hist(group_level_params_psycho)+facet_wrap(~variable, nrow = 2)+theme_minimal()
  
  # full histogram on group level means parameters
  means_group_level_params_psycho = group_level_params_psycho %>% filter(grepl("mu",variable))


  # stats for histograms:
  datas = data.frame()
  
  parameters  = c("mu_threshold", "mu_slope","mu_lapse","mu_rt_int", "mu_rt_beta", "mu_rt_sd")
  
  for (parameter in parameters){
    
    ranks = data.frame(means_group_level_params_psycho %>% filter(variable == parameter)) %>% .$rank
    
    observed <- table(factor(ranks, levels = 1:399))  # Force levels to include all 1:399
    
    # Calculate expected frequencies (uniform distribution)
    n_ranks <- 399  # Total number of unique ranks (1:399)
    total_count <- length(ranks)      # Total number of observations
    expected <- rep(total_count / n_ranks, n_ranks)
    
    qq = chisq.test(x = observed, p = expected / total_count)
    print(paste0("parameter = ", parameter, " p = ",qq$p.value))
    
    dd = data.frame(variable = parameter, pvalue = round(qq$p.value,2), chisq = round(qq$statistic,2))
    
    datas = rbind(dd,datas)
  }

  
  datas$variable = c("μ[RT_sd]","μ[RT_slope]","μ[RT_int]","μ[lapse]","μ[slope]","μ[threshold]")
  
  
  variables  = c("μ[threshold]","μ[slope]","μ[lapse]","μ[RT_int]","μ[RT_slope]","μ[RT_sd]")
  
  #only group means
  means_level_params_psycho = group_level_params_psycho %>% 
    filter(grepl("mu",variable)) %>% 
    mutate(variable = ifelse(variable == "mu_rt_beta", "mu_rt_slope",variable))
  
  plot_data_psycho = means_level_params_psycho %>%
    mutate(variable = str_replace(variable, "^mu_", "μ["),
           variable = paste0(variable, "]"),
           variable = str_replace(variable, "rt", "RT"))%>% 
    mutate(variable = factor(variable, levels = variables))
  
  facet_labels <- as_labeller(setNames(paste0("italic('", plot_data_psycho$variable, "')"), plot_data_psycho$variable), default = label_parsed)
  
  
  #scatter plot:
  pp_recov_psycho = plot_sim_estimated(plot_data_psycho%>% filter(variable != "μ[lapse]"), alpha = 0.05)+
    facet_wrap(~variable, nrow = 1, scales = "free", labeller = label_parsed)+theme_classic()+
    theme+text+
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          legend.position = "center",
          legend.box = "horizontal",
          legend.direction="horizontal",
          legend.justification='center',
          legend.key.height = unit(0.5, 'cm'),
          legend.key.width = unit(0.5, 'cm'),
          legend.spacing.x = unit(0.35, "cm"),
          plot.title = element_text(hjust = 0.5,size = font_size), #change legend key height
          legend.title = element_blank(), #change legend title font size
          legend.text = element_text(size=font_size_small),
          axis.text=element_text(size=font_size_small),
          axis.title=element_text(size=font_size),
          axis.line=element_line(size=axis_width),
          axis.ticks=element_line(size=axis_width),
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))+
    facetted_pos_scales(
      y = list(
        variable == "μ[threshold]" ~ scale_y_continuous(limits = c(-15-0.25,15), breaks = c(-15,0,15), labels = c("-15","0","15")),
        variable == "μ[slope]" ~ scale_y_continuous(limits = c(-3-0.25,0), breaks = c(-3,-1.5,0), labels = c("-3","-1.5","0")),
        variable == "μ[lapse]" ~ scale_y_continuous(limits = c(-6-0.25,0), breaks = c(-6,-3,0), labels = c("-6","-3","0")),
        variable == "μ[RT_int]" ~ scale_y_continuous(limits = c(-1.5-0.25,0), breaks = c(-1.5,-0.75,0), labels = c("-1.6","-0.8","0")),
        variable == "μ[RT_slope]" ~ scale_y_continuous(limits = c(0-0.25,3), breaks = c(0, 1.5,3), labels = c("0","0.5","3")),
        variable == "μ[RT_sd]" ~ scale_y_continuous(limits = c(-1.5-0.25,0), breaks = c(-1.5,-0.75,0), labels = c("-1.6","-0.8","0"))
      )
    ) +
    facetted_pos_scales(
      x = list(
        variable == "μ[threshold]" ~ scale_x_continuous(limits = c(-15-0.25,15), breaks = c(-15,0,15), labels = c("-15","0","15")),
        variable == "μ[slope]" ~ scale_x_continuous(limits = c(-3-0.25,0), breaks = c(-3,-1.5,0), labels = c("-3","-1.5","0")),
        variable == "μ[lapse]" ~ scale_x_continuous(limits = c(-6-0.25,0), breaks = c(-6,-3,0), labels = c("-6","-3","0")),
        variable == "μ[RT_int]" ~ scale_x_continuous(limits = c(-1.5-0.25,0), breaks = c(-1.5,-0.75,0), labels = c("-1.6","-0.8","0")),
        variable == "μ[RT_slope]" ~ scale_x_continuous(limits = c(0-0.25,3), breaks = c(0, 1.5,3), labels = c("0","1.5","3")),
        variable == "μ[RT_sd]" ~ scale_x_continuous(limits = c(-1.5-0.25,0), breaks = c(-1.5,-0.75,0), labels = c("-1.6","-0.8","0"))
      )
    )+
    scale_x_continuous("Simulated value",breaks = scales::pretty_breaks(n = 3))+
    scale_y_continuous("Estimated",breaks = scales::pretty_breaks(n = 3))+
    ylab("Estimated [q5 ; q95]")+
    xlab("Simulated value")
  pp_recov_psycho
  
  
  # %>% filter(variable != "μ[lapse]")
  #histograms
  histograms_psycho = plot_rank_hist(plot_data_psycho%>% filter(variable != "μ[lapse]"))+
    facet_wrap(~variable, nrow = 1, labeller = label_parsed)+theme_classic()+
    # geom_text(data = datas%>% filter(variable != "μ[lapse]") %>%  mutate(variable = factor(variable, levels = variables)),
    #           aes(x = 200, y = 22, label = paste0("X^2 = ", chisq, "\np = ", pvalue)), size = font_size_small/3.5)+
    ggtitle("Psychophysical Paradigm")+
    scale_x_continuous(breaks = c(0,200,400), labels = c(0,200,400))+
    scale_y_continuous("Count",breaks = c(0,10,20), labels = c(0,10,20))+
    coord_cartesian(ylim = c(0,21))+
    theme+text+
    theme(strip.background = element_blank(),
          strip.text = element_text(font_size),
          legend.position = "center",
          legend.box = "horizontal",
          legend.direction="horizontal",
          legend.justification='center',
          legend.key.height = unit(0.5, 'cm'),
          legend.key.width = unit(0.5, 'cm'),
          legend.spacing.x = unit(0.35, "cm"),
          plot.title = element_text(hjust = 0.5,size = font_size), #change legend key height
          legend.title = element_blank(), #change legend title font size
          legend.text = element_text(size=font_size_small),
          axis.text=element_text(size=font_size_small),
          axis.title=element_text(size=font_size),
          axis.line=element_line(size=axis_width ),
          axis.ticks=element_line(size=axis_width),
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))
  
  histograms_psycho
  
  comb_psycho = (histograms_psycho) / pp_recov_psycho
  comb_psycho

  
  histograms_learning_rw+histograms_psycho+plot_layout(axis_titles = "collect")
  
  pp_recov_psycho+pp_recov_learning+plot_layout(axis_titles = "collect")
  
  comb_learning = (histograms_learning_rw / pp_recov_learning)
  
  combined = comb_psycho | comb_learning
  
  
  # combined = combined+plot_layout(tag_level = 'new') +
  #   plot_annotation(tag_levels = list(c("A)", "",'B)')))
  
  combined
    
  ggsave(here::here("Figures","plot_2.tiff"),combined, dpi = 600,
         height = 6,width = 12)
  
  # ggsave(here::here("Figures","plot_2_psy.tiff"),comb_psycho, dpi = 600,
  #        width = 12, height = 6)
  # 
  # ggsave(here::here("Figures","plot_2_rw.tiff"),comb_learning, dpi = 600,
  #        width = 12, height = 6)
  # 
  return(combined)
  
}

plot_figure3 = function(){
    
  load(here::here("Simulations","Psychometric","model_recovery_estimation","results","1000_fitter_ddm_1000_samples1.RData"))
  
  #estimation time (DDM)
  map_dfr(results, 3)%>% filter(max_div == 0, max_tree == 0) %>% summarize(median = median(estimation_time),IQR = IQR(estimation_time))
      
  ddm = map_dfr(results, 1)
  
  ddms = ddm %>% 
    mutate(sim_id = rep(1:1000, each = 2)) %>% 
    filter(elpd_diff != 0,div == 0, tree == 0) %>% 
    mutate(elpd_diff = ifelse(models == "ddm", -elpd_diff, elpd_diff)) %>% 
    mutate(fitted_model = "ddm")
  
  
  ddm_model_recov = ddms %>% select(fitted_model,simulated_model,models,elpd_diff,se_diff) %>% 
    mutate(signi = ifelse((simulated_model != models) & (0 < -elpd_diff - 2* se_diff), 1, 0)) %>% 
    group_by(simulated_model, fitted_model) %>%   summarize(
      n_success = sum(signi),  # Count of signi = 1
      n_total = n(),  # Total count per simulated_model
      lower_CI = qbeta(0.025, n_success + 1, n_total - n_success + 1),  # 95% CI Lower Bound
      upper_CI = qbeta(0.975, n_success + 1, n_total - n_success + 1)   # 95% CI Upper Bound
    )
  ddm_model_recov
  
  
  # Load second dataset
  load(here::here("Simulations","Psychometric","model_recovery_estimation","results","1000_fitter_cop_rt_1000_samples1.RData"))
  
  results <- Filter(is.list, results)
  
  # estimation times:
  map_dfr(results, 3) %>% filter(max_div == 0, max_tree == 0) %>% summarize(median = median(estimation_time),IQR = IQR(estimation_time))
  map_dfr(results, 4) %>% filter(max_div == 0, max_tree == 0) %>% summarize(median = median(estimation_time),IQR = IQR(estimation_time))
  

  rts = map_dfr(results, 1)
  
  rtss = rts %>% 
    mutate(sim_id = rep(1001:2000, each = 2)) %>% 
    filter(elpd_diff != 0,div == 0, tree == 0) %>% 
    mutate(elpd_diff = ifelse(models == "ddm", -elpd_diff, elpd_diff)) %>% 
    mutate(fitted_model = "rts")
  
  
  rts_model_recov = rtss %>% select(fitted_model,simulated_model,models,elpd_diff,se_diff) %>% 
    mutate(signi = ifelse((simulated_model != models) & (0 < elpd_diff - 2* se_diff), 1, 0)) %>% 
    group_by(simulated_model, fitted_model) %>%   summarize(
      n_success = sum(signi),  # Count of signi = 1
      n_total = n(),  # Total count per simulated_model
      lower_CI = qbeta(0.025, n_success + 1, n_total - n_success + 1),  # 95% CI Lower Bound
      upper_CI = qbeta(0.975, n_success + 1, n_total - n_success + 1)   # 95% CI Upper Bound
    )
  rts_model_recov
  
  
  # Combine datasets first
  dd = rbind(ddms, rtss)
  
  
  # Arrange by elpd_diff globally and ensure sim_id is a factor
  dd = dd %>%  
    arrange(desc(elpd_diff)) %>% 
    mutate(sim_id = factor(sim_id, levels = unique(sim_id))) %>% 
    mutate(sim_id = as.numeric(sim_id))
  
  
  #Define x-axis breaks and labels
  x_breaks <- c(-100, 0,400)
  x_labels <- c("100\nFavors DDM", "0", "400\nFavors CBM")
  

  # Plot with correctly ordered sim_id
  ddm_vs_rt_psycho = dd  %>% 
    mutate(fitted_model = ifelse(fitted_model == "ddm","DDM","CBM"),
           paradigm = "Psychophysical Paradigm") %>% 
    rename(`Simulated model` = fitted_model) %>% 
    ggplot(aes(y = sim_id, x = elpd_diff, col = `Simulated model`)) +
    geom_pointrange(aes(xmin = elpd_diff - 2 * se_diff, xmax = elpd_diff + 2 * se_diff), alpha = 0.1, size = 0.1) +
    geom_point(size = 0.1) +
    geom_vline(aes(xintercept = 0), linetype = 2) +
    scale_color_manual(values = c("blue","red"))+
    scale_x_continuous(" ",breaks = x_breaks, labels = x_labels, limits = c(-150,600)) +
    theme_classic()+
    ggtitle("Psychophysical Paradigm")+
    ylab("Ordered simulations") +   
    theme+text+
    theme(strip.background = element_blank(),
          strip.text = element_text(size=font_size),
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
  
  
  
  ddm_vs_rt_psycho
  

  
  ## norts vs rts
  load(here::here("Simulations","Psychometric","model_recovery_estimation","results","1000_fitter_cop_rt_1000_samples1.RData"))
  results <- Filter(is.list, results)
  
  rts = map_dfr(results,2)
  
  
  rtss = rts %>% mutate(sim_id = rep(1:(1000),each = 2)) %>% 
    filter(elpd_diff != 0, div == 0) %>% 
    mutate(elpd_diff = ifelse(models == "nort",-elpd_diff,elpd_diff)) %>% 
    mutate(fitted_model = "rts")
  
  # number of times the Response time model significantly wins
  rt_binary_psycho = rtss %>% select(fitted_model,simulated_model,models,elpd_diff,se_diff) %>% 
    mutate(signi = ifelse((models == "nort") & (0 < elpd_diff- 2* se_diff), 1, 0)) %>% 
    group_by(simulated_model, fitted_model) %>%   summarize(
      n_success = sum(signi),  # Count of signi = 1
      n_total = n(),  # Total count per simulated_model
      lower_CI = qbeta(0.025, n_success + 1, n_total - n_success + 1),  # 95% CI Lower Bound
      upper_CI = qbeta(0.975, n_success + 1, n_total - n_success + 1)   # 95% CI Upper Bound
    )
  rt_binary_psycho
  
  
  rt_binary_psycho_loss = rtss %>% select(fitted_model,simulated_model,models,elpd_diff,se_diff) %>% 
    mutate(signi = ifelse((models == "rt") & (0 > elpd_diff), 1, 0)) %>% 
    group_by(simulated_model, fitted_model) %>%   summarize(
      n_total = n(),  # Total count per simulated_model
      n_success = sum(signi),  # Count of signi = 1
      lower_CI = qbeta(0.025, n_success + 1, n_total - n_success + 1),  # 95% CI Lower Bound
      upper_CI = qbeta(0.975, n_success + 1, n_total - n_success + 1)   # 95% CI Upper Bound
    )
  
  rt_binary_psycho_loss
  
  
  # losses significant
  rt_binary_psycho_loss = rtss %>% select(fitted_model,simulated_model,models,elpd_diff,se_diff) %>% 
    mutate(signi = ifelse((models == "rt") & (0 > elpd_diff + 2*se_diff), 1, 0)) %>% 
    group_by(simulated_model, fitted_model) %>%   summarize(
      n_total = n(),  # Total count per simulated_model
      n_success = sum(signi),  # Count of signi = 1
      lower_CI = qbeta(0.025, n_success + 1, n_total - n_success + 1),  # 95% CI Lower Bound
      upper_CI = qbeta(0.975, n_success + 1, n_total - n_success + 1)   # 95% CI Upper Bound
    )
  
  rt_binary_psycho_loss
  
  
  #Define x-axis breaks and labels
  x_breaks <- c(-12, 0,12)
  x_labels <- c("-12\nFavors with only choices", "0", "12\nFavors CBM")
  
  
  data_rt_vs_nort_psycho = rtss %>% 
    mutate(paradigm = "Psychometric Paradigm") %>% 
    arrange(desc(elpd_diff)) %>% 
    mutate(sim_id = factor(sim_id, levels = unique(sim_id)))
  
  # Plot with reordered sim_id
  rt_vs_nort_psycho_df = rtss %>% 
    mutate(paradigm = "Psychometric Paradigm") %>% 
    arrange(desc(elpd_diff)) %>% drop_na() %>% 
    mutate(sim_id = factor(sim_id, levels = unique(sim_id))) %>% 
    mutate(sim_id = as.numeric(as.factor(sim_id)))
  
  rt_vs_nort_psycho = rtss %>% 
    mutate(paradigm = "Psychometric Paradigm") %>% 
    arrange(desc(elpd_diff)) %>% drop_na() %>% 
    mutate(sim_id = factor(sim_id, levels = unique(sim_id))) %>% 
    mutate(sim_id = as.numeric(as.factor(sim_id))) %>% 
    ggplot(aes(y = sim_id, x = elpd_diff)) +
    geom_pointrange(aes(xmin = elpd_diff - 2 * se_diff, xmax = elpd_diff + 2 * se_diff), col = "blue", alpha = 0.1, size = 0.1)+
    geom_point(col = "blue",alpha = 0.5, size = 0.1)+
    geom_vline(aes(xintercept = 0), linetype = 2)+
    scale_x_continuous("ELPD Difference",breaks = x_breaks, labels = x_labels, limits = c(-15,19)) +
    coord_cartesian(ylim = c(20,nrow(rtss)-20))+
    theme_classic()+
    ylab("Ordered simulations") +
    theme+text+
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
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
  
  
  rt_vs_nort_psycho
  
  ddm_vs_rt_psycho/rt_vs_nort_psycho
  


  
  ##############################################################################
  
  ## Reinforcement learning
  
  load(here::here("Simulations","Learning","model_recovery_estimation","results","1000_rw_fitter_ddm_1000_samples1.RData"))
  
  ddm = map_dfr(results, 1)
  
  # estimation times:
  map_dfr(results, 3) %>% filter(max_div == 0, max_tree == 0) %>% summarize(median = median(estimation_time),IQR = IQR(estimation_time))
  
  
  ddms_rw = ddm %>% 
    mutate(sim_id = rep(1:1000, each = 2)) %>% 
    filter(elpd_diff != 0,div == 0, tree == 0) %>% 
    mutate(elpd_diff = ifelse(models == "ddm", -elpd_diff, elpd_diff)) %>% 
    mutate(fitted_model = "ddm")
  
  ddm_model_recov_rw = ddms_rw %>% select(fitted_model,simulated_model,models,elpd_diff,se_diff) %>% 
    mutate(signi = ifelse((simulated_model != models) & (0 < -elpd_diff - 2* se_diff), 1, 0)) %>% 
    group_by(simulated_model, fitted_model) %>%   summarize(
      n_success = sum(signi),  # Count of signi = 1
      n_total = n(),  # Total count per simulated_model
      lower_CI = qbeta(0.025, n_success + 1, n_total - n_success + 1),  # 95% CI Lower Bound
      upper_CI = qbeta(0.975, n_success + 1, n_total - n_success + 1)   # 95% CI Upper Bound
    )
  ddm_model_recov_rw
  
  
  
  ddm_model_recov_rw_significant_loss = ddms_rw %>% select(fitted_model,simulated_model,models,elpd_diff,se_diff) %>% 
    mutate(signi = ifelse((simulated_model != models) & (0 < -elpd_diff), 1, 0)) %>% 
    group_by(simulated_model, fitted_model) %>%   summarize(
      n_total = n(),  # Total count per simulated_model
      n_success = n_total-sum(signi),  # Count of signi = 1
      lower_CI = qbeta(0.025, n_success + 1, n_total - n_success + 1),  # 95% CI Lower Bound
      upper_CI = qbeta(0.975, n_success + 1, n_total - n_success + 1)   # 95% CI Upper Bound
    )
  
  ddm_model_recov_rw_significant_loss
  
  
  
  # Load second dataset
  load(here::here("Simulations","Learning","model_recovery_estimation","results","1000_rw_fitter_rt_1000_samples1.RData"))
  results <- Filter(is.list, results)
  
  # estimation times:
  map_dfr(results, 3) %>% filter(max_div == 0, max_tree == 0) %>% summarize(median = median(estimation_time),IQR = IQR(estimation_time))
  map_dfr(results, 4) %>% filter(max_div == 0, max_tree == 0) %>% summarize(median = median(estimation_time),IQR = IQR(estimation_time))
  
  
  rts = map_dfr(results, 1)
  
  rtss_rw = rts %>% 
    mutate(sim_id = rep(1001:2000, each = 2)) %>% 
    filter(elpd_diff != 0,div == 0, tree == 0) %>% 
    mutate(elpd_diff = ifelse(models == "ddm", -elpd_diff, elpd_diff)) %>% 
    mutate(fitted_model = "rts")
  
  rt_model_recov_rw = rtss_rw %>% select(fitted_model,simulated_model,models,elpd_diff,se_diff) %>% 
    mutate(signi = ifelse((simulated_model != models) & (0 < elpd_diff - 2* se_diff), 1, 0)) %>% 
    group_by(simulated_model, fitted_model) %>%   summarize(
      n_success = sum(signi),  # Count of signi = 1
      n_total = n(),  # Total count per simulated_model
      lower_CI = qbeta(0.025, n_success + 1, n_total - n_success + 1),  # 95% CI Lower Bound
      upper_CI = qbeta(0.975, n_success + 1, n_total - n_success + 1)   # 95% CI Upper Bound
    )
  rt_model_recov_rw
  
  
  
  
  # Combine datasets first
  dd_rw = rbind(ddms_rw, rtss_rw)
  
  
  # Arrange by elpd_diff globally and ensure sim_id is a factor
  dd_rw = dd_rw %>%  
    arrange(desc(elpd_diff)) %>% 
    mutate(sim_id = factor(sim_id, levels = unique(sim_id))) %>% 
    mutate(sim_id = as.numeric(sim_id))
  

  
  
  
  #Define x-axis breaks and labels
  x_breaks <- c(-200, 0,400)
  x_labels <- c("200\nFavors DDM", "0", "400\nFavors RT")
  
  # Plot with correctly ordered sim_id
  ddm_vs_rt_rw = dd_rw  %>% 
    mutate(fitted_model = ifelse(fitted_model == "ddm","DDM","CBM"),
           paradigm = "Learning Paradigm") %>% 
    rename(`Simulated model` = fitted_model) %>% 
    ggplot(aes(y = sim_id, x = elpd_diff, col = `Simulated model`)) +
    geom_pointrange(aes(xmin = elpd_diff - 2 * se_diff, xmax = elpd_diff + 2 * se_diff), alpha = 0.1, size = 0.1, show.legend = F) +
    geom_point(size = 0.1, show.legend = F) +
    scale_color_manual(values = c("blue","red"))+
    geom_vline(aes(xintercept = 0), linetype = 2) +
    ggtitle("Learning Paradigm")+
    scale_x_continuous(" ",breaks = x_breaks, labels = x_labels, limits = c(-200,600)) +
    coord_cartesian(xlim= c(-200,600))+
    theme_classic()+
    ylab(" ") + 
    theme+text+
    theme(strip.background = element_blank(),
          strip.text = element_text(size=font_size),
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
  
  ddm_vs_rt_rw
  
  
  ## norts vs rts
  
  load(here::here("Simulations","Learning","model_recovery_estimation","results","1000_rw_fitter_rt_1000_samples1.RData"))
  results <- Filter(is.list, results)
  
  rts_rw = map_dfr(results,2)
  
  rtss_rw = rts_rw %>% mutate(sim_id = rep(1:(1000),each = 2)) %>% 
    filter(elpd_diff != 0, div == 0, tree == 0) %>% 
    mutate(elpd_diff = ifelse(models == "nort",-elpd_diff,elpd_diff)) %>% 
    mutate(fitted_model = "rts")
  
  
  
  # number of times the Response time model significantly wins
  rt_model_recov_rw = rtss_rw %>% select(fitted_model,simulated_model,models,elpd_diff,se_diff) %>% 
    mutate(signi = ifelse((models == "nort") & (0 < elpd_diff- 2* se_diff), 1, 0)) %>% 
    group_by(simulated_model, fitted_model) %>%   summarize(
      n_success = sum(signi),  # Count of signi = 1
      n_total = n(),  # Total count per simulated_model
      lower_CI = qbeta(0.025, n_success + 1, n_total - n_success + 1),  # 95% CI Lower Bound
      upper_CI = qbeta(0.975, n_success + 1, n_total - n_success + 1)   # 95% CI Upper Bound
    )
  rt_model_recov_rw
  
  # losses non-significant
  rt_binary_rw_loss = rtss_rw %>% select(fitted_model,simulated_model,models,elpd_diff,se_diff) %>% 
    mutate(signi = ifelse((models == "rt") & (0 > elpd_diff), 1, 0)) %>% 
    group_by(simulated_model, fitted_model) %>%   summarize(
      n_total = n(),  # Total count per simulated_model
      n_success = sum(signi),  # Count of signi = 1
      lower_CI = qbeta(0.025, n_success + 1, n_total - n_success + 1),  # 95% CI Lower Bound
      upper_CI = qbeta(0.975, n_success + 1, n_total - n_success + 1)   # 95% CI Upper Bound
    )
  
  rt_binary_rw_loss

  # losses significant
  rt_binary_rw_loss = rtss_rw %>% select(fitted_model,simulated_model,models,elpd_diff,se_diff) %>% 
    mutate(signi = ifelse((models == "rt") & (0 > elpd_diff + 2*se_diff), 1, 0)) %>% 
    group_by(simulated_model, fitted_model) %>%   summarize(
      n_total = n(),  # Total count per simulated_model
      n_success = sum(signi),  # Count of signi = 1
      lower_CI = qbeta(0.025, n_success + 1, n_total - n_success + 1),  # 95% CI Lower Bound
      upper_CI = qbeta(0.975, n_success + 1, n_total - n_success + 1)   # 95% CI Upper Bound
    )
  
  rt_binary_rw_loss
  
    
  
  #Define x-axis breaks and labels
  x_breaks <- c(-12, 0,12)
  x_labels <- c("-12\nFavors with only choices", "0", "12\nFavors CBM")
  
  
  data_rt_vs_nort_rw = rtss_rw %>% 
    mutate(paradigm = "Learning Paradigm") %>% 
    arrange(desc(elpd_diff)) %>% mutate(sim_id = factor(sim_id, levels = unique(sim_id)))
  
  # Plot with reordered sim_id
  rt_vs_nort_rw = rtss_rw  %>% 
    mutate(paradigm = "Learning Paradigm") %>% 
    arrange(desc(elpd_diff)) %>% mutate(sim_id = factor(sim_id, levels = unique(sim_id))) %>% 
    ggplot(aes(y = sim_id, x = elpd_diff)) +
    facet_wrap(~paradigm)+
    geom_pointrange(aes(xmin = elpd_diff - 2 * se_diff, xmax = elpd_diff + 2 * se_diff), col = "blue", alpha = 0.1, size = 0.1)+
    geom_point(col = "blue",alpha = 0.5, size = 0.1)+
    geom_vline(aes(xintercept = 0), linetype = 2)+
    scale_x_continuous("ELPD Difference",breaks = x_breaks, labels = x_labels, limits = c(-13,13)) +
    coord_cartesian(ylim = c(20,1000))+
    theme_classic()+
    ylab(" ") + 
    theme+text+
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
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
  
  
  
  
  rt_vs_nort_rw
  
  #((ddm_vs_rt_psycho/posterior_sd_psycho) | (ddm_vs_rt_rw/posterior_sd_rw)) + plot_layout(guides = "collect")
  
  ddm_vs_rt = ((ddm_vs_rt_psycho) | (ddm_vs_rt_rw))+
    plot_layout(axis_titles = "collect")

  rt_vs_nort = ((rt_vs_nort_psycho) | (rt_vs_nort_rw))+
    plot_layout(axis_titles = "collect")
  
  
  combined_plot = (ddm_vs_rt / rt_vs_nort)+plot_layout(tag_level = 'new', guides = "collect", axis_titles = "collect")&
    theme(legend.position = "top",
          legend.text = element_text(size = font_size),  # Increase legend text size
          legend.title = element_text(size = font_size, face = "bold"),  # Increase title size if needed
          legend.key.size = unit(1, "cm"))& 
    guides(
      color = guide_legend(override.aes = list(
        linewidth = c(5, 5)
      ))
    )
  
  combined_plot
  

  ggsave(here::here("Figures","plot_3.tiff"),combined_plot,
         width = 12,
         height = 6,
         dpi = 600)
  
  
  return(combined_plot)
  
  
}

plot_figure4 = function(){
  # plot 4
  load(here::here("Simulations","Psychometric","model_recovery_estimation","results","1000_fitter_cop_rt_1000_samples1.RData"))
  results <- Filter(is.list, results)
  
  
  point_size = 2
  # sds plots:
  qqq = map_dfr(results,3)%>% 
    mutate(sim_id = rep(1:1000,each = 92)) %>% 
    filter(max_div == 0, max_tree == 0) %>% 
    select(sim_id,sd,variable, simulated, mean,subj_id) %>% 
    mutate(bias = simulated - mean) %>% select(sim_id,sd,variable, bias,subj_id) 
  
  qqq2 = map_dfr(results,4)%>% 
    mutate(sim_id = rep(1:1000,each = 36)) %>% 
    filter(max_div == 0, max_tree == 0) %>% 
    select(sim_id,sd,variable, simulated, mean, subj_id) %>% 
    mutate(bias_NoRT = simulated - mean,
           sd_NoRT = sd) %>% select(sim_id,variable, bias_NoRT,sd_NoRT,subj_id) 
  
  
  # inner_join(qqq,qqq2) %>% 
  #   ggplot(aes(x = bias, y = bias_NoRT))+
  #   geom_point(alpha = 0.1)+
  #   facet_wrap(~variable, scales = "free")+geom_abline()
  # 
  # inner_join(qqq,qqq2) %>% 
  #   ggplot(aes(x = sd, y = sd_NoRT))+
  #   geom_point(alpha = 0.1)+
  #   facet_wrap(~variable, scales = "free")+geom_abline()
  
  
  variables  = c("mu_threshold", "mu_slope","mu_lapse", "threshold", "slope","lapse")
  
  
  myfun <- function(x) x
  
  
  posterior_sd_psycho = inner_join(qqq,qqq2) %>% filter(variable %in% variables) %>% 
    mutate(variable = factor(variable, levels = variables),
           category = ifelse(grepl("^mu_", variable), "mu", "non-mu")) %>%
    ggplot()+
    geom_point(aes(x = sd, y = sd_NoRT,alpha = category), size = point_size)+
    facet_wrap(~variable, scales = "free")+
    geom_abline(linetype = 2)+
    geom_ribbon(stat = 'function', fun = myfun,
                mapping = aes(ymin = after_stat(y), ymax = Inf, fill = "Smaller SD with CBM"),
                alpha = 0.2) +
    geom_ribbon(stat = 'function', fun = myfun,
                mapping = aes(ymin = after_stat(y), ymax = -Inf, fill = "Smaller SD with only choices"),
                alpha = 0.2) +
    scale_fill_manual(values = c("Smaller SD with only choices" = "lightblue", "Smaller SD with CBM" = "orange"),
                      name = " ")+
    scale_alpha_manual(values = c("mu" = 0.15, "non-mu" = 0.015),  # Adjust alpha levels
                       name = "Variable Type")+
    guides(alpha = "none")+
    xlab("Estimated  SD with CBM")+
    ylab("Estimated  SD with only choices")+
    theme_classic()
  
  posterior_sd_psycho
  
  
  ### sds plots
  
  load(here::here("Simulations","Learning","model_recovery_estimation","results","1000_rw_fitter_rt_1000_samples1.RData"))
  results <- Filter(is.list, results)
  
  
  qqq_rw = map_dfr(results,3)%>% 
    mutate(sim_id = rep(1:1000,each = 80)) %>% 
    filter(max_div == 0, max_tree == 0) %>% 
    select(sim_id,sd,variable, simulated, mean,subj_id) %>% 
    mutate(bias = simulated - mean) %>% 
    select(sim_id,sd,variable, bias,subj_id) 
  
  qqq2_rw = map_dfr(results,4)%>% 
    mutate(sim_id = rep(1:1000,each = 24)) %>% 
    filter(max_div == 0, max_tree == 0) %>% 
    select(sim_id,sd,variable, simulated, mean, subj_id) %>% 
    mutate(bias_NoRT = simulated - mean,
           sd_NoRT = sd) %>% select(sim_id, variable, bias_NoRT,sd_NoRT,subj_id) 
  
  
  # inner_join(qqq,qqq2) %>% 
  #   ggplot(aes(x = bias, y = bias_NoRT))+
  #   geom_point(alpha = 0.1)+
  #   facet_wrap(~variable, scales = "free")+geom_abline()
  
  
  
  
  variables_rw  = c("mu_learningrate", "mu_e0", "learningrate", "e0")
  
  
  posterior_sd_rw = inner_join(qqq_rw,qqq2_rw) %>% filter(variable %in% variables_rw) %>% 
    mutate(variable = factor(variable, levels = c("mu_learningrate", "mu_e0", "learningrate", "e0")),
           category = ifelse(grepl("^mu_", variable), "mu", "non-mu")) %>%
    ggplot()+
    geom_point(aes(x = sd, y = sd_NoRT,alpha = category), size = point_size)+
    facet_wrap(~variable, scales = "free")+
    geom_abline(linetype = 2)+
    geom_ribbon(stat = 'function', fun = myfun,
                mapping = aes(ymin = after_stat(y), ymax = Inf, fill = "Smaller SD with CBM"),
                alpha = 0.2) +
    geom_ribbon(stat = 'function', fun = myfun,
                mapping = aes(ymin = after_stat(y), ymax = -Inf, fill = "Smaller SD with only choices"),
                alpha = 0.2) +
    scale_fill_manual(values = c("Smaller SD with only choices" = "lightblue", "Smaller SD with CBM" = "orange"),
                      name = " ")+
    scale_alpha_manual(values = c("mu" = 0.15, "non-mu" = 0.015),  # Adjust alpha levels
                       name = "Variable Type")+guides(alpha = "none")+theme_classic()+
    theme(legend.position = "top")+
    xlab("Estimated  SD with CBM")+
    ylab("Estimated  SD with only choices")+
    theme_classic()
  
  
  
  
  
  
  
  
  
  ## combined sds
  sds_rw = inner_join(qqq_rw,qqq2_rw) %>% filter(variable %in% variables_rw) %>% 
    mutate(variable = factor(variable, levels = c("mu_learningrate", "mu_e0", "learningrate", "e0")),
           category = ifelse(grepl("^mu_", variable), "mu", "non-mu"))
  
  sds_psycho = inner_join(qqq,qqq2) %>% filter(variable %in% variables) %>% 
    mutate(variable = factor(variable, levels = variables),
           category = ifelse(grepl("^mu_", variable), "mu", "non-mu"))
  
  
  
  variables  = c("μ[threshold]","μ[slope]","μ[lapse]","μ[learning]","μ[e0]",
                 "threshold","slope","lapse", "learning", "e0")
  
  plotdata = rbind(sds_psycho,sds_rw) %>% mutate(variable = as.character(variable))
  
  #only group means
  plotdata = plotdata %>% 
    mutate(variable = ifelse(variable == "mu_rt_beta", "mu_rt_slope",variable),
           variable = ifelse(variable == "rt_beta", "rt_slope",variable),
           variable = ifelse(variable == "mu_learningrate","mu_learning",variable),
           variable = ifelse(variable == "learningrate","learning",variable))
  
  plotdata_means = plotdata %>%
    filter(grepl("mu_",variable)) %>% 
    mutate(variable = str_replace(variable, "^mu_", "μ["),
           variable = paste0(variable, "]"),
           variable = str_replace(variable, "rt", "RT"))
  
  subject = plotdata %>% 
    filter(!grepl("mu_",variable))
   
  plot_data = rbind(plotdata_means,subject) %>% mutate(variable = factor(variable, levels = variables))
  
  facet_labels <- as_labeller(setNames(paste0("italic('", plot_data$variable, "')"), plot_data$variable), default = label_parsed)
  
  
  full_sds_comb =  plot_data %>% 
    ggplot()+
    geom_point(aes(x = sd, y = sd_NoRT,alpha = category), size = point_size)+
    facet_wrap(~variable, nrow = 2, scales = "free", labeller = label_parsed) +
    geom_abline(linetype = 2)+
    geom_ribbon(stat = 'function', fun = myfun,
                mapping = aes(ymin = after_stat(y), ymax = Inf, fill = "Smaller SD with RT"),
                alpha = 0.2) +
    geom_ribbon(stat = 'function', fun = myfun,
                mapping = aes(ymin = after_stat(y), ymax = -Inf, fill = "Smaller SD with only choices"),
                alpha = 0.2) +
    scale_fill_manual(values = c("Smaller SD with only choices" = "lightblue", "Smaller SD with RT" = "orange"),
                      name = " ")+
    scale_alpha_manual(values = c("mu" = 0.15, "non-mu" = 0.015),  # Adjust alpha levels
                       name = "Variable Type")+guides(alpha = "none")+theme_classic()+
    theme(legend.position = "top")+
    xlab("Estimated  SD with RT")+
    ylab("Estimated  SD with only choices")+
    theme_classic()+
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+ 
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
  
  
  full_sds_comb
  
  
  
  
  pyscho_sds_comb =  plot_data %>% filter(variable %in% c("threshold","slope",
                                                          #"μ[lapse]","lapse",
                                                          "μ[threshold]","μ[slope]")) %>% 
    mutate(alphas = as.character(as.numeric(as.factor(variable)))) %>% 
    
    ggplot()+
    geom_point(aes(x = sd, y = sd_NoRT,alpha = alphas), size = point_size)+
    facet_wrap(~variable, nrow = 2, scales = "free", labeller = label_parsed) +
    geom_abline(linetype = 2)+
    geom_ribbon(stat = 'function', fun = myfun,
                mapping = aes(ymin = after_stat(y), ymax = Inf, fill = "Smaller SD with CBM"),
                alpha = 0.2) +
    geom_ribbon(stat = 'function', fun = myfun,
                mapping = aes(ymin = after_stat(y), ymax = -Inf, fill = "Smaller SD with only choices"),
                alpha = 0.2) +
    scale_fill_manual(values = c("Smaller SD with only choices" = "lightblue", "Smaller SD with CBM" = "orange"),
                      name = " ")+
    scale_alpha_manual(values = c("1" = 0.15, "2" = 0.05,"6" = 0.025,"7" = 0.025),  # Adjust alpha levels
                       name = "Variable Type")+
    guides(alpha = "none")+theme_classic()+
    theme(legend.position = "top")+
    xlab("Estimated  SD with CBM")+
    ylab("Estimated  SD with only choices")+
    ggtitle("Psychometric Paradigm")+
    theme_classic()+
    facetted_pos_scales(
      y = list(
        variable == "μ[threshold]" ~ scale_y_continuous(limits = c(0,3), breaks = c(0,1.5,3)),
        variable == "μ[slope]" ~ scale_y_continuous(limits = c(0,0.6), breaks = c(0,0.3,0.6)),
        variable == "Threshold" ~ scale_y_continuous(limits = c(0,6), breaks = c(0,3,6)),
        variable == "slope" ~ scale_y_continuous(limits = c(0,1), breaks = c(0, 0.5,1))
      )
    ) +
    facetted_pos_scales(
      x = list(
        variable == "μ[threshold]" ~ scale_x_continuous(limits = c(0,3), breaks = c(0,1.5,3)),
        variable == "μ[slope]" ~ scale_x_continuous(limits = c(0,0.6), breaks = c(0,0.3,0.6)),
        variable == "Threshold" ~ scale_x_continuous(limits = c(0,6), breaks = c(0,3,6)),
        variable == "slope" ~ scale_x_continuous(limits = c(0,0.8), breaks = c(0, 0.4,0.8))
      )
    ) +
    theme+text+
    theme(strip.background = element_blank(),
          legend.position = "bottom",
          legend.box = "horizontal",
          legend.direction="horizontal",
          legend.justification='center',
          legend.key.height = unit(0.5, 'cm'),
          legend.key.width = unit(0.5, 'cm'),
          plot.title = element_text(hjust = 0.5, size = font_size), #change legend key height
          legend.text = element_text(size=font_size_small),
          axis.text=element_text(size=font_size_small),
          axis.title=element_text(size=font_size),
          axis.line=element_line(size=axis_width),
          axis.ticks=element_line(size=axis_width),
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))
  
  pyscho_sds_comb
  
  
  
  variables  = c("μ[threshold]","μ[slope]","μ[lapse]","μ[learning]","μ[e0]",
                 "threshold","slope","lapse", "learning", "e0")
  
  
  
  rw_sds_comb =  plot_data %>% filter(variable %in% c("learning","e0","μ[learning]","μ[e0]")) %>% 
    mutate(alphas = as.character(as.numeric(as.factor(variable)))) %>% 
    ggplot()+
    geom_point(aes(x = sd, y = sd_NoRT,alpha = alphas), size = point_size)+
    facet_wrap(~variable, nrow = 2, scales = "free", labeller = label_parsed) +
    geom_abline(linetype = 2)+
    geom_ribbon(stat = 'function', fun = myfun,
                mapping = aes(ymin = after_stat(y), ymax = Inf, fill = "Smaller SD with CBM"),
                alpha = 0.2, show.legend = F) +
    geom_ribbon(stat = 'function', fun = myfun,
                mapping = aes(ymin = after_stat(y), ymax = -Inf, fill = "Smaller SD with only choices"),
                alpha = 0.2, show.legend = F) +
    scale_fill_manual(values = c("Smaller SD with only choices" = "lightblue", "Smaller SD with CBM" = "orange"),
                      name = " ")+
    scale_alpha_manual(values = c("4" = 0.04, "5" = 0.15, "9" = 0.015, "10" = 0.015),  # Adjust alpha levels
                       name = "Variable Type")+guides(alpha = "none")+theme_classic()+
    theme(legend.position = "top")+
    xlab("Estimated  SD with CBM")+
    ylab("Estimated  SD with only choices")+
    ggtitle("Learning Paradigm")+
    theme_classic()+
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+ 
    theme+text+
    facetted_pos_scales(
      y = list(
        variable == "μ[learning]" ~ scale_y_continuous(limits = c(0.1,0.4), breaks = c(0.1,0.25,0.4)),
        variable == "μ[e0]" ~ scale_y_continuous(limits = c(0.12,0.2), breaks = c(0.12,0.16,0.2)),
        variable == "learning" ~ scale_y_continuous(limits = c(0,0.1), breaks = c(0,0.05,0.1)),
        variable == "e0" ~ scale_y_continuous(limits = c(0.05,0.15), breaks = c(0.05, 0.1,0.15))
      )
    ) +
    facetted_pos_scales(
      x = list(
        variable == "μ[learning]" ~ scale_x_continuous(limits = c(0.1,0.4), breaks = c(0.1,0.25,0.4)),
        variable == "μ[e0]" ~ scale_x_continuous(limits = c(0.12,0.2), breaks = c(0.12,0.16,0.2)),
        variable == "learning" ~ scale_x_continuous(limits = c(0,0.1), breaks = c(0,0.05,0.1)),
        variable == "e0" ~ scale_x_continuous(limits = c(0.05,0.15), breaks = c(0.05, 0.1,0.15))
      )
    ) +
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
          axis.line=element_line(size=axis_width),
          axis.ticks=element_line(size=axis_width),
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))
  
  rw_sds_comb
  
  
  plot4 = (pyscho_sds_comb | rw_sds_comb)+plot_layout(tag_level = 'new', guides = "collect", axis_titles = "collect")&
    theme(legend.position = "top",
          legend.text = element_text(size = font_size),  # Increase legend text size
          legend.title = element_text(size = font_size, face = "bold"),  # Increase title size if needed
          legend.key.size = unit(1, "cm"))& 
    guides(
      color = guide_legend(override.aes = list(
        linewidth = c(1.5, 1.5)
      ))
    )
  

  ggsave(here::here("Figures","plot_4.tiff"),plot4, dpi = 600,
         height = 6,width = 12)
  
  return(plot4)

}



plot_figure5 = function(rerun = F){

  
  if(rerun){
    
    
  ################### 
  
  # reinforcement learning:
  df_rw <- read_csv(here::here("Real data","Learning","data","Data_as_csv.csv"))%>% 
    group_by(Subject) %>% 
    mutate(trial = 1:n(),
           rts  = exp(rts) / 1000) %>% drop_na()
  
  df_filtered <- df_rw %>%
    group_by(Subject) %>%
    filter(rts > min(rts))
  
  df_rw = df_filtered
  
  source(here::here("Real data","Learning","utility.R"))
  # models:
  rw_ddm_nolow <- readRDS(here::here("Real data","Learning","models","rw_ddm_nolow.rds"))
  
  ## ddm
  
  #draws to plot
  draw_id = sample(1:4000,500)
  
  big_df_ddm = data.frame()
  
  for(s in unique(df_rw$Subject)){
    print(s)
    us = df_rw %>% filter(Subject == s) %>% .$u
    
    minRT = min(df_rw %>% filter(Subject == s) %>% .$rts)
    
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
  #summary of each draw (marginal mean and median)
  qq_summar_ddm = big_df_ddm %>%  ungroup() %>%  mutate(rts = (rts),bin_resp = resp) %>% drop_na() %>%
    group_by(trial,draw) %>% 
    summarize(resp = mean(bin_resp), rt = median(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n()))
  #the mean / median of that
  qq_summar_sum_ddm = qq_summar_ddm %>%   group_by(trial) %>% 
    summarize(resps = mean(resp), rts = median(rt), se_resp =  (mean(resp) * (1- mean(resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n()))
  
  ## check the plot
  # df_rw %>%   ungroup() %>%  mutate(rts = (rts),bin_resp = bin_resp) %>% drop_na() %>%
  #   group_by(trial) %>% 
  #   summarize(resp = mean(bin_resp), rt = median(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>%
  #   ggplot(aes(x = trial))+
  #   geom_line(data = qq_summar_ddm, aes(x = trial, y = rt, group = draw), alpha = 0.01, col = "orange")+
  #   geom_pointrange(aes(y = rt, ymin = rt-2*se_rts, ymax = rt+2*se_rts))+
  #   geom_line(data = qq_summar_sum_ddm, aes(x = trial, y = rts), col = "red")+
  #   scale_y_continuous(limits = c(0,1))+
  #   theme_classic()
  
  ## check the rts
  # df_rw %>%   ungroup() %>%  mutate(rts = (rts),bin_resp = bin_resp) %>% drop_na() %>%
  #   group_by(trial) %>% 
  #   summarize(resp = mean(bin_resp), rt = median(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>%
  #   ggplot(aes(x = trial))+
  #   geom_line(data = qq_summar_ddm, aes(x = trial, y = resp, group = draw), alpha = 0.01, col = "orange")+
  #   geom_line(data = qq_summar_sum_ddm, aes(x = trial, y = resps), col = "red")+
  #   geom_pointrange(aes(y = resp, ymin = resp-2*se_resp, ymax = resp+2*se_resp))+
  #   theme_classic()

  
  ############ Rts
  
  rw_rt_nolow <- readRDS(here::here("Real data","Learning","models","rw_rt_nolow.rds"))
  
  
  # subject level 
  big_df_rt = data.frame()
  
  for(s in unique(df_rw$Subject)){
    print(s)
    us = df_rw %>% filter(Subject == s) %>% .$u
    
    minRT = min(df_rw %>% filter(Subject == s) %>% .$rts)
    
    parameters = paste0(c("learningrate","e0","rt_int","rt_beta","rt_ndt","rt_sd","rho"), "[",s,"]")
    
    dfq = as_draws_df(rw_rt_nolow$draws(parameters))%>% select(-contains(".")) %>% mutate(draw = 1:n()) %>% 
      filter(draw %in% draw_id) %>% 
      rename_with(~c("learningrate","e0","rt_int","rt_beta","rt_ndt","rt_sd","rho","draw")) %>% 
      group_by(draw) %>% 
      summarize(list(generate_trialwise_data_rt_copula(u = us, learningrate = learningrate,e0 = e0,participant = s,
                                                       rt_int = rt_int,rt_beta = rt_beta,rt_sd = rt_sd,rt_ndt = rt_ndt,rho = rho))) %>% unnest()
    
    big_df_rt = rbind(big_df_rt,dfq)
  }
  
  # summary of each draw (marginal mean and median)
  qq_summar_rt = big_df_rt %>%  ungroup() %>%  mutate(rts = (rts),bin_resp = resp) %>% drop_na() %>%
    group_by(trial,draw) %>% 
    summarize(resp = mean(bin_resp), rt = median(rts))
  #mean / median of that
  qq_summar_sum_rt = qq_summar_rt %>% 
    group_by(trial) %>% 
    summarize(resps = mean(resp), rts = median(rt))
  
  # # plot to ensure correctness
  # df_rw %>%   ungroup() %>%  mutate(rts = (rts),bin_resp = bin_resp) %>% drop_na() %>%
  #   group_by(trial) %>% 
  #   summarize(resp = mean(bin_resp), rt = median(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>%
  #   ggplot(aes(x = trial))+
  #   geom_line(data = qq_summar_rt, aes(x = trial, y = rt, group = draw), alpha = 0.01, col = "orange")+
  #   geom_pointrange(aes(y = rt, ymin = rt-2*se_rts, ymax = rt+2*se_rts))+
  #   geom_line(data = qq_summar_sum_rt, aes(x = trial, y = rts), col = "red")+
  #   scale_y_continuous(limits = c(0,1))+
  #   theme_classic()

  # # plot to ensure correctness
  # df_rw %>%   ungroup() %>%  mutate(rts = (rts),bin_resp = bin_resp) %>% drop_na() %>%
  #   group_by(trial) %>% 
  #   summarize(resp = mean(bin_resp), rt = median(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>%
  #   ggplot(aes(x = trial))+
  #   geom_line(data = qq_summar_rt, aes(x = trial, y = resp, group = draw), alpha = 0.01, col = "orange")+
  #   geom_pointrange(aes(y = resp, ymin = resp-2*se_resp, ymax = resp+2*se_resp))+
  #   geom_line(data = qq_summar_sum_rt, aes(x = trial, y = resps), col = "red")+
  #   theme_classic()
  
  
  
  ## subject level  no rts
  rw_rt_nolow <- readRDS(here::here("Real data","Learning","models","rw_nort_nolow.rds"))
  
  
  big_df_nort = data.frame()
  
  for(s in unique(df_rw$Subject)){
    print(s)
    us = df_rw %>% filter(Subject == s) %>% .$u
    
    parameters = paste0(c("learningrate","e0"), "[",s,"]")
    
    dfq = as_draws_df(rw_nort_nolow$draws(parameters))%>% select(-contains(".")) %>% mutate(draw = 1:n()) %>% 
      filter(draw %in% draw_id) %>% 
      rename_with(~c("learningrate","e0","draw")) %>% 
      group_by(draw) %>% 
      summarize(list(generate_learn(u = us, learningrate = learningrate,e0 = e0,participant = s))) %>% unnest()
    
    big_df_nort = rbind(big_df_nort,dfq)
  }
  
  # median / mean
  qq_summar_nort = big_df_nort %>% mutate(bin_resp = resp) %>% 
    ungroup() %>% drop_na() %>% 
    group_by(trial,draw) %>% 
    summarize(resp = mean(bin_resp))
  #median mean
  
  qq_summar_sum_nort = qq_summar_nort %>%   
    group_by(trial) %>% 
    summarize(resps = mean(resp))
  
  
  # plot for correctness:
  # df_rw %>%   ungroup() %>% drop_na() %>%
  #   group_by(trial) %>% 
  #   summarize(resp = mean(bin_resp), rt = median(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>%
  #   ggplot(aes(x = trial))+
  #   geom_line(data = qq_summar_nort, aes(x = trial, y = resp, group = draw), alpha = 0.01, col = "#00C853")+
  #   geom_line(data = qq_summar_sum_nort, aes(x = trial, y = resps), col = "green")+
  #   geom_pointrange(aes(y = resp, ymin = resp-se_resp, ymax = resp+se_resp))+
  #   theme_classic()
  
  
  
  
  ## combine
  

  alpha = 0.3
  
  ## Combine
  rw_rw = df_rw %>%   ungroup() %>% drop_na()  %>%
    group_by(trial)  %>% 
    summarize(resp = mean(bin_resp), rt = median(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>%
    ggplot(aes(x = trial))+
    geom_line(data = qq_summar_rt %>% filter(draw %in% draw_id), aes(x = trial, y = resp, group = draw), alpha = alpha, col = "#6CEEF8")+
    # geom_line(data = qq_summar_nort, aes(x = Difficulty_bin, y = resp, group = draw), alpha = 0.01, col = "black")+
    # geom_line(data = qq_summar_sum_nort, aes(x = trial, y = resps, col = "No RT"))+
    geom_line(data = qq_summar_ddm%>% filter(draw %in% draw_id), aes(x = trial, y = resp, group = draw), alpha = alpha, col = "orange")+
    geom_pointrange(aes(y = resp, ymin = resp-2*se_resp, ymax = resp+2*se_resp, col = "Data"))+
    geom_line(data = qq_summar_sum_ddm , aes(x = trial, y = resps, col = "DDM"), linewidth = 1.1)+
    # geom_line(data = qq_summar_sum_nort, aes(x = trial, y = resps, col = "No RT"), linewidth = 1.1)+
    geom_line(data = qq_summar_sum_rt, aes(x = trial, y = resps, col = "CBM"), linewidth = 1.1)+
    theme_classic()+
    scale_color_manual(name = "Model Predictions & Data", 
                       values = c("DDM" = "red", 
                                  "CBM" = "blue", 
                                  "Data" = "black"),
                       breaks = c("Data", "CBM", "DDM"))+
    scale_x_continuous(" ", breaks = scales::pretty_breaks(n = 3))+
    scale_y_continuous(" ", breaks = scales::pretty_breaks(n = 3))+ 
    ggtitle("Learning Paradigm")+
    theme+text+
    theme(strip.background = element_blank(),
          strip.text = element_text(size=font_size),
          legend.box = "horizontal",
          legend.position = "top",
          legend.direction="horizontal",
          legend.justification='center',
          legend.key.height = unit(0.5, 'cm'),
          legend.key.width = unit(0.5, 'cm'),
          legend.spacing.x = unit(0.35, "cm"),
          plot.title = element_text(hjust = 0.5, size = font_size), #change legend key height
          legend.text = element_text(size=font_size_small),
          axis.text=element_text(size=font_size_small),
          axis.title=element_text(size=font_size),
          axis.line=element_line(size=axis_width),
          axis.ticks=element_line(size=axis_width),
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))
  
  
  rw_rw
  
  # rts
  
  
  rts_rw = df_rw %>%   ungroup() %>% drop_na() %>%
    group_by(trial) %>% 
    summarize(resp = mean(bin_resp), rt = median(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>%
    ggplot(aes(x = trial))+
    geom_line(data = qq_summar_rt%>% filter(draw %in% draw_id), aes(x = trial, y = rt, group = draw), alpha = alpha, col = "#6CEEF8")+
    geom_line(data = qq_summar_ddm%>% filter(draw %in% draw_id), aes(x = trial, y = rt, group = draw), alpha = alpha, col = "orange")+
    geom_pointrange(aes(y = rt, ymin = rt-2*se_rts, ymax = rt+2*se_rts, col = "Data"), show.legend = F)+
    geom_line(data = qq_summar_sum_rt, aes(x = trial, y = rts,col = "CBM"), linewidth = 1.1, show.legend = F)+
    geom_line(data = qq_summar_sum_ddm, aes(x = trial, y = rts, col = "DDM"), linewidth = 1.1, show.legend = F)+
    scale_color_manual(name = "Model Predictions & Data", 
                       values = c("DDM" = "red", 
                                  "CBM" = "blue", 
                                  "Data" = "black"),
                       breaks = c("Data", "CBM", "DDM"))+
    theme_classic()+
    scale_x_continuous("Trial", breaks = scales::pretty_breaks(n = 3))+
    scale_y_continuous(" ", limits = c(0.275,0.625), breaks = c(0.3,0.4,0.5,0.6), labels = c(0.3,0.4,0.5,0.6))+
    theme+text+
    theme(strip.background = element_blank(),
          strip.text = element_text(size=font_size),
          legend.position = "top",
          legend.box = "horizontal",
          legend.direction="horizontal",
          legend.justification='center',
          legend.key.height = unit(0.5, 'cm'),
          legend.key.width = unit(0.5, 'cm'),
          legend.spacing.x = unit(0.35, "cm"),
          plot.title = element_text(hjust = 0.5, size = font_size), #change legend key height
          legend.text = element_text(size=font_size_small),
          axis.text=element_text(size=font_size_small),
          axis.title=element_text(size=font_size),
          axis.line=element_line(size=axis_width),
          axis.ticks=element_line(size=axis_width),
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))
  
  
  rts_rw
  
  
  
  ############ 
  # psychometric:
  df_psy <- read_csv(here::here("Real data","Psychometric","data","data_Bang_2019_Exp1.csv"))
  
  df_psy = df_psy %>% filter(Condition == 1, Day == 1) %>% filter(SN != max(SN) & SN != min(SN)) %>%
    mutate(ACC = ifelse(Stimulus == Response,1,0),
           Difficulty  = ifelse(Stimulus == 2, SN , -SN),
           Response = Response-1,
           Confidence = Confidence/max(Confidence,na.rm = T)) %>% rename(participant = Subj_idx) %>% 
    group_by(participant) %>% mutate(trial = 1:n())
  
  
  df_filtered <- df_psy %>%
    group_by(participant) %>%
    filter(RT_dec > min(RT_dec))
  
  df_psy = df_filtered
  
  draw_id = sample(1:4000,500)
  
  # ddm
  psy_ddm <- readRDS(here::here("Real data","Psychometric","models","psy_ddm.rds"))
  
  source(here::here("Real data","Psychometric","utility.R"))
  
  # subject level 
  big_df_ddm_psy = data.frame()
  
  for(s in unique(df_psy$participant)){
    print(s)
    xs = seq(-0.3,0.3,by = 0.01)
    
    minRT = min(df_psy %>% filter(participant == s) %>% .$RT_dec)
    
    parameters = paste0(c("threshold","slope","lapse","alpha","beta","delta","rt_ndt"), "[",s,"]")
    
    dfq = as_draws_df(psy_ddm$draws(parameters))%>% select(-contains(".")) %>% mutate(draw = 1:n()) %>% 
      filter(draw %in% draw_id) %>% 
      rename_with(~c("threshold","slope","lapse","alpha","beta","delta","rt_ndt","draw")) %>% 
      group_by(draw) %>% 
      summarize(list(generate_trial_ddm(x = xs, threshold = threshold,slope = slope, lapse = lapse,participant = s,
                                        alpha = alpha,beta = beta,delta = delta,rt_ndt = rt_ndt))) %>% unnest() %>% 
      rename(RT_dec = q) %>% 
      mutate(resp = ifelse(resp == "upper",1,0))
    
  
    big_df_ddm_psy = rbind(big_df_ddm_psy,dfq)
    
  }

  
  
  
  # summary of the median and mean
  qq_summar_psy = big_df_ddm_psy %>%  ungroup() %>%  mutate(rts = (RT_dec),bin_resp = resp, Difficulty = x) %>% drop_na() %>%
    mutate(Difficulty_bin = cut((Difficulty), breaks = 10, labels = FALSE),
           Difficulty = abs(Difficulty),
           Difficulty_bin = (Difficulty_bin - min(Difficulty_bin)) / (max(Difficulty_bin) - min(Difficulty_bin)) * (0.3 - (-0.3)) + (-0.3)
    ) %>%
    group_by(Difficulty_bin,draw) %>% 
    summarize(resp = mean(bin_resp), rt = median(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n()))
  
  # same
  qq_summar_sum_psy = qq_summar_psy %>%   group_by(Difficulty_bin) %>% 
    summarize(resps = mean(resp), rts = median(rt), se_resp =  (mean(resp) * (1- mean(resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n()))
  
  
  # # check the results:
  # df %>%   ungroup() %>%  mutate(rts = (RT_dec),bin_resp = Response) %>% drop_na() %>%
  #   mutate(Difficulty_bin = cut((Difficulty), breaks = 5, labels = FALSE),
  #          Difficulty = abs(Difficulty),
  #          Difficulty_bin = (Difficulty_bin - min(Difficulty_bin)) / (max(Difficulty_bin) - min(Difficulty_bin)) * (0.3 - (-0.3)) + (-0.3)
  #   ) %>%
  #   group_by(Difficulty_bin) %>% 
  #   summarize(resp = mean(bin_resp), rt = median(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>%
  #   ggplot(aes(x = Difficulty_bin))+
  #   geom_line(data = qq_summar, aes(x = Difficulty_bin, y = rt, group = draw), alpha = 0.01, col = "orange")+
  #   geom_line(data = qq_summar_sum, aes(x = Difficulty_bin, y = rts), col = "red")+
  #   geom_pointrange(aes(y = rt, ymin = rt-2*se_rts, ymax = rt+2*se_rts))+
  #   theme_classic()


  # df %>%   ungroup() %>%  mutate(rts = (RT_dec),bin_resp = Response) %>% drop_na() %>%
  #   mutate(Difficulty_bin = cut((Difficulty), breaks = 5, labels = FALSE),
  #          Difficulty = abs(Difficulty),
  #          Difficulty_bin = (Difficulty_bin - min(Difficulty_bin)) / (max(Difficulty_bin) - min(Difficulty_bin)) * (0.3 - (-0.3)) + (-0.3)
  #   ) %>%
  #   group_by(Difficulty_bin) %>% 
  #   summarize(resp = mean(bin_resp), rt = median(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>%
  #   ggplot(aes(x = Difficulty_bin))+
  #   geom_line(data = qq_summar, aes(x = Difficulty_bin, y = resp, group = draw), alpha = 0.01, col = "orange")+
  #   geom_line(data = qq_summar_sum, aes(x = Difficulty_bin, y = resps), col = "red")+
  #   geom_pointrange(aes(y = resp, ymin = resp-se_resp, ymax = resp+se_resp))+
  #   theme_classic()
  
  
  
  ### rts

  psy_rt <- readRDS(here::here("Real data","Psychometric","models","psy_rt.rds"))
  
  # subject level   
  big_df_rt_psy = data.frame()
  for(s in unique(df_psy$participant)){
    print(s)
    xs = seq(-0.3,0.3,by = 0.01)
    
    minRT = min(df_psy %>% filter(participant == s) %>% .$RT_dec)
    
    parameters = paste0(c("threshold","slope","lapse","rt_int","rt_beta","rt_ndt","rt_sd","rho"), "[",s,"]")
    
    dfq = as_draws_df(psy_rt$draws(parameters))%>% select(-contains(".")) %>% mutate(draw = 1:n()) %>% 
      filter(draw %in% draw_id) %>% 
      rename_with(~c("threshold","slope","lapse","rt_int","rt_beta","rt_ndt","rt_sd","rho","draw")) %>% 
      group_by(draw) %>% 
      summarize(list(generate_trial_shannon_entropy_rt(x = xs, threshold = threshold,slope = slope, lapse = lapse,participant = s,
                                                       rt_int = rt_int,rt_beta = rt_beta,rt_ndt = rt_ndt,rt_sd = rt_sd,rho = rho))) %>% unnest()
    big_df_rt_psy = rbind(big_df_rt_psy,dfq)
    
  }
  
  # summaries
  
  qq_summar_rt_psy = big_df_rt_psy %>% rename(Response = expectation,RT_dec = rts) %>%   
    ungroup() %>%  
    mutate(rts = (RT_dec),bin_resp = resp, Difficulty = x) %>% drop_na() %>%
    mutate(Difficulty_bin = cut((Difficulty), breaks = 10, labels = FALSE),
           Difficulty = abs(Difficulty),
           Difficulty_bin = (Difficulty_bin - min(Difficulty_bin)) / (max(Difficulty_bin) - min(Difficulty_bin)) * (0.3 - (-0.3)) + (-0.3)
    ) %>%
    group_by(Difficulty_bin,draw) %>% 
    summarize(resp = mean(bin_resp), rt = median(rts))
  
  
  
  qq_summar_sum_rt_psy = qq_summar_rt_psy %>%   group_by(Difficulty_bin) %>% 
    summarize(resps = mean(resp), rts = median(rt))
  
  
  # # sanity plots:
  # df_psy %>%   ungroup() %>%  mutate(rts = (RT_dec),bin_resp = Response) %>% drop_na() %>%
  #   mutate(Difficulty_bin = cut((Difficulty), breaks = 5, labels = FALSE),
  #          Difficulty = abs(Difficulty),
  #          Difficulty_bin = (Difficulty_bin - min(Difficulty_bin)) / (max(Difficulty_bin) - min(Difficulty_bin)) * (0.3 - (-0.3)) + (-0.3)
  #   ) %>%
  #   group_by(Difficulty_bin) %>% 
  #   summarize(resp = mean(bin_resp), rt = median(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>%
  #   ggplot(aes(x = Difficulty_bin))+
  #   geom_line(data = qq_summar_rt, aes(x = Difficulty_bin, y = rt, group = draw), alpha = 0.01, col = "#6CEEF8")+
  #   geom_line(data = qq_summar_sum_rt, aes(x = Difficulty_bin, y = rts), col = "blue")+
  #   geom_pointrange(aes(y = rt, ymin = rt-se_rts, ymax = rt+se_rts))+
  #   theme_classic()
  
  
  # df_psy %>%   ungroup() %>%  mutate(rts = (RT_dec),bin_resp = Response) %>% drop_na() %>%
  #   mutate(Difficulty_bin = cut((Difficulty), breaks = 5, labels = FALSE),
  #          Difficulty = abs(Difficulty),
  #          Difficulty_bin = (Difficulty_bin - min(Difficulty_bin)) / (max(Difficulty_bin) - min(Difficulty_bin)) * (0.3 - (-0.3)) + (-0.3)
  #   ) %>%
  #   group_by(Difficulty_bin) %>% 
  #   summarize(resp = mean(bin_resp), rt = median(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>%
  #   ggplot(aes(x = Difficulty_bin))+
  #   geom_line(data = qq_summar_rt, aes(x = Difficulty_bin, y = resp, group = draw), alpha = 0.01, col = "#6CEEF8")+
  #   geom_line(data = qq_summar_sum_rt, aes(x = Difficulty_bin, y = resps), col = "blue")+
  #   geom_pointrange(aes(y = resp, ymin = resp-se_resp, ymax = resp+se_resp))+
  #   theme_classic()
  
  
  
  # no rts:
  psy_nort <- readRDS(here::here("Real data","Psychometric","models","psy_nort.rds"))
  
  big_df_nort_psy = data.frame()
  for(s in unique(df_psy$participant)){
    print(s)
    xs = seq(-0.3,0.3,by = 0.01)
    
    
    parameters = paste0(c("threshold","slope","lapse"), "[",s,"]")
    
    dfq = as_draws_df(psy_nort$draws(parameters))%>% select(-contains(".")) %>% mutate(draw = 1:n()) %>% 
      filter(draw %in% draw_id) %>% 
      rename_with(~c("threshold","slope","lapse","draw")) %>% 
      group_by(draw) %>% 
      summarize(list(generate_psycho(x = xs, threshold = threshold,slope = slope, lapse = lapse,participant = s))) %>% unnest()
    
    big_df_nort_psy = rbind(big_df_nort_psy,dfq)
    
  }
  
  #summary
  qq_summar_nort_psy = big_df_nort_psy %>% 
    ungroup() %>%  mutate(bin_resp = resp, Difficulty = x) %>% drop_na() %>%
    mutate(Difficulty_bin = cut((Difficulty), breaks = 10, labels = FALSE),
           Difficulty = abs(Difficulty),
           Difficulty_bin = (Difficulty_bin - min(Difficulty_bin)) / (max(Difficulty_bin) - min(Difficulty_bin)) * (0.3 - (-0.3)) + (-0.3)
    ) %>%
    group_by(Difficulty_bin,draw) %>% 
    summarize(resp = mean(bin_resp))
  
  
  qq_summar_sum_nort_psy = qq_summar_nort_psy %>%   
    group_by(Difficulty_bin) %>% 
    summarize(resps = mean(resp))
  
  
  # #plot
  # df %>%   ungroup() %>%  mutate(rts = (RT_dec),bin_resp = Response) %>% drop_na() %>%
  #   mutate(Difficulty_bin = cut((Difficulty), breaks = 5, labels = FALSE),
  #          Difficulty = abs(Difficulty),
  #          Difficulty_bin = (Difficulty_bin - min(Difficulty_bin)) / (max(Difficulty_bin) - min(Difficulty_bin)) * (0.3 - (-0.3)) + (-0.3)
  #   ) %>%
  #   group_by(Difficulty_bin) %>% 
  #   summarize(resp = mean(bin_resp), rt = median(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>%
  #   ggplot(aes(x = Difficulty_bin))+
  #   geom_line(data = qq_summar_nort, aes(x = Difficulty_bin, y = resp, group = draw), alpha = 0.01, col = "#00C853")+
  #   geom_line(data = qq_summar_sum_nort, aes(x = Difficulty_bin, y = resps), col = "green")+
  #   geom_pointrange(aes(y = resp, ymin = resp-se_resp, ymax = resp+se_resp))+
  #   theme_classic()
  
  
  
  
  
  
  ## psychometric fits
  alpha = 0.3
  
  psycho_psycho_psy = df_psy %>%   ungroup() %>%  mutate(rts = (RT_dec),bin_resp = Response) %>% drop_na() %>%
    mutate(Difficulty_bin = cut((Difficulty), breaks = 10, labels = FALSE),
           Difficulty = abs(Difficulty),
           Difficulty_bin = (Difficulty_bin - min(Difficulty_bin)) / (max(Difficulty_bin) - min(Difficulty_bin)) * (0.3 - (-0.3)) + (-0.3)
    ) %>%
    group_by(Difficulty_bin) %>% 
    summarize(resp = mean(bin_resp), rt = median(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>%
    ggplot(aes(x = Difficulty_bin))+
    geom_line(data = qq_summar_rt_psy %>% filter(draw %in% draw_id), aes(x = Difficulty_bin, y = resp, group = draw), alpha = alpha, col = "#6CEEF8")+
    
    # geom_line(data = qq_summar_nort, aes(x = Difficulty_bin, y = resp, group = draw), alpha = 0.01, col = "black")+
    
    geom_line(data = qq_summar_psy %>% filter(draw %in% draw_id), aes(x = Difficulty_bin, y = resp, group = draw), alpha = alpha, col = "orange")+
    geom_line(data = qq_summar_sum_psy, aes(x = Difficulty_bin, y = resps, col = "DDM"), linewidth = 1.1, show.legend = F)+
    geom_line(data = qq_summar_sum_rt_psy, aes(x = Difficulty_bin, y = resps, col = "CBM"), linewidth = 1.1, show.legend = F)+
    # geom_line(data = qq_summar_sum_nort, aes(x = Difficulty_bin, y = resps, col = "No RT"), linewidth = 1.1)+
    geom_pointrange(aes(y = resp, ymin = resp-2*se_resp, ymax = resp+2*se_resp, col = "Data"), show.legend = F)+
    theme_classic()+
    scale_color_manual(name = "Model Predictions & Data", 
                       values = c("DDM" = "red", 
                                  "CBM" = "blue", 
                                  "Data" = "black"),
                       breaks = c("Data", "CBM", "DDM"))+
    scale_x_continuous(" ", breaks = scales::pretty_breaks(n = 3))+
    scale_y_continuous("P(Response == 1)", breaks = scales::pretty_breaks(n = 3))+  
    ggtitle("Psychophysical Paradigm")+
    theme+text+
    theme(strip.background = element_blank(),
          strip.text = element_text(size=font_size),
          legend.position = "top",
          legend.box = "horizontal",
          legend.direction="horizontal",
          legend.key.height = unit(0.5, 'cm'),
          legend.key.width = unit(0.5, 'cm'),
          legend.spacing.x = unit(0.35, "cm"),
          plot.title = element_text(hjust = 0.5, size = font_size), #change legend key height
          legend.text = element_text(size=font_size_small),
          axis.text=element_text(size=font_size_small),
          axis.title=element_text(size=font_size),
          axis.line=element_line(size=axis_width),
          axis.ticks=element_line(size=axis_width),
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))
  
  
  
  psycho_psycho_psy
  ## rts
  
  rts_psycho_psy = df_psy %>%   ungroup() %>%  mutate(rts = (RT_dec),bin_resp = Response) %>% drop_na() %>%
    mutate(Difficulty_bin = cut((Difficulty), breaks = 10, labels = FALSE),
           Difficulty = abs(Difficulty),
           Difficulty_bin = (Difficulty_bin - min(Difficulty_bin)) / (max(Difficulty_bin) - min(Difficulty_bin)) * (0.3 - (-0.3)) + (-0.3)
    ) %>%
    group_by(Difficulty_bin) %>% 
    summarize(resp = mean(bin_resp), rt = median(rts), se_resp =  (mean(bin_resp) * (1- mean(bin_resp))) / sqrt(n()), se_rts = sd(rts)/sqrt(n())) %>%
    ggplot(aes(x = Difficulty_bin))+
    geom_line(data = qq_summar_rt_psy%>% filter(draw %in% draw_id), aes(x = Difficulty_bin, y = rt, group = draw), alpha = alpha, col = "#6CEEF8")+
    geom_line(data = qq_summar_psy%>% filter(draw %in% draw_id), aes(x = Difficulty_bin, y = rt, group = draw), alpha = alpha, col = "orange")+
    geom_line(data = qq_summar_sum_rt_psy, aes(x = Difficulty_bin, y = rts,col = "CBM"), show.legend = F, linewidth = 1.1)+
    geom_line(data = qq_summar_sum_psy, aes(x = Difficulty_bin, y = rts, col = "DDM"), show.legend = F, linewidth = 1.1)+
    geom_pointrange(aes(y = rt, ymin = rt-2*se_rts, ymax = rt+2*se_rts, col = "Data"), show.legend = F)+
    scale_color_manual(name = "Model Predictions & Data", 
                       values = c("DDM" = "red", 
                                  "CBM" = "blue", 
                                  "Data" = "black"),
                       breaks = c("Data", "CBM", "DDM"))+
    theme_classic()+
    scale_x_continuous("Binned stimulus intensity", breaks = scales::pretty_breaks(n = 3))+
    scale_y_continuous("Response time (s)", breaks = c(0.4,0.7,1.0,1.3), labels = c("0.4","0.7","1.0","1.3"))+
    theme+text+
    theme(strip.background = element_blank(),
          strip.text = element_text(size=font_size),
          legend.position = "top",
          legend.box = "horizontal",
          legend.direction="horizontal",
          legend.justification='center',
          legend.key.height = unit(0.5, 'cm'),
          legend.key.width = unit(0.5, 'cm'),
          legend.spacing.x = unit(0.35, "cm"),
          plot.title = element_text(hjust = 0.5, size = font_size), #change legend key height
          legend.text = element_text(size=font_size_small),
          axis.text=element_text(size=font_size_small),
          axis.title=element_text(size=font_size),
          axis.line=element_line(size=axis_width),
          axis.ticks=element_line(size=axis_width),
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))
  
  rts_psycho_psy
  
  
  
  combined_plot = (psycho_psycho_psy | rw_rw) /  (rts_psycho_psy | rts_rw) +plot_layout(tag_level = 'new', guides = "collect")&
    theme(legend.position = "top",
          legend.direction="horizontal",
          legend.text = element_text(size = font_size),  # Increase legend text size
          legend.title = element_blank(),  # Increase title size if needed
          legend.key.size = unit(1, "cm"))& 
    guides(
      color = guide_legend(override.aes = list(
        linewidth = c(0.5, 2, 2)
      ))
    )
  
  combined_plot
  
  
  ggsave(here::here("Figures","plot_5.tiff"),combined_plot, dpi = 600,
         height = 6,width = 12)
  
  return(plot)
  
}
 
  load(here::here("Plots","workspace","Figure_5.RData"))
  
  combined_plot = (psycho_psycho_psy | rw_rw) /  (rts_psycho_psy | rts_rw) +plot_layout(tag_level = 'new', guides = "collect")&
    theme(legend.position = "top",
          legend.direction="horizontal",
          legend.text = element_text(size = font_size),  # Increase legend text size
          legend.title = element_blank(),  # Increase title size if needed
          legend.key.size = unit(1, "cm"))& 
    guides(
      color = guide_legend(override.aes = list(
        linewidth = c(0.5, 2, 2)
      ))
    )
  
  combined_plot
  
  
  ggsave(here::here("Figures","plot_5.tiff"),combined_plot, dpi = 600,
         height = 6,width = 12)
  
  
  return(plot)
  
  
  
  }

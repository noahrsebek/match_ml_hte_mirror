library(dplyr)
library(ggplot2)
# digging into graduation ML output
source('mltx_globals.R')


# -- reading in fullsample models
fullsample_models <- readr::read_rds("all_fullsample_cf_models.Rds")

# -- reading in splitsample models
grf_rds_files <- list.files()[stringr::str_detect(list.files(), 'all_splitsample_models_final')]


all_grf_splitsamples <- NULL

for (grffile in grf_rds_files){
  temp_grffile <- readRDS(grffile)
  all_grf_splitsamples <- all_grf_splitsamples %>% append(temp_grffile)
}


# initiate with the first iteration, and then loop through 2:all grf models
splitsample_models <- all_grf_splitsamples[[1]]$causal_forest_models

for (i in 2:length(all_grf_splitsamples)){
  current_splitsample_iteration <- all_grf_splitsamples[[i]]$causal_forest_models
  n_outcomes <- length(current_splitsample_iteration)
  
  # for each outcome...
  for (j in 1:n_outcomes){
    # get each dataframe
    n_tables <- length(current_splitsample_iteration[[j]])
    
    for (k in 1:n_tables){
      # rbind the kth table to the same table in the combined model
      splitsample_models[[j]][[k]] <- bind_rows(splitsample_models[[j]][[k]],
                                             current_splitsample_iteration[[j]][[k]])
    }
  }
  
}

splitsample_models_full <- splitsample_models

for (j in 1:n_outcomes){
  for (k in 1:n_tables){
    splitsample_models[[j]][[k]] <-
      splitsample_models[[j]][[k]] %>%
      group_by_at(1) %>%
      summarise_all(list(mean = mean#,
                         #sd = sd
      )) %>%
      ungroup()
  }
}


# making mean version of the coefficient plot
plot_data <- NULL

for (outcome in outcomes_of_interest_labels){
  
  # get fullsample calibration test
  plot_data <- plot_data %>%
    bind_rows(
      fullsample_models[[outcome]]$calibration_test %>%
        filter(term %in% 'differential.forest.prediction') %>% 
        dplyr::select(est=estimate, lo=lower_CI_full, hi=upper_CI_full) %>% 
        mutate(type = 'fullsample', outcome=outcome)
    )
  
  # get splitsample calibration test
  plot_data <- plot_data %>%
    bind_rows(
      splitsample_models[[outcome]]$calibration_test %>%
        filter(term %in% 'differential.forest.prediction') %>% 
        dplyr::select(est=estimate_mean, lo=lower_CI_mean, hi=upper_CI_mean) %>% 
        mutate(type = 'splitsample', outcome=outcome)
    )
}


full_het_plot <- plot_data %>%
  #filter(outcome %in% names(outcome_and_label_list)[which(outcome_and_label_list %in% main_outcomes)]) %>% 
  ggplot(aes(y=est,
             x=factor(outcome, levels = unique(outcome)),
             color=outcome,
             group=factor(type, levels=unique(type)))) +
  geom_point(position = position_dodge(0.4)) +
  geom_errorbar(aes(ymin=lo, ymax=hi),
                width=0.6,
                position = position_dodge(0.4)) +
  geom_hline(yintercept = 1, linetype='dashed', alpha=0.5) + 
  geom_hline(yintercept = 0, linetype='dashed', alpha=0.5) +
  scale_x_discrete(labels = scales::wrap_format(10)) + 
  ggtitle(label = 'Calibration Test: Detecting Heterogeneity',
          subtitle = paste0("For each outcome, the estimate and confidence intervals for the differential",
                            " forest prediction coefficient are shown for the full-sample (left) and split-sample procedure (right).",
                            " Horizontal lines at 1 and 0.")) +
  theme_bw() +   theme(legend.position='none') + 
  xlab('Outcome of interest') + ylab('Coefficient estimate')

print(full_het_plot)

# density plot of estimates + out fullsample estimate
fullsample_grad_estimate <- fullsample_models$`Graduated on-time`$calibration_test %>% 
  filter(term %in% 'differential.forest.prediction') %>% pull(estimate)

splitsample_estimates <- splitsample_models_full$`Graduated on-time`$calibration_test %>%
  filter(term %in% 'differential.forest.prediction')


splitsample_estimates %>% 
  ggplot(aes(x=estimate)) + geom_density() +
  geom_vline(xintercept = fullsample_grad_estimate, linetype='dashed') +
  ggtitle('Calibration test estimates from N=1000 runs of splitsample GRF',
          subtitle = paste0('Outcome: graduated on-time; vertical line at ',
                            round(fullsample_grad_estimate, 3), " (our main fullsample estimate)"))

mean(abs(splitsample_estimates$estimate) > abs(fullsample_grad_estimate))

mean(splitsample_estimates$estimate < fullsample_grad_estimate)

mean(splitsample_estimates$estimate > 0)

# rerun the fullsample thing many times for graduated on-time

all_fullsample_grad_models <- readr::read_rds('fullsample_graduation_1k_runs.Rds')



combine_models <- function(file){
  temp_file <- readr::read_rds(file)
  combined_models_df <- temp_file[[1]]
  for (i in 2:length(temp_file)){
    combined_models_df <- combined_models_df %>% bind_rows(temp_file[[i]])
  }
  combined_models_df
}


all_fullsample_grad_models <- combine_models('fullsample_graduation_1k_runs.Rds')



all_fullsample_grad_models %>%
  filter(term %in% 'differential.forest.prediction') %>% 
  ggplot(aes(x=estimate)) + geom_density() +
  geom_vline(xintercept = fullsample_grad_estimate, linetype='dashed') +
  ggtitle('Calibration test estimates from N=1000 runs of fullsample GRF',
          subtitle = paste0('Outcome: graduated on-time; vertical line at ',
                            round(fullsample_grad_estimate, 3), " (our main estimate)"))


mean(abs(combined_fullsample_models_df %>%
           filter(term %in% 'differential.forest.prediction') %>%
           pull(estimate)) > abs(fullsample_grad_estimate))

mean((combined_fullsample_models_df %>%
        filter(term %in% 'differential.forest.prediction') %>%
        pull(estimate)) < fullsample_grad_estimate)



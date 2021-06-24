# making calibration plots for an arbitrary outcome

# outcome_of_interest <- names(final_forests_missingness)[1]


make_calibration_plot <- function(forest, outcome_of_interest,
                                  master=master_pool){
  # -> get tau's
  tau_df <- forest$tau_df
  
  # -> add ventiles
  qcut <- function(x, n) {
    cut(x, quantile(x, seq(0, 1, length = n + 1)), labels = seq_len(n),
        include.lowest = TRUE)
  }
  
  n_calibration_quantiles <- 20
  
  calibration_df <- tau_df %>% mutate(calibration_quantiles = grf_tau_hat %>% qcut(n_calibration_quantiles))
  
  
  
  master <- master %>% left_join(calibration_df, by='sid') %>% mutate(blocknum = as_factor(blocknum))
  
  
  itt_vars_calibration <- c("d13andunder","d14","d15","d16","d17andover","dlearningdisabled","dfreelunch",
                            "dblack","dhispanic","dother","dgrade9","dgrade10","gpa_pre_zeros",
                            "numAs_pre","numBs_pre","numCs_pre","numDs_pre","numFs_pre","missing_gpa_pre",
                            "days_absent_pre_zeros","missing_attend_pre","mathxil_z_pre_np_zeros","readxil_z_pre_np_zeros",
                            "mathxil_z_pre_missing","readxil_z_pre_missing","oss_dis_pre_zeros","incidents_pre_zeros",
                            "any_arrests_pre","violent_pre","property_pre","drug_pre", "blocknum")
  
  # add 
  quantile_itts <- lm(paste0(outcome_of_interest, " ~ ",
                             "dmatch:calibration_quantiles + calibration_quantiles +",
                             paste(itt_vars_calibration, collapse='+')),
                      data=master) %>% broom::tidy() %>% filter(term %>% startsWith('dmatch:')) %>% 
    mutate( term = term %>% stringr::str_replace('dmatch:calibration_quantiles', "")) %>% 
    rename(calibration_quantiles = term)
  
  # get avg pte for each bucket
  avg_pte_df <- calibration_df %>% group_by(calibration_quantiles) %>% summarise(pte = mean(grf_tau_hat))
  
  
  
  
  # get the statistic from the linear test; this will be the plot subtitle
  calibration_test_statistic <- lm(paste0(outcome_of_interest, " ~ ",
                                          "dmatch + grf_tau_hat + grf_tau_hat*dmatch +",
                                          paste0(itt_vars_calibration, collapse="+")),
                                   data=master) %>% broom::tidy() %>% filter(term %in% "dmatch:grf_tau_hat")
  
  
  calibration_test_note <- paste0("PTE x Treatment interaction estimate is ",
                                  calibration_test_statistic$estimate %>% round(3), 
                                  " (", calibration_test_statistic$std.error %>% round(3), "), with a p-value of ",
                                  calibration_test_statistic$p.value %>% round(3))
  
  
  quantile_calibration_plot <- quantile_itts %>% left_join(avg_pte_df) %>%
    ggplot(aes(x=pte, y=estimate)) + 
    geom_point() +
    geom_smooth(method='lm', se=F, color='black', linetype = 'dashed') +
    geom_abline(intercept = 0, slope = 1) + 
    xlab("Average Quantile Individual PTE") + 
    ylab("ITT Estimate for Quantile") + 
    ggtitle(paste0("Calibration Plot: ", outcome_of_interest),
            subtitle = calibration_test_note)
  
  
  
  
  
  return(quantile_calibration_plot)
  
}


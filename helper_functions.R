

makeDummies <- function(dataframe,
                        flipped=F){
  dataframe$d_positive <- ifelse(dataframe$grf_tau_hat > 0, 1, 0)
  dataframe$d_negative <- ifelse(dataframe$grf_tau_hat <= 0, 1, 0)
  
  dataframe$d_positive <- ifelse(dataframe$grf_tau_hat > 0, 1, 0)
  dataframe$d_negative <- ifelse(dataframe$grf_tau_hat <= 0, 1, 0)
  
  #qcut function taken from Hadley Wickham answer on stackoverflow
  qcut <- function(x, n) {
    cut(x, quantile(x, seq(0, 1, length = n + 1)), labels = seq_len(n),
        include.lowest = TRUE)
  }
  
  dataframe$tau_quartile <- qcut(dataframe$grf_tau_hat, 4) %>% as.numeric()
  
  if (flipped ==T){
    dataframe$tau_quartile <- abs(4 - dataframe$tau_quartile) + 1
  }
  
  dummy_df <- data.frame(dummies::dummy(dataframe$tau_quartile, sep="_"))
  
  #note: dummy function has issues with column names; overwrite them manually here
  colnames(dummy_df) <- c("d_quart1", "d_quart2", "d_quart3", "d_quart4")
  dataframe <- cbind(dataframe, dummy_df)
  return(dataframe)
  
}



scalar1 <- function(x) {x / sqrt(sum(x^2, na.rm=T))}

add_pairwise_interactions <- function(data, list_of_controls){
  
  temp_df <- data
  
  for (control_1 in list_of_controls){
    for (control_2 in list_of_controls){
      # if control1==control2, pass
      if (control_1 == control_2) next
      temp_df[,paste0(control_1, "_x_", control_2)] <- scalar1(temp_df[,control_1]) * scalar1(temp_df[,control_2])
      
    }
  }
  
  return(temp_df)
  
}


get_group_avg_tau_conf_interval <- function(group1, group2){
  group1[["estimate"]] - group2[["estimate"]] +
    c(-1, 1) * qnorm(0.975) * sqrt(group1[["std.err"]]^2 + group2[["std.err"]]^2)
}

clean_estimate_w_95_conf_intervals <- function(est, se){
  paste0(round( est , 3), " +/- ", round( qnorm(0.975) * se , 3))
}


make_quartile_baseline_table <- function(input_augmented_df,
                                         subsample_tau_avgs,
                                         baselines = table_baselines,
                                         baselines_labels = table_baselines_labels){
  quartile_baseline_table <- list()
  
  tau_avgs <- subsample_tau_avgs
  tau_avg_df <- tribble(~tau_quartile, ~avg,
                        1, tau_avgs$quartile_1[['estimate']],
                        2, tau_avgs$quartile_2[['estimate']],
                        3, tau_avgs$quartile_3[['estimate']],
                        4, tau_avgs$quartile_4[['estimate']])
  
  for (quartile in 1:4){
    
    # get N in quartile and get avg tau hat 
    n_quart <- input_augmented_df %>% filter(tau_quartile %in% quartile) %>% nrow()
    avg_tau_hat <- input_augmented_df %>% filter(tau_quartile %in% quartile) %>% pull(grf_tau_hat) %>% mean(na.rm=T)
    #avg_tau_hat <- tau_avg_df %>% filter(tau_quartile %in% quartile) %>% pull(avg)
    
    # make column name
    col_name <- paste0("$\\hat{\\tau}$ Quartile ", quartile)
    temp_quart_column <- c(avg_tau_hat, n_quart)
    
    for (bl_var in baselines){
      mean_val <- input_augmented_df %>% filter(tau_quartile %in% quartile) %>% pull(bl_var) %>% mean(na.rm=T)
      temp_quart_column <- c(temp_quart_column, mean_val)
    }
    quartile_baseline_table[[col_name]] <- temp_quart_column
  }
  
  quartile_baseline_table <- quartile_baseline_table %>% as_tibble()
  quartile_baseline_table <- bind_cols(
    tibble(Baseline = c("\\textit{Mean} $\\hat{\\tau}$",
                        "\\textit{N}",
                        baselines_labels)),
    quartile_baseline_table)
  
  return(quartile_baseline_table)
}



make_ate_summary_table <- function(subsample_tau_avgs,
                                   group_tau_avgs){
  
  subsample_ATEs_mean_se <- subsample_tau_avgs
  
  quartile_95_ci_col <- c()
  quartile_se_col <- c()
  for (i in 4:1){
    current_quartile <- paste0('quartile_', i)
    current_ATE <- subsample_ATEs_mean_se[[current_quartile]]
    quartile_95_ci_col <- c(quartile_95_ci_col,
                            current_ATE[1])
    quartile_se_col <- c(quartile_se_col, current_ATE[[2]])
  }
  
  
  
  
  # Quartile
  ate_95ci_table_pte_quartiles <- tibble('Sample'  = paste("Individual PTE Quartile", 4:1),
                                         "Avg. Treatment Effect" = quartile_95_ci_col,
                                         "Standard Error" = quartile_se_col)
  
  # Overall
  ate_95ci_table_overall <- tibble("Sample" = "Whole Sample",
                                   "Avg. Treatment Effect" = group_tau_avgs['overall.estimate'],
                                   "Standard Error" = group_tau_avgs['overall.std.err'])
  
  # Other subgroups
  other_subgroups <- c("above_median", "below_median", "bottom_three_quartiles")
  other_subgroup_labels <- c("Top 2 PTE Quartiles", "Bottom 2 PTE Quartiles", "Bottom 3 PTE Quartiles")
  
  other_subgroup_95_ci_col <- c()
  other_subgroup_se_col <- c()
  
  for (group in other_subgroups){
    current_ATE <- subsample_ATEs_mean_se[[group]]
    other_subgroup_95_ci_col <- c(other_subgroup_95_ci_col,
                                  current_ATE[1])
    other_subgroup_se_col <- c(other_subgroup_se_col, current_ATE[2])
  }
  
  ate_95ci_table_subgroup <- tibble('Sample'  = other_subgroup_labels,
                                    "Avg. Treatment Effect" = other_subgroup_95_ci_col,
                                    "Standard Error" = other_subgroup_se_col)
  
  
  
  bind_rows(ate_95ci_table_overall,
            ate_95ci_table_pte_quartiles,
            ate_95ci_table_subgroup)
}


make_subsample_difference_table <- function(subsample_tau_avgs){
  
  highest_quartile_vs_rest_ci <- get_group_avg_tau_conf_interval(
    subsample_tau_avgs$highest_quartile,
    subsample_tau_avgs$bottom_three_quartiles)
  
  above_median_vs_below_ci <- get_group_avg_tau_conf_interval(
    subsample_tau_avgs$above_median, 
    subsample_tau_avgs$below_median)
  
  top_vs_bottom_ci <- get_group_avg_tau_conf_interval(
    subsample_tau_avgs$quartile_4,
    subsample_tau_avgs$quartile_1)
  
  
  
  tribble(
    ~"Group 1", ~"Group 2", ~"95\\% Confidence Interval of Difference (low)", ~"95\\% CI (high)",
    "Quartile 4", "Bottom 3 Quartiles", highest_quartile_vs_rest_ci[1], highest_quartile_vs_rest_ci[2],
    "Top 2 Quartiles (3 \\& 4)", "Bottom 2 Quartiles (1 \\& 2)", above_median_vs_below_ci[1], above_median_vs_below_ci[2],
    "Quartile 4", "Quartile 1", top_vs_bottom_ci[1], top_vs_bottom_ci[2]
  )
}

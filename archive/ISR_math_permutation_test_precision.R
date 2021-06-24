


# overall df: changing to question-by-student level dataset
question_by_student_level_df <- overall_df %>% 
  pivot_longer(cols = starts_with('Q_'), names_to ="var", values_to = 'response') %>%
  left_join(pplus_quartile_df)

p_plus_values <- sort(unique(question_by_student_level_df$p_plus))


#input_df <- question_by_student_level_df

  
get_TxC_ITT_difference_by_difficulty <- function(input_df){
  #creating formulas with the question
  itt_vars <- c("dmatch","d13andunder","d14","d15","d16","d17andover","dlearningdisabled","dfreelunch",
                "dblack","dhispanic","dother","dgrade9","dgrade10","gpa_pre_zeros",
                "numAs_pre","numBs_pre","numCs_pre","numDs_pre","numFs_pre","missing_gpa_pre",
                "days_absent_pre_zeros","missing_attend_pre","mathxil_z_pre_np_zeros","readxil_z_pre_np_zeros",
                "mathxil_z_pre_missing","readxil_z_pre_missing","oss_dis_pre_zeros","incidents_pre_zeros",
                "any_arrests_pre","violent_pre","property_pre","drug_pre")
  
  
  
  questionlevel_interaction_df <- tibble()
  
  
  
  for (difficulty in p_plus_values){
    
    # filter to just this level of difficulty
    temp_df <- input_df %>% filter(p_plus %in% difficulty) %>% drop_na(response)
    
    itt_formula <- as.formula(paste0("response ~ ", paste(c('above_median', 'dmatch_above_median', 'blocknum', itt_vars),
                                                            collapse=" + ")))
    
    
    # RUNNING ITT
    #print(paste0('running ITT ', counter))
    itt_model <- lm(itt_formula, data=temp_df, weights = isr_weight_14)
    itt_model_robust <- coeftest(itt_model, vcov=vcovHC(itt_model, type="HC1"))
    # itt_model_cluster <- summary(
    #   miceadds::lm.cluster(itt_formula, data=temp_df, cluster='sid',
    #         weights=temp_df$isr_weight_14))
    
    temp_row <- as.data.frame(matrix(ncol=4))
    colnames(temp_row) <- c("var",
                            # "n", "ctrl_mean",
                            # "itt", "itt_se", "itt_pval",
                            # "above_median", "above_median_se", "above_median_pval",
                            "treat_x_abovemed", "treat_x_abovemed_se", "treat_x_abovemed_pval")
    
    #adding question into
    temp_row$var <- difficulty
    temp_row$treat_x_abovemed <- itt_model_robust['dmatch_above_median', 'Estimate']
    temp_row$treat_x_abovemed_se <- itt_model_robust['dmatch_above_median', 'Std. Error']
    temp_row$treat_x_abovemed_pval <- itt_model_robust['dmatch_above_median', "Pr(>|t|)"]
    
    questionlevel_interaction_df <- bind_rows(questionlevel_interaction_df, temp_row)
    
  }
  
  
  individual_question_itt_difference_df <-
    questionlevel_interaction_df %>%
    mutate(hi_lo_abs_diff = abs(treat_x_abovemed)) %>% 
    dplyr::select(p_plus = var, hi_lo_abs_diff)
  
  return(individual_question_itt_difference_df)
  
}

get_TxC_ITT_difference_by_difficulty(question_by_student_level_df)




real_ITT_diffs_precision <- get_TxC_ITT_difference_by_difficulty(question_by_student_level_df)



# permute the above/below dummy (and recalc the interaction term)
run_permuted_TxC_ITT_diff_precision <- function(i){
  permuted_df <- suppressMessages(question_by_student_level_df %>%
                                    group_by(blocknum) %>%
                                    mutate(above_median = sample(above_median)) %>% 
                                    ungroup %>% 
                                    mutate(dmatch_above_median = dmatch * above_median))
  
  
  permuted_output <- suppressMessages(get_TxC_ITT_difference_by_difficulty(permuted_df))
  
  return(permuted_output)
}




set.seed(20210115)
t_pre <- Sys.time()
permuted_ITT_diffs <- lapply(1:10000, run_permuted_TxC_ITT_diff_precision)
t_post <- Sys.time()
t_post - t_pre
permuted_ITT_diffs %>% readr::write_rds('isr_permutation_test/itt_diffs_permuted.rds')

get_statistic <- function(x){x %>% pull(hi_lo_abs_diff) %>% mean}

overall_true_stat <- real_ITT_diffs_precision %>% get_statistic()
permutation_ensemble <- sapply(permuted_ITT_diffs, get_statistic)

mean(permutation_ensemble > overall_true_stat)


# using weighted avgs ---

p_plus_weights <- overall$data %>%
  dplyr::select(var, p_plus) %>%
  group_by(p_plus) %>%
  summarise(weight=n()) %>%
  ungroup


get_weighted_avg_diff <- function(x){
  temp_x <- x %>% left_join(p_plus_weights, by='p_plus')
  weighted.mean(temp_x$hi_lo_abs_diff, temp_x$weight)
}

# get overall statistic for 'true' and permuted
true_summary_stat <- get_weighted_avg_diff(real_ITT_diffs_precision)
permuted_summary_stats <- sapply(permuted_ITT_diffs, get_weighted_avg_diff)

mean(permuted_summary_stats >= true_summary_stat)

# extra r code for ISR math test: permutation test
get_TxC_ITT_difference_by_question <- function(overall_df){
  #creating formulas with the question
  itt_vars <- c("dmatch","d13andunder","d14","d15","d16","d17andover","dlearningdisabled","dfreelunch",
                "dblack","dhispanic","dother","dgrade9","dgrade10","gpa_pre_zeros",
                "numAs_pre","numBs_pre","numCs_pre","numDs_pre","numFs_pre","missing_gpa_pre",
                "days_absent_pre_zeros","missing_attend_pre","mathxil_z_pre_np_zeros","readxil_z_pre_np_zeros",
                "mathxil_z_pre_missing","readxil_z_pre_missing","oss_dis_pre_zeros","incidents_pre_zeros",
                "any_arrests_pre","violent_pre","property_pre","drug_pre")
  
  
  
  questionlevel_interaction_df <- tibble()
  for (question in mathitem_cols){
    
    itt_formula <- as.formula(paste0(question, " ~ ", paste(c('above_median', 'dmatch_above_median', 'blocknum', itt_vars),
                                                            collapse=" + ")))
    
    
    # RUNNING ITT
    #print(paste0('running ITT ', counter))
    itt_model <- lm(itt_formula, data=overall_df[which(!is.na(overall_df[,question])),], weights = isr_weight_14)
    
    itt_model_robust <- coeftest(itt_model, vcov=vcovHC(itt_model, type="HC1"))
    
    temp_row <- as.data.frame(matrix(ncol=4))
    colnames(temp_row) <- c("var",
                            # "n", "ctrl_mean",
                            # "itt", "itt_se", "itt_pval",
                            # "above_median", "above_median_se", "above_median_pval",
                            "treat_x_abovemed", "treat_x_abovemed_se", "treat_x_abovemed_pval")
    
    #adding question into
    temp_row$var <- question
    # temp_row$n <- sum(!is.na(overall_df[,question]))
    # temp_row$ctrl_mean <- mean(overall_df %>% filter(dmatch==0) %>% pull(question),na.rm=T)
    # temp_row$itt      <- itt_model_robust['dmatch','Estimate']
    # temp_row$itt_se   <- itt_model_robust['dmatch','Std. Error']
    # temp_row$itt_pval <- itt_model_robust['dmatch',"Pr(>|t|)"]
    # temp_row$above_median <- itt_model_robust['above_median', 'Estimate']
    # temp_row$above_median_se <- itt_model_robust['above_median', 'Std. Error']
    # temp_row$above_median_pval <- itt_model_robust['above_median', "Pr(>|t|)"]
    temp_row$treat_x_abovemed <- itt_model_robust['dmatch_above_median', 'Estimate']
    temp_row$treat_x_abovemed_se <- itt_model_robust['dmatch_above_median', 'Std. Error']
    temp_row$treat_x_abovemed_pval <- itt_model_robust['dmatch_above_median', "Pr(>|t|)"]
    
    questionlevel_interaction_df <- bind_rows(questionlevel_interaction_df, temp_row)
    
  }
  
  
  individual_question_itt_difference_df <-
    questionlevel_interaction_df %>%
    mutate(hi_lo_abs_diff = abs(treat_x_abovemed)) %>% 
    dplyr::select(var, hi_lo_abs_diff) %>% 
    left_join(overall$data %>% dplyr::select(var, p_plus)) %>% 
    group_by(p_plus) %>% summarise(hi_lo_abs_diff = mean(hi_lo_abs_diff))
  
  return(individual_question_itt_difference_df)
    
}



# getting real (average absolute) differences for each P-plus value
real_ITT_diffs <- get_TxC_ITT_difference_by_question(overall_df)



# permute the above/below dummy (and recalc the interaction term)
run_permuted_TxC_ITT_diff <- function(i){
  permuted_df <- suppressMessages(overall_df %>%
    group_by(blocknum) %>%
    mutate(above_median = sample(above_median)) %>% 
    ungroup %>% 
    mutate(dmatch_above_median = dmatch * above_median))
  
  
  permuted_output <- suppressMessages(get_TxC_ITT_difference_by_question(permuted_df))
  
  return(permuted_output)
}
  



set.seed(20210104)
t_pre <- Sys.time()
permuted_ITT_diffs <- lapply(1:10000, run_permuted_TxC_ITT_diff)
t_post <- Sys.time()
t_post - t_pre


xil_95 <- function(x){quantile(x, 0.95)}
xil_05 <- function(x){quantile(x, 0.05)}

permutation_summary_df <- 
  permuted_ITT_diffs %>% bind_rows() %>% 
  group_by(p_plus) %>%
  summarize_all(.funs = list('mean'=mean, 'hi'=xil_95, 'lo'=xil_05))

# making p-plus plot (distribution of perm test against 'true' differences at each p-plus value)
permutation_summary_df %>% 
ggplot(aes(x=p_plus)) + geom_point(aes(y=mean), alpha=0.5) +
  geom_errorbar(aes(ymin=lo, ymax=hi), alpha=0.5) + 
  geom_point(aes(y=hi_lo_abs_diff), data = real_ITT_diffs, color='red', size=2)



# making weighted average for the overall statistic
# -> we have the df at the p-plus level; we want an overall *question-level* statistic
# -> here the 'weights' are just the counts of questions at a given  p-plus value

p_plus_weights <- overall$data %>%
  select(var, p_plus) %>%
  group_by(p_plus) %>%
  summarise(weight=n()) %>%
  ungroup


get_weighted_avg_diff <- function(x){
  temp_x <- x %>% left_join(p_plus_weights, by='p_plus')
  weighted.mean(temp_x$hi_lo_abs_diff, temp_x$weight)
}

# get overall statistic for 'true' and permuted
true_summary_stat <- get_weighted_avg_diff(real_ITT_diffs)
permuted_summary_stats <- sapply(permuted_ITT_diffs, get_weighted_avg_diff)

# distribution of permutation estimates
ggplot(data=tibble(statistic = permuted_summary_stats), aes(x=statistic)) +
  geom_density() + geom_vline(xintercept = true_summary_stat)


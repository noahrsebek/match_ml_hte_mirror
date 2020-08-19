itt_tot_analysis <- function(input_dataset,
                             outcome_list,
                             outcome_labels=NULL,
                             block_or_school='block',
                             extra_baselines = NULL,
                             robust=T,
                             cluster=NA){
  library(lmtest)
  library(sandwich)
  library(ivpack)
  
  itt_vars <- c("dmatch","d13andunder","d14","d15","d16","d17andover","dlearningdisabled","dfreelunch",
                "dblack","dhispanic","dother","dgrade9","dgrade10","gpa_pre_zeros",
                "numAs_pre","numBs_pre","numCs_pre","numDs_pre","numFs_pre","missing_gpa_pre",
                "days_absent_pre_zeros","missing_attend_pre","mathxil_z_pre_np_zeros","readxil_z_pre_np_zeros",
                "mathxil_z_pre_missing","readxil_z_pre_missing","oss_dis_pre_zeros","incidents_pre_zeros",
                "any_arrests_pre","violent_pre","property_pre","drug_pre")
  
  tot_vars_1 <-c("treat_post1","d13andunder","d14","d15","d16","dlearningdisabled",
                 "dfreelunch","dblack","dhispanic","dgrade9","gpa_pre_zeros","numAs_pre","numBs_pre",
                 "numCs_pre","numDs_pre","numFs_pre","missing_gpa_pre","days_absent_pre_zeros",
                 "missing_attend_pre","mathxil_z_pre_np_zeros","readxil_z_pre_np_zeros",
                 "mathxil_z_pre_missing","readxil_z_pre_missing","oss_dis_pre_zeros","incidents_pre_zeros",
                 "any_arrests_pre","violent_pre","property_pre","drug_pre","other_pre")
  
  tot_vars_2 <- c("dmatch","d13andunder","d14","d15","d16","dlearningdisabled","dfreelunch",
                  "dblack","dhispanic","dgrade9","gpa_pre_zeros","numAs_pre","numBs_pre","numCs_pre",
                  "numDs_pre","numFs_pre","missing_gpa_pre","days_absent_pre_zeros","missing_attend_pre",
                  "mathxil_z_pre_np_zeros","readxil_z_pre_np_zeros","mathxil_z_pre_missing",
                  "readxil_z_pre_missing","oss_dis_pre_zeros","incidents_pre_zeros","any_arrests_pre",
                  "violent_pre","property_pre","drug_pre","other_pre")
  
  if (block_or_school == 'block'){
    itt_vars <- c(itt_vars, 'blocknum')
    tot_vars_1 <- c(tot_vars_1, 'blocknum')
    tot_vars_2 <- c(tot_vars_2, 'blocknum')
  }
  if (block_or_school == 'school'){
    itt_vars <- c(itt_vars, 'schlid')
    tot_vars_1 <- c(tot_vars_1, 'schlid')
    tot_vars_2 <- c(tot_vars_2, 'schlid')
  }
  
  
  if (!is.null(extra_baselines)){
    itt_vars <- c(itt_vars, extra_baselines)
    tot_vars_1 <- c(tot_vars_1, extra_baselines)
    tot_vars_2 <- c(tot_vars_2, extra_baselines)
  }
  
  master = input_dataset
  
  master$blocknum <- as.factor(master$blocknum)
  master$schlid <- as.factor(master$schlid)
  
  # defining CCM within this namespace
  ccm <- function(outcome) {
    pA <- master %>% filter(dmatch==0) %>% pull(treat_post1) %>% mean(na.rm=T)
    pB <- master %>% filter(dmatch==1) %>% pull(treat_post1) %>% mean(na.rm=T)
    pN <- 1-pB
    pC <- 1-pA-pN
    
    EYAC <- master %>% filter(dmatch==1, treat_post1==1) %>% pull(outcome) %>% mean(na.rm=T)
    EYA <- master %>% filter(dmatch==0, treat_post1==1) %>% pull(outcome) %>% mean(na.rm=T)
    EYA[is.na(EYA)] <- 0
    EY <- (EYAC*(1-pN)-EYA*pA)/pC
    EY
  }
  
  
  out_df <- as.data.frame(matrix(ncol=10, nrow=0))
  colnames(out_df) <- c("var", "n", "ctrl_mean", "itt", "itt_se", "itt_pval", "tot", "tot_se", "tot_pval", "ccm")
  
  final_outcome_list <- c()
  final_outcome_labels <- c()
  i = 0
  # running analysis on each outcome of interest
  for (outcome in outcome_list){ # for each outcome...
    
    i = i + 1
    
    temp_row <- as.data.frame(matrix(ncol=10))
    colnames(temp_row) <- c("var", "n", "ctrl_mean", "itt", "itt_se", "itt_pval", "tot", "tot_se", "tot_pval", "ccm")
    
    #adding outcome into
    temp_row$var <- outcome
    temp_row$n <- sum(!is.na(master[,outcome]))
    temp_row$ctrl_mean <- master %>% filter(dmatch==0) %>% pull(outcome) %>% mean(na.rm=T)
    
    
    #creating formulas with the outcome
    itt_formula <- as.formula(paste0(outcome, " ~ ", paste(itt_vars, collapse=" + ")))
    tot_formula <- as.formula(paste0(outcome, " ~ ", 
                                     paste(tot_vars_1, collapse=" + "), "|",
                                     paste(tot_vars_2, collapse=" + ")))
    
    # RUNNING ITT
    if (nrow(master[which(!is.na(master[,outcome])),]) > 0){
      
      if (is.na(cluster)){
        itt_model <- lm(itt_formula, data=master[which(!is.na(master[,outcome])),])
      } else {
        cluster_var <- cluster
        itt_model <- miceadds::lm.cluster(itt_formula,
                                          data=master[which(!is.na(master[,outcome])),],
                                          cluster=cluster_var)
        
      }
      
    } else {
      
      # remove 
      
      next}
    
    # add outcome and label to the final list if we've made it this far
    final_outcome_list <- c(final_outcome_list, outcome)
    final_outcome_labels <- c(final_outcome_labels, outcome_labels[i])
    
    
    if (robust==TRUE & is.na(cluster)){
      itt_model_robust <- coeftest(itt_model, vcov=vcovHC(itt_model, type="HC1"))
      temp_row$itt      <- itt_model_robust['dmatch','Estimate']
      temp_row$itt_se   <- itt_model_robust['dmatch','Std. Error']
      temp_row$itt_pval <- itt_model_robust['dmatch',"Pr(>|t|)"]
    }
    
    if (robust==FALSE & is.na(cluster)){
      temp_row$itt      <- coef(summary(itt_model))['dmatch','Estimate']
      temp_row$itt_se   <- coef(summary(itt_model))['dmatch','Std. Error']
      temp_row$itt_pval <- coef(summary(itt_model))['dmatch',"Pr(>|t|)"]
    }
    
    if (!is.na(cluster)){
      temp_row$itt      <- summary(itt_model)['dmatch','Estimate']
      temp_row$itt_se   <- summary(itt_model)['dmatch','Std. Error']
      temp_row$itt_pval <- summary(itt_model)['dmatch',"Pr(>|t|)"]
    }
    
    # RUNNING TOT
    #print(paste0('running TOT ', counter))
    tot_model <- AER::ivreg(tot_formula, data=master[which(!is.na(master[,outcome])),])
    
    
    if (!is.na(cluster)){
      cluster_var <- unlist(master[which(!is.na(master[,outcome])),cluster])
      tot_model_cluster <- cluster.robust.se(tot_model, cluster_var)
      
      temp_row$tot <- tot_model_cluster['treat_post1','Estimate']
      temp_row$tot_se <- tot_model_cluster['treat_post1','Std. Error']
      temp_row$tot_pval <- tot_model_cluster['treat_post1','Pr(>|t|)'] 
    }
    
    
    else if (robust==TRUE){
      tot_model_robust <- robust.se(tot_model)
      temp_row$tot <- tot_model_robust['treat_post1','Estimate']
      temp_row$tot_se <- tot_model_robust['treat_post1','Std. Error']
      temp_row$tot_pval <- tot_model_robust['treat_post1','Pr(>|t|)']    
    }
    
    else if (robust==FALSE){
      temp_row$tot <- coef(summary(tot_model))['treat_post1','Estimate']
      temp_row$tot_se <- coef(summary(tot_model))['treat_post1','Std. Error']
      temp_row$tot_pval <- coef(summary(tot_model))['treat_post1','Pr(>|t|)']     
    }
    
    # adding CCM using output from TOT
    #  note: output from CCM
    temp_row$ccm <- ccm(outcome) - temp_row$tot
    
    
    # adding to final dataframe
    out_df <- rbind(out_df, temp_row)
  }
  
  #print("done with main loop")
  row.names(out_df) <- NULL
  
  
  
  
  
  # LABELLING COLUMNS
  colnames(out_df) <- c("Outcome", "N", "Ctrl Mean", "ITT", "ITT SE", "ITT P-val", "TOT", "TOT SE", "TOT P-val", "CCM")
  
  if (!is.null(outcome_labels)){
    out_df[,"Outcome"] <- final_outcome_labels
  }
  
  
  # # ADDING STARS
  # out_df <- out_df %>% mutate(ITT = round(ITT,3),
  #                             TOT = round(TOT, 3),
  #                             ITT_stars = case_when(`ITT P-val` >= 0.1 ~ "",
  #                                                   `ITT P-val` < 0.1 & `ITT P-val` >= 0.05 ~ "*",
  #                                                   `ITT P-val` < 0.05 & `ITT P-val` >= 0.01 ~ "**",
  #                                                   `ITT P-val` < 0.01 ~ "***"),
  #                             TOT_stars = case_when(`TOT P-val` >= 0.1 ~ "",
  #                                                   `TOT P-val` < 0.1 & `TOT P-val` >= 0.05 ~ "*",
  #                                                   `TOT P-val` < 0.05 & `TOT P-val` >= 0.01 ~ "**",
  #                                                   `TOT P-val` < 0.01 ~ "***"),
  #                             ITT = paste0(ITT, ITT_stars),
  #                             TOT = paste0(TOT, TOT_stars)) %>% dplyr::select(-ITT_stars, -TOT_stars)
  
  return(out_df)
}

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
  
  
  # -> get ITT (and avg PTE) for each bucket
  
  quantile_calibration_plot_df <- NULL
  
  for (i in 1:n_calibration_quantiles){
    # get students in that quantile
    subset_of_youth <- calibration_df %>% filter(calibration_quantiles==i) %>% pull(sid)
    
    # get their average PTE
    quantile_pte_avg <- calibration_df %>% filter(calibration_quantiles==i) %>% pull(grf_tau_hat) %>% mean()
    
    
    # run robust ITT (without covariates), get estimate + SE
    temp_master_dataset <- master %>% filter(sid %in% subset_of_youth)
    
    new_quantile_row <- temp_master_dataset %>% itt_tot_analysis(outcome_and_label_list[[outcome_of_interest]]) %>%
      select(N, ITT, `ITT SE`) %>% 
      mutate("Avg PTE" = quantile_pte_avg,
             "Quantile" = i)
    
    quantile_calibration_plot_df <- bind_rows(quantile_calibration_plot_df, new_quantile_row)
    
    
  }
  
  
  quantile_calibration_plot <- quantile_calibration_plot_df %>%
    ggplot(aes(x=`Avg PTE`, y=ITT)) +
    geom_point() + geom_smooth(method='lm', se=F, color='black', linetype = 'dashed') +
    geom_abline(intercept = 0, slope = 1) + 
    xlab("Average Quantile Individual PTE") + 
    ylab("ITT Estimate for Quantile") + 
    ggtitle(paste0("Calibration Plot: ", outcome_of_interest))
  
  
  return(quantile_calibration_plot)
  
}


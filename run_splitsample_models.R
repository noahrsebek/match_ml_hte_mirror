
# OVERALL DATA WORK ----
# - adding inv probability weights
# - rand block mean imputation
# - 
block_treatment_probability_df <- master_pool %>% group_by(blocknum, study) %>% summarize(p_block = mean(dmatch)) %>% ungroup()

singleton_blocks <- master_pool %>% group_by(blocknum, study) %>% filter(n()==1) %>% ungroup() %>% pull(blocknum)

master_pool <- master_pool %>%
  left_join(block_treatment_probability_df, by=c('study', 'blocknum')) %>% 
  mutate(inv_prob_weight = dmatch/p_block + (1-dmatch)/(1-p_block)) %>% 
  filter(! blocknum %in% singleton_blocks)


# # we want each core to run ALL outcomes
# for each iteration:
#   randomly split sample (block/treatment stratfication)
#   for each method:
#     for each outcome:
#       run method(outcome)  
  


seeds <- round(1000000*runif(n_splitsample_iterations))

# deduping and prepping data for split
input_dataset <- master_pool
input_controls <- controls_sans_missingness_dummies


duplicate_sids <- input_dataset %>% 
  group_by(sid) %>% filter(n()>1) %>% 
  distinct(sid) %>% pull(sid)

input_dataset <- input_dataset %>% filter(!sid %in% duplicate_sids)

input_dataset <- input_dataset %>%  drop_na(inv_prob_weight) %>%  mutate(row_id = row_number())


# running models ----

seeds <- round(1000000*runif(n_splitsample_iterations))


for (i in seeds){
  
  # set seed
  set.seed(i)
  
  estimation_df <- input_dataset %>% group_by(blocknum, dmatch) %>% sample_frac(splitsample_fraction)
  holdout_df <- input_dataset %>% anti_join(estimation_df, by='row_id') 
  
  # run cf for each outcome
  causal_forest <- make_causal_forest_for_each_outcome(estimation_df,
                                                       splitsample_df = holdout_df,
                                                       input_controls,
                                                       outcomes = outcomes_of_interest[1:2],
                                                       outcomes_of_interest_labels[1:2],
                                                       flipped_outcomes)
  
  # run xrf for each outcome
  
  
  # run BART method(s)
  
  
}






# generate list of random seeds

# loading packages and data ----

library('dplyr')
library('readr')
library('grf')
library('tidyr')
library('magrittr')
library('ggplot2')
library('furrr')
library('purrr')


n_splitsample_iterations <- 125
splitsample_fraction <- 0.5
default_n_trees <- 5000
to_cluster <- F
drop_duplicates <- T



bivariate_pte_comparisons <- tribble(
  ~outcome1, ~outcome2,
  "mathxil_z_post1_np", "mathgpa_post1",
  "mathxil_z_post1_np", "mathfail_post1",
  "mathxil_z_post1_np", "mathfailpercent_post1",
  "mathxil_z_post1_np", "nonmathgpa_post1",
  "mathgpa_post1", "nonmathgpa_post1",
  "mathfail_post1", "nonmathgpa_post1",
  "mathfailpercent_post1", "nonmathgpa_post1"
)


# load data ----
setwd("/export/projects/migrated-BAM/Match_AnalysisFiles_May2016/Results/2020/grf")
source(file = 'helper_functions.R')
source(file = 'grf_globals.R')

# overriding flipped outcomes (for the time being)
flipped_outcomes <- c()



master_pool <- load_master_dataset()


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




# deduping and prepping data for split
duplicate_sids <- master_pool %>% 
  group_by(sid) %>% filter(n()>1) %>% 
  distinct(sid) %>% pull(sid)

master_pool <- master_pool %>% filter(!sid %in% duplicate_sids)

master_pool <- master_pool %>%  drop_na(inv_prob_weight) %>%  mutate(row_id = row_number())

# editing our code to it runs with BART ----
library(bcf)

dataset <- master_pool
controls <- controls_sans_missingness_dummies
outcome = outcomes_of_interest[1]
n_burn = 2000
n_sim = 2000


# keep observations complete on outcome and that has an inv prob weight
working_df <- dataset %>% drop_na(all_of(outcome)) #%>%  mutate(row_id = row_number())

working_df <- working_df[c(controls, 'inv_prob_weight', outcome, 'dmatch')] %>% drop_na()

n.obs <- nrow(working_df)

X <- working_df[,c(controls, 'inv_prob_weight')] %>% as.matrix()
W <- working_df %>% pull(dmatch)
Y <- working_df %>% pull(tidyselect::all_of(outcome))
sample_weights <- working_df %>% pull(inv_prob_weight)


pi_hat <- lm(dmatch ~ ., data=working_df[c(controls, 'dmatch')]) %>% predict()

bcf_fit <- bcf(Y, W, X, X, pi_hat, nburn=n_burn, nsim=n_sim)
tau.hat <- bcf_fit$tau %>% colMeans()


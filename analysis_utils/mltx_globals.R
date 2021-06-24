

load_master_dataset <- function(addvars=T,
                                load_cached = T,
                                new_grade_vars=T){
  require(dplyr)
  
  if (load_cached == T){
    master_pool <- readr::read_csv("/export/projects/migrated-BAM/Match_AnalysisFiles_May2016/Results/2020/grf/master_dataset_cached.csv") %>% 
      mutate(baseline_mathxil_inschool_decile_zeros = ifelse(is.na(baseline_mathxil_inschool_decile), 0, baseline_mathxil_inschool_decile),
             missing_baseline_withinschool_mathtest = ifelse(is.na(baseline_mathxil_inschool_decile), 1, 0))
  }  else {
    master_pool <- readr::read_csv("/export/projects/migrated-BAM/Match_AnalysisFiles_May2016/SourceFiles/2017/analysisdata_unique_20180111_SS_alldupes.csv")
    master_pool <- master_pool %>%
      filter(pilotblocks %in% c(0, NA)) %>%
      distinct(sid, dmatch, study, randomized_date, .keep_all=T) %>% 
      mutate(mathfailpercent_post1 = mathfail_post1 / math_courses_post1)
    
    
    if (addvars==T){
      
      # add graduation outcomes ----
      graduation_df <- readr::read_csv('/export/projects/migrated-BAM/Match_AnalysisFiles_May2016/Results/2020/graduation_df.csv') %>% distinct()
      grad_outcomes <- c('graduated_ontime', 'graduated_ever',
                         'd_dropout', 'd_corrections', 'd_transfer', 'd_unknown',
                         'd_noexit', 'd_deceased')
      grad_outcome_labels <- c("Graduated On-Time", "Ever Graduated",
                               "Leave code: Dropout", "Leave code: Corrections", "Leave code: Transfer", "Leave code: Unknown",
                               "Leave code: No exit", "Leave code: Deceased")
      master_pool <- master_pool %>% left_join(graduation_df[,c('sid', grad_outcomes)], by='sid')
      
      # making sure 'treated in year 2 [of study 1]' is NA for 
      master_pool <- master_pool %>% mutate(treat_post2 = ifelse(study2==1, NA, treat_post2))
      
      
      # if (drop_duplicates==T){
      #   # get duplicate sids
      #   duplicate_sids <- master_pool %>% 
      #     group_by(sid) %>% filter(n()>1) %>% 
      #     distinct(sid) %>% pull(sid)
      #   
      #   master_pool <- master_pool %>% filter(!sid %in% duplicate_sids)
      # }
      source('/export/projects/migrated-BAM/Match_AnalysisFiles_May2016/Results/2020/get_11th_grade_outcomes.R')
      
      master_w_11th_grade_outcomes <- get_11th_grade_outcomes() %>% select(sid, gpa11_math, eleventh_grade_math_z, study) %>% distinct()
      
      
      master_pool <- master_pool %>% left_join(master_w_11th_grade_outcomes, by=c('sid', 'study'))}
   
    
    if (new_grade_vars==T){
      new_grade_vars <- read_csv("/export/projects/migrated-BAM/Match_AnalysisFiles_May2016/Results/2020/new_grade_variables.csv")
      new_grade_vars <- new_grade_vars
      master_pool <- master_pool %>% left_join(new_grade_vars, by=c('sid', 'study'))
    }
     
  }
  return(master_pool)
  
}

# --> outcomes of interest (and labels)
outcomes_of_interest <- c('mathxil_z_post1_np',
                          "math_gpa_full",
                          "math_failures_full",
                          "math_failures_full_percent",
                          #"mathfailpercent_post1",
                          "nonmathcore_gpa_all",
                          "nonmathcore_fails_all",
                          #"nonmathcore_fail_all_percent",
                          "nonmathcore_fails_percent_all",
                          "nonmathcore_gpa_bytopic_highgrade",
                          "nonmathcore_gpa_bytopic_lowgrade",
                          "nonmathcore_gpa_topthree_eachsem",
                          "nonmathcore_gpa_topsix",
                          "graduated_ontime",
                          "graduated_ever",
                          "gpa11_math",
                          "eleventh_grade_math_z")

outcomes_of_interest_labels <- c('Math test score (Z)',
                                 'Math GPA',
                                 "Math Course Failures",
                                 "Percent of Math Courses Failed",
                                 "Non-Math Core GPA (all non-math core courses)",
                                 "Non-Math Core Failures (all non-math core courses)",
                                 "Non-Math Core Percent Failed (all non-math core courses)",
                                 "Non-Math Core GPA (highest grade in each nonmath core topic each sem)",
                                 "Non-Math Core GPA (lowest grade in each nonmath core topic each sem)",
                                 "Non-Math Core GPA (top 3 nonmath core classes each sem)",
                                 "Non-Math Core GPA (top 6 nonmath core classes during year)",
                                 "Graduated on-time",
                                 "Graduated ever",
                                 "11th Grade Math GPA",
                                 "11th Grade Math Test Score (Z)")


additional_table_outcomes <- c("incidents_post1", "days_absent_post1", "oss_dis_post1")
additional_table_outcome_labels <- c("Disciplinary Incidents", "Days Absent", "Out-of-School Suspensions")
additional_table_outcome_label_list <- additional_table_outcomes
names(additional_table_outcome_label_list) <- additional_table_outcome_labels



# outcomes_of_interest <- c('mathxil_z_post1_np',
#                           "math_gpa_full",
#                           "math_failures_full",
#                           #"mathfailpercent_post1",
#                           "nonmathcore_gpa_all",
#                           "nonmathcore_fails_all",
#                           "nonmathcore_gpa_bytopic_highgrade",
#                           "nonmathcore_gpa_bytopic_lowgrade",
#                           "nonmathcore_gpa_topthree_eachsem",
#                           "nonmathcore_gpa_topsix")
# 
# outcomes_of_interest_labels <- c('Math test score (Z)',
#                                  'Math GPA',
#                                  "Math Course Failures",
#                                  "Non-Math Core GPA (all non-math core courses)",
#                                  "Non-Math Core Failures (all non-math core courses)",
#                                  "Non-Math Core GPA (highest grade in each nonmath core topic each sem)",
#                                  "Non-Math Core GPA (lowest grade in each nonmath core topic each sem)",
#                                  "Non-Math Core GPA (top 3 nonmath core classes each sem)",
#                                  "Non-Math Core GPA (top 6 nonmath core classes during year)")


# these are the outcomes we choose to display in the tables; we still run everything with all outcomes listed above
#   so we can make our coefficient plot
main_outcomes <- c('mathxil_z_post1_np',
                   "math_gpa_full",
                   "math_failures_full",
                   "nonmathcore_gpa_all",
                   "nonmathcore_fails_all",
                   "nonmathcore_gpa_bytopic_highgrade",
                   "nonmathcore_gpa_bytopic_lowgrade"#,
                   # "graduated_ontime",
                   # "graduated_ever"#,
                   # "nonmathcore_gpa_topthree_eachsem",
                   # "nonmathcore_gpa_topsix"
                   )



outcome_and_label_list <- outcomes_of_interest
names(outcome_and_label_list) <- outcomes_of_interest_labels


# list of controls/covariates/features to use in the model
# - these are all the 'standard' covariates we include in the ITT/TOT models
# - ... buuut without imputing 0's for the missing values
# - we still include missing dummies, however
controls_with_missingness_dummies <- c(
  #"age_pre",
  "d13andunder","d14","d15","d16","d17andover",
  "dlearningdisabled","dfreelunch",
  "dblack","dhispanic","dother",
  "dgrade9","dgrade10",
  "gpa_pre_zeros", "numAs_pre","numBs_pre","numCs_pre","numDs_pre","numFs_pre","missing_gpa_pre",
  "dfemale",
  "days_absent_pre_zeros","missing_attend_pre",
  "mathxil_z_pre_np_zeros","readxil_z_pre_np_zeros",
  "mathxil_z_pre_missing","readxil_z_pre_missing",
  "oss_dis_pre_zeros","incidents_pre_zeros",
  "any_arrests_pre","violent_pre","property_pre","drug_pre",
  'baseline_mathxil_inschool_decile_zeros',
  'missing_baseline_withinschool_mathtest')





# --> baseline variables to show in the quartile table
table_baselines <- c(#"d13andunder", "d14", "d15", "d16", "d17andover",
  "age_pre",
  'dfemale',
  "dlearningdisabled",
  "dfreelunch",
  "dblack",
  "dhispanic",
  "dother",
  "dgrade9",
  "dgrade10",
  "gpa_pre",
  "numAs_pre",
  "numBs_pre",
  "numCs_pre",
  "numDs_pre",
  "numFs_pre",
  "missing_gpa_pre",
  "days_absent_pre",
  "missing_attend_pre",
  "mathxil_z_pre_np", "readxil_z_pre_np",
  "mathxil_z_pre_missing", "readxil_z_pre_missing",
  "oss_dis_pre", "incidents_pre",
  "any_arrests_pre", "violent_pre", "property_pre", "drug_pre",
  'baseline_mathxil_inschool_decile',
  'study2',
  'treat_post1')

table_baselines_labels <- c("Age",
                            'Female',
                            "Has IEP",
                            "Has Free/Reduced Lunch",
                            "Black",
                            "Hispanic",
                            "Other Race",
                            "In 9th Grade",
                            "In 10th Grade",
                            "Baseline GPA",
                            "Num. A's",
                            "Num. B's",
                            "Num. C's",
                            "Num. D's",
                            "Num. F's",
                            "Missing Baseline GPA/Grades",
                            "Days Absent",
                            "Missing Attendance Data",
                            "Math Test Score (Z)", "Reading Test Score (Z)",
                            "Missing Math Test", "Missing Reading Test",
                            "Out-of-School Suspensions", "Disciplinary Incidents",
                            "Any Arrests at Baseline", "Arrests: Violent Crime", "Arrests: Property Crime", "Arrests: Drug Crime",
                            'Math Score - Decile in Previous School',
                            'In Study 2',
                            'Participated in Year 1 of Study')
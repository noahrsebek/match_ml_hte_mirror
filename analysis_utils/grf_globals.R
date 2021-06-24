

load_master_dataset <- function(addvars=T,
                                load_cached = T){
  require(dplyr)
  
  if (load_cached == T){
    master_pool <- readr::read_csv("/export/projects/migrated-BAM/Match_AnalysisFiles_May2016/Results/2020/grf/master_dataset_cached.csv")
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
    
  }
  return(master_pool)
  
}

# --> outcomes of interest (and labels)
outcomes_of_interest <- c('mathxil_z_post1_np',
                          "mathgpa_post1",
                          "mathfail_post1",
                          "mathfailpercent_post1",
                          "nonmathgpa_post1",
                          "readxil_z_post1_np",
                          "graduated_ontime",
                          "graduated_ever",
                          #'treat_post2',
                          "gpa11_math",
                          "eleventh_grade_math_z"
)

flipped_outcomes <- c('mathfail_post1',
                      'mathfailpercent_post1')

outcomes_of_interest_labels <- c('Math test score (Z)',
                                 'Math GPA',
                                 "Math Course Failures",
                                 "Percent of Math Courses Failed",
                                 "Non-math GPA",
                                 "Reading test score (Z)",
                                 "Graduated on-time",
                                 "Graduated ever",
                                 #'Participated in Study 1 Year 2',
                                 "11th Grade Math GPA",
                                 "11th Grade Math Test Score (Z)")


outcome_and_label_list <- outcomes_of_interest
names(outcome_and_label_list) <- outcomes_of_interest_labels


# list of controls/covariates/features to use in the model
# - these are all the 'standard' covariates we include in the ITT/TOT models
# - ... buuut without imputing 0's for the missing values
# - we still include missing dummies, however
controls <- c("d13andunder", "d14", "d15", "d16", "d17andover", "dlearningdisabled", "dfreelunch", "dblack", "dhispanic", "dother",
              "dgrade9", "dgrade10", "gpa_pre", "numAs_pre", "numBs_pre", "numCs_pre", "numDs_pre", "numFs_pre", "missing_gpa_pre",
              'dfemale',
              "days_absent_pre", "missing_attend_pre",
              "mathxil_z_pre_np", "readxil_z_pre_np",
              "mathxil_z_pre_missing",
              "readxil_z_pre_missing",
              "oss_dis_pre", "incidents_pre", "any_arrests_pre", "violent_pre", "property_pre", "drug_pre",
              'baseline_mathxil_inschool_decile')



controls_sans_missingness_dummies <- c("d13andunder", "d14", "d15", "d16", "d17andover",
                                       "dlearningdisabled", "dfreelunch", "dblack", "dhispanic", "dother",
                                       "dgrade9", "dgrade10", "gpa_pre",
                                       "numAs_pre", "numBs_pre", "numCs_pre", "numDs_pre", "numFs_pre",
                                       #"missing_gpa_pre",
                                       'dfemale',
                                       "days_absent_pre",
                                       #"missing_attend_pre",
                                       "mathxil_z_pre_np", "readxil_z_pre_np",
                                       #"mathxil_z_pre_missing",
                                       #"readxil_z_pre_missing",
                                       "oss_dis_pre", "incidents_pre", "any_arrests_pre",
                                       "violent_pre", "property_pre", "drug_pre",
                                       'baseline_mathxil_inschool_decile')




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
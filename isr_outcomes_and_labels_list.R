
isr_outcomes_2014 <- c("snw1_post1_recode1",
                       "snw2_post1_zscore_match",
                       "snw2a_post1_zscore_match",
                       "snw2b_post1_zscore_match",
                       "snw2c_post1_zscore_match",
                       'snw3_post1_recode1',
                       "snw4_post1_recode1",
                       
                       #Academic
                       "yed1_post1_zscore_match",
                       "yed2_post1_zscore_match",
                       "yed3_post1_zscore_match",
                       
                       "yed4_post1_recode1",
                       "yed5_post1_recode1",
                       "yed7_post1_zscore_match",
                       "yed8_post1_zscore_match",
                       "yed9a_post1_zscore_match",
                       "yed9b_post1_zscore_match",
                       "yed9c_post1_zscore_match",
                       "yed9d_post1_zscore_match",
                       "yed9e_post1_zscore_match",
                       "olk2_post1_zscore_match",
                       
                       #Adults
                       "yas1_post1",
                       "yas2_post1",
                       "yas5_post1_recode",
                       
                       #Thinking/impulse
                       "aut1_post1_zscore_match",
                       "aut2_post1_zscore_match",
                       "rdm1_post1_recode1",
                       
                       #Emotions
                       "apa1_post1_inverted_zscore_match",
                       "olk1_post1_inverted_zscore_match",
                       "apa2_post1_recode1",
                       "apa2a_post1_inverted_zscore_match",
                       
                       #agency
                       "lcs1_post1_zscore_match",
                       "lcs2_post1_zscore_match",
                       "lcs3_post1_zscore_match",
                       
                       "lcs4_post1_zscore_match",
                       "lcs5_post1_inverted_zscore_match",
                       
                       #workethic
                       "ncg2_post1_inverted_zscore_match",
                       "ncg4_post1_inverted_zscore_match",
                       "ncg6_post1_zscore_match",
                       "ncg8_post1_inverted_zscore_match",
                       "ncg7_post1_inverted_zscore_match",
                       'ncc2_post1_inverted_zscore_match',
                       
                       #indices
                       "aut_index",
                       "aut_index2",
                       'aut_index3',
                       #"sel_index",
                       "apa_index",
                       "lcs_index",
                       "olk_index",
                       "ncg_index",
                       'ncc_index',
                       
                       # "ncp1_index",
                       # "ncp2_index",
                       # "ncp3_index",
                       # "ncp4_index",
                       # "ncp5_index",
                       # "ncp6_index",
                       # 
                       "snw_index",
                       "yed_index",
                       "yas_index",
                       
                       "ncc1_post1_inverted_zscore_match",
                       "ncc2_post1_inverted_zscore_match",
                       "ncc3_post1_zscore_match")

# -> isr labels
isr_outcome_labels_2014 <- c("Reports No Close Friends (Dummy)",
                             "Friends think it is important to attend classes regularly (Z)",
                             "Friends think it is important to get good grades (Z)",
                             "Friends think it is important to study (Z)",
                             "Friends think it is important to continue education to college (Z)",
                             "Have stopped hanging around with someone (Recoded Dummy)",
                             "Have started hanging around with someone (Recoded Dummy)",
                             
                             "Disruptions by others get in the way of my learning (Z)",
                             "Weekly time spent on homework (Z)",
                             "How much of your homework do you finish? (Z)",
                             "Expect to graduate high school and continue education (Dummy)",
                             "Want to graduate high school and continue education (Dummy)",
                             "How important are good grades to you? (Z)",
                             "Feeling safe at school (Z)",
                             "Like math (Z)",
                             "Good grades in math (Z)",
                             "Like school in general (Z)",
                             "Like reading (Z)",
                             "Good grades in reading (Z)",
                             "Agree: I am as smart as others my age (Z)",
                             
                             "Number of adults to talk to (No Change)",
                             "Number of adults who care (No Change)",
                             "Would talk to adults at school (Dummy)",
                             
                             "Agree: It is important to think before act (Z)",
                             "Agree: It is important to think when making decisions (Z)",
                             "Chooses to postpone cash (Dummy)",
                             
                             "How often do you have arguments? (Z)",
                             "How happy are you? (High = Very Happy) (Z)",
                             "Have you ever fought another person? (Dummy)",
                             "How often do you have fights? (Z)",
                             
                             "Agree: I have control over direction of life (Z)",
                             "Disagree: Every time I try to get ahead, something or somebody stops me (Z)",
                             "Disagree: Luck is more important than hard work (Z)",
                             "Disagree: My plans never work out, planning makes me unhappy (Z)",
                             "Agree: I can make plans work (Z)",
                             
                             "Agree: Setbacks don't discourage me (Z)",
                             "Agree: I am a hard worker (Z)",
                             "Disagree: I have difficulty maintaining focus (Z)",
                             "Agree: I am diligent (Z)",
                             "Agree: I finish what I begin (Z)",
                             "Agree: I can continue until everything is perfect (Z)",
                             
                             "Automaticity",
                             "Automaticity with SEL3",
                             "Automaticity with SEL3 and SEL4",
                             #"Self Distancing (-)",
                             "Awareness of Past Action",
                             "Locus of Control",
                             "Outlook",
                             "Grit",
                             "Conscientiousness",
                             # "Peer Conflict Vignette 1 (-)",
                             # "Peer Conflict Vignette 2 (-)",
                             # "Peer Conflict Vignette 3 (-)",
                             # "Peer Conflict Vignette 4 (-)",
                             # "Peer Conflict Vignette 5 (-)",
                             # "Peer Conflict Vignette 6 (-)",
                             "Social Networks",
                             "Education and Schooling",
                             "Adult Supports",
                             
                             "Agree: I am always prepared (Z)",
                             "Agree: I continue until everything is perfect (Z)",
                             "Agree: I leave a mess in my room (Z)")



isr_outcomes_2015 <- c("snw1_post2_recode1",
                       "snw2_post2_zscore_match",
                       "snw2a_post2_zscore_match",
                       "snw2b_post2_zscore_match",
                       "snw2c_post2_zscore_match",
                       "snw3_post2_recode1",
                       "snw4_post2_recode1",
                       "yed6_post2_zscore_match",
                       "yed2_post2_recode1",
                       "yed1_post2_recode1",
                       "yed7_post2_zscore_match",
                       "yed8_post2_zscore_match",
                       "yed9c_post2_zscore_match",
                       "yed9d_post2_zscore_match",
                       "yed10_post2_zscore_match",
                       "yed9a_post2_zscore_match",
                       "yed9b_post2_zscore_match",
                       
                       "yas1_post2",
                       "yas2_post2",
                       
                       #Impulse
                       "aut1_post2_zscore_match",
                       "aut2_post2_zscore_match",
                       "rdm10_post2_recode1",
                       
                       #Emotions
                       "olk1_post2_inverted_zscore_match",
                       "hit1_post2_inverted_zscore_match",
                       "hit2_post2_inverted_zscore_match",
                       "hit4_post2_inverted_zscore_match",
                       "hit6_post2_inverted_zscore_match",
                       
                       #Ethic
                       "ncg2_post2_inverted_zscore_match",
                       "ncg4_post2_inverted_zscore_match",
                       "ncg6_post2_zscore_match",
                       "ncg8_post2_inverted_zscore_match",
                       
                       #Belief
                       "gmt1_post2_inverted_zscore_match",
                       "gmt2_post2_inverted_zscore_match",
                       "gmt3_post2_inverted_zscore_match",
                       "gmt4_post2_inverted_zscore_match",
                       "gmt5_post2_inverted_zscore_match",
                       "gmt6_post2_inverted_zscore_match",
                       "olk2_post2_zscore_match",
                       
                       #Behavior
                       "sep1_post2_zscore_match",
                       "sep2_post2_inverted_zscore_match",
                       "sep3_post2_inverted_zscore_match",
                       "sep4_post2_zscore_match",
                       
                       "rdb1_post2_zscore_match",
                       "rdb2_post2_zscore_match",
                       "rdb3_post2_zscore_match",
                       "rdb4_post2_zscore_match",
                       "rdb5_post2_zscore_match",
                       
                       "rdb6_post2_recode1",
                       "rdb7_post2_recode1",
                       "rdb8a_post2_recode1",
                       "rdb8b_post2_recode1",
                       
                       "rdb10_post2_inverted_zscore_match",
                       "rdb12_post2_inverted_zscore_match",
                       
                       "rdb13_post2_zscore_match",
                       "rdb14_post2_zscore_match",
                       "rdb15_post2_zscore_match",
                       "rdb16_post2_zscore_match",
                       "rdb17_post2_zscore_match",
                       "rdb18_post2_zscore_match",
                       "rdb19_post2_zscore_match",
                       "rdb20_post2_zscore_match",
                       
                       #Violence
                       "cvp1_post2_zscore_match",
                       "cvp2_post2_zscore_match",
                       "cvp3_post2_zscore_match",
                       "cvp4_post2_zscore_match",
                       
                       #indices
                       "aut_index",
                       "sel_index",
                       "ncg_index",
                       "ses_index",
                       # "ncp1_index",
                       # "ncp2_index",
                       # 
                       # "ncp3_index",
                       # "ncp4_index",
                       "hit_index",
                       #"con_index",
                       "yed_index",
                       "snw_index",
                       
                       "yas_index",
                       "gmt_index",
                       "sep_index",
                       "olk_index",
                       "meh_index",
                       #"cvp_index",
                       
                       "rdb_index",
                       
                       "con1_post2_inverted_zscore_match",
                       "con2_post2_zscore_match",
                       "con3_post2_inverted_zscore_match",
                       "con4_post2_zscore_match",
                       "con5_post2_zscore_match",
                       "con6_post2_inverted_zscore_match",
                       "con7_post2_inverted_zscore_match",
                       "con8_post2_inverted_zscore_match",
                       "con9_post2_zscore_match")


isr_outcome_labels_2015 <- c("Reports No Close Friends (Dummy)",
                             "Friends think it is important to attend classes regularly (Z)",
                             "Friends think it is important to get good grades (Z)",
                             "Friends think it is important to study (Z)",
                             "Friends think it is important to continue education to college (Z)",
                             "Have stopped hanging around with someone (Recoded Dummy)",
                             "Have started hanging around with someone (Recoded Dummy)",
                             
                             "Weekly time spent on homework (Z)",
                             "Expect to graduate high school and continue education (Dummy)",
                             "Want to graduate high school and continue education (Dummy)",
                             "How important are good grades to you? (Z)",
                             "Feeling safe at school (Z)",
                             "Like math (Z)",
                             "Good grades in math (Z)",
                             "Like school in general (Z)",
                             "Like reading (Z)",
                             "Good grades in reading (Z)",
                             
                             "Number of adults to talk to (No Change)",
                             "Number of adults who care (No Change)",
                             
                             "Agree: It is important to think before act (Z)",
                             "Agree: It is important to think when making decisions (Z)",
                             "Chooses to postpone cash (Dummy)",
                             
                             "How happy are you? (High = Very Happy) (Z)",
                             "Disagree: When I get mad, I don't care who I hurt",
                             "Disagree: When I lose my temper, it's because people try to make me mad",
                             "Disagree: Only a coward would ever walk away from a fight",
                             "Disagree: It is no use trying to stay out of fights",
                             
                             "Agree: Setbacks don't discourage me (Z)",
                             "Agree: I am a hard worker (Z)",
                             "Disagree: I have difficulty maintaining focus (Z)",
                             "Agree: I am diligent (Z)",
                             
                             #Intelligence Belief
                             "Disagree: You can learn but not change intelligence (Z)",
                             "Disagree: Intelligence can't be changed much (Z)",
                             "Disagree: You have a certain amount of intelligence, can't change (Z)",
                             "Disagree: Moral character can't be changed (Z)",
                             "Disagree: Responsibility can't be changed (Z)",
                             "Disagree: Moral traits can't be changed (Z)",
                             "Agree: I am as smart as other students (Z)",
                             
                             #Risk Behavior
                             "Percent chance student in regular school in 1 year (Z)",
                             "Percent chance you will get someone pregnant within 1 year (Z)",
                             "Percent chance of arrest within 1 year (Z)",
                             "Percent chance of graduating high school by 20 (Z)",
                             
                             "During your life, how many days have you had at least one drink of alcohol? (Z)",
                             "During the past 30 days, on how many days did you have at least one drink of alcohol? (Z)",
                             "During your life, how many times have you used marijuana? (Z)",
                             "During the past 30 days, how many times did you use marijuana? (Z)",
                             "During your life, how many times have you tried any other sort of illegal drug/inhalant/prescription drug? (Z)",
                             
                             "Do any of your brothers, sisters, cousins, or friends belong to a gang? (Dummy)",
                             "Do you belong to a gang? (Dummy)",
                             "Have you ever sold marijuana or any other drug to your friends? (Dummy)",
                             "Have you ever sold marijuana or any other drug to people you didn't know? (Dummy)",
                             
                             "During the past 3 months with how many people did you have sexual intercourse? (Z)",
                             "How many times have you gotten someone pregnant? (Z)",
                             
                             "In the past year, how many times did you get in a physical fight in which you were so badly injured that you were treated by a doctor or a nurse? (Z)",
                             "In the past year, how often did you hurt someone badly enough in a physical fight that he or she needed to be treated by a doctor or nurse? (Z)",
                             "During the past 30 days, on how many days did you carry a weapon -- such as a gun, knife, or club -- to school? (Z)",
                             "In the past year, how often did you paint graffiti or signs on someone else's property or in a public place? (Z)",
                             "In the past year, how often did you deliberately damage property that didn't belong to you? (Z)",
                             "In the past year, how often did you take something from a store without paying for it? (Z)",
                             "In the past year, how often did you drive a car without owner's permission? (Z)",
                             "In the past year, how often did you break into someone's home in order to steal? (Z)",
                             
                             #Violence
                             "In the past year, how often did someone pulle a gun/knife on you? (Z)",
                             "In the past year, how often did you get into a physical fight? (Z)",
                             "In the past year, how often did you get jumped? (Z)",
                             "In the past year, how often did you get beaten up and something was stolen from you? (Z)",
                             
                             "Automaticity",
                             "Self Distancing",
                             "Grit",
                             "Sensation Seeking",
                             # "Peer Conflict Vignette 1 (-)",
                             # "Peer Conflict Vignette 2 (-)",
                             # 
                             # "Peer Conflict Vignette 3 (-)",
                             # "Peer Conflict Vignette 4 (-)",
                             "How I Think",
                             #"Conscientiousness (-)",
                             "Education and Schooling",
                             "Social Networks",
                             
                             "Adult Supports",
                             "Growth Mindset",
                             "Subjective Expectations",
                             "Outlook",
                             "Mental Health",
                             #"Crime Victimization (-)",
                             
                             "Risky Behavior",
                             
                             # Conscientiousness
                             "I see myself as someone who does things carefully and completely.",
                             "I see myself as someone who can be somewhat careless.",
                             "I see myself as someone who is a reliable worker.",
                             "I see myself as someone who tends to be disorganized.",
                             "I see myself as someone who tends to be lazy.",
                             "I see myself as someone who keeps working until things are done.",
                             "I see myself as someone who does things efficiently (quickly and correctly).",
                             "I see myself as someone who makes plans and sticks to them.",
                             "I see myself as someone who is easily distracted; has trouble paying attention.")
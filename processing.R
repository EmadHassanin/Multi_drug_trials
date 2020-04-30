#Processing with the focal vs. whole data
#Contains all the processing, loading, statistics (survival)

#### Loading packages
library(tidyverse)
library(DescTools)
library(data.table)
library(timetools)
library(eeptools)
library(DiagrammeR)
library(knitr)
library(survival)
library(survminer)
library(bookdown)
library(kableExtra)
library(pander)
library(readxl)
library(lubridate)
library(scales)
library(diagram)
library(captioner)
library(officer)

#### Function used later on

#To load tables

### Load tables
# Caching?


load_table <- function(table_name){
  file_name <- paste0(tpref, table_name, allsuf)
  read_excel(file.path(dir, file_name ), guess_max = 10000) }

#Function to remove duplicates of AED combination
uniquecomb <- function(aed_final) {
  aed_mono <- aed_final %>% 
    split(aed_final$S_FKEY) #split according to subjects
  for(l in 1:length(aed_mono)){
    aed_mono[[l]] <- as.data.frame(aed_mono[[l]]) %>% 
      rowwise %>%  #grouping by row
      mutate(temp = paste(sort(c(D120_AED_GENERIC_ABBREV.x, D120_AED_GENERIC_ABBREV.y)), collapse= "+")) %>%  #sort names of each drug combination(each row)
      distinct(temp, .keep_all=T) #%>%  #filter unique sorted combination
    #select(-temp) 
  }
  aed_mono <- do.call("rbind",aed_mono)
  return(aed_mono)
}

#Gives censor data based on mode of action for survival analysis in GA and SCB models
censor_aed <- function(aed_input, mode_action){
  relevel_key <- paste(mode_action, mode_action, sep="+")
  aed_censorm <- aed_input %>% filter(str_detect(key, mode_action)) %>% filter(int_days < 9500)
  aed_censorm$key <- factor(aed_censorm$key)
  aed_censorm$key <- relevel(aed_censorm$key, relevel_key)
  return(aed_censorm)
}

#Kaplan meyer for GA and SCB models

kaplan_meyer <- function(aed_input, mode_action){
  # mode_action one of "SCB", "GA" (which table)
  # Error message
  if(!mode_action %in% c("SCB", "GA","SV","MM")) {
    stop(paste("Mode of action ", mode_action, " not defined."))
  }
  aed_censorr <- censor_aed(aed_input, mode_action)
  # surv_object <-Surv(time = aed_censorr$int_days, event = aed_censorr$censor)
  fit <- survfit(Surv(time = int_days, event = censor) ~ key , data = aed_censorr)
  names(fit$strata) <- str_remove(names(fit$strata),"key=")
  ggsurvplot(fit,data= aed_censorr, pval = TRUE,
             #xlim = c(0,6000),
             scale_colour_brewer = "Dark2")+
    ylab("Persistence probability") +
    xlab ("Time (days)") -> survplot
  return(list(fit, survplot))
}

# Cox proportional hazard regression function for SCB and GA models

cox_model <- function(aed_input, mode_action){
  
  # fit2 <- survfit(surv_object_sc ~ key , data = aed_censor_sc)
  # names(fit2$strata) <- str_remove(names(fit2$strata),"key=")
  # surv_object_sc <-Surv(time = aed_censor_sc$int_days, event = aed_censor_sc$censor) 
  aed_censor_m <- censor_aed(aed_input, mode_action)
  surv_object <- Surv(time = aed_censor_m$int_days, event = aed_censor_m$censor)
  fit.coxph <- coxph(surv_object  ~ key , data = aed_censor_m)
  return(fit.coxph)
}

#Summary from cox analysis

cox_summary <- function(fit_coxph,combination){
  
  summary(fit_coxph)  -> sum_cox1
  sum_cox1$coefficients -> prop1
  rownames(prop1) <-  str_remove(rownames(prop1),"key")
  prop1 <- setDT(as.data.frame(prop1), keep.rownames = TRUE)
  colnames(prop1) <-  c(combination,"coef","coeff" ,"se","z_score","P")
  sum_cox1$conf.int -> cox1
  rownames(cox1) <-  str_remove(rownames(cox1),"key")
  cox1 <- setDT(as.data.frame(cox1), keep.rownames = TRUE)
  colnames(cox1) <-  c(combination,"HR","coeff" ,"lower_CI","upper_CI")
  
  return(list(cox1,prop1))
  
}

#For estimation of hazard ratios based on the p value

p_value_output <- function(p_value){
  p_value <- case_when(p_value < 0.001 ~ "< 0.001",
                       TRUE ~ sprintf("%.3f", p_value))
  return(p_value)
}




#### Loading the data

#source("load.epipgx.tables.R")
source("manuscript_setup/config.R")

### Load tables

table_names <- c("ADR", "AED", "INVESTIGATION", "PREGN_DRUGS", "PREGNANCY", "SEIZ_REM",
                 "SUBJECT")
# e <- map(table.names, load.table) 
for (tbl_name in table_names){
  tdf <- load_table(tbl_name)
  assign(tbl_name, tdf)
}

## Change to date formats for all fields in AED and SUBJECT
# Change the formate to numeric first then to date, otherwise the excel file shows it as characters
# TODO: Rework date formats for all fields
# AED$D108_START_date <-  as.Date(as.numeric(AED$D108_START_date), origin = "1900-01-01") 
# AED$D113_END_date <-  as.Date(as.numeric(AED$D113_END_date), origin = "1900-01-01") 

#cols_AED <- c( "D108_START_date", "D109_START_yyyymm", "D113_END_date", "D114_END_yyyymm")
cols_AED <- c( 14, 16, 21, 22)
AED[,cols_AED] <- lapply(AED[,cols_AED], FUN= as.numeric)
AED[,cols_AED] <- lapply(AED[,cols_AED], function(x) as.Date(x,origin="1900-01-01") )

# cols_subject <- c(10, 22, 46, 52, 58, 65, 140, 180)
# SUBJECT[,cols_subject] <- lapply(SUBJECT[,cols_subject], FUN= as.numeric)
# SUBJECT[,cols_subject] <- lapply(SUBJECT[,cols_subject], function(x) as.Date(x,origin="1900-01-01") )


# Assign AED to its MOA

aed_moa <- tribble(
  ~D120_AED_GENERIC_ABBREV, ~MOA, ~Desc,
  "LEV" , "SV" , "synaptic vesicle potein 2A",
  "LCM" , "SCB" , "sodium channel blocker",
  "LTG" , "SCB" , "sodium channel blocker",
  "OXC" , "SCB" , "sodium channel blocker",
  "PHT" , "SCB" , "sodium channel blocker",
  "ETN" , "SCB" , "sodium channel blocker",
  "ESL" , "SCB" , "sodium channel blocker",
  "CBZ" , "SCB" , "sodium channel blocker",
  "CLB" , "GA" , "gaba analog",
  "CNZ" , "GA" , "gaba analog",
  "PB" , "GA" , "gaba analog",
  "GBP" , "GA" , "gaba analog",
  "VGB" , "GA" , "gaba analog",
  "PGB" , "GA" , "gaba analog",
  "PRM" , "GA" , "gaba analog",
  "TGB" , "GA" , "gaba analog",
  "DZP" , "GA" , "gaba analog",
  "LZP" , "GA" , "gaba analog",
  "NZP" , "GA" , "gaba analog",
  "ZNS" , "MM" , "multiple mechanisms",
  "VPA" , "MM" , "multiple mechanisms",
  "FBM" , "MM" , "multiple mechanisms",
  "TPM" , "MM" , "multiple mechanisms"
)

#### Preparing the data

subject_clean <- SUBJECT %>%               
  filter(!str_detect(S102_SUBJ_ID, "^S120")) %>% 
  mutate(S158_LATEST_VISIT_date =  as.Date(as.numeric(S158_LATEST_VISIT_date), origin = "1900-01-01"),
         S108_DOB_DATE = as.Date(as.numeric(S108_DOB_DATE), origin = "1900-01-01"),
         dob = S108_DOB_DATE)
          


## change End_date  with Last visit date if D113_END_date= NA and  D115_END_NA = "AED tx ongoing"
AED %>%  filter(!str_detect(D101_SUBJ_ID, "^S120")) %>%  
  left_join(aed_moa , by =c("D120_AED_GENERIC_ABBREV"="D120_AED_GENERIC_ABBREV")) %>%
  left_join(subject_clean %>% 
              filter(S158_LATEST_VISIT_date > 0) %>% filter(!is.na(S108_DOB_DATE)) %>% 
              select(S_PKEY, S158_LATEST_VISIT_date, S103_GENDER, dob ), 
            by = c("S_FKEY"="S_PKEY")) %>% 
  
  mutate(end_date = case_when(
    is.na(D113_END_date) & str_detect(D115_END_NA, "ongoing") ~ S158_LATEST_VISIT_date, 
    TRUE ~ D113_END_date)) %>% 
  
  select( S_FKEY, D_PKEY, 
         D113_END_date,
         D120_AED_GENERIC_ABBREV,
         MOA,
         Desc,
         D108_START_date,
         end_date,
         D115_END_NA,
         D140_OUTCOME,
         dob,
         S103_GENDER) %>%
         { nrow(.) ->> n_aed_start;. } %>%
         { length(unique(.$S_FKEY)) ->> n_patients_start;. } -> aed_tmp  #serves to make a summary for patients with 2 and more aed trials

#patients have more than one trial

cases_more_than_2_AED <-aed_tmp %>% 
  group_by(S_FKEY) %>% 
  summarize(aed_count = n()) %>%
  filter(aed_count > 1)

#### Extracting the data

#joining to obtain drug combinations 
# TODO: check why the aed_final number of trials differes with age calculation and last filtering step when done one or the other way
aed_tmp %>% 
  inner_join(aed_tmp [,-c(11,12)],
             by="S_FKEY" ) %>% 
  filter(D108_START_date.x <= D108_START_date.y & 
           D120_AED_GENERIC_ABBREV.x != D120_AED_GENERIC_ABBREV.y)%>% #find the intersect between two drugs in each row
  #mutate(interval = D108_START_date.y - end_date.x) -> aed_tmp
  #calculate the interval of each ttt combination in days, choose which more than zero
  mutate(int_days = end_date.x - D108_START_date.y, units = "day") %>% 
  filter(int_days > 0 ) %>% 
  #calculte interval of final ttt combination
  mutate(ovr_interval = lubridate::intersect(lubridate::interval(D108_START_date.x,end_date.x),
                                             lubridate::interval(D108_START_date.y,end_date.y))) %>%
  mutate(start_ovr = int_start(ovr_interval)) %>%
  mutate(end_ovr   = int_end(ovr_interval)) %>%
  mutate(int_days  = difftime(int_end(ovr_interval),
                              int_start(ovr_interval), units = "day"))%>%
  select(-ovr_interval) %>%
  { length(unique(.$S_FKEY)) ->> n_patients_concomitant;. } %>%
  { nrow(.)->>n_aed_concomitant;. } %>%
  
  #filtering for concomitant trials longer than 90 days
  filter(int_days > 90) %>%
  {length(unique(.$S_FKEY)) ->> n_patients_concomitant_90_days;. }%>%
  {nrow(.)->>n_aed_concomitant_90_days;. } %>%
  
  #excluding cases without moa information
  mutate(MOA_comb_same = (MOA.x == MOA.y)) %>% 
  filter(!is.na(MOA.x),!is.na(MOA.y)) %>%
  rowwise %>%
  mutate(key = paste(sort(c(MOA.x, MOA.y)), collapse="+")) %>%
  
  #filtering for unique moa combinations
  uniquecomb() %>%
  # 
  # aedtmp <- aedtmp%>%     
  # left_join(subject_clean %>% 
  #             filter(!is.na(S108_DOB_DATE)) %>% 
  #             select(S_PKEY, S108_DOB_DATE,S103_GENDER, dob), by = c("S_FKEY"="S_PKEY")) %>% 
  # filter(dob.y < D108_START_date.x) %>% 
  # mutate(start_date = D108_START_date.x)  %>% 
  # mutate(age = age_calc(dob.y, enddate = start_date, units = "years")) %>% #end date needs to be determined 
  # filter(!is.na(S103_GENDER.y)) 

  #calculating age and removing cases without gender info
  filter(dob < D108_START_date.x) %>% 
  # mutate(start_date = D108_START_date.x)  %>% 
  mutate(age = age_calc(dob, enddate = D108_START_date.x, units = "years")) %>% #end date needs to be determined 
  filter(!is.na(S103_GENDER)) -> aed_final



# Focal cases (final with the investigation)

aed_final_focal <- aed_final %>% 
  left_join(subject_clean %>%
              select(S_PKEY,S171_EPI_diagnosis_1) , by = c("S_FKEY" = "S_PKEY")) %>% 
  filter(str_detect(S171_EPI_diagnosis_1 , regex('local', ignore_case=TRUE))) 

#number of trials
aed_final_focal %>% group_by(S_FKEY) %>% summarize(n())  %>% nrow 

#Results
#Study sample

#### Monotrials frequencies reporting
# the whole data
moa_monotrials <- AED %>% 
  left_join(aed_moa , by = c("D120_AED_GENERIC_ABBREV" = "D120_AED_GENERIC_ABBREV")) %>% 
  filter(!is.na(MOA)) %>% 
  count(MOA) %>% mutate(per = n/sum(n)) 

aed_monotrials <- AED %>% 
  left_join(aed_moa , by = c("D120_AED_GENERIC_ABBREV" = "D120_AED_GENERIC_ABBREV")) %>% 
  filter(!is.na(MOA)) %>% 
  group_by( D120_AED_GENERIC_NAME, MOA) %>% 
  summarise(trials =  length(D120_AED_GENERIC_ABBREV), 
            per= trials/ length(.$D120_AED_GENERIC_ABBREV))

#Monotrial frequencies in focal cases reporting
# aed with investigation
moa_monotrials_focal <- AED %>% 
  left_join(subject_clean %>% select(S_PKEY,S171_EPI_diagnosis_1) ,
            by = c("S_FKEY" = "S_PKEY")) %>% 
  filter(str_detect(S171_EPI_diagnosis_1 , regex('local', ignore_case=TRUE)))  %>% 
  left_join(aed_moa , by = c("D120_AED_GENERIC_ABBREV" = "D120_AED_GENERIC_ABBREV")) %>% 
  filter(!is.na(MOA)) %>% 
  count(MOA) %>% 
  mutate(per = n/sum(n))


aed_monotrials_focal <- AED %>% 
  left_join(subject_clean %>% select(S_PKEY,S171_EPI_diagnosis_1) ,
            by = c("S_FKEY" = "S_PKEY")) %>% 
  filter(str_detect(S171_EPI_diagnosis_1 , regex('local', ignore_case=TRUE)))  %>% 
  left_join(aed_moa , by = c("D120_AED_GENERIC_ABBREV" = "D120_AED_GENERIC_ABBREV")) %>% 
  filter(!is.na(MOA)) %>% 
  group_by( D120_AED_GENERIC_NAME, MOA) %>% 
  summarise(trials =  length(D120_AED_GENERIC_ABBREV), 
            per= trials/ length(.$D120_AED_GENERIC_ABBREV))

#### Demographic distribution over the cohort 

# The number of patients after data filtering

# All cases
aed_demographics <- aed_final %>% 
  group_by(key) %>%
  summarise(age_mean = mean(age) , age_median = median(age), 
            male = sum(S103_GENDER == 1)/n() , female = sum(S103_GENDER == 2)/n())
# The subgroup of patients included in the study

clinical_details <- aed_final %>%  
  group_by(key) %>%
  summarise(aed_trials = n() , trials_per = aed_trials/nrow(aed_final),
            patient_count= n_distinct(S_FKEY), patient_per = patient_count/n_patients_concomitant_90_days) 

# A dataframe with summarized values from previous one per grouped combinations

#TODO: recheck numbers in total aed_trials vs total aed_trials number in aed_final
moaX <- aed_final %>% select(S_FKEY, D_PKEY.x, MOA.x) 
names(moaX) <- c("S_FKEY", "D_PKEY", "MOA")

moaY <- aed_final %>% select(S_FKEY, D_PKEY.y, MOA.y)
names(moaY) <- c("S_FKEY", "D_PKEY", "MOA")
moaXmoaY <- rbind(moaX, moaY)

moaXmoaY%>%
  group_by(MOA) %>% 
  summarise(aed_trials = n_distinct(D_PKEY) , trials_per = aed_trials/nrow(aed_final), #length(unique(moaXmoaY$D_PKEY)),
            patient_count= n_distinct(S_FKEY), patient_per = patient_count/n_patients_concomitant_90_days) -> combination_detailsMOA
names(combination_detailsMOA) <- c("key", "aed_trials", "trials_per", "patient_count", "patient_per")


AEDclinical_details <- rbind(clinical_details, combination_detailsMOA)

#Individual combinations
aed_final %>%
group_by(key, add = TRUE) %>%
 count(temp) %>%
group_nest(key) %>%
 deframe %>%
 map(arrange,desc(n)) %>%
 map(function(x) setNames(x, c("combination","count"))) %>%
 pander()

# Focal cases
aed_focal_demographics <- aed_final_focal %>%  #patient characteristics
  group_by(key) %>%
  summarise(age_mean = mean(age) , age_median = median(age), 
            male = sum(S103_GENDER == 1)/n() , female = sum(S103_GENDER == 2)/n()) 

# The subgroup of patients included in the study
aed_focal_summarized <- aed_final_focal %>%  
  group_by(key) %>%
  summarise(aed_trials = n() , trials_per = aed_trials/nrow(aed_final),
            patient_count = n_distinct(S_FKEY), patient_per = patient_count/n_patients_concomitant_90_days) 
#A dataframe with summarized values from previous one per grouped combinations

moaXfocal <- aed_final %>% select(S_FKEY, D_PKEY.x, MOA.x) 
names(moaXfocal) <- c("S_FKEY", "D_PKEY", "MOA")

moaYfocal <- aed_final %>% select(S_FKEY, D_PKEY.y, MOA.y)
names(moaYfocal) <- c("S_FKEY", "D_PKEY", "MOA")
moaXmoaYfocal <- rbind(moaXfocal, moaYfocal)

moaXmoaYfocal%>%
  group_by(MOA) %>% 
  summarise(aed_trials = n_distinct(D_PKEY) , trials_per = aed_trials/nrow(aed_final), #length(unique(moaXmoaY$D_PKEY)),
            patient_count= n_distinct(S_FKEY), patient_per = patient_count/n_patients_concomitant_90_days) -> combination_detailsMOA_focal
names(combination_detailsMOA_focal) <- c("key", "aed_trials", "trials_per", "patient_count", "patient_per")

focal_clinical_details <- rbind(aed_focal_summarized, combination_detailsMOA_focal)

# Individualcombinations - focal cases
aed_final_focal %>%
group_by(key) %>%
count(temp) %>%
group_nest(key) %>%
deframe %>%
map(arrange,desc(n)) %>%
map(function(x) setNames(x, c("combination","count"))) %>%
pander()


#### Drug persistance (duration of the AED, by MOA) 
# AED in all cases
comb_per <- aed_final %>% 
  group_by(key)  %>% 
  summarise(per_mean= mean(int_days) , per_median= median(int_days) ) 

aed_final %>% 
  group_by(key)  %>% 
  summarise(per_mean= mean(int_days) , per_median= median(int_days) ) %>%
  gather(key = "persistence" , value = "value" , - key) %>% 
  filter(key != "MM+MM")

# Focal AED

comb_per_focal <- aed_final_focal %>%  #persistence by MOA
  group_by(key)  %>%
  summarise(per_mean= mean(int_days) , per_median= median(int_days) ) 


aed_final_focal %>%  #persistence by MOA
  group_by(key)  %>% summarise(per_mean= mean(int_days) , per_median= median(int_days) ) %>%
  gather(key = "persistence" , value = "value" , - key) %>% 
  filter(key != "MM+MM") 


#### Survival analysis for all cases

aed_censor <- aed_final %>% 
  mutate(censor = case_when(
    str_detect(D115_END_NA.x, "ongoing")  & str_detect(D115_END_NA.y, "ongoing") ~ 1,TRUE ~ 2))


#Kaplan Meyer survival curves of persistence comparisons among AEDs combinations for GA combination.
#surv_object <-Surv(time = aed_censor$int_days, event = aed_censor$censor)

#model1: GA trials and GA+GA is the reference
kaplan_meyer(aed_censor, "GA") ->fit_GA

## Kaplan Meyer for survival curves of persistence comparisons among AEDs combinations for SCB combination.
#Table include number of trials and number of events and median persistence time. 
#model2: SCB trials and SCB+SCB is the reference
kaplan_meyer(aed_censor, "SCB")->fit_SCB


## Cox proportional hazard
#model1 > GA trials and GA+GA is the reference

aed_censor_GA <- censor_aed(aed_censor, "GA")
fit.coxph_GA <- cox_model(aed_censor,"GA")
ggforest(fit.coxph_GA, data = aed_censor_GA)

cox1 <- cox_summary(fit.coxph_GA,"GA_combination")


#Hazard Ratios of model1 GA combination therapy discontinuation.
cox_ga <- cox1[[1]] %>% select(-coeff) %>% left_join(cox1[[2]] %>% select(GA_combination,P) %>% 
                                        mutate( P = p_value_output(P),
                                                reference = "GA+GA"), by= c("GA_combination"="GA_combination")) %>% 
  select(GA_combination,reference,HR ,lower_CI,upper_CI, P) 


#model2 > SCB trials and SCB+SCB is the reference
aed_censor_SC <- censor_aed(aed_censor, "SCB")
fit.coxph_SC <- cox_model(aed_censor,"SCB")
ggforest(fit.coxph_SC, data = aed_censor_SC)

cox2 <- cox_summary(fit.coxph_SC,"SCB_combination")

#Hazard Ratios of model2 SCB combination therapy discontinuation.
cox_scb <- cox2[[1]]  %>% select(-coeff) %>% left_join(cox2[[2]] %>% select(SCB_combination,P) %>% 
                                         mutate( P = p_value_output(P),
                                                 reference = "SCB+SCB"), 
                                       by= c("SCB_combination"="SCB_combination")) %>% 
  select(SCB_combination,reference,HR ,lower_CI,upper_CI, P) 



#### Outcome of the combinations for all cases

out_comb <- aed_final %>%
  filter(D140_OUTCOME.y == D140_OUTCOME.x) %>%
  count( key, D140_OUTCOME.x) %>% 
  spread(D140_OUTCOME.x, n)

sum <- out_comb %>% 
  select(-key) %>% 
  rowSums() 
cbind(out_comb,sum) -> out_comb
colnames(out_comb) <- c("AED_combination","response","failure","unclassified","unknown","total") 

out_comb %>% 
  mutate(per = response/total) %>% 
  select(AED_combination,response,per) 


#### Survival analysis for focal cases

### Kaplan Meier survival curves of persistence comparisons among AEDs combinations for GA combination (focal)

aed_censor_focal <- aed_final_focal %>% 
  mutate(censor = case_when(
    str_detect(D115_END_NA.x, "ongoing")  & str_detect(D115_END_NA.y, "ongoing") ~ 1,TRUE ~ 2))

#Kaplan Meyer
#surv_object <-Surv(time = aed_censor$int_days, event = aed_censor$censor)
#model1: GA trials and GA+GA is the reference
kaplan_meyer(aed_censor_focal, "GA") -> fit_focalGA

#model2: SCB trials and SCB+SCB is the reference
kaplan_meyer(aed_censor_focal, "SCB")-> fit_focalSCB

# Cox proportional hazard (focal)
#model1 > GA trials and GA+GA is the reference
focal_censor_GA <- censor_aed(aed_censor_focal, "GA")
fit.coxph_GA_focal <- cox_model(aed_censor_focal,"GA")
ggforest(fit.coxph_GA_focal, data = focal_censor_GA)

cox1focal <- cox_summary(fit.coxph_GA_focal,"GA_combination")


cox1focal[[1]] %>% select(-coeff) %>% left_join(cox1focal[[2]] %>% select(GA_combination,P) %>% 
                                        mutate( P = p_value_output(P),
                                                reference = "GA+GA"), by= c("GA_combination"="GA_combination")) %>% 
  select(GA_combination,reference,HR ,lower_CI,upper_CI, P) -> cox_ga_focal

#model2 > SCB trials and SCB+SCB is the reference
focal_censor_SC <- censor_aed(aed_censor_focal, "SCB")
fit.coxph_SC_focal <- cox_model(aed_censor_focal,"SCB")
ggforest(fit.coxph_SC_focal, data = focal_censor_SC)

cox2focal <- cox_summary(fit.coxph_SC_focal,"SCB_combination")

cox2focal[[1]] %>% select(-coeff) %>% left_join(cox2focal[[2]] %>% select(SCB_combination,P) %>% 
                                         mutate( P = p_value_output(P),
                                                 reference = "SCB+SCB"), 
                                       by= c("SCB_combination"="SCB_combination")) -> cox_scb_focal

#### The outcome of the combinations in focal cases
aed_final_focal %>%
  filter(D140_OUTCOME.y == D140_OUTCOME.x) %>%
  count( key, D140_OUTCOME.x) %>% 
  spread(D140_OUTCOME.x, n)  -> out_comb_focal

out_comb_focal %>% select(-key) %>% rowSums() -> sum
cbind(out_comb_focal,sum) -> out_comb_focal
colnames(out_comb_focal) <- c("AED_combination","response","failure","unclassified","unknown","total") 
out_comb_focal %>% 
  mutate(per = response/total) %>% select(AED_combination,response,per) -> out_comb_focal

#### To compare drug persistance between same-MOA and different-MOA combinations

## All cases - AEDs same vs. different-moa based 
aed_censorMOA <- aed_censor %>% filter(int_days< 9500)
aed_censorMOA$MOA_comb_same <- as.factor(aed_censorMOA$MOA_comb_same)
aed_censorMOA$MOA_comb_same <- relevel(aed_censorMOA$MOA_comb_same, "TRUE")

surv_object_sameMOA <-Surv(time = aed_censorMOA$int_days, event = aed_censorMOA$censor)
fit_same_moa <- survfit(surv_object_sameMOA  ~ MOA_comb_same , data = aed_censorMOA)
ggsurvplot(fit_same_moa,data= aed_censorMOA, pval = TRUE, legend.labs = c( "Same-MOA based combinations", "Different-MOA based combinations"),
                           #xlim = c(0,6000),
                             scale_colour_brewer = "Dark2")+
    ylab("Persistence probability") +
       xlab ("Time (days)")

#Focal cases - AEDs same vs. different-moa based 
aed_censor_focalMOA <- aed_censor_focal %>% filter(int_days< 9500)
aed_censor_focalMOA$MOA_comb_same <- as.factor(aed_censor_focalMOA$MOA_comb_same)
aed_censor_focalMOA$MOA_comb_same <- relevel(aed_censor_focalMOA$MOA_comb_same, "TRUE")

surv_object_same_focals <-Surv(time = aed_censor_focalMOA$int_days, event = aed_censor_focalMOA$censor)
fit_focal_same_moa <- survfit(surv_object_same_focals  ~ MOA_comb_same , data = aed_censor_focalMOA)
ggsurvplot(fit_focal_same_moa,data= aed_censor_focalMOA, pval = TRUE, legend.labs = c("Different-MOA based combinations", "Same-MOA based combinations"),
                        #xlim = c(0,6000),
                           scale_colour_brewer = "Dark2")+
      ylab("Persistence probability") +
       xlab ("Time (days)")


#### Cox proportional hazard for all cases same vs diff moa
aed_censorMOA <- within(aed_censorMOA, MOA_comb_same <- factor(MOA_comb_same, labels = c("Same-MOA", "Different-MOA")) )


fit.coxph_same_diff <- coxph(surv_object_sameMOA  ~ MOA_comb_same , data = aed_censorMOA)
ggforest(fit.coxph_same_diff, data = aed_censorMOA)

summary(fit.coxph_same_diff)  -> sum_cox1
sum_cox1$coefficients -> prop1
rownames(prop1) = str_remove(rownames(prop1),"MOA_comb_same")
prop1 <- setDT(as.data.frame(prop1), keep.rownames =TRUE)
colnames(prop1) = c("combination","coef","coeff" ,"se","z_score","P")
sum_cox1$conf.int -> cox_same_diff
rownames(cox_same_diff) = str_remove(rownames(cox_same_diff),"MOA_comb_same")
cox_same_diff <- setDT(as.data.frame(cox_same_diff), keep.rownames = TRUE)
colnames(cox_same_diff) = c("combination","HR","coeff" ,"lower_CI","upper_CI")


cox_same_diff %>% select(-coeff) %>% left_join(prop1 %>% select(combination,P) %>% 
                                                 mutate( P = p_value_output(P), reference = "Same-MOA based combination"), 
                                               by= c("combination"="combination")) %>%
  select(combination,reference,HR ,lower_CI,upper_CI, P) -> cox_sameMOA

#### Cox proportional hazard for focal epilepsy cases same vs diff moa
aed_censor_focalMOA <- within(aed_censor_focalMOA, MOA_comb_same <- factor(MOA_comb_same, labels = c("Same-MOA", "Different-MOA")) )

fit.coxph_focal_same_diff <- coxph(surv_object_same_focals  ~ MOA_comb_same , data = aed_censor_focalMOA)
ggforest(fit.coxph_focal_same_diff, data = aed_censor_focalMOA)

summary(fit.coxph_focal_same_diff)  -> sum_cox1
sum_cox1$coefficients -> prop1
rownames(prop1) = str_remove(rownames(prop1),"MOA_comb_same")
prop1 <- setDT(as.data.frame(prop1), keep.rownames = TRUE)
colnames(prop1) = c("combination","coef","coeff" ,"se","z_score","P")
sum_cox1$conf.int -> cox_focal_same_diff
rownames(cox_focal_same_diff) = str_remove(rownames(cox_focal_same_diff),"MOA_comb_same")
cox_focal_same_diff <- setDT(as.data.frame(cox_focal_same_diff), keep.rownames = TRUE)
colnames(cox_focal_same_diff) = c("combination","HR","coeff" ,"lower_CI","upper_CI")


cox_focal_same_diff %>% select(-coeff) %>% left_join(prop1 %>% select(combination,P) %>% 
                                                 mutate( P = p_value_output(P), reference = "Same-MOA based combination"), 
                                               by= c("combination"="combination")) %>%
  select(combination,reference,HR ,lower_CI,upper_CI, P) -> cox_focal_sameMOA


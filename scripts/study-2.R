#Load packages
library(qmethod) #for Q analysis
library(psych) #for CFAs, reliability analyses, correlations,e tc.
library(RGenData) #for selecting factors
library(tidyverse) #for data wrangling
library(glmnet) #for regressions
library("parallel")
library("doParallel")
library("doRNG") #[DFP] Used to make dopar loops reproducible 

#variables
ready_for_analysis <- T #change this to TRUE once all data is collected and ready for analysis


######
#DEVIATIONS FROM PREREGISTRATION
######
#Deviations from preregistration are flagged in comments with [DFP]
#Reasons are provided

######
#PER THE PREREGISTRATION
######
#Some lines of code could not be written without knowledge of the data (e.g., classification of fast responders)
#Intent or rules for calculation are specified in preregistration
#These lines are flagged with [PTP]

############
#Data wrangling
############

############
#Anonymise data
############
rawDir <- "../data-raw/study-2.csv"
dataDir <- gsub('-raw', '', rawDir)
idDir <- "../data-raw/study-2-id.csv"
#Check if 'raw' data exists, if so, place an anonymised version of raw data in /data
#Public repository will not contain the raw file
if(file.exists(rawDir)){
  dat.raw <- read.csv(rawDir, stringsAsFactors = F, header = T)
  #Check if data is new by accessing list of IDs associated with existing data set
  #List of IDs is stored in a private directory, contains sensitive information
  if(file.exists(idDir)){
    exist.id <- read.csv(file = idDir, stringsAsFactors = F)
  }
  #Add unique ID
  set.seed(4583147) #set seed
  #Ignore first two rows (additional variable information)
  #Randomly assign a ID to each participant
  #completion_order assigns a ranking to each participant, indicating order of data collection (value 1 is 1st)
  dat.id <- dat.raw %>%
    slice(-c(1,2)) %>%
    mutate(completion_order= 1:nrow(.)) %>%
    sample_frac(1, replace = F) %>%
    mutate(id = 1:nrow(.)) %>%
    select(ResponseId, id, completion_order)
  #Map internal IDs with Qualtrics IDs, in case any participants require replacement
  #This will be only be visible in private repository
  dat.id %>%
    write.csv(file = idDir, row.names = F)
  #Join Id with ResponseId
  dat.raw <- dat.raw %>%
    left_join(dat.id, by = "ResponseId")
  #Remove variables which identify Qualtrics accounts
  #Qualtrics adds participant data to files
  #If sensitive variables are named differently than expected
  #the following lines of codes will require changing AFTER PREREGISTRATION
  #to adequately deidentify participants
    #Remove sensitive data
    sensitive_names <- c("IPAddress", "Status", "RecordedDate", "ResponseId", "RecipientLastName", 
                         "RecipientFirstName", "RecipientEmail", "ExternalReference", "LocationLatitude",	
                         "LocationLongitude",	"DistributionChannel",	"UserLanguage", "opp", "rid",	"RISN", "V")
    dat.raw <- dat.raw[ , !(colnames(dat.raw) %in% sensitive_names)] 
    #Export data
    dat.raw  %>%
      arrange(id) %>%
      arrange_at(vars(id), funs(desc(is.na(.)))) %>%
      write.csv(file = dataDir, row.names = F)
}
###########
#Clean data
##########
set.seed(54865) #set seed for reproducable code
#Qualtrics csv output is (m)essy
dat.m <- read.csv(dataDir, stringsAsFactors = F)
#1. Remove redundant rows
dat.m <- dat.m[-c(1:2), ] %>%
    as.tibble()
#2. Embedded data and other variables are used by Qualtrics as background variables
# Not required for analysis
embedVars <- c("UniqueID", "code", "Progress", "Finished", "gc", "Q_TotalDuration")
ev <- embedVars[{embedVars %in% colnames(dat.m)}]

dat.m <- dat.m %>%
    select(-ev)

#3. Rename variables
#Create function to assist
rename_v <- function(tb, orig, new){
  #Take tibble tb, replaces part of a variable name (orig) with a new part (new), for all vars
  tb %>%
    rename_at(vars(contains(!!orig)), 
              funs(str_replace(., orig, new))) %>%
    return()
}
#Rename variables
#Firstly, for Qualtrics administration of Qsort
dat.clean <- dat.m %>% 
  rename(time_start = StartDate) %>%
  rename(time_end = EndDate) %>%
  rename(time_to_complete = Duration..in.seconds.) %>%
  rename_v("Q2_", "pif") %>%
  rename(pcf = Q3) %>%
  rename_v("Q5_", "pcf_") %>%
  rename(age = Q7) %>%
  rename(gender = Q8) %>%
  rename_v("Q9_", "demographics_") %>%
  rename_v("Q11", "code_entered") %>%
  rename_v("Q12", "qualtrics_qsort") %>%
  rename_v("Q14", "qualtrics_suvery_info") %>%
  rename_v("CharID", "qualtrics_qsort_id") %>%
  rename_v("Q15_", "SS_") %>%
  rename_v("Q16_", "SS_") %>%
  rename(PA = Q17) %>%
  rename_v("Q18_", "PA_") %>%
  rename_v("Q19_", "NCS_6_") %>%
  rename_v("Q20_", "NCS_6_") %>%
  rename(MMS_human = Q21) %>%
  rename_v("Q22_", "MMS_human_") %>%
  rename_v("Q23_", "MMS_cause_") %>%
  rename_v("Q24_", "MMS_cause_") %>%
  rename_v("Q25_", "MMS_consequence_") %>%
  rename_v("Q26_", "MMS_consequence_") %>%
  rename_v("Q27_", "MMS_mitigation_") %>%
  rename_v("Q28_", "MMS_mitigation_") %>%
  rename_v("Q29_", "SVSS_") %>%
  rename_v("Q30_", "SVSS_") %>%
  rename(KV = Q31) %>%
  rename_v("Q32_", "KV_") %>%
  rename_v("Q33_", "SJ_") %>%
  rename_v("Q34_", "SJ_") %>%
  rename_v("Q35_", "CI_") %>%
  rename_v("Q36_", "CI_") %>%
  rename(W = Q37) %>%
  rename_v("Q38_", "W_") %>%
  rename_v("Q39_", "CFC_S_") %>%
  rename_v("Q40_", "CFC_S_") %>%
  rename_v("Q41_", "BFI_10_") %>%
  rename_v("Q42_", "BFI_10_") %>%
  rename_v("Q43_", "EWS_") %>%
  rename_v("Q44_", "EWS_") %>%
  rename_v("Q999_", "specs_")
#4. Import Q sort data
dat.m_qsort <- read.csv("../data-raw/study-2-qsort.csv", stringsAsFactors = F, header = F)
#Assign data column names, based on schema used to create MySQL db, see ../study2/data/schema.txt
colnames(dat.m_qsort) <- c(
  paste0("pref_sta_", 1:30),
  paste0("reason_sta_", 1:30),
  paste0("sort_sta_", 1:30),
  paste0("time_", c("instruct", "pilot", "pref", "reason", "sort")),
  paste0("count_", c("noroom", "incon")),
  paste0("supress_", c("noroom", "incon")),
  paste0("pilot_response", 1:2),
  paste0("place_", c("grid", "avail", "spare", "spare_max")),
  "qs",
  "pass",
  "time_start_qsort"
)
dat.m_qsort <- dat.m_qsort %>% 
  as.tibble()
#Convert passcodes and query strings into CHARACTER
passcodes_to_character <- function(c){
  #Turns code (c, integer) into a 7-digit string
  c <- as.character(c)
  if(is.na(c)){return(NA)}
  short_by <- 7-nchar(c)
  if(short_by > 0){
    rep(0, short_by) %>%
      paste(collapse = "") %>%
      paste0(c) %>%
      return()
  } else{
    return(c)
  }
}
dat.m_qsort <- dat.m_qsort %>%
  rowwise %>%
  mutate_at(vars(qs, pass), passcodes_to_character)
dat.clean <- dat.clean %>%
  rowwise() %>%
  mutate_at(vars(qualtrics_qsort_id, code_entered), passcodes_to_character)

#5. Match Qualtrics issued IDs to Q sort information
match_qsort_ids <- function(q, s, t1, t2, t3, t4, t5){
  #Function to match Q sort data to a unique participant ID (id var in dat.clean)
  #Uses timing information (start time (s) and phase times (t1:t5)) and issued query string (q)
  #Returns a participant ID
  #1. Identify matching query strings
  match <- dat.clean %>%
    filter(qualtrics_qsort_id == q)
  if(nrow(match) != 1){
    #Duplicates found
    #Refine further using timing data
    #Participant must have started survey before starting Q sort
    match <- match %>%
      filter(as.POSIXct(time_start, format = '%Y-%m-%d %H:%M:%S') <= as.POSIXct(s, format = "%Y-%m-%d %H:%M:%S")) #[DFP] line was changed due to formating error
    #Participant must have completed Q sort before completing survey
    match <- match %>%
      filter(as.POSIXct(time_end, format = '%Y-%m-%d %H:%M:%S') >= (as.POSIXct(s, format = "%Y-%m-%d %H:%M:%S")+t1+t2+t3+t4+t5)) #[DFP] line was changed due to formating error
  }
  if(nrow(match) != 1){
    #Then data cannot be successfully matched
    #Assign an ID of -1
    return(-1)
  }
  else{
    #Data has been successfully matched
    #Return matched id
    match %>%
      pull(id) %>%
      return
  }
}
dat.m_qsort <- dat.m_qsort %>%
  rowwise %>%
  mutate(id = match_qsort_ids(qs, time_start_qsort, time_instruct, time_pilot, time_pref, time_reason, time_sort)) %>%
  group_by(id) %>%
  mutate(id_is_diagnostic = !(n() > 1)) %>% 
  rowwise %>% 
  mutate(id = ifelse(id_is_diagnostic, id, -1))  %>% 
  select(-id_is_diagnostic)
#Check for any duplicate IDs
#If any exist, process has not been diagnostic and participants should be removed/replaced



#6. Join Q sort data with Qualtrics data
dat.clean <- dat.clean %>%
  full_join(dat.m_qsort, by = "id")
#7. Remove any excess Q sort data
dat.clean <- dat.clean %>%
  mutate(pcf = as.integer(pcf)) %>%
  filter(pcf == 1)
#8. Export a clean data file
#Export to /out
cleanDir <- "../out/study-2-clean.csv"
dat.clean %>%
  write_csv(cleanDir)


#######
#Exclusion criteria
#######

#1. Exclude fast responders
#[PTP] classification of fast responders
dat.clean <- dat.clean %>%
  rowwise() %>%
  mutate(fast_responder = as.integer(time_to_complete) < 872.5)
dat.clean <- dat.clean
#Throw out data from all irrelevant participants
dat.clean <- dat.clean %>%
  filter(fast_responder == F | is.na(fast_responder)) %>%
  filter(id > 0)
#2. Check number of responses for each part of study
n_total <- dat.clean %>%
  nrow()

#####################
#Survey scale scoring
#####################
rscore <- function(score, max, min){
  #Reverse scores data, given score and a maximum and minimum of total scale
  if(missing(min)){min <- 1}
  return(max - score + min)
}
dat.clean <- dat.clean %>%
  mutate_at(vars(SS_1:latin), as.numeric) %>%
  rowwise %>%
  mutate(BFI_10_E = mean(c(rscore(BFI_10_1, 5), BFI_10_6))) %>%
  mutate(BFI_10_A = mean(c(rscore(BFI_10_7, 5), BFI_10_2))) %>%
  mutate(BFI_10_C = mean(c(rscore(BFI_10_3, 5), BFI_10_8))) %>%
  mutate(BFI_10_N = mean(c(rscore(BFI_10_4, 5), BFI_10_9))) %>%
  mutate(BFI_10_O = mean(c(rscore(BFI_10_5, 5), BFI_10_10))) %>%
  mutate(CFC_S_F = mean(c(CFC_S_1, CFC_S_2, CFC_S_5, CFC_S_6))) %>%
  mutate(CFC_S_I = mean(c(CFC_S_3, CFC_S_4, CFC_S_7, CFC_S_8, CFC_S_9))) %>%
  mutate(CI = mean(c(CI_1, CI_2, CI_3, CI_4, CI_5, CI_6))) %>%
  mutate(EWS_E = mean(c(EWS_1, EWS_3, EWS_5, EWS_6, EWS_7, EWS_10))) %>%
  mutate(EWS_D = mean(c(EWS_2, EWS_4, EWS_8, EWS_9, EWS_11, EWS_12))) %>%
  mutate(NCS_6 = mean(c(NCS_6_1, NCS_6_2, rscore(NCS_6_3, 5), rscore(NCS_6_4, 5), NCS_6_5, NCS_6_6))) %>%
  mutate(PA_r = PA - 4) %>%
  mutate(SJ = mean(c(SJ_1, SJ_2, rscore(SJ_3, 9), SJ_4, SJ_5, SJ_6, rscore(SJ_7, 9), SJ_8))) %>%
  mutate(SS_E = mean(c(SS_1, SS_2, SS_4, SS_5, SS_8, SS_9, SS_13, SS_15))) %>%
  mutate(SS_R = mean(c(SS_3, SS_6, SS_7, SS_10, SS_11, SS_12, SS_14))) %>%
  mutate(SVSS_C = .82+.05*SVSS_1+.06*SVSS_2-.04*SVSS_3-.09*SVSS_4-.18*SVSS_5-.16*SVSS_6+.03*SVSS_7+.16*SVSS_8+.18*SVSS_9+.11*SVSS_10) %>%
  mutate(SVSS_ST= -0.60-.19*SVSS_1-.14*SVSS_2-.09*SVSS_3-.11*SVSS_4+.01*SVSS_5+.10*SVSS_6+.13*SVSS_7+.07*SVSS_8+.06*SVSS_9+.02*SVSS_10) %>%
  mutate(MMS_human_r = MMS_human - 4) %>%
  mutate(MMS_consequence_P = mean(c(MMS_consequence_2, MMS_consequence_6, MMS_consequence_11))-4) %>%
  mutate(MMS_consequence_S = mean(c(MMS_consequence_1, MMS_consequence_3, MMS_consequence_4, MMS_consequence_5, MMS_consequence_7, MMS_consequence_8, MMS_consequence_9, MMS_consequence_10))-4) %>%
  mutate(MMS_mitigation_C = mean(c(MMS_mitigation_1, MMS_mitigation_7, MMS_mitigation_8))-4) %>%
  mutate(MMS_mitigation_E = mean(c(MMS_mitigation_2, MMS_mitigation_5, MMS_mitigation_10))-4) %>%
  mutate(MMS_mitigation_G = mean(c(MMS_mitigation_3, MMS_mitigation_4, MMS_mitigation_6, MMS_mitigation_9, MMS_mitigation_11))-4)

###Factor analysis of MMS_cause
mms_cause_composite_scales_completed <- T #change to TRUE when composites have been calculated
if(!mms_cause_composite_scales_completed){
  #Number of factors to extract
  mms_cause_factor_num <- dat.clean %>% 
    filter(!is.na(MMS_cause_1)) %>% 
    select(MMS_cause_1:MMS_cause_13) %>% 
    as.data.frame() %>%
    {capture.output(EFACompData(., f.max = 10))} %>%
    substr(nchar(.)-1, nchar(.)-1) %>%
    as.integer()
  #Calculate factor structure
  mms_cause_factor_structure <- dat.clean %>% 
    filter(!is.na(MMS_cause_1)) %>% 
    select(MMS_cause_1:MMS_cause_13) %>% 
    as.data.frame() %>%
    principal(mms_cause_factor_num, rotate = "varimax")
  stop("Look at the following output and determine composite scores")
}else{
  mms_cause_factor_structure <- dat.clean %>% 
    select(MMS_cause_1:MMS_cause_13) %>% 
    as.data.frame() %>%
    principal(3, rotate = "varimax")
  #When study is analysed (FOLLOWING PREREGISTRATION), this code will need to be uncommented & ammended to calculate composite scores
  dat.clean <- dat.clean %>%
    rowwise %>%
    mutate(MMS_cause_C = mean(c(MMS_cause_1, MMS_cause_2, MMS_cause_3, MMS_cause_4, MMS_cause_7, MMS_cause_8, MMS_cause_11))) %>%
    mutate(MMS_cause_E = mean(c(MMS_cause_5, MMS_cause_6, MMS_cause_9, MMS_cause_10))) %>%
    mutate(MMS_cause_N = mean(c(MMS_cause_12, MMS_cause_13)))
}


############
#RQ1 (Q analysis)
############
if(ready_for_analysis & n_total >= 400){
#Convert statement positions to integers
#Not this code only works when < 10 columns
if(is.character(pull(select(dat.clean, sort_sta_1)))){
dat.clean <- dat.clean %>%
  mutate_at(vars(starts_with('sort_sta_')), function(x){return(as.integer(substr(x, 2, 2)) - 5)}) %>%
  mutate_at(vars(starts_with('sort_sta_')), as.integer)
}

#Contains information on statements, where each row is a seperate statement, columns are individual participants and statement text
#[PTP] updated code to run with newer version of tidyverse
statement_information <- dat.clean %>% 
  filter(!is.na(sort_sta_1)) %>% 
  select(id,starts_with('sort_sta_') ) %>% 
  gather(key = statement, vars, -id) %>%
  mutate(id = paste0('id_', id)) %>%
  spread(id, vars)

statement_information <- statement_information %>%
  rowwise %>%
  mutate(order = as.integer(substr(statement, 10, nchar(statement))))

statement_text <- read.csv('../study2/study/q-statements-ordered.csv', stringsAsFactors = F, header = T) %>%
  select(-id) %>%
  rename(text = statement)

statement_information <- statement_information %>%
  full_join(statement_text, by = 'order') %>%
  ungroup %>% 
  arrange(order)
#Need to convert to specific format for QMETHOD package
q_sorts <- statement_information %>% 
  select(-order,-text,-statement) %>% 
  as.data.frame
rownames(q_sorts) <- paste0('sta_', rownames(q_sorts))


#Determine factors to retain for CFA, using technique from Ruscio and Roche (2012)
#[DFP] EFACompData algorithm unstable, decision to retain factors based on scree plot
q_factor_num <- 1 #[DFP] this is derived from the code below. Comment out this line to confirm
if(!exists('q_factor_num')){ #[DFP] added for functionality
q_factor_num <- capture.output(EFACompData(q_sorts, f.max = 10)) %>%
  substr(nchar(.)-1, nchar(.)-1) %>%
  as.integer()
} #[DFP] added for functionality


#Calculate loading matrix
q_fa <- principal(q_sorts, q_factor_num, rotate = "varimax")
q_loa <- q_fa$loadings %>%
  unclass

#Flag Q sorts for factor with highest loading
#[DFP] calculate loads for bipolar factors
#Following Brown (1980, p. 253)
q_flagged <- matrix(F, ncol(q_sorts), q_factor_num*2)
n <- 1 #iterate this over the loop
tol <- 1e-8 #tolerance value for comparing doubles
set.seed(1746985) #set seed
while(n <= ncol(q_sorts)){
  #[DFP] added variable track flagged factor and flagged pole
  #Track factor which largest loading
  flag_f <- q_loa[n, ] %>%
    abs %>%
    {which(. > max(.) - tol)}
  #Track sign of loading: T for +ve loading, F for -ve loading
  flag_p <- q_loa[n, flag_f] >= 0-tol
  #check in case tiebreaker (commented this out, as I think it is irrelevant)
  # if(length(which(q_flagged[n, ])) > 1){
  #   #Then randomly select one TRUE value to remain TRUE, others substituted with FALSE
  #   q_flagged[n, sample(which(q_flagged[n, ]))[-1]] <- F
  # }
  #Assign flag
  if(abs(q_loa[n, flag_f]) >= 1.96/sqrt(30) - tol){
  q_flagged[n, ifelse(flag_p, flag_f, flag_f+q_factor_num)] <- T
  }
  #Incriment n
  n <- n + 1
}

#[DFP] split loadings into bipolar reflections
q_loa <- q_loa %>%
  cbind(-q_loa)

colnames(q_flagged) <- paste0("flag_f", 1:(q_factor_num*2)) #to match q_loa
row.names(q_flagged) <- row.names(q_loa) #to match q_loa
colnames(q_loa) <- paste0("f", 1:(q_factor_num*2)) #to match q_flagged

#Run analysis
q_results <- qzscores(q_sorts, q_factor_num*2, q_loa, q_flagged) #runs full Q analysis
q_results$qdc <- qdc(q_sorts, q_factor_num*2, q_results$zsc, q_results$f_char$sd_dif) #calculates consensus/distinguishing statements

#[ATP] export some of the Q sort results:
df_qres <- q_results$zsc %>%
  cbind(q_results$zsc_n) %>%
  rownames_to_column(var = "sta_id")
df_qres %>%
  write.csv('../out/study-2-factor-scores.csv', row.names = F)

#########
#Construct crib sheet
#########

#Initialise crib sheet data
cribble <- tibble(
  n = 1:nrow(q_sorts),
  id = paste0('sta_', n),
  text = as.character(statement_information$text),
  dist.and.cons = q_results$qdc$dist.and.cons
) %>%
  right_join(data.frame(as.data.frame(q_results$zsc_n), id = rownames(q_results$zsc_n)), by = 'id')

#Determine consensus statements
cribble <- cribble %>%
  rowwise %>%
  mutate(cons = (dist.and.cons == "Consensus"))

#Determine which statements distinguish a factor over all others (distinguishing statements)
if((q_factor_num*2) <= 2){
  cribble <- cribble %>% 
    mutate(dist.and.cons = ifelse(dist.and.cons == "Distinguishing", "Distinguishing all", dist.and.cons))
}
cribble <- cribble %>%
  mutate(dist = ifelse(str_detect(dist.and.cons, "all"), list(paste0("f", 1:(q_factor_num*2))),  str_extract_all(dist.and.cons, "f\\d*"))) %>%
  unnest %>%
  select(n, dist) %>%
  mutate(dist = paste0("dist_", dist)) %>%
  mutate(d = T) %>%
  spread(dist, d, fill = F) %>%
  right_join(cribble, by = "n") %>%
  mutate_at(vars(starts_with("dist_")), ~replace_na(., F)) %>%
  select(-dist.and.cons)


#For Q factor interpretation, create empty comment variables for ea. factor
cribble[,paste0("com_f", rep(1:(q_factor_num*2)))] <- ""


#Recode this to a factor
#Different crib sheets are created for each factor
#Items are displayed on the crib sheet if they are salient features which facilitate holistic interpretation of a factor
#To begin, these are considered to be extremely high/low factor scores, or...
#...Rankings higher than other factor arrays
crib_factor_levels = c("Not on Crib Sheet",
                       "Items Ranked at Lowest Extreme",
                       "Items Ranked Lower in this Factor Array than in Other Factor Arrays",
                       "Items Ranked Higher in this Factor Array than in Other Factor Arrays",
                       "Items Ranked at Highest Extreme"
)
#Create variable to indicate which items to display for each factors' crib sheet
cribble[,paste0("crib_f", rep(1:(q_factor_num*2)))] <- crib_factor_levels[1]
cribble <- cribble %>%
  mutate_at(vars(starts_with("crib_f")), funs(factor(., levels = crib_factor_levels)))



#Create functions to convert cribble from wide to long (and long to wide) forms
wide_to_long_cribble <- function(c){
  #Take input cribble, and converts it to long form (each row is a statement, for each factor)
  c %>%
    gather(key, value, -n,-id,-text,-cons) %>%
    separate(key, into = c("var", "f"), sep = "_f") %>%
    spread(key = var, value = value) %>%
    mutate_at(vars(f, fsc), funs(as.integer)) %>%
    mutate_at(vars(dist), funs(as.logical)) %>%
    mutate(crib = factor(crib, levels = crib_factor_levels)) %>%
    return()
}
long_to_wide_cribble <- function(c){
  #Takes input cribble, and converts to wide form
  c %>%
    gather(key = var, value = value, -n, -id, -text, -cons, -f) %>%
    unite("key", var, f, sep = "_f") %>%
    spread(key, value) %>%
    mutate_at(vars(starts_with("dist")), funs(as.logical)) %>%
    mutate_at(vars(starts_with("fsc")), funs(as.integer)) %>%
    mutate_at(vars(starts_with("crib_f")), funs(factor(., levels = crib_factor_levels))) %>%
    return()
}

#Create a long form cribble
cribble_long <- cribble %>%
  wide_to_long_cribble

#Create a function to extract the factor scores for a statement (item)
#Does not include factor scores belonging to the factor in question
other_fsc <- function(item, factor){
  #Takes an item (n) and factor (f) from cribble_long
  pull(filter(cribble_long, f != factor, n == item), fsc)
}

#For each factor, determine the highest & lowest factor score for statements (not of that factor)
cribble_long <- cribble_long %>%
  rowwise %>%
  mutate(othermax = max(other_fsc(n, f))) %>%
  mutate(othermin = min(other_fsc(n, f)))

#Determine the most extreme points of the Q sort distribution
q_extremes <- cribble_long %>%
  pull(fsc) %>%
  {c(min(.), max(.))}

#Determine whether item will be initalised as displaying on the crib sheet
cribble_long <- cribble_long %>%
  mutate(crib = as.character(crib)) %>%
  mutate(crib = ifelse(fsc < othermin, 3, 1)) %>%
  mutate(crib = ifelse(fsc > othermax, 4, crib)) %>%
  mutate(crib = ifelse(fsc == q_extremes[1], 2, crib)) %>%
  mutate(crib = ifelse(fsc == q_extremes[2], 5, crib)) %>%
  mutate(crib = factor(crib, levels = 1:5, labels = crib_factor_levels)) %>%
  select(-othermax, -othermin)


#Convert cribble_long back into wide format, to save to an external file
cribble <- cribble_long %>%
  long_to_wide_cribble()

#Export cribble, for use by ShinyApp to interpret factors
cribbleDir <- '../out/study-2-cribsheet.csv'
if(!file.exists(cribbleDir)){
  print("Crib sheets will be created (.csv files)")
  
  #Export factor array
  cribble %>%
    write.csv(file = cribbleDir, row.names = F)
  
  #Export a crib sheet for each factor
  sapply(1:(q_factor_num*2), function(x){
    #Initalise names of variables for factor x
    crib_name <- paste0('crib_f', x)
    fsc_name <- paste0('fsc_f', x)
    com_name <- paste0('com_f', x)
    dist_name <- paste0('dist_f', x)
    
    #Create path for exported crib sheet
    cribsheetDir <- cribbleDir %>% 
      substr(1, nchar(.) - 4) %>%
      paste0('-f', x ,'.csv')
    
    
    #Create a crib sheet for factor x
    cribble %>%
      select(UQ(as.name(crib_name)), n, text, UQ(as.name(fsc_name)), cons,
             UQ(as.name(dist_name)), UQ(as.name(com_name))
             ) %>%
      arrange(desc(UQ(as.name(crib_name))), 
              desc(UQ(as.name(fsc_name)))
      ) %>%
      write.csv(file = cribsheetDir, row.names = F)
    
    #Output message
    sprintf('Crib sheet for factor %d was exported', x) %>%
      return()
  })
  
}

#Add Q sort results to dat.clean file

q_participant_loa <- q_results$loa %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'id') %>% 
  as.tibble %>%
  mutate(id = substr(id, 4, nchar(id))) %>%
  mutate(id = as.integer(id))

q_participant_flagged <- q_results$flagged %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'id') %>% 
  as.tibble %>%
  mutate(id = substr(id, 4, nchar(id))) %>%
  mutate(id = as.integer(id))

#Create variable 'flag' to denote which Q factor participant was flagged as
#This variable is the DV in the following ridge regressions
q_participant_flagged <- q_participant_flagged %>% 
  gather(key = key, value = value, -id) %>% 
  filter(value) %>% 
  mutate(flag = substr(key, 6, nchar(key))) %>% 
  select(id, flag) %>% 
  right_join(q_participant_flagged, by = 'id') %>% 
  mutate(flag = as.factor(flag))

dat.clean <- dat.clean %>%
  left_join(q_participant_loa, by = 'id') %>%
  left_join(q_participant_flagged, by = 'id')

#Assign labels
dat.clean <- dat.clean %>%
  mutate(flag = ifelse(is.na(flag), 3, flag)) %>%
  mutate(flag = factor(flag, levels = 1:3, labels = c("f1", "f2", "f3")))
#f1 = Acceptors
#f2 = Sceptics
#f3 = Fencesitters

###Check comments is intended to be used with crib sheet
check_comments <- function(f, sta, pos){
  #For a given factor (f), find the qualitative comments which justify the extreme positions (pos) of statements
  #User can limit these to a particular statement (sta), indicated by an integer
  #Output is a tibble with participant comments who load onto factor f, ordered according to loadings
  f_varname <- paste0('f', f) #[DFP] changed this to work as intended
  flag_varname <- paste0('flag_f', f)
  if(missing(pos)){pos <- 4}
  if(missing(sta)){sta <- 1:30}
  
  #Output tibble
  o <- dat.clean %>% 
    filter(flag == f_varname) %>% #[DFP] static variable used, changed to function as intended 
    select(id, !!f_varname, sort_sta_1:sort_sta_30, reason_sta_1:reason_sta_30) %>%
    gather(key, value, -!!f_varname, -id) %>%
    mutate(order = as.integer(sub('^.*sta_', "", key))) %>%
    mutate(key = sub('_sta_.*$', "", key)) %>%
    spread(key, value) %>%
    mutate(sort = as.integer(sort)) %>%
    filter(sort %in% pos) %>%
    mutate(reason = sub("&rsquo;", "\'", reason)) %>%
    filter(!is.na(reason)) %>%
    filter(order %in% sta)
  
  #Join with statement text
   o %>%
     left_join(statement_text, by = 'order') %>%
     arrange(desc(UQ(as.name(f_varname)))) %>%
     rename(sta_text = text) %>%
     return()
   
   #If no matches, output is an empty tibble
} 

###
#RQ 2
###

##Ridge regression

#Specify predictors
#[PTP] added causes subscales
#[DFP] scaled predictors for better interpretability
#[DFP] changed sign of mitigation scales for better interpretability
#[DFP] Scale Mitigation and Consequence scales for better interpretability; though, this will not effect GLM estimates of interest

######
#Scaling variables:
preds <- dat.clean %>%
  select(MMS_human,
         MMS_cause_C,
         MMS_cause_E,
         MMS_cause_N,
         MMS_consequence_P,
         MMS_consequence_S,
         MMS_mitigation_G,
         MMS_mitigation_C,
         MMS_mitigation_E) %>%
  ungroup %>%
  mutate_at(vars(starts_with('MMS_mitigation')), .funs = list(~ -.+4)) %>%
  mutate_at(vars(starts_with('MMS_consequence_')), .funs = list(~ .+4)) %>%
  mutate_all(scale, center = TRUE, scale = TRUE) %>%
  as.matrix



#To be predicted by the following ridge regressions
flags <- pull(dat.clean, flag)

glmres <- glmnet(preds, flags,
                 family="multinomial", alpha = 0)
cvres <- cv.glmnet(preds, flags,
                   family="multinomial", alpha = 0)
print('Ridge regression for all predictors')
print(summary(cvres$glmnet.fit))
print(coef(glmres, cvres$lambda.min))
print(coef(glmres, cvres$lambda.1se))
glmres_mm_scaled <- glmres
cvres_mm_scaled <- cvres


###Estimating confidence intervals with bootstrap procedure
if(!file.exists('../out/study-2-glm_mm_scaled.csv')){
  cl <- makeCluster(detectCores()-1) # create a cluster with 2 cores
  registerDoParallel(cl)
  
  nboots <- 1000
  #set.seed(00134) [DFP] This does not function as intended. Does not make analysis reproducable. Instead, use a dorng loop below
  tt <- {}
  tt <- foreach (brep=1:nboots,
                 .packages = c("tidyverse","glmnet"),
                 .options.RNG=00134) %dorng% {
                   write(brep, "../out/parlog.txt") #take current bootstrap
                   
                   #generate bootstrap
                   #select data 
                   dat.boot <- dat.clean %>%
                     sample_frac(1, replace = T)
                   #extract predictors
                   #[PTP] added causes scales
                   #[DFP] added Worry and Knowledge Volume scales: initially accidently ommited from code, despite declaration in online registration
                   preds.boot <- dat.boot %>%
                     select(MMS_human,
                            MMS_cause_C,
                            MMS_cause_E,
                            MMS_cause_N,
                            MMS_consequence_P,
                            MMS_consequence_S,
                            MMS_mitigation_G,
                            MMS_mitigation_C,
                            MMS_mitigation_E) %>%
                     ungroup %>%
                     mutate_all(.funs = list(~. - mean(.))) %>%
                     mutate_at(vars(starts_with('MMS_mitigation')), .funs = list(~ -.+4)) %>%
                     mutate_at(vars(starts_with('MMS_consequence_')), .funs = list(~ .+4)) %>%
                     mutate_all(scale, center = TRUE, scale = TRUE) %>%
                     as.matrix
                   #extract flags
                   flags.boot <- pull(dat.boot, flag) #[DFP] incorrectly specified variable 'flags'
                   #Run ridge regression
                   tglm <- glmnet(preds.boot, flags.boot,
                                  family="multinomial", alpha=0)
                   tcvres <- cv.glmnet(preds.boot, flags.boot,
                                       family="multinomial", alpha=0)
                   return(coef(tglm, tcvres$lambda.min))
                 }
  stopCluster(cl)
  
  glmCI <- {}
  
  for (distr in 1:length(levels(flags))){
    kk <- lapply(tt,function(x){x[[distr]]})
    tCI <- mat.or.vec(dim(preds)[2],2)
    colnames(tCI) <- paste0(c("CI_.025_", "CI_.975_"), "f", distr)
    for (vv in 1:dim(preds)[2]){
      vest <- lapply(kk,function(x) x[vv+1])
      tCI[vv,] <- quantile(unlist(vest),probs = c(.025,.975))
    }
    glmCI <- cbind(glmCI,tCI)
  }
  row.names(glmCI) <- colnames(preds)
  glmCI_mm_scaled <- glmCI
  #[DFP] added output file to save time on repeated runs of script
  glmCI_mm_scaled %>% 
    as.data.frame() %>% 
    rownames_to_column('pred') %>%
    write.csv(file = '../out/study-2-glm_mm_scaled.csv', row.names = F)
  print(glmCI_mm_scaled)
}else{
  glmCI_mm_scaled <- read.csv(file = '../out/study-2-glm_mm_scaled.csv', stringsAsFactors = F)
}


###
#RQ 3
###

#Specify predictors
#[PTP] added causes subscales
#[DFP] added Worry and Knowledge volume. Despite being declared in preregistration, variables were accidently omitted
#[DFP] scaled predictors for better interpretability
#[DFP] changed sign of mitigation scales for better interpretability

preds <- dat.clean %>%
  select(MMS_human,
         MMS_cause_C,
         MMS_cause_E,
         MMS_cause_N,
         MMS_consequence_P,
         MMS_consequence_S,
         MMS_mitigation_G,
         MMS_mitigation_C,
         MMS_mitigation_E,
         BFI_10_E:SVSS_ST, W, KV) %>%
  ungroup %>%
  mutate_at(vars(starts_with('MMS_mitigation')), .funs = list(~ -.+4)) %>%
  mutate_at(vars(starts_with('MMS_consequence_')), .funs = list(~ .+4)) %>%
  mutate(PA_r = PA_r + 4) %>%
  mutate_all(scale, center = TRUE, scale = TRUE) %>%
  as.matrix


#To be predicted by the following ridge regressions
flags <- pull(dat.clean, flag)

glmres <- glmnet(preds, flags,
                 family="multinomial", alpha = 0)
cvres <- cv.glmnet(preds, flags,
                   family="multinomial", alpha = 0)
print('Ridge regression for all predictors')
print(summary(cvres$glmnet.fit))
print(coef(glmres, cvres$lambda.min))
print(coef(glmres, cvres$lambda.1se))
glmres_all_scaled <- glmres
cvres_all_scaled <- cvres


###Estimating confidence intervals with bootstrap procedure
if(!file.exists('../out/study-2-glm_all_scaled.csv')){
  cl <- makeCluster(detectCores()-1) # create a cluster with 2 cores
  registerDoParallel(cl)
  
  nboots <- 1000
 #set.seed(00134) [DFP] This does not function as intended. Does not make analysis reproducable. Instead, use a dorng loop below
  tt <- {}
  tt <- foreach (brep=1:nboots,
                 .packages = c("tidyverse","glmnet"),
                 .options.RNG=00134) %dorng% {
                   write(brep, "../out/parlog.txt") #take current bootstrap
                   
                   #generate bootstrap
                   #select data 
                   dat.boot <- dat.clean %>%
                     sample_frac(1, replace = T)
                   #extract predictors
                   #[PTP] added causes scales
                   #[DFP] added Worry and Knowledge Volume scales: initially accidently ommited from code, despite declaration in online registration
                   preds.boot <- dat.boot %>%
                     select(MMS_human,
                            MMS_cause_C,
                            MMS_cause_E,
                            MMS_cause_N,
                            MMS_consequence_P,
                            MMS_consequence_S,
                            MMS_mitigation_G,
                            MMS_mitigation_C,
                            MMS_mitigation_E,
                            BFI_10_E:SVSS_ST, W, KV) %>%
                     ungroup %>%
                     mutate_at(vars(starts_with('MMS_mitigation')), .funs = list(~ -.+4)) %>%
                     mutate_at(vars(starts_with('MMS_consequence_')), .funs = list(~ .+4)) %>%
                     mutate(PA_r = PA_r + 4) %>%
                     mutate_all(scale, center = TRUE, scale = TRUE) %>%
                     as.matrix
                   #extract flags
                   flags.boot <- pull(dat.boot, flag) #[DFP] incorrectly specified variable 'flags'
                   #Run ridge regression
                   tglm <- glmnet(preds.boot, flags.boot,
                                  family="multinomial", alpha=0)
                   tcvres <- cv.glmnet(preds.boot, flags.boot,
                                       family="multinomial", alpha=0)
                   return(coef(tglm, tcvres$lambda.min))
                 }
  stopCluster(cl)
  
  glmCI <- {}
  
  for (distr in 1:length(levels(flags))){
    kk <- lapply(tt,function(x){x[[distr]]})
    tCI <- mat.or.vec(dim(preds)[2],2)
    colnames(tCI) <- paste0(c("CI_.025_", "CI_.975_"), "f", distr)
    for (vv in 1:dim(preds)[2]){
      vest <- lapply(kk,function(x) x[vv+1])
      tCI[vv,] <- quantile(unlist(vest),probs = c(.025,.975))
    }
    glmCI <- cbind(glmCI,tCI)
  }
  row.names(glmCI) <- colnames(preds)
  glmCI_all_scaled <- glmCI
  #[DFP] added output file to save time on repeated runs of script
  glmCI_all_scaled %>% 
    as.data.frame() %>% 
    rownames_to_column('pred') %>%
    write.csv(file = '../out/study-2-glm_all_scaled.csv', row.names = F)
  print(glmCI_all_scaled)
}else{
  glmCI_all_scaled <- read.csv(file = '../out/study-2-glm_all_scaled.csv', stringsAsFactors = F)
}


#################################
#####[ATP]Figures and Tables#####
#################################
index_fig <- 2 #to index figures (e.g., Figure 2)
index_tbl <- 2 #to index tables (e.g., Table 2)
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k)) #to round numbers

####Scree plot of Q sort components
apatheme=theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        text=element_text(family='sans', color='black'),
        legend.title=element_blank(),
        legend.position=c(.7,.8),
        axis.line.x = element_line(color='black'),
        axis.line.y = element_line(color='black'),
        axis.text = element_text(colour = 'black'))
factor_range <- 1:15
tibble(num = factor_range, eigen = head(q_fa$values, tail(factor_range, 1))) %>%
  ggplot(aes(x = num, y = eigen)) +
  geom_line()+
  geom_point(size = 4)+
  scale_y_continuous(name='Eigenvalue')+
  scale_x_continuous(name='Factor number', breaks = factor_range)+
  geom_vline(xintercept = 1.5, linetype = 'dashed')+
  apatheme
ggsave(sprintf('../out/figures/paper-figure-%s.png', index_fig), width=5, height=5, unit='in', dpi=300)
index_fig <- index_fig + 1 #increment counter

####Correlations of psychological characteristics
cormat <- dat.clean %>%
  select(MMS_human,
         MMS_cause_C,
         MMS_cause_E,
         MMS_cause_N,
         MMS_consequence_P,
         MMS_consequence_S,
         MMS_mitigation_G,
         MMS_mitigation_C,
         MMS_mitigation_E,
         BFI_10_E:SVSS_ST, W, KV
  ) %>%
  mutate_at(vars(starts_with('MMS_mitigation')), .funs = list(~ -.)) %>%
  as.matrix %>%
  cor

construct_terms <- tibble(
  full = rownames(cormat),
  abbrev = c(paste0('MMS_', c('HUM', 'CAU_C', 'CAU_E', 'CAU_N', 'CON_P', 'CON_S', 'MIT_G', 'MIT_C', 'MIT_E')),
             paste0('BFI_', c('E', 'A', 'C', 'N', 'O')),
             'CFC_F', 'CFC_I',
             'CI',
             'EWS_E', 'EWS_D',
             'NCS', 'PI', 'SJ',
             'SS_E', 'SS_R',
             'SVSS_C', 'SVSS_ST',
             'W', 'KV'))

rownames(cormat) <- construct_terms$abbrev
colnames(cormat) <- construct_terms$abbrev


df <- read_csv('../paper/assets/scales.csv', col_types = cols(
  abbrev = col_character(),
  measure = col_character(),
  example = col_character(),
  reference = col_character(),
  range = col_character(),
  num = col_integer()
))

construct_terms  <- construct_terms %>%
  arrange(abbrev) %>%
  left_join(df, by = 'abbrev')


df <- dat.clean %>%
  select(id, flag, construct_terms$full) %>%
  pivot_longer(BFI_10_A:W) %>%
  mutate(value = ifelse(grepl('MMS_consequence_', name), value+4, value)) %>%
  mutate(value = ifelse(grepl('MMS_mitig', name), -value+4, value)) %>%
  mutate(value = ifelse(name == 'PA_r', value+4, value))

#for grand mean
construct_terms <- df %>%
  group_by(name) %>%
  summarise(m = mean(value), sd = sd(value)) %>%
  mutate(m = specify_decimal(m, 2)) %>%
  mutate(sd = specify_decimal(sd, 2)) %>%
  mutate(des = paste0(m, ' (', sd, ')')) %>%
  select(-m, -sd) %>%
  rename(full = name) %>%
  left_join(construct_terms, by = 'full')

#for flag specific
construct_terms <- df %>%
  group_by(name, flag) %>%
  summarise(m = mean(value), sd = sd(value)) %>%
  mutate(m = specify_decimal(m, 2)) %>%
  mutate(sd = specify_decimal(sd, 2))  %>%
  group_by(name, flag) %>%
  summarise(des = paste0(m, ' (', sd, ')')) %>%
  ungroup %>%
  rename(full = name) %>%
  pivot_wider(names_from = flag, values_from = des) %>%
  left_join(construct_terms, by = 'full')


cronbach_alpha <- function(var, items, key){
  #Takes tibble (dat) for scale (var), with item numbers (items), and reverse keyed items (key)
  #Returns standardised alpha
  #E.g., cronbach_alpha(dat.clean, "CI_", 1:6)), cronbach_alpha("BFI_10_", c(1,6), 1))
  item_names <- paste0(var, items)
  subscale <- dat.clean %>%
    select(!!item_names)
  if(!missing(key)){
    rev_item_names <- paste0(var, key)
    a <- subscale %>%
      psych::alpha(keys = rev_item_names)
  } else {
    a <- subscale %>%
      psych::alpha()
  }
  a$total$raw_alpha %>%
    return()
}
reliability_of_key <- function(k){
  #Requires key k (e.g., SJ)
  rel <- 0
  if(k == "BFI_10_A"){
    rel <- cronbach_alpha("BFI_10_", c(2,7), 7)
  }
  if(k == "BFI_10_C"){
    rel <- cronbach_alpha("BFI_10_", c(3,8), 3)
  }
  if(k == "BFI_10_E"){
    rel <- cronbach_alpha("BFI_10_", c(1,6), 1)
  }
  if(k == "BFI_10_N"){
    rel <- cronbach_alpha("BFI_10_", c(4,9), 4)
  }
  if(k == "BFI_10_O"){
    rel <- cronbach_alpha("BFI_10_", c(5,10), 10)
  }
  if(k == "CFC_S_F"){
    rel <- cronbach_alpha("CFC_S_", c(1,2,5,6))
  }
  if(k == "CFC_S_I"){
    rel <- cronbach_alpha("CFC_S_", c(3,4,7,8,9))
  }
  if(k == "EWS_E"){
    rel <- cronbach_alpha("EWS_", c(1,3,5,6,7,10))
  }
  if(k == "EWS_D"){
    rel <- cronbach_alpha("EWS_", c(2,4,8,9,11,12))
  }
  if(k == "CI"){
    rel <- cronbach_alpha("CI_", 1:6)
  }
  if(k == "NCS_6"){
    rel <- cronbach_alpha("NCS_6_", 1:6, c(3,4))
  }
  if(k == "SJ"){
    rel <- cronbach_alpha("SJ_", 1:8, c(3,7))
  }
  if(k == "SS_E"){
    rel <- cronbach_alpha("SS_", c(1,2,4,5,8,9,13,15))
  }
  if(k == "SS_R"){
    rel <- cronbach_alpha("SS_", c(3,6,7,10,11,12,14))
  }
  if(k == "MMS_cause_C"){
    rel <- cronbach_alpha("MMS_cause_", c(1,2,3,4,7,8,11))
  }
  if(k == "MMS_cause_E"){
    rel <- cronbach_alpha("MMS_cause_", c(5,6,9,10))
  }
  if(k == "MMS_cause_N"){
    rel <- cronbach_alpha("MMS_cause_", c(12,13))
  }
  if(k == "MMS_consequence_P"){
    rel <- cronbach_alpha("MMS_consequence_", c(2,6,11))
  }
  if(k == "MMS_consequence_S"){
    rel <- cronbach_alpha("MMS_consequence_", c(1,3,4,5,7,8,9,10))
  }
  if(k == "MMS_mitigation_C"){
    rel <- cronbach_alpha("MMS_mitigation_", c(1,7,8))
  }
  if(k == "MMS_mitigation_E"){
    rel <- cronbach_alpha("MMS_mitigation_", c(2,5,10))
  }
  if(k == "MMS_mitigation_G"){
    rel <- cronbach_alpha("MMS_mitigation_", c(3,4,6,9,11))
  }
  if(k == 'SVSS_ST'){
    dat_ST <- dat.clean %>%
      select(SVSS_1:SVSS_10) %>%
      mutate(SVSS_1 = -.19*SVSS_1) %>%
      mutate(SVSS_2 = -.14*SVSS_2) %>%
      mutate(SVSS_3 = -.09*SVSS_3) %>%
      mutate(SVSS_4 = -.11*SVSS_4) %>%
      mutate(SVSS_5 = .01*SVSS_5) %>%
      mutate(SVSS_6 = .10*SVSS_6) %>%
      mutate(SVSS_7 = .13*SVSS_7) %>%
      mutate(SVSS_8 = .07*SVSS_8) %>%
      mutate(SVSS_9 = .06*SVSS_9) %>%
      mutate(SVSS_10 = .02*SVSS_10) %>%
      mutate(ST_r = -.6+SVSS_1+SVSS_2+SVSS_3+SVSS_4+SVSS_5+SVSS_6+SVSS_7+SVSS_8+SVSS_9+SVSS_10)
    rel <- psych::alpha(dat_ST, delete=F, warnings=F)$total$raw_alpha
  }
  if(k == 'SVSS_C'){
    #Calculates alpha for C value scale
    dat_C <- dat.clean %>%
      select(SVSS_1:SVSS_10) %>%
      mutate(SVSS_1 = .05*SVSS_1) %>%
      mutate(SVSS_2 = .06*SVSS_2) %>%
      mutate(SVSS_3 = -.04*SVSS_3) %>%
      mutate(SVSS_4 = -.09*SVSS_4) %>%
      mutate(SVSS_5 = -.18*SVSS_5) %>%
      mutate(SVSS_6 = -.16*SVSS_6) %>%
      mutate(SVSS_7 = .03*SVSS_7) %>%
      mutate(SVSS_8 = .16*SVSS_8) %>%
      mutate(SVSS_9 = .18*SVSS_9) %>%
      mutate(SVSS_10 = .11*SVSS_10) %>%
      mutate(ST_r = .82+SVSS_1+SVSS_2+SVSS_3+SVSS_4+SVSS_5+SVSS_6+SVSS_7+SVSS_8+SVSS_9+SVSS_10)
    rel <- psych::alpha(dat_C, delete=F, warnings=F)$total$raw_alpha
  }
  return(rel)
}

construct_terms <- construct_terms %>%
  rowwise %>%
  mutate(alp = ifelse(num > 1, reliability_of_key(full), NA)) %>%
  mutate(alp = ifelse(is.na(alp), NA, specify_decimal(alp, 2))) %>%
  ungroup

#determine theoretical min-max of svss
#SVSS_ST
st_min <- -.6+-.19*8-.14*8-.09*8-.11*8
st_max <- -.6+.01*8+.10*8+.13*8+.07*8+.06*8+.02*8
st_range <- paste0(st_min,'-',st_max)

#SVSS_C
c_min <- .82-.04*8-.09*8-.18*8-.16*8
c_max <- .82+.05*8+.06*8+.03*8+.16*8+.18*8+.11*8
c_range <- paste0(c_min,'-',c_max)

construct_terms <- construct_terms %>%
  mutate(range = ifelse(is.na(range), ifelse(full == 'SVSS_C', c_range, st_range), range))

#make into nice table, inc. rearrange columents and dropping full
construct_terms <- construct_terms %>%
  mutate_all(.funs = list(~ ifelse(is.na(.), 'na', .)))

df <- construct_terms %>%
  select(-f1,-f2,-f3,-des) %>%
  select(abbrev:measure, num, alp, range, example, reference) 

dd <- as.dist((1-cormat)/2)
hc <- hclust(dd)
cormat <-cormat[hc$order, hc$order]
cormat[upper.tri(cormat)] <- NA #only one triangle is required

max_cor <- cormat %>%
  as.data.frame %>%
  rownames_to_column('var1') %>%
  gather(var2, value, -var1) %>%
  filter(var1 != var2) %>%
  filter(!is.na(value)) %>%
  mutate(abs_value = abs(value)) %>%
  arrange(desc(abs_value)) %>%
  slice(1)

cor_specific <- cormat %>%
  as.data.frame() %>%
  rownames_to_column('var1') %>%
  gather(var2, value, -var1) %>%
  mutate(var2 = fct_relevel(var2, colnames(cormat))) %>%
  mutate(var1 = fct_relevel(var1, colnames(cormat))) %>%
  filter(var1 == 'MMS_HUM' | var2 == 'MMS_HUM') %>%
  filter(!is.na(value)) %>%
  filter(var1 != var2) %>%
  filter(abs(value) > .30) %>%
  summarise(med = median(abs(value)), min = min(abs(value)), max = max(abs(value))) %>%
  mutate_all(.funs = list(~ specify_decimal(., 2)))

cor_mm <- cormat %>%
  as.data.frame() %>%
  rownames_to_column('var1') %>%
  gather(var2, value, -var1) %>%
  mutate(var2 = fct_relevel(var2, colnames(cormat))) %>%
  mutate(var1 = fct_relevel(var1, colnames(cormat))) %>%
  filter(grepl('MMS_', var1), grepl('MMS_', var2)) %>%
  filter(!is.na(value)) %>%
  filter(var1 != var2) %>%
  mutate(condition = ifelse(grepl('_MIT_', var1) & grepl('_MIT_', var2), 3, 
                            ifelse(grepl('_MIT_', var1) | grepl('_MIT_', var2),2,1)
  )) %>%
  group_by(condition) %>%
  summarise(med = median(abs(value)), min = min(abs(value)), max = max(abs(value))) %>%
  mutate_all(.funs = list(~ specify_decimal(., 2)))

cormat %>%
  as.data.frame() %>%
  rownames_to_column('var1') %>%
  gather(var2, value, -var1) %>%
  mutate(var2 = fct_relevel(var2, colnames(cormat))) %>%
  mutate(var1 = fct_relevel(var1, colnames(cormat))) %>%
  #mutate(value = ifelse(var2 == var1, NA, value)) %>%
  ggplot(aes(var2, var1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    mid = "#FBFEF9",
    low = "#0C6291",
    high = "#A63446" ,
    midpoint = 0,
    limit = c(-1, 1),
    space = "Lab",
    name = "Pearson\nCorrelation",
    na.value = 'white'
  ) +
  theme_classic() + # minimal theme
  theme(axis.text.x = element_text(
    angle = 45,
    vjust = 1,
    hjust = 0
  )) +
  scale_x_discrete(position = "top") +
  coord_fixed() +
  xlab("") +
  ylab("") +
  theme(
    legend.title = element_blank(),
    axis.text=element_text(family="sans", size=5, color = '#000000'),
    legend.text = element_text(family="sans", size=5, hjust = 1)
  )+
  guides(fill= guide_colorbar(barheight=4, barwidth = .5,
                              ticks.colour = "black",
                              frame.colour = "black"))
ggsave(sprintf('../out/figures/paper-figure-%s.png', index_fig), width=5, height=5, unit='in', dpi=300)
index_fig <- index_fig + 1 #increment counter

####Regression model figure (paper)
results_glm <- cvres_all_scaled
results_bootstrap <- glmCI_all_scaled
#results_classification <- 'mm'

tmp_coeffs <- coef(results_glm, s = "lambda.min")
tmp_names <- names(tmp_coeffs)
#retrieve results
glm_coefs <- do.call(rbind, tmp_coeffs) %>%
  summary %>%
  as.data.frame %>%
  as.tibble %>%
  mutate(pred = rep(rownames(tmp_coeffs[[1]]), times = length(tmp_names))) %>%
  mutate(outcome = rep(tmp_names, each = nrow(tmp_coeffs[[1]]))) %>%
  select(-i, -j) %>%
  rename(coef = x)
#now add CIs
if('pred' %in% colnames(results_bootstrap)){
  glm_cis <- results_bootstrap %>%
    as.data.frame() %>%
    gather(var, value, -pred) %>%
    mutate(outcome = gsub('^CI_.*_', '', var)) %>%
    mutate(ci = ifelse(grepl('.*\\.025.*', var), 'lower', 'upper')) %>%
    select(-var) %>%
    spread(ci, value, sep = '_') %>%
    as.tibble 
}else{
  glm_cis <- results_bootstrap %>%
    as.data.frame() %>%
    rownames_to_column('pred') %>%
    gather(var, value, -pred) %>%
    mutate(outcome = gsub('^CI_.*_', '', var)) %>%
    mutate(ci = ifelse(grepl('.*\\.025.*', var), 'lower', 'upper')) %>%
    select(-var) %>%
    spread(ci, value, sep = '_') %>%
    as.tibble 
}



#combine coefs and CIs for graphical/table display
glm_comb <- glm_cis %>%
  left_join(glm_coefs, by = c('pred', 'outcome'))
df_pred_order <- glm_comb %>%
  left_join(rename(construct_terms, pred = full), by = 'pred') %>%
  mutate(pred = abbrev) %>%
  mutate(intersects_zero = sign(ci_lower) != sign(ci_upper)) %>%
  group_by(pred) %>%
  summarise(m = max(abs(coef)), all_intersects_zero = all(intersects_zero)) %>%
  mutate(construct = gsub('_.*$', '', pred)) %>%
  group_by(m) %>%
  arrange(construct, .by_group = TRUE)

pred_order_levels <- pull(df_pred_order, pred)

pd <- position_dodge(width = 0.5)
segment_palette <- c('#AA3377', '#CCBB44', "#4477AA") #sceptic, opaque, acceptor

#num of constructs where all CIs intersect zero
zero_constructs <- glm_comb %>%
  group_by(pred) %>%
  summarise(all_intersect_zero = all(sign(ci_lower) != sign(ci_upper))) %>%
  filter(all_intersect_zero) %>%
  nrow()
p <- glm_comb %>%
  left_join(rename(construct_terms, pred = full), by = 'pred') %>%
  mutate(pred = abbrev) %>%
  mutate(pred = fct_relevel(pred, pred_order_levels)) %>%
  mutate(segment = ifelse(
    outcome == 'f1',
    'Acceptor',
    ifelse(outcome == 'f2', 'Sceptic', 'Opaque')
  )) %>%
  mutate(segment = fct_relevel(segment, c('Sceptic', 'Opaque', 'Acceptor'))) %>%
  ggplot(aes(pred, coef, color = segment)) +
  geom_point(stat = 'identity', position = pd) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), position = pd) +
  geom_hline(yintercept = 0, linetype = "dashed") 

#Red indicates where all coeffecients cross zero
df_annotes <- df_pred_order %>%
  ungroup %>%
  rownames_to_column('index') %>%
  mutate(index = as.integer(index)) %>%
  filter(all_intersects_zero) %>%
  mutate(xmax = index + 0.5) %>%
  mutate(xmin = index - 0.5) %>%
  mutate(xmin = ifelse(xmin < lag(xmax), lag(xmax), xmin)) %>%
  mutate(xmin = ifelse(index == 1, 0, xmin))

p <- p +
  annotate("rect", xmin=df_annotes$xmin, xmax=df_annotes$xmax, ymin=-Inf, ymax=Inf, alpha=0.12, fill="#303030")+
  coord_flip() +
  xlab('Predictors') +
  ylab('Regression Coefficient') +
  scale_color_manual(values = segment_palette) +
  theme_classic() +
  labs(color = "Segment")+
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    axis.text= element_text(family="sans", size=8, color = '#000000'),
    axis.title = element_text(family="sans", size=9, color = '#000000')
  )
p
ggsave(sprintf('../out/figures/paper-figure-%s.png', index_fig), width=5, height=5, unit='in', dpi=300)
index_fig <- index_fig + 1 #increment counter

####Table 2 of paper
tab_fsc <- q_results$zsc_n %>%
  mutate(order = 1:nrow(.)) %>%
  left_join(statement_text, by = 'order') %>%
  mutate(theme = c(1, 3, 5, 2, 5, 4, 3, 2, 5, 4, 3, 2, 2, 4, 3, 3, 2, 4, 5, 5, 3, 4, 1, 1, 1, 5, 1, 2, 4, 1))
#List statements from highest z score difference to least
tab_fsc_diff_order <- q_results$zsc_n %>% 
  mutate(order = 1:nrow(.)) %>%
  mutate(diff = abs(fsc_f1-fsc_f2)) %>%
  mutate(nonsig = (q_results$qdc$sig_f1_f2=="")) %>%
  mutate(max_rank = max(fsc_f1, fsc_f2)) %>%
  group_by(nonsig, -diff) %>%
  arrange(desc(max_rank), .by_group = TRUE) %>%
  ungroup %>%
  select(order, nonsig)
tab_fsc <- tab_fsc_diff_order %>%
  left_join(tab_fsc, by = 'order') %>%
  slice(1:26, 29, 27, 28, 30)
tab_fsc <- tab_fsc %>% 
  mutate(id = 1:nrow(.)) %>%
  rowwise %>%
  mutate(text = paste0(c(id, ".", " ", text), sep = "", collapse = "")) %>%
  ungroup
tab_fsc %>%
  select(text, fsc_f1, fsc_f2) %>%
  rename(Statement = text) %>%
  rename(Acceptor = fsc_f1) %>%
  rename(Sceptic = fsc_f2) %>%
  write_csv(sprintf('../out/tables/paper-table-%s.csv', index_tbl))
#Only one statement does not distinguish Acceptors and Sceptics,
q_results$qdc

index_tbl <- index_tbl + 1


####Table of survey scale means/sd
tbl_scales <- construct_terms %>%
  select(abbrev, des, f1, f3, f2)
colnames(tbl_scales) <- c('Scale', 'Overall', 'Acceptor', 'Fencesitter', 'Sceptic')
write_csv(tbl_scales, sprintf('../out/tables/paper-table-%s.csv', index_tbl))


##### Supplementary Figures
index_fig <- 1

results_glm <- cvres_mm_scaled
results_bootstrap <- glmCI_mm_scaled

tmp_coeffs <- coef(results_glm, s = "lambda.min")
tmp_names <- names(tmp_coeffs)
#retrieve results
glm_coefs <- do.call(rbind, tmp_coeffs) %>%
  summary %>%
  as.data.frame %>%
  as.tibble %>%
  mutate(pred = rep(rownames(tmp_coeffs[[1]]), times = length(tmp_names))) %>%
  mutate(outcome = rep(tmp_names, each = nrow(tmp_coeffs[[1]]))) %>%
  select(-i, -j) %>%
  rename(coef = x)
#now add CIs
if('pred' %in% colnames(results_bootstrap)){
  glm_cis <- results_bootstrap %>%
    as.data.frame() %>%
    gather(var, value, -pred) %>%
    mutate(outcome = gsub('^CI_.*_', '', var)) %>%
    mutate(ci = ifelse(grepl('.*\\.025.*', var), 'lower', 'upper')) %>%
    select(-var) %>%
    spread(ci, value, sep = '_') %>%
    as.tibble 
}else{
  glm_cis <- results_bootstrap %>%
    as.data.frame() %>%
    rownames_to_column('pred') %>%
    gather(var, value, -pred) %>%
    mutate(outcome = gsub('^CI_.*_', '', var)) %>%
    mutate(ci = ifelse(grepl('.*\\.025.*', var), 'lower', 'upper')) %>%
    select(-var) %>%
    spread(ci, value, sep = '_') %>%
    as.tibble 
}
#combine coefs and CIs for graphical/table display
glm_comb <- glm_cis %>%
  left_join(glm_coefs, by = c('pred', 'outcome'))

construct_terms <- tibble(
  full = sort(unique(glm_comb$pred)),
  abbrev = c(paste0('MMS_', c( 'CAU_C', 'CAU_E', 'CAU_N', 'CON_P', 'CON_S', 'HUM', 'MIT_C', 'MIT_E', 'MIT_G')))
)
#For all...
#Arrange graph such in order of differences in predictor coefs
df_pred_order <- glm_comb %>%
  left_join(rename(construct_terms, pred = full), by = 'pred') %>%
  mutate(pred = abbrev) %>%
  mutate(intersects_zero = sign(ci_lower) != sign(ci_upper)) %>%
  group_by(pred) %>%
  summarise(m = max(abs(coef)), all_intersects_zero = all(intersects_zero)) %>%
  mutate(construct = gsub('_.*$', '', pred)) %>%
  group_by(m) %>%
  arrange(construct, .by_group = TRUE)

pred_order_levels <- pull(df_pred_order, pred)

pd <- position_dodge(width = 0.5)
segment_palette <- c('#AA3377', '#CCBB44', "#4477AA") #sceptic, opaque, acceptor

#num of constructs where all CIs intersect zero
zero_constructs <- glm_comb %>%
  group_by(pred) %>%
  summarise(all_intersect_zero = all(sign(ci_lower) != sign(ci_upper))) %>%
  filter(all_intersect_zero) %>%
  nrow()
#create annotation rect
max_point <- glm_comb %>%
  pivot_longer(ci_lower:ci_upper, names_to = 'pt', values_to = 'val') %>%
  pull(val) %>%
  abs() %>%
  max %>%
  {.*10} %/% 2 %>%
  {(.+1)*2} %>%
  {./10}



p <- glm_comb %>%
  left_join(rename(construct_terms, pred = full), by = 'pred') %>%
  mutate(pred = abbrev) %>%
  mutate(pred = fct_relevel(pred, pred_order_levels)) %>%
  mutate(segment = ifelse(
    outcome == 'f1',
    'Acceptor',
    ifelse(outcome == 'f2', 'Sceptic', 'Opaque')
  )) %>%
  mutate(segment = fct_relevel(segment, c('Sceptic', 'Opaque', 'Acceptor'))) %>%
  ggplot(aes(pred, coef, color = segment)) +
  geom_point(stat = 'identity', position = pd) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), position = pd) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(limits = c(-max_point, max_point), breaks = seq(-max_point*10, max_point*10, 2 - (max_point*10 %% 2))/10)

#Red indicates where all coeffecients cross zero
df_annotes <- df_pred_order %>%
  ungroup %>%
  rownames_to_column('index') %>%
  mutate(index = as.integer(index)) %>%
  filter(all_intersects_zero) %>%
  mutate(xmax = index + 0.5) %>%
  mutate(xmin = index - 0.5) %>%
  mutate(xmin = ifelse(xmin < lag(xmax), lag(xmax), xmin)) %>%
  mutate(xmin = ifelse(index == 1, 0, xmin))

p <- p +
  annotate("rect", xmin=df_annotes$xmin, xmax=df_annotes$xmax, ymin=-Inf, ymax=Inf, alpha=0.12, fill="#303030")+
  coord_flip() +
  xlab('Predictors') +
  ylab('Regression Coefficient') +
  scale_color_manual(values = segment_palette) +
  theme_classic() +
  labs(color = "Segment")+
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    axis.text= element_text(family="sans", size=8, color = '#000000'),
    axis.title = element_text(family="sans", size=9, color = '#000000')
  )
p
ggsave(sprintf('../out/figures/SM-figure-%s.png', index_fig), width=5, height=5, unit='in', dpi=300)

###Scree plot of MMS causes scale
index_fig <- 6 #increment counter
apatheme=theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        text=element_text(family='sans', color='black'),
        legend.title=element_blank(),
        legend.position=c(.7,.8),
        axis.line.x = element_line(color='black'),
        axis.line.y = element_line(color='black'),
        axis.text = element_text(colour = 'black'))
eig <- dat.clean %>% 
  filter(!is.na(MMS_cause_1)) %>% 
  select(MMS_cause_1:MMS_cause_13) %>%
  as.matrix %>%
  cor %>%
  eigen(only.values = TRUE)
tibble(num = 1:length(eig$values), eigen = eig$values) %>%
  ggplot(aes(x = num, y = eigen)) +
  geom_line()+
  geom_point(size = 4)+
  scale_y_continuous(name='Eigenvalue')+
  scale_x_continuous(name='Factor number', breaks = 1:length(eig$values))+
  geom_vline(xintercept = 3.5, linetype = 'dashed')+
  apatheme+
  theme(axis.text = element_text(family="sans", size=14, color = '#000000'),
        axis.title.x = element_text(family="sans", size=16, color = '#000000'),
        axis.title.y = element_text(family="sans", size=16, color = '#000000'))
ggsave(sprintf('../out/figures/SM-figure-%s.png', index_fig), width=5, height=5, unit='in', dpi=300)

#Factor structure of Mental Model Scale (cause)
df_scale <- read_csv('../study2/study/surveyscales-mms-with-subscales.csv', col_names = TRUE, col_types = cols(
  item.original = col_character(),
  scale = col_character(),
  item.text = col_character(),
  order.new = col_double(),
  item.id = col_character(),
  subscale = col_character()
))
tab <- df_scale %>%
  filter(scale == 'cause') %>%
  select(item.text, order.new) %>%
  rename(order = order.new) %>%
  left_join(mutate(as.data.frame(loadings(mms_cause_factor_structure)[]), order = 1:n()), by = 'order')
#order according to loading
tab <- tab %>%
  slice(c(7,4,11,8,3,2,1,9,10,6,5,13,12))
#bold largest loadings
tab <- tab %>%
  pivot_longer(RC1:RC2, names_to = 'factor', values_to = 'loading') %>%
  mutate(loading = specify_decimal(loading, 2)) %>%
  group_by(item.text) %>%
  pivot_wider(names_from = 'factor', values_from = 'loading')
#descriptives
tab_des <- dat.clean %>%
  select(MMS_cause_1:MMS_cause_13) %>%
  ungroup %>%
  pivot_longer(MMS_cause_1:MMS_cause_13, names_to = 'item', values_to = 'score') %>%
  group_by(item) %>%
  summarise(m = mean(score), se = sd(score)/sqrt(n()))
tab <- tab %>%
  mutate(item = paste0('MMS_cause_', order)) %>%
  left_join(tab_des, by = 'item') %>%
  select(item.text,  m, se, RC1:RC2) %>%
  mutate_at(vars(m, se), .funs = list(~ specify_decimal(., 2)))
index_tbl <- 1
write_csv(tab, sprintf('../out/tables/SM-table-%s.csv', index_tbl))


}
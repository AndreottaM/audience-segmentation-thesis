#Load packages
library(qmethod) #for Q analysis
library(psych) #for CFAs, reliability analyses, correlations,e tc.
library(tidyverse) #for data wrangling
library('lme4')
library("parallel")
library("doParallel")

############
#Data wrangling
############

############
#Anonymise data
############
rawDir <- "../data-raw/study-3.csv"
dataDir <- gsub('-raw', '', rawDir)
idDir <- "../data-raw/study-3-id.csv"
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
  set.seed(50741) #set seed
  #Ignore first two rows (additional variable information)
  #Randomly assign a ID to each participant
  #completion_order assigns a ranking to each participant, indicating order of data collection (value 1 is 1st)
  dat.id <- dat.raw %>%
    slice(-1) %>%
    mutate(completion_order= 1:nrow(.)) %>%
    sample_frac(1, replace = F) %>%
    mutate(id = 1:nrow(.)) %>%
    select(QSEDResponseID, id, completion_order)
  #Map internal IDs with Qualtrics IDs, in case any participants require replacement
  #This will be only be visible in private repository
  dat.id %>%
    write.csv(file = idDir, row.names = F)
  #Join Id with ResponseId
  dat.raw <- dat.raw %>%
    left_join(dat.id, by = "QSEDResponseID")
  #Remove variables which identify Qualtrics accounts
  #Qualtrics adds participant data to files
  #If sensitive variables are named differently than expected
  #the following lines of codes will require changing AFTER PREREGISTRATION
  #to adequately deidentify participants
  #Remove sensitive data
  sensitive_names <- c("status", "X_recordId",	"ipAddress", "recordedDate",	"_recordId",	"recipientLastName",	"recipientFirstName",	"recipientEmail",	"externalDataReference",	"locationLatitude",	"locationLongitude",	"distributionChannel",	"userLanguage", "opp", "rid",	"RISN", "V", "QSEDResponseID")
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
set.seed(09940) #set seed for reproducable code
#Qualtrics csv output is (m)essy
dat.m <- read.csv(dataDir, stringsAsFactors = F)
#1. Remove redundant rows
dat.m <- dat.m[-1, ] %>%
  as.tibble()
#2. Embedded data and other variables are used by Qualtrics as background variables
# Not required for analysis
embedVars <- c("UniqueID", "code", "progress",	"finished", "gc", "Q_TotalDuration", "LS", "term")
ev <- embedVars[{embedVars %in% colnames(dat.m)}]

dat.m <- dat.m %>%
  select(-ev)

#3. Rename variables
#Create function to assist
rename_v <- function(tb, new, orig){
  #Take tibble tb, replaces part of a variable name (orig) with a new part (new), for all vars
  tb %>%
    rename_at(vars(contains(!!orig)), 
              funs(str_replace(., orig, new))) %>%
    return()
}
#Rename variables
dat.clean <- dat.m %>%
  rename(time_start = startDate) %>%
  rename(time_end = endDate) %>%
  rename(time_to_complete = duration) %>%
  rename_v("specs_", "QID99_") %>%
  rename_v("pif_", "QID103_") %>%
  rename(pcf = QID101) %>%
  rename_v("pcf_", "QID103") %>%
  rename(age = QID105_TEXT) %>%
  rename(gender = QID106) %>%
  rename_v("demographics_", "QID107_") %>%
  rename(code_entered = QID109_TEXT) %>%
  rename_v("qualtrics_qsort_", "QID110_") %>%
  rename_v("qualtrics_suvery_info_", "QID112_") %>%
  rename(qualtrics_qsort_id = CharID) %>%
  rename(trust_sci = QID49) %>%
  rename(trust_cat = QID119) %>%
  rename_v("trust_", "QID203_") %>%
  rename_v("cause_emissions_pre_", "QID271_") %>%
  rename_v("cause_emissions_pre_", "QID278_") %>%
  rename_v("cause_emissions_post_", "QID272_") %>%
  rename_v("cause_emissions_post_", "QID279_") %>%
  rename_v("cause_activity_pre_", "QID276_") %>%
  rename_v("cause_activity_post_", "QID277_") %>%
  rename_v("cause_activity_post_", "QID283_") %>%
  rename_v("conseq_pre_", "QID290_") %>%
  rename_v("conseq_pre_", "QID287_") %>%
  rename_v("conseq_post_", "QID286_") %>%
  rename_v("conseq_post_", "QID288_") %>%
  rename_v("cause_gen_pre_", "QID274_") %>%
  rename_v("cause_gen_pre_", "QID280_") %>%
  rename_v("cause_gen_post_", "QID275_") %>%
  rename_v("cause_gen_post_", "QID281_") %>%
  rename_v("mitig_approve_pre_", "QID260_") %>%
  rename("mitig_success_pre" = QID190_TEXT) %>%
  rename("mitig_direction_pre" = QID261.1_1) %>%
  rename("mitig_amount_pre" = QID261.2_1_1) %>%
  rename_v("mitig_pre", "QID193") %>%
  rename_v("mitig_approve_post_", "QID294_") %>%
  rename("mitig_success_post" = QID295_TEXT) %>%
  rename("mitig_direction_1_post" = QID259.1_1) %>%
  rename("mitig_amount_1_post" = QID259.2_1_1) %>%
  rename("mitig_direction_2_post" = QID300.1_1) %>%
  rename("mitig_amount_2_post" = QID300.2_1_1) %>%
  rename_v("mitig_post", "QID213") %>%
  rename_v("conseq_sentiment_", "QID313_") %>%
  rename_v("conseq_sentiment_", "QID317_") %>%
  rename("mitig_sentiment" = QID320) %>%
  rename_v("mitig_sentiment_", "QID315_") %>%
  rename(ls = int)
#4. Import Q sort data
dat.m_qsort <- read.csv("../data-raw/study-3-qsort.csv", stringsAsFactors = F, header = F)
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
cleanDir <- "../out/study-3-clean.csv"
dat.clean %>%
  write_csv(cleanDir)

####
#Recoding
####
dat.clean <- dat.clean %>% 
  mutate_at(vars(trust_sci:mitig_sentiment, -ends_with("_CLICK"),  -ends_with("_COUNT"),  -ends_with("_SUBMIT")), as.numeric) %>%
  mutate_at(vars(age,time_to_complete), as.numeric)
#Remove the _1
dat.clean <- dat.clean %>%
  rename_at(.vars = vars(starts_with("cause_"), starts_with("conseq_"), -contains("_sentiment_")),
            .funs = funs(sub("[_]1$", "", .)))
dat.clean <- dat.clean %>%
  rename(conseq_pre_2 = conseq_pre_17) %>%
  rename(conseq_post_2 = conseq_post_17) %>%
  rename(conseq_pre_3 = conseq_pre_18) %>%
  rename(conseq_post_3 = conseq_post_18) %>%
  rename(conseq_pre_4 = conseq_pre_39) %>%
  rename(conseq_post_4 = conseq_post_39) %>%
  rename(conseq_pre_5 = conseq_pre_56) %>%
  rename(conseq_post_5 = conseq_post_56) %>%
  rename(conseq_pre_6 = conseq_pre_40) %>%
  rename(conseq_post_6 = conseq_post_40) %>%
  rename(conseq_pre_7 = conseq_pre_41) %>%
  rename(conseq_post_7 = conseq_post_41) %>%
  rename(conseq_pre_8 = conseq_pre_42) %>%
  rename(conseq_post_8 = conseq_post_42) %>%
  rename(conseq_pre_9 = conseq_pre_43) %>%
  rename(conseq_post_9 = conseq_post_43) %>%
  rename(conseq_sentiment_2 = conseq_sentiment_14) %>%
  rename(conseq_sentiment_3 = conseq_sentiment_15) %>%
  rename(conseq_sentiment_4 = conseq_sentiment_16) %>%
  rename(conseq_sentiment_5 = conseq_sentiment_17) %>%
  rename(conseq_sentiment_6 = conseq_sentiment_18) %>%
  rename(conseq_sentiment_7 = conseq_sentiment_19) %>%
  rename(conseq_sentiment_8 = conseq_sentiment_20) %>%
  rename(conseq_sentiment_9 = conseq_sentiment_21)

dat.clean <- dat.clean %>%
  mutate(mitig_direction_post = ifelse(is.na(mitig_direction_1_post), mitig_direction_2_post, mitig_direction_1_post)) %>%
  mutate(mitig_amount_post = ifelse(is.na(mitig_amount_1_post), mitig_amount_2_post, mitig_amount_1_post))

dat.clean <- dat.clean %>%
  mutate(mitig_change_pre = ifelse(mitig_direction_pre == 1, 100+mitig_amount_pre,
                                      ifelse(mitig_direction_pre == 3, 100-mitig_amount_pre, 100))) %>%
  mutate(mitig_change_post = ifelse(mitig_direction_post == 1, 100+mitig_amount_post,
                                  ifelse(mitig_direction_post == 3, 100-mitig_amount_post, 100))) 
#Recode close-ended scale questions
#Approval questions on scale of -3 to +3 (0 is neutral)
#sentiment questions on scale -2 to +2 (0 is neutral)
dat.clean <- dat.clean %>%
  mutate_at(vars(starts_with('mitig_approve')), ~ .- 7) %>%
  mutate_at(vars(contains('_sentiment'), -ends_with("_CLICK"),  -ends_with("_COUNT"),  -ends_with("_SUBMIT")), ~ .- 3)



#######
#Exclusion criteria
#######

#1. Exclude fast responders
dat.clean <- dat.clean %>%
  rowwise() %>%
  mutate(fast_responder = as.integer(time_to_complete) < 664)
#Throw out data from all irrelevant participants
dat.clean <- dat.clean %>%
  filter(fast_responder == F | is.na(fast_responder)) %>%
  filter(id > 0)
#2. Check number of responses for each part of study
n_total <- dat.clean %>%
  nrow()

####
#Q sort analysis
####

dat.clean <- dat.clean %>%
  mutate_at(vars(starts_with('sort_sta_')), function(x){return(as.integer(substr(x, 2, 2)) - 5)}) %>%
  mutate_at(vars(starts_with('sort_sta_')), as.integer)

#Contains information on statements, where each row is a seperate statement, columns are individual participants and statement text
statement_information <- dat.clean %>% 
  filter(!is.na(sort_sta_1)) %>% 
  select(id,starts_with('sort_sta_')) %>%
  gather(key = statement, rank, -id) %>% 
  mutate(id = paste0("id_", id)) %>% 
  spread(id, rank)

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

#Eigenvalues:
q_sorts %>% 
  cor() %>% 
  eigen() %>% 
  .$values %>% 
  .[1:10] %>% 
  plot()

q_factor_num <- 1

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
cribbleDir <- '../out/study-3-cribsheet.csv'
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

dat.clean <- dat.clean %>% 
  mutate(flag = ifelse(is.na(flag), 3, flag))  %>%
  mutate(flag = factor(flag, levels = 1:3, labels = c("f1", "f2", "f3")))

# model_demo <- lmer(update ~ 1 + segment + (1|id), data = d)
# summary(model_demo)

####Check stratified sampling procedure:####
df_strats <- tibble(
  bracket_name = c("18-24", "25-34", "35-44", "45-54", "55-64", "65+", "male", "female"),
  bracket_code = c(1:6,1,2),
  strat = c(rep("age",6), rep("gender", 2)),
  percent = c(12.3,19.3,17.2,16.8,14.8,19.6, 49.3, 50.7) #per Australian 2016 census
)

find_extreme_proportions <- function(n, z, t, low){
  #z is z-score for confidence interval (1.96 corresponds to 95% CIs)
  #n is sample size
  #t is true proportion within population
  #if low is true, return min proprtion (otherwise max)
  const <- z^2/n
  a <- 1+const
  b <- -(const + 2*t)
  c <- t^2 
  dir <- ifelse(low, -1, 1)
  p <- (-b + dir*(sqrt(b^2 - 4*a*c)))/(2*a)
  return(p)
}


df_strats <- df_strats %>%
  rowwise %>%
  mutate(min = find_extreme_proportions(n_total, 1.96, percent/100, T)*n_total) %>%
  mutate(max = find_extreme_proportions(n_total, 1.96, percent/100, F)*n_total) %>%
  mutate(min = ceiling(min)) %>%
  mutate(max = floor(max))

df_strats <- dat.clean %>% 
  select(age) %>% 
  mutate(age = ifelse(age > 1900, 2020-age, age)) %>%
  mutate(bracket_code = findInterval(age, c(18, 25, 35, 45, 55, 65))) %>%
  count(bracket_code) %>%
  add_column(strat = 'age') %>%
  right_join(df_strats, by = c('bracket_code', 'strat'))

df_strats <- dat.clean %>% 
  select(gender) %>%
  mutate(bracket_code = as.integer(gender)) %>% 
  count(bracket_code) %>%
  add_column(strat = 'gender') %>%
  right_join(df_strats, by = c('bracket_code', 'strat'))

df_strats <- df_strats %>%
  mutate(n = ifelse(is.na(n.x), n.y, n.x)) %>%
  select(-n.x, -n.y) %>%
  mutate(quota_met = n >= min & n <= max)
df_strats

#Extract stimuli data
df_but <- read.csv('../study3/stimuli/stimuli.csv', stringsAsFactors = F) %>%
  as.tibble %>% 
  rename(statement = Statement) %>% 
  rename(block = Block) %>% 
  select(block, order, statement, p)
#Turn block into factor
df_but <- df_but %>% 
  mutate(block = as.factor(block)) %>% 
  mutate(block = as.integer(block))
#Add row for mitigation question
df_but <- df_but %>%
  add_row(block = 5, order = 1, statement = "Under current policies in place (including the Emissions Reduction Fund), what will be the change in Australia's carbon dioxide emissions by the year 2030 (compared to 2005 levels)?", p = "Increase by 8%")
#Clarify p
df_but <- df_but %>%
  mutate(p_text = p) %>%
  mutate(p = ifelse(p_text == "Near 100", "100",
                    ifelse(p_text == "Near 0", "0",
                           ifelse(p_text == "Increase by 8%", "108", p_text)))) %>%
  mutate(p = as.integer(p))

######
#Q sort: do profiles replicate?
#####
#From Study 2
dir_results_study2 <- '../out/study-2-factor-scores.csv'

if(file.exists(dir_results_study2)){
df_qres2 <- read_csv(
  dir_results_study2,
  col_types = cols(
    sta_id = col_character(),
    zsc_f1 = col_double(),
    zsc_f2 = col_double(),
    fsc_f1 = col_integer(),
    fsc_f2 = col_integer()
  )
)

#From Study 3
zsc_s3 <- q_results$zsc %>%
  rownames_to_column('sta_id') %>%
  as.tibble %>%
  gather(zsc_f, z_s3, -sta_id)

fsc_s3 <- q_results$zsc_n %>%
  rownames_to_column('sta_id') %>%
  as.tibble %>%
  gather(fsc_f, f_s3, -sta_id)

#Factor z scores
df_qcompar_zsc <- df_qres2 %>%
  select(-fsc_f1, -fsc_f2) %>%
  gather(zsc_f, z_s2, -sta_id) %>%
  right_join(zsc_s3, by = c('sta_id', 'zsc_f')) %>%
  group_by(zsc_f) %>%
  mutate(f = strsplit(zsc_f, 'zsc_f', fixed=T)[[1]][2]) %>%
  ungroup %>%
  mutate(f = as.integer(f)) %>%
  select(-zsc_f)

df_qcompar_fsc <- df_qres2 %>%
  select(-zsc_f1, -zsc_f2) %>%
  gather(fsc_f, f_s2, -sta_id) %>%
  right_join(fsc_s3, by = c('sta_id', 'fsc_f')) %>%
  group_by(fsc_f) %>%
  mutate(f = strsplit(fsc_f, 'fsc_f', fixed=T)[[1]][2]) %>%
  ungroup %>%
  mutate(f = as.integer(f)) %>%
  select(-fsc_f)

#Combine for easy management
df_qcompar <- df_qcompar_zsc %>%
  right_join(df_qcompar_fsc, by = c('sta_id', 'f'))
}



  
####
#Belief updating
####

#Relative update is dependent variable of interest:

df_update <- dat.clean %>% 
  select(id,
         cause_emissions_pre_1:cause_emissions_pre_6, cause_emissions_post_1:cause_emissions_post_6,
         cause_activity_pre_1:cause_activity_pre_6, cause_activity_post_1:cause_activity_post_6,
         cause_gen_pre_1:cause_gen_pre_2, cause_gen_post_1:cause_gen_post_2,
         conseq_pre_1:conseq_pre_9, conseq_post_1:conseq_post_9,
         mitig_change_pre, mitig_change_post) %>%
  gather(key, belief, -id)

#Create session variable (i.e., 1 = pre or 2 =  post)
#Create block variable (one of 5)
#Create item variable (which event/driver is being updated)
df_update <- df_update %>%
  group_by(key) %>%
  mutate(session = ifelse(length(grep('_pre', key)) > 0, 1, 2)) %>%
  mutate(block = ifelse(session == 1, sub('_pre.*', '', key), sub('_post.*', '', key))) %>%
  mutate(order = ifelse(length(grep('mitig_change_.*', key)) > 0, "1",
                       ifelse(session == 1, sub('^.*_pre_', '', key), sub('^.*_post_', '', key)))) %>%
  ungroup %>%
  mutate(order = as.integer(order))

#Add true probabilities
df_update <- df_but %>%
  group_by(block) %>%
  mutate(block_num = block) %>%
  group_by(block_num) %>%
  mutate(block = as.character(block)) %>%
  mutate(block = c("cause_activity", "cause_gen", "cause_emissions", "conseq", "mitig_change")[block_num]) %>%
  ungroup %>%
  select(block, order, p) %>%
  left_join(df_update, by = c('block', 'order'))

#Spread based on session
#Calculate update towards p
df_update <- df_update %>%
  select(-key) %>%
  spread(session, belief, sep = "_") %>%
  mutate(update = ifelse(session_1 < p, session_2 - session_1,
                    ifelse(session_1 > p, session_1 - session_2, 
                           ifelse(session_1 != session_2, -abs(session_2 - session_1), NA))))

#Add in Q sort flags
df_update <- dat.clean %>%
  select(id, flag) %>%
  left_join(df_update, by = 'id') %>%
  ungroup



#Extract trust and estimation error information
df_update <- dat.clean %>%
  select(id, trust_sci, trust_cat) %>%
  right_join(df_update, by = c('id')) %>%
  ungroup %>%
  mutate(trust = ifelse(block == 'mitig_change', trust_cat, trust_sci)) %>%
  mutate(estimation_error = abs(session_1 - p)) %>%
  ungroup

#Add sentiment to update information
df_update <- dat.clean %>%
  select(id, contains('_sentiment'), -ends_with("_CLICK"),  -ends_with("_COUNT"),  -ends_with("_SUBMIT")) %>%
  gather(key, sentiment, -id) %>%
  mutate(block = ifelse(key == 'mitig_sentiment', 'mitig_change', 'conseq')) %>%
  mutate(order = ifelse(block == 'mitig_change', 1, sub('^.*_sentiment_', '', key))) %>%
  ungroup %>%
  mutate(order = as.integer(order)) %>%
  select(-key) %>%
  right_join(df_update, by = c('id', 'block', 'order'))

#Add ordering information to update tibble
df_update <- dat.clean %>%
  select(id, ls) %>%
  right_join(df_update, by = 'id')

#Determine whether news is good/bad for individuals.
#good_news is TRUE when participant has recieved good news, FALSE when bad news, and NA when neutral news or no sentiment information
#calculation based on consequences, then reversed for mitigation (as relationship with p and sentiment is reversed)
df_update <- df_update %>%
  rowwise %>%
  mutate(good_news = ifelse(sentiment == 0 | is.na(sentiment) | estimation_error == 0, NA, 
                            ifelse(sentiment < 0 & p < session_1, TRUE , #event is negative, and P has discovered its probability is lower than estimated
                                   ifelse(sentiment > 0 & p > session_1, TRUE, FALSE)))) %>% #event is positive, and P has discovered probability is higher than expected
  mutate(good_news = ifelse(block == 'mitig_change', !good_news, good_news))

#Collapse all cause blocks into a 'cause' domain
df_update <- df_update %>%
  rowwise %>%
  mutate(domain = sub('_.*$', '', block)) %>%
  ungroup %>% 
  mutate(domain = as.factor(domain))
#Create another DV, relative update
df_update <- df_update %>%
  mutate(update_relative = ifelse(estimation_error != 0, update/estimation_error, NA))

#Remove participant 292. Answered mitigation question strangely
#Large answer on single-item domain prevents models from converging
df_update <- df_update %>%
  filter(id != 292)

#Add item ID as grouping for random effect
df_update <- df_update %>% 
  mutate(item_id = paste0(block, '_', order)) %>% 
  mutate(item_id = as.factor(item_id))

#Show relevant means
df_update %>%
  group_by(domain, flag) %>%
  summarise_at(vars(update, update_relative), mean, na.rm = T)

####
#Model 1
#DV: (Unscaled) update per domain
#IV: Flag
###

df_update_r <- df_update %>%
  mutate(flag = relevel(flag, 'f2'))

#d1 = cause
ms_d1_v1 <- lmer(update_relative ~ 1 + (1|id) + (1|item_id) + (1|ls), filter(df_update, domain == 'cause'), REML = F)
ms_d1_v2 <- lmer(update_relative ~ 1 + flag + (1|id) + (1|item_id) + (1|ls), filter(df_update, domain == 'cause'), REML = F)
#lowest AIC is v2, see below:
AIC(ms_d1_v1, ms_d1_v2) %>% 
  .[ , 2] %>% 
  which.min
#run contrasts on best model:
ms_d1_v2_c <- lmer(update_relative ~ 1 + flag + (1|id) + (1|item_id) + (1|ls), filter(df_update, domain == 'cause'), REML = F, contrasts = list(flag = contr.sum(3)))
ms_d1_v2_r <- lmer(update_relative ~ 1 + flag + (1|id) + (1|item_id) + (1|ls), filter(df_update_r, domain == 'cause'), REML = F)



#d2 = conseq
ms_d2_v1 <- lmer(update_relative ~ 1 + (1|id) + (1|item_id) + (1|ls), filter(df_update, domain == 'conseq'), REML = F)
ms_d2_v2 <- lmer(update_relative ~ 1 + flag + (1|id) + (1|item_id) + (1|ls), filter(df_update, domain == 'conseq'), REML = F)
#lowest AIC is v2, see below:
AIC(ms_d2_v1, ms_d2_v2) %>% 
  .[ , 2] %>% 
  which.min
#run contrasts on best model:
ms_d2_v2_c <- lmer(update_relative ~ 1 + flag + (1|id) + (1|item_id) + (1|ls), filter(df_update, domain == 'conseq'), REML = F, contrasts = list(flag = contr.sum(3)))
ms_d2_v2_r <- lmer(update_relative ~ 1 + flag + (1|id) + (1|item_id) + (1|ls), filter(df_update_r, domain == 'conseq'), REML = F)

#d3 = mitig
ms_d3_v1 <- lmer(update_relative ~ 1 + (1|ls), filter(df_update, domain == 'mitig'), REML = FALSE)
ms_d3_v2 <- lmer(update_relative ~ 1 + flag + (1|ls), filter(df_update, domain == 'mitig'), REML = FALSE)
#lowest AIC is v2, see below:
AIC(ms_d3_v1, ms_d3_v2) %>% 
  .[ , 2] %>% 
  which.min
ms_d3_v2_c <- lmer(update_relative ~ 1 + flag + (1|ls), filter(df_update, domain == 'mitig'), REML = F, contrasts = list(flag = contr.sum(3)))
ms_d3_v2_r <- lmer(update_relative ~ 1 + flag + (1|ls), filter(df_update_r, domain == 'mitig'), REML = F)


######
#Trust*Flag per domain
######
#Center trust values, such that trust_c == 0 is equivalent to average trust
df_update <- df_update %>%
  group_by(domain, ls) %>%
  mutate(trust_c = trust - mean(trust))

#Specify the fixed effects of models of interest
trust_models <- tibble(
  model_id = 1:5,
  fixed_effects = c(
    "",
    "flag",
    "trust_c",
    "flag + trust_c",
    "flag + flag:trust_c"
  )
)
trust_models

#Seperate models for each domain
data_filter <- c('cause', 'conseq', 'mitig')

for(i in 1:length(data_filter)){
  var_AIC <- paste0(data_filter[i], '_AIC')
  trust_models <- trust_models %>%
    add_column(AIC = double(5)) #renamed to var_AIC following completion of loop
  for(j in 1:5){
    fixed_effects <- trust_models %>%
      filter(model_id == j) %>%
      select(fixed_effects) %>%
      slice(1)
    if(data_filter[i] != 'mitig'){
      random_effects <- paste('(1|id)', '(1|item_id)', '(1|ls)', sep = ' + ')
    }else{
      random_effects <- '(1|ls)'
    }
    f <- ifelse(nchar(fixed_effects) == 0, "update_relative~1",
                  paste("update_relative", '~', paste('1', fixed_effects, sep = ' + ')))
    m <- lm(f, filter(df_update, domain == data_filter[i])) #lmer model
    trust_models <- tibble(model_id = j, AIC = AIC(m)) %>%
      right_join(trust_models, by = c('model_id')) %>%
      mutate(AIC = ifelse(is.na(AIC.x), AIC.y, AIC.x)) %>%
      select(-AIC.x, -AIC.y)
  }
  #Rename AIC
  trust_models <- trust_models %>%
    mutate(!!var_AIC := AIC) %>%
    select(-AIC)
}
#Best models:
trust_models %>%
  gather(domain, AIC, -model_id, -fixed_effects) %>%
  group_by(domain) %>%
  top_n(-1, AIC)
#best models:
ms_d1_trust <- lmer(update_relative ~ 1 + flag+flag:trust_c + (1|id) + (1|item_id) + (1|ls), filter(df_update, domain == 'cause'), REML = FALSE)
ms_d2_trust <- lmer(update_relative ~ 1 + flag+trust_c + (1|id) + (1|item_id) + (1|ls), filter(df_update, domain == 'conseq'), REML = FALSE)
#interaction models:
ms_d1_trust_int <- lmer(update_relative ~ 1 + flag+flag:trust_c + (1|id) + (1|item_id) + (1|ls), filter(df_update, domain == 'cause'), REML = FALSE)
ms_d2_trust_int <- lmer(update_relative ~ 1 + flag+flag:trust_c + (1|id) + (1|item_id) + (1|ls), filter(df_update, domain == 'conseq'), REML = FALSE)
ms_d3_trust_int <- lmer(update_relative ~ 1 + flag+flag:trust_c +  (1|ls), filter(df_update, domain == 'mitig'), REML = FALSE)



#Good news/Bad news:
news_models <- trust_models %>%
  select(model_id, fixed_effects) %>%
  mutate(fixed_effects = gsub('trust_c', 'good_news', fixed_effects))
news_models


#Ignore cause domain
for(i in 2:length(data_filter)){
  var_AIC <- paste0(data_filter[i], '_AIC')
  news_models <- news_models %>%
    add_column(AIC = double(5)) #renamed to var_AIC following completion of loop
  for(j in 1:5){
    fixed_effects <- news_models %>%
      filter(model_id == j) %>%
      select(fixed_effects) %>%
      slice(1)
    if(data_filter[i] != 'mitig'){
      random_effects <- paste('(1|id)', '(1|item_id)', '(1|ls)', sep = ' + ')
    }else{
      random_effects <- '(1|ls)'
    }
      f <- paste("update_relative", '~', paste('1', fixed_effects, random_effects, sep = ' + ')) #lmer formula
      m <- lmer(f, filter(df_update, domain == data_filter[i] & !is.na(good_news)), REML = FALSE) #lmer model, only include valid good_news rows to ensure each model is modelling the same ID
    news_models <- tibble(model_id = j, AIC = AIC(m)) %>%
      right_join(news_models, by = c('model_id')) %>%
      mutate(AIC = ifelse(is.na(AIC.x), AIC.y, AIC.x)) %>%
      select(-AIC.x, -AIC.y)
  }
  #Rename AIC
  news_models <- news_models %>%
    mutate(!!var_AIC := AIC) %>%
    select(-AIC)
}
#Best models:
news_models %>%
  gather(domain, AIC, -model_id, -fixed_effects) %>%
  group_by(domain) %>%
  top_n(-1, AIC)




######Does update predict approval change in mitigation?
df_change <- dat.clean %>% 
  mutate(change_1 = mitig_approve_post_1-mitig_approve_pre_1) %>%
  mutate(change_2 = mitig_approve_post_2-mitig_approve_pre_2) %>%
  mutate(change_success = mitig_success_post-mitig_success_pre) %>% 
  select(id, change_1, change_2, change_success) %>%
  filter(id != 292) %>%
  left_join(df_update, by = 'id') %>%
  filter(domain == 'mitig') %>%
  mutate(belief_change = session_1 - session_2)#+ve = policy is better than initially thought, -ve = policy is worse than initially thought

#center changes on ls
df_change <- df_change %>%
  group_by(ls) %>%
  mutate(change_success = change_success - mean(change_success)) %>%
  mutate(belief_change = belief_change - mean(belief_change))

#For both change_1 & change_2
#Three IVs of interest: flag, belief_change, and change_success
#Mark these by A, B, and C
iv <- c('flag','belief_change','change_success')
#main effects
fixed_effects <- c(iv,
                   apply(combn(iv, 2), 2, paste0, sep = "", collapse = "+"),
                   paste0(iv, collapse = "+", sep = ""))
#interaction effects
interactions <- lapply(fixed_effects, function(x){
  if(grepl('[+]', x)){
    #then there are at least two main effects
    #return all interaction models
    main <- strsplit(x, '[+]')[[1]] #main effects
    main_missing <- iv[!(iv %in% main)] #missing main effects
    #All twoway interacts
    m_int <- apply(combn(main,2), 2, function(y){c(y[1], paste0(y, sep = "", collapse = ":"))}) %>%
       apply(2, paste0, sep = "", collapse = "+")
    if(length(main) == 3){
      #Add a three way interaction
      m_int <- c(m_int, paste0(c(main[1], paste0(main, sep = "", collapse = ":")), sep = "", collapse = "+"))
    }
    if(length(main_missing) > 0){
      #For models missing a main effect, create a model with the main effect
      #For three main effects, only one effect can be missing for interaction models
      m_int <- paste0(c(m_int, main_missing), sep = "", collapse = "+") %>%
        c(m_int)
    }
    return(m_int)
  }
}) %>%
  unlist
#All fixed effects:
fixed_effects <- c("", fixed_effects, interactions) %>%
  unique
random_effects <- "(1|ls)"
#Create tibble to house models
policy_models <- tibble(
  model_id = rep(1:length(fixed_effects), times = 2),
  fixed_effects = rep(fixed_effects, times = 2),
  dv = rep(c('change_1', 'change_2'), each = length(fixed_effects)/2),
  random_effects = random_effects
)
#Create model equations
policy_models <- policy_models %>%
  rowwise %>%
  mutate(intercept_term = ifelse(fixed_effects == "", "1", "1+")) %>%
  mutate(lmer_eq = paste0(c(dv, '~', intercept_term, fixed_effects, '+', random_effects), sep = "", collapse = ""))
#Record AICs
policy_models <- policy_models %>%
  mutate(AIC = AIC(lmer(lmer_eq, df_change, REML = FALSE))) %>%
  ungroup %>%
  group_by(dv) %>%
  mutate(AIC_diff = AIC - min(AIC))
#Three models with best fit
policy_models %>%
  filter(AIC_diff < 2 + 1e-05)
change_v <- policy_models %>%
  filter(AIC_diff < 4 + 1e-05) %>%
  pull(lmer_eq)
m_change_v1 <- lmer(change_v[1], df_change, REML = FALSE)
m_change_v2 <- lmer(change_v[2], df_change, REML = FALSE)
m_change_v3 <- lmer(change_v[3], df_change, REML = FALSE)
m_change_v4 <- lmer(change_v[4], df_change, REML = FALSE)






##########################
####Tables and Figures####
##########################

specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k)) #for rounding

#Figure 5 from Paper

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

segment_palette <- c('#AA3377', '#CCBB44', "#4477AA") #sceptic, opaque, acceptor
flaglabels <- c("Acceptor", "Opaque", "Sceptic")
segment_palette <- rev(segment_palette)
pd <- position_dodge(.9) 

df_update %>% 
  group_by(id, domain, flag) %>%
  summarise(mean_per_n = mean(update_relative, na.rm = TRUE)) %>%
  group_by(domain, flag) %>%
  summarise(mean = mean(mean_per_n, na.rm = TRUE), sd = sd(mean_per_n, na.rm = TRUE), n = n()) %>%
  mutate(se = sd/sqrt(n)) %>%
  mutate(flag = fct_relevel(flag, "f2", after = 2)) %>%
  ggplot(aes(domain, mean, fill=flag)) +
  geom_bar(position=pd, stat="identity") +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = .3, position = pd)+
  ylab("Mean relative belief update")+
  xlab("Domain")+
  ylim(-0.5, 1)+
  apatheme +
  theme(legend.position = "none")+
  scale_x_discrete(labels= c('Cause', 'Consequence', 'Mitigation'))+
  scale_fill_manual(values = segment_palette)
ggsave(sprintf('../out/figures/paper-figure-%s.png', 5), width=5, height=5, unit='in', dpi=300)


#Figure 6 from Paper

library(ggpubr)

dat_news <- df_update %>% 
  group_by(id, domain, flag, good_news) %>%
  summarise(mean_per_n = mean(update_relative, na.rm = TRUE)) %>%
  group_by(domain, flag, good_news) %>%
  summarise(mean = mean(mean_per_n, na.rm = TRUE), sd = sd(mean_per_n, na.rm = TRUE), n = n()) %>%
  mutate(se = sd/sqrt(n)) %>%
  filter(!is.na(good_news)) %>%
  ungroup %>%
  mutate(flag = fct_relevel(flag, "f2", after = 2))
p <- lapply(c('conseq', 'mitig'), function(x){
  #for each domain,
  dat_news %>%
    filter(domain == x) %>%
    ggplot(aes(good_news, mean, fill=flag)) +
    geom_bar(position=pd, stat="identity") +
    geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = .3, position = pd)+
    ylab("Mean relative belief update")+
    xlab("Type of News")+
    ylim(-0.5, 1)+
    apatheme +
    theme(legend.position = "none")+
    scale_x_discrete(labels= c('Bad', 'Good'))+
    scale_fill_manual(values = segment_palette)
})

ggarrange(p[[1]], p[[2]], 
          labels = c("(a)", "(b)"),
          ncol = 2, nrow = 1, hjust = 0)
ggsave(sprintf('../out/figures/paper-figure-%s.png', 6), width=5, height=5, unit='in', dpi=300)


#Figure 7 from Supplementary Material

factor_range <- 1:15
tibble(num = factor_range, eigen = head(q_fa$values, tail(factor_range, 1))) %>%
  ggplot(aes(x = num, y = eigen)) +
  geom_line()+
  geom_point(size = 4)+
  scale_y_continuous(name='Eigenvalue', limits = c(0, 150))+
  scale_x_continuous(name='Factor number', breaks = factor_range)+
  geom_vline(xintercept = 1.5, linetype = 'dashed')+
  apatheme
ggsave(sprintf('../out/figures/SM-figure-%s.png', 7), width=5, height=5, unit='in', dpi=300)


#Table 6 from Supplementary Material
tab_ind <- 6
tab_fsc <- df_qcompar %>%
  select(-z_s2, -z_s3) %>%
  pivot_longer(f_s2:f_s3, names_to = 'study', values_to = 'f_score') %>%
  mutate(f = ifelse(f == 1, 'Acceptor', 'Sceptic')) %>%
  mutate(study = gsub('f_s', '', study)) %>%
  mutate(var = paste0('Study_', study, '_', f)) %>%
  mutate(order = as.integer(gsub('sta_', '', sta_id))) %>%
  left_join(statement_text, by = 'order') %>%
  select(-sta_id, -f, -study) %>%
  pivot_wider(names_from = 'var', values_from = 'f_score')
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
  slice(1:26, 29, 27, 28, 30) %>%
  mutate(id = 1:nrow(.)) %>%
  rowwise %>%
  mutate(text = paste0(c(id, ".", " ", text), sep = "", collapse = "")) %>%
  ungroup
tab_fsc <- tab_fsc %>%
  select(-order, -nonsig, -id) %>%
  select(text, Study_2_Acceptor, Study_2_Sceptic, Study_3_Acceptor, Study_3_Sceptic) %>%
  rename(Statement = text)
write_csv(tab_fsc, sprintf('../out/tables/SM-table-%s.csv', tab_ind))
tab_ind <- tab_ind + 1

#Table 7 from Supplementary Materials
tab_trust <- trust_models %>%
  mutate_at(vars(ends_with('AIC')), .funs = list(~ . - min(.))) %>%
  mutate_at(vars(ends_with('AIC')), .funs = list(~ specify_decimal(., 2))) %>%
  mutate(model = c('Intercept only', 'Main effect of segment membership', 'Main effect of trust', 'Main effect of segment membership, main effect of trust', 'Main effect of segment membership, main effect of trust, interaction between segment membership and trust')) %>%
  select(model, cause_AIC:mitig_AIC)
write_csv(tab_trust, sprintf('../out/tables/SM-table-%s.csv', tab_ind))
tab_ind <- tab_ind + 1

#Table 8 from Supplementary Materials
tab_news <- news_models %>%
  mutate_at(vars(ends_with('AIC')), .funs = list(~ . - min(.))) %>%
  mutate_at(vars(ends_with('AIC')), .funs = list(~ specify_decimal(., 2))) %>%
  mutate(model = c('Intercept only', 'Main effect of segment membership', 'Main effect of news type', 'Main effect of segment membership, main effect of news type', 'Main effect of segment membership, main effect of news type, interaction between segment membership and news type')) %>%
  select(model, conseq_AIC:mitig_AIC)
write_csv(tab_news, sprintf('../out/tables/SM-table-%s.csv', tab_ind))
tab_ind <- tab_ind + 1

#Table 9 from Supplementary Materials

tab_policy_v1 <- policy_models %>%
  group_by(dv) %>%
  mutate(AIC_diff = AIC_diff - min(AIC_diff)) %>%
  mutate(AIC_diff = specify_decimal(AIC_diff, 2)) %>%
  mutate(model = c('Intercept only.',
                   'Main effect of S.', 'Main effect of BC1.', 'Main effect of BC2.',
                   'Main effect of S and BC1.', 'Main effect of S and BC2.', 'Main effect of BC1 and BC2.', 'Main effect of S, BC1, and BC2.',
                   'Main effect of S, BC1, and BC2, Interaction of S with BC1.', 'Main effect of S and BC1, Interaction of S with BC1.',
                   'Main effect of S, BC1, and BC2, Interaction of S with BC2.', 'Main effect of S and BC2, Interaction of S with BC2.',
                   'Main effect of S, BC1, and BC2, Interaction of BC1 with BC2.', 'Main effect of BC1 and BC2, Interaction of BC1 with BC2',
                   'Main effect of S, BC1, and BC2, Interaction of S with BC1 with BC2.'
  )) %>%
  mutate(order = c(1:8, 10, 9, 12, 11, 14, 13, 15)) %>%
  ungroup %>%
  select(dv, model, AIC_diff, order) %>%
  pivot_wider(names_from = 'dv', values_from = 'AIC_diff') %>%
  arrange(order) %>%
  select(-order)
write_csv(tab_policy_v1, sprintf('../out/tables/SM-table-%s.csv', tab_ind))
tab_ind <- tab_ind + 1


#Table 10 from Supplementary Materials
tab_policy_v2 <- policy_models %>%
  group_by(dv) %>%
  mutate(AIC_diff = AIC_diff - min(AIC_diff)) %>%
  mutate(AIC_diff = specify_decimal(AIC_diff, 2)) %>%
  mutate(model = c('Intercept only.',
                   'Main effect of S.', 'Main effect of BC1.', 'Main effect of BC2.',
                   'Main effect of S and BC1.', 'Main effect of S and BC2.', 'Main effect of BC1 and BC2.', 'Main effect of S, BC1, and BC2.',
                   'Main effect of S, BC1, and BC2, Interaction of S with BC1.', 'Main effect of S and BC1, Interaction of S with BC1.',
                   'Main effect of S, BC1, and BC2, Interaction of S with BC2.', 'Main effect of S and BC2, Interaction of S with BC2.',
                   'Main effect of S, BC1, and BC2, Interaction of BC1 with BC2.', 'Main effect of BC1 and BC2, Interaction of BC1 with BC2',
                   'Main effect of S, BC1, and BC2, Interaction of S with BC1 with BC2.'
  )) %>%
  ungroup %>%
  select(dv, model, AIC_diff)

tab_policy_v2 <- tab_policy_v2 %>%
  mutate(order = rep(rank(as.double(tab_policy_v2$AIC_diff[1:15])), 2)) %>%
  pivot_wider(names_from = 'dv', values_from = 'AIC_diff') %>%
  arrange(order) %>%
  select(-order)
write_csv(tab_policy_v2, sprintf('../out/tables/SM-table-%s.csv', tab_ind))
tab_ind <- tab_ind + 1





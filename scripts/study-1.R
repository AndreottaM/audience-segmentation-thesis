###Libraries required
library(tidyverse) #For data wrangling
###Variables
annotation_complete <- TRUE #TRUE if '../out/study-1-list.csv' has been annotated by a human coder.


############
#Anonymise data
############
rawDir <- "../data-raw/study-1.csv"
dataDir <- "../data/study-1.csv"
idDir <- "../data-raw/study-1-id.csv"
#Check if 'raw' data exists, if so, place an anonymised version of raw data in /data
#Public repository will not contain the raw file
if(file.exists(rawDir)){
  dat.raw <- read.csv(rawDir, stringsAsFactors = F)
  #Check if data is new by accessing list of IDs associated with existing data set
  #List of IDs is stored in a private directory, contains sensitive information
  if(file.exists(idDir)){
    exist.id <- read.csv(file = idDir, stringsAsFactors = F)
  }
  #Add unique ID
  set.seed(193522) #set seed
  #Ignore first two rows (additional variable information)
  #Randomly assign a ID to each participant
  #completion_order assigns a ranking to each participant, indicating order of data collection (value 1 is 1st)
  dat.id <- dat.raw %>%
    slice(-c(1,2)) %>%
    mutate(completion_order = 1:nrow(.)) %>%
    sample_frac(1, replace = F) %>%
    mutate(id = 1:nrow(.)) %>%
    select(ResponseId, id, completion_order)
  #Map internal IDs with Qualtrics IDs, in case any participants require replacement (i.e., non-flagged data set < 70 participants)
  #This will be hidden in public repository
  dat.id %>%
    write.csv(file = idDir, row.names = F)
  #Join Id with ResponseId
  dat.raw <- dat.raw %>%
    left_join(dat.id)
  #Identify new data (e.g., data acquired after soft launch, data replaced)
  dat.raw <- dat.raw %>%
    rowwise %>%
    mutate(new = if_else(exists("exist.id"), !(ResponseId %in% exist.id$ResponseId), F))
  #First 25 participants were in the soft launch, where data quality was inspected
  dat.raw <- dat.raw %>%
    mutate(in.soft.launch = completion_order <= 22)
  #Remove variables which identify Qualtrics accounts
  dat.raw[ , -c(1:5, 7:17, (dim(dat.raw)[2]-10):(dim(dat.raw)[2] - 4))] %>%
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
#1. Remove extraneous rows
dat.m <- dat.m[-c(1:2), ] %>%
  as.tibble()
#2. Remove embedded data
# Embedded data is used by Qualtrics as background variables
# Not required for analysis
embedVars <- c("resp1", "resp2", "firstlist", "otherprompt", "qID1", "qID2", "qID3", "preamble", "freelist", "listitem", "scalenumber", "itemrefer", "showprompt")
dat.m <- dat.m %>%
  select(-embedVars)
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
dat.m <- dat.m %>% 
  rename(time = Duration..in.seconds.) %>%
  rename_v("Q54", "pif") %>%
  rename_v("Q48", "pcf") %>%
  rename_v("Q55", "pcf") %>%
  rename_v("Q42", "age") %>%
  rename_v("Q44", "gender") %>%
  rename_v("Q45", "education") %>%
  rename_v("Q46", "willing") %>%
  rename_v("Q53", "demographics") %>%
  rename_v("Q10", "define") %>%
  rename_v("Q50", "define") %>%
  rename_v("Q11", "happening") %>%
  rename_v("Q51", "happening") %>%  
  rename_v("Q13", "possible") %>%
  rename_v("Q52", "possible") %>%
  rename_v("Q1", "list.test") %>%
  rename_v("Q56", "list.test") %>%
  rename_v("Q57", "list.test.check") %>%
  rename_v("Q56", "list.test") %>%
  rename_v("X1_Q20", "list.cause") %>%
  rename_v("X1_Q24", "list.cause.comment") %>%
  rename_v("X1_Q61", "list.cause") %>%
  rename_v("X2_Q20", "list.conseq") %>%
  rename_v("X2_Q24", "list.conseq.comment") %>%
  rename_v("X2_Q61", "list.conseq") %>%
  rename_v("X3_Q20", "list.mitig") %>%
  rename_v("X3_Q24", "list.mitig.comment") %>%
  rename_v("X3_Q61", "list.mitig") %>%
  rename_v("X3_Q25", "can.mitigate") %>%
  rename_v("X3_Q60", "can.mitigate") %>%
  rename_v("Q28", "mms.cause") %>%
  rename_v("Q37", "mms.cause.comment") %>%
  rename_v("Q64", "mms.cause") %>%
  rename_v("Q32", "mms.conseq") %>%
  rename_v("Q39", "mms.conseq.comment") %>%
  rename_v("Q65", "mms.conseq") %>%
  rename_v("Q35", "mms.mitig") %>%
  rename_v("Q40", "mms.mitig.comment") %>%
  rename_v("Q66", "mms.mitig")
#4. Remove unseen questions
dat.m <- dat.m %>%
  select(-c(X1_Q25:X1_Q60_Click.Count, X2_Q25:X2_Q60_Click.Count))
#Export to /out
cleanDir <- "../out/study-1-clean.csv"
dat.m %>%
  write_csv(cleanDir)


######
#Create an csv of generated free list items
#####
#Import clean data
dat.full <- read.csv(cleanDir, stringsAsFactors = F)
#Dtermine cutoff (calculated at end of soft-launch to prevent slow responders from being recorded)
cutoff <- dat.full %>% 
  filter(in.soft.launch) %>%
  #Three responses of poor quality were removed from the sample, times were used to calculate cutoff
  bind_rows(tibble(time = c(233, 386, 232))) %>%
  summarise(cutoff = trunc(median(time)/3)) %>% 
  pull
#Flag any responses where time taken to complete is below cutoff
#Cutoff determined at soft-launch 
dat.full <- dat.full %>%
  mutate(flagged = if_else(time < cutoff, 1, 0))
#Create a tibble for free list data
tb.list <- dat.full %>%
  select(id,condition,can.mitigate,list.cause,list.cause.comment,list.conseq,list.conseq.comment,list.mitig,list.mitig.comment,flagged,new)
#Create long format for list contents and comment contents
tb.list <- tb.list %>%
   gather(list_key, list_content, list.cause, list.conseq, list.mitig) %>%
   gather(comment_key, comment_content, list.cause.comment, list.conseq.comment, list.mitig.comment)
#Currently tb.list contains incompatible duplicates (e.g., mitigation comments on cause lists)
#Remove incompatible pairings
tb.list <- tb.list %>%
  filter(paste0(list_key, ".comment") == comment_key)
#list_content requires unpacking (each free list item is stored as a newline)
tb.list <- tb.list %>%
  rowwise %>%
  mutate(list_content = strsplit(list_content, "\n")) %>%
  unnest(list_content) %>%
  group_by(id, list_key, add = F) %>%
  mutate(rank = row_number())
#only some lists are directly compatible with the Bostrom et al. (2012) scales
#Indicate these lists as a seperate variable
tb.list <- tb.list %>%
  rowwise %>%
  mutate(mms_compat = ifelse(condition == 3, 0, 
                             ifelse(list_key == "list.mitig" & can.mitigate == 3, 0, 1)
                             ))
#allow human coder to classify individual responses
tb.list <- tb.list %>%
  mutate(human_code = -1) %>%
  mutate(human_comment = "")
listDir <- "../out/study-1-list.csv" #output file
if(file.exists(listDir)){
  #Import file, make human_comment readable
  dat.list <- read.csv(listDir, stringsAsFactors = F) %>%
    mutate(human_comment = as.character(human_comment)) %>%
    mutate(human_comment = if_else(is.na(human_comment), "", human_comment))
  #Check if any new id
  #  union(tb.list)
}
#Arrange rows, such that compatible items are addressed first
#Easier to examine each participant lists
tb.list <- tb.list %>%
  arrange(list_key, flagged, desc(mms_compat), id, rank)
#If listDir already exists, then append new data to list
#Otherwise write all data to list
if(!file.exists(listDir)){
  #Export this tibble as a csv file
  #Remove 'new' variable (as to not bias data collection) and change arrangement of variables to faciliate readability
  tb.list %>%
    select(-new) %>%
    arrange(list_key, flagged, desc(mms_compat), id, rank) %>%
    write_csv(listDir, na = "")
} else{
  #Read old.list
  old.list <- read_csv(listDir, col_types = list(id = col_integer(), condition = col_integer(), can.mitigate = col_integer(), flagged = col_number(), list_key = col_character(), comment_key = col_character(), comment_content = col_character(), list_content = col_character(), rank = col_integer(), mms_compat = col_number(), human_code = col_character(), human_comment = col_character()))
  #Add new IDs to old list
  old.list <- old.list %>%
    #First, match old IDs to completion_order
    left_join(exist.id) %>%
    #Remove old IDs
    select(-id, -ResponseId) %>%
    #Add new IDs
    left_join(select(dat.full, c(id, completion_order))) %>%
    #Remove completion_order
    select(-completion_order)
  #Append new data
  tb.list %>%
    mutate(human_code = as.character(human_code)) %>%
    filter(new) %>%
    select(-new) %>%
    bind_rows(old.list) %>%
    arrange(list_key, flagged, desc(mms_compat), id, rank) %>%
    write_csv(listDir, na = "")
}

###############
#Analyse results
##############

if(annotation_complete){
#Give each code its own line (some participant responses may contain multiple items per line)
tb.res <- read.csv(listDir, stringsAsFactors = F) %>%
  mutate(human_code = strsplit(human_code, ",")) %>%
  unnest(human_code) %>%
  mutate(human_code = as.integer(human_code))
if(file.exists(idDir)){
  #Extract flagged IDs
  flagged_ids <-  tb.res %>%
    group_by(id) %>%
    filter(any(flagged == 1)) %>%
    pull(id) %>%
    unique
  #Identify Qualtrics IDs which correspond to these IDs
  dat.id %>%
    filter(id %in% flagged_ids) %>%
    select(ResponseId) %>%
    write_csv("../out/study-1-flagged-ids.csv")
}
tb.res <- tb.res %>%
  group_by(id) %>%
  filter(all(flagged == 0))
#Examine compatible lists
tb.res_com <- tb.res %>%
  filter(mms_compat == 1) %>%
  filter(human_code > -1)
#Calculate relative proportions
tb.res_com <- tb.res_com %>%
  group_by(list_key, human_code) %>%
  mutate(n = length(unique(id))) %>%
  group_by(list_key) %>%
  mutate(freq = n/length(unique(id)))
#Add completion order information
tb.res_com <- tb.res_com %>%
  ungroup %>%
  left_join(select(dat.full, c('id', 'completion_order')))


###
#Examine mms_compat lists
###
library("parallel")
library("doParallel")
library("doRNG") #Used to make dopar loops reproducible 


#Collect data of interest
#Remove duplicated items from list
keys <- c('list.cause', 'list.conseq', 'list.mitig')
tb_mms <- tb.res_com %>%
  group_by(list_key) %>%
  select(id, list_key, human_code, n, freq) %>%
  distinct

tb_mms_random <- tb_mms %>%
  group_by(list_key) %>%
  mutate(unique_concepts = length(unique(human_code))) %>%
  group_by(id, list_key) %>%
  summarise(list_length = length(human_code), unique_concepts = unique(unique_concepts))

#Determine probability that concept would appear in a given list, if concepts were randomly selected from a the total set of unique concepts
tb_mms_random <- tb_mms_random %>%
  mutate(all_concept_comb = factorial(unique_concepts)/(factorial(unique_concepts - list_length)*factorial(list_length))) %>%
  mutate(specific_concept_comb = factorial(unique_concepts - 1)/(factorial((unique_concepts-1)-(list_length-1))*factorial(list_length-1))) %>%
  mutate(specific_concept_prob = specific_concept_comb/all_concept_comb)

sample_resultsDir <- '../out/study-1-sample_results.csv'

if(!file.exists(sample_resultsDir)){
  #Determine if frequency of each item differs from randomly generated lists (bootstraps)
  cl <- makeCluster(detectCores()-1) # create a cluster with 2 cores
  registerDoParallel(cl)
  
  nboots <- 1000
  tt <- {}
  tt <- foreach (brep=1:nboots,
                 .packages = "tidyverse",
                 .options.RNG=48320) %dorng% {
                   write(brep, "../out/parlog.txt") #take current bootstrap
                   #loop through each list, and compute a subsample of N-5
                   list_boot <- lapply(keys, function(x){
                     #Extract relevant data
                     dat.boot <- tb_mms %>%
                       select(-n, -freq) %>%
                       filter(list_key == x)
                     #Extract relevant concepts
                     list_concepts <- dat.boot %>%
                       filter(list_key == x) %>%
                       pull(human_code) %>%
                       unique
                     #Create a randomly generated list for each participant
                     dat.boot <- dat.boot %>%
                       group_by(id) %>%
                       mutate(human_code_boot = sample(list_concepts, n(), replace = FALSE))
                     #Calculate frequency of items
                     dat.boot <- dat.boot %>%
                       gather(var, value, -id, -list_key) %>%
                       group_by(var) %>%
                       mutate(num_of_participants = length(unique(id))) %>%
                       select(-id) %>%
                       group_by(var, value) %>%
                       count(value, list_key, num_of_participants) %>%
                       mutate(freq = n/num_of_participants)                   
                     return(dat.boot)
                   })
                   return(do.call("rbind", list_boot))
                 }
  stopCluster(cl)
  
  #Contains frequency data
  df_freq <- tb_mms %>%
    ungroup %>%
    select(list_key, human_code, freq) %>%
    distinct
  
  kk <- do.call('rbind', tt) #contains all bootstrap outcomes
  
  #calculate bootstrap means and CIs
  df_freq <- kk %>%
    ungroup %>%
    filter(var == 'human_code_boot') %>%
    group_by(list_key, value) %>%
    summarise(mean_boot = mean(freq),
              ci_low_boot = quantile(freq,probs = .025),
              ci_high_boot = quantile(freq, probs = .975)) %>%
    rename(human_code = value) %>%
    right_join(df_freq, by = c('list_key', 'human_code')) %>%
    ungroup
  
  df_freq_sample <- df_freq
  
  #Add codebook to frequency data
  #Update results tibble to indicate salient items (before cutoff)
  #and peripheral items (after cutoff)
  #all items stated by 1 person are below threshold
  df_freq_sample <- df_freq_sample %>%
    rowwise %>%
    mutate(salient = freq > (ci_high_boot - 1e-8)) %>%
    ungroup
  
  
  df_freq_sample %>%
    arrange(list_key, desc(freq)) %>%
    write_csv(sample_resultsDir)
}else{
  df_freq_sample <- read_csv(sample_resultsDir, cols(
    list_key = col_character(),
    human_code = col_integer(),
    freq = col_double(),
    boot_mean = col_double(),
    boot_ci_low = col_double(),
    boot_ci_high = col_double(),
    salient = col_logical()), col_names = TRUE)
}

#######
#Calculate if salient items differ in a smaller subsample
#####


subsample_resultsDir <- '../out/study-1-subsample_results.csv'

if(!file.exists(subsample_resultsDir)){
cl <- makeCluster(detectCores()-1) # create a cluster with 2 cores
registerDoParallel(cl)

start_time <- Sys.time()
nboots <- 1000
tt <- {}
tt <- foreach (brep=1:nboots,
               .packages = "tidyverse",
               .options.RNG=48320) %dorng% {
                 write(brep, "../out/parlog.txt") #take current bootstrap
                 #loop through each list, and compute a subsample of N-5
                 list_boot <- lapply(keys, function(x){
                   #Extract relevant data
                   dat.boot <- tb_mms %>%
                     select(-n, -freq) %>%
                     filter(list_key == x) %>%
                     ungroup %>%
                     mutate(ss = ifelse(id %in% sample(unique(id), 5, replace = FALSE), 0, 1)) %>%
                     filter(ss == 1) %>%
                     select(-ss)
                     
                   #Extract relevant concepts
                   list_concepts <- dat.boot %>%
                     filter(list_key == x) %>%
                     pull(human_code) %>%
                     unique
                   #Create a randomly generated list for each participant
                   dat.boot <- dat.boot %>%
                     group_by(id) %>%
                     mutate(human_code_boot = sample(list_concepts, n(), replace = FALSE))
                   #Calculate frequency of items
                   dat.boot <- dat.boot %>%
                     gather(var, value, -id, -list_key) %>%
                     group_by(var) %>%
                     mutate(num_of_participants = length(unique(id))) %>%
                     select(-id) %>%
                     group_by(var, value) %>%
                     count(value, list_key, num_of_participants) %>%
                     mutate(freq = n/num_of_participants)                   
                   return(dat.boot)
                 })
                 return(do.call("rbind", list_boot))
               }
stopCluster(cl)
end_time <- Sys.time()
#Time to run: 23.56778 secs

kk <- do.call('rbind', tt) #combine bootstrap results into single tibble

#Calculate freqency data for subsamples
df_freq <- kk %>%
  filter(var == 'human_code') %>%
  group_by(list_key, value) %>%
  summarise(freq_ss = mean(freq)) %>%
  rename(human_code = value)


#calculate bootstrap means and CIs
df_freq <- kk %>%
  ungroup %>%
  filter(var == 'human_code_boot') %>%
  group_by(list_key, value) %>%
  summarise(mean_boot_ss = mean(freq),
            ci_low_boot_ss = quantile(freq,probs = .025),
            ci_high_boot_ss = quantile(freq, probs = .975)) %>%
  rename(human_code = value) %>%
  right_join(df_freq, by = c('list_key', 'human_code')) %>%
  ungroup

df_freq_subsample <- df_freq

# Calculate which are salient
df_freq_subsample <- df_freq_subsample %>%
  rowwise %>%
  mutate(salient_ss = freq_ss > (ci_high_boot_ss - 1e-8)) %>%
  ungroup

df_freq_subsample %>%
  arrange(list_key, desc(freq_ss)) %>%
  write_csv(subsample_resultsDir)
}else{
  df_freq_subsample <- read_csv(subsample_resultsDir, col_names = TRUE, col_types = cols(
    list_key = col_character(),
    human_code = col_integer(),
    mean_boot_ss = col_double(),
    ci_low_boot_ss = col_double(),
    ci_high_boot_ss = col_double(),
    freq_ss = col_double(),
    salient_ss = col_logical()
  ))
}

#Determine which concepts are salient in sample XOR subsample
df_freq_sample %>%
  left_join(df_freq_subsample, by = c('list_key', 'human_code')) %>%
  filter(salient != salient_ss) 
#One borderline-salient concepts were not deemed salient in the N-5 subsamples


##########################
####Figures and Tables####
##########################

#codebook of concepts
codebook <-lapply(c('cause', 'consequence', 'mitigation'), function(x) {
  paste0(c('../codebook/study-1/', x, '.csv'), collapse = '', sep = '') %>%
    read_csv(col_names = TRUE, col_types =  cols(
      code = col_double(),
      Item = col_character(),
      mms_compat = col_factor(),
      item.id = col_character()
    )) %>%
    mutate(list_key = ifelse(x == 'cause', keys[1], ifelse(x == 'consequence', keys[2], keys[3])))
}) %>%
{do.call('rbind', .)} %>%
  rename(human_code = code) %>%
  rename(item_name = Item)
#APA theme
#Source: Sakaluk, J. K., & Short, S. D. (2016). A Methodological Review of Exploratory Factor Analysis in Sexuality Research: Used Practices, Best Practices, and Data Analysis Resources. Journal of Sex Research.
#Made own alterations
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

fill_mms_compat <- c('#686868', '#181818')



#Generate figures for Supplementary Material
p <- lapply(keys, function(x){
  df<-df_freq_sample %>%
    filter(list_key == x) %>%
    arrange(desc(freq)) %>%
    head(12) %>%
    left_join(select(filter(codebook, list_key == x), human_code, item_name, mms_compat), by = 'human_code') %>%
    mutate(index = fct_reorder(as.factor(row_number()), -freq))
  
  img <- df %>%
    ggplot(aes(x = index, y=freq, group = 1, fill = mms_compat)) +
    geom_bar(stat = 'identity')+
    scale_fill_manual(values = fill_mms_compat) +
    geom_bar(aes(y=mean_boot), stat = 'identity', fill = '#E8E8E8', alpha = 0.70) + # '#606060'
    geom_errorbar(aes(ymin=ci_low_boot, ymax = ci_high_boot), colour="#505050", width=.5, position = position_dodge(0.70), size = .70) +
    scale_y_continuous(name = "Proportion", limits = c(0, ceiling(max((df$freq))*10)/10))+
    scale_x_discrete(name='Concept')+
    scale_shape_manual(values=1)+
    apatheme+
    theme(
      legend.title = element_blank(),
      legend.position = "none"
    )
  
  return(list(dat = df, plot = img))
})

#Output figures and legend for supplementary material
lapply(keys, function(x){
  #Index number
  ind <- which(x == keys)
  #Figure number
  num <- ind + 2
  #Figure 3, 4, 5
  ggsave(sprintf('../out/figures/SM-figure-%s.png', num), plot = p[[ind]][[2]], width=4, height=5, unit='in', dpi=300)
  #Concept labels for Figure 3, 4, 5
  p[[ind]][[1]] %>%
    select(index, item_name) %>%
    write_csv(sprintf('../out/figures/SM-figure-%s-legend.csv', num))
  return('Saved in figure directory')
})

}
library(tidyverse) #for ...everything

####
#Purpose of script
####
#To compile stimuli into a single file
#Randomise order where appropriate

#Import seperate data files
dat.con <- read.csv('consequence.csv', stringsAsFactors = F) %>%
  as.tibble %>%
  rename(p = Probability....) %>%
  mutate(p = as.character(p)) %>%
  mutate(Block = "consequence")
dat.cause <- read.csv('cause.csv', stringsAsFactors = F)  %>%
  as.tibble %>%
  mutate(Include = T) %>%
  rename(p = Proportion) %>%
  mutate(Block = paste0("cause-",Block))

#Create a stimuli tibble
stimuli <- dat.con %>%
  as.tibble %>%
  bind_rows(dat.cause)

#Filter out irrelevant stimuli
stimuli <- stimuli %>%
  filter(Include) %>%
  arrange(Block)

#For all other blocks than cause-GHG-eq, randomly assign an order of administration
#Administer block cause-GHG-eq as is. Current order facilitates better reading
set.seed(361948)
stimuli <- stimuli %>%
  filter(Block != "cause-GHG-eq") %>%
  group_by(Block) %>%
  sample_frac(1, replace = F) %>%
  bind_rows(filter(stimuli, Block == "cause-GHG-eq")) %>%
  mutate(order = 1:length(Block)) %>%
  ungroup %>%
  arrange(Block)

#Export file
stimuli %>%
  write.csv(file = 'stimuli.csv', row.names = F)

#Blocks
#Four seperate sets of questions exist: Mitigation (mitig), Consequence (conseq), Cause specific (cause-s: "cause-GHG-eq" & "cause-RF-emitted-compound") and cause-general (cause-g: "cause-RF-anthropogenic"). Two seperate blocks are subsumed under cause-s as there are dependencies between the two. The description for "cause-GHG-eq" includes information useful for "cause-RF-emitted-compound" and therefore should always be administered after it.
#As per previous study, a diagram-balanced latin square design will be used
create_dbls <- function(n, as.letter = F){
  #Takes number of conditions (n) and returns the diagram-balanced Latin square(s)
  #If as.letter = T, return squares in letter format, if false/null, use numbers
  #If n > 26, letters cannot be used
  if(as.letter){
    if(n > 26){
      print("number of conditions are greater than number of available letters (26), using numbers instead")
      as.letter <- F
    }
  }
  #Firstly, create one empty cyclic Latin square
  cls <- matrix(1:n, n, n, byrow = T)
  #Create two empty diagram-balanced Latin squares
  dbls <- matrix(0, n, n)
  for (i in 1:n){
    #Cycle previous row to create a cyclic square, unnecessary for first row
    if(i > 1){
      old.seq <- cls[i - 1, ] #old sequence, to be shifted one to the left
      cls[i, ] <- old.seq %% n + 1 #shifts sequence to the left
    }
    #Interleave sequence from cyclic latin square
    il.seq <- suppressWarnings(as.vector(rbind(cls[i, ], rev(cls[i, ]))))
    dbls[i, ] <- il.seq[1:n]
  }
  #Create a warning to ensure users know which squares ensure counterbalance
  msg <- if_else((n %% 2) > 0, "Use both Diagram-Balanced Latin Squares to ensure adequate counterbalance",
                 "Use either Diagram-Balanced Latin Square to ensure adequate counterbalance")
  print(msg)
  #If necessary, conver to letters
  if(as.letter){
    dbls <- apply(dbls, c(1,2), function(x){return(LETTERS[x])})
  }
  #Return the diagram-balanced Latin squares
  list(dbls, dbls[ , n:1]) %>%
    return()
}
#Set seed
set.seed(891739)
set_of_names <- c("mitig", "conseq", "cause-s", "cause-m")
#For even n, two dbls exist, we select one
dbls <- create_dbls(4) %>%
  {.[[sample(1:2,1)]]}
dbls <- dbls %>%
  apply(1:2, function(x){return(set_of_names[x])})
#Write dbls to .csv
dbls %>%
  write.table(file="latin-square.txt", row.names=FALSE, col.names=FALSE)
#Display latin square
dbls













  



library(tidyverse) #load tidyverse
#retrieve Q sort statements
dat.q <- read_csv("q-statements.csv", col_names = T, 
                  cols(id = col_integer(), 
                       statement = col_character()))
#randomise statement order
set.seed(10243218) #set seed
df.r <- dat.q %>%
  sample_n(nrow(.), replace = F) %>%
  rowid_to_column("order")
#export order
df.r %>%
  write_csv("q-statements-ordered.csv")
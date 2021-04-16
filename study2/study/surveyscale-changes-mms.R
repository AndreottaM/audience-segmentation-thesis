library(tidyverse) #for data wrangling

mms <- '../../study1/mental-model-scales-study1.csv' %>%
  read_csv() %>%
  mutate(scale = factor(scale)) %>%
  select(-starts_with('order.item'), -item.id) %>%
  rename(item.text = item.study1)
  
head(mms)

#Add new item
item_other <- 'Other natural causes (e.g., fluctuations in the sun, changes in the Earth's axis, the Earth's magnetic field, etc.)'
mms <- mms %>%
  add_row(item.original = NA,
          scale = 'cause',
          item.text = item_other)

#Add variable to fix order of cause scale
#Fixed ordering maintains the ordering between items, to facilitate comparisons
mms <- mms %>%
  mutate(fixed_order = if_else(item.text == 'Volcanic eruptions', 1, 
                               if_else(item.text == item_other, 2, 0
                                       )
                               )
         )

#Randomly arrange items (whilst maintaing any fixed order)
set.seed(116311414) #set seed for reproducible sampling
mms <- mms %>%
  group_by(scale, fixed_order) %>%
  sample_frac(1, replace = F)

#Assign an identifier for each item
mms <- mms %>%
  group_by(scale) %>%
  mutate(order.new = row_number()) %>%
  mutate(item.id = paste0(scale, '-', order.new))

head(mms)

#Export scale to CSV
mmsDir <- 'surveyscales-mms.csv'
mms %>%
  select(-fixed_order) %>%
  arrange(scale) %>%
  write.csv(mmsDir, row.names = F)
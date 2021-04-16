#requires library(tidyverse)
#converts a query string of 7 digits into another 7 digit sequence
qsconvert <- function(qs){
lastd <- function(num){
  #returns last digit of character or integer
  num %>%
    substr(nchar(.), nchar(.)) %>%
    as.integer() %>%
    return()
}

swap <- function(vec, x,y){
  #swaps two elements (x,y) in a vector vec
  b <- vec[x]
  vec[x] <- vec[y]
  vec[y] <- b
  return(vec)
}
uid <- qs #fetch string query
#convert string query to vector of integers
id <- strsplit(uid, split = '') %>%
  unlist %>%
  as.integer()
#used for validation of Q sort completion
#id is a vector of length 7, containing each digit of user's unique ID
#swap first four elements around
id <- id %>%
  swap(1, 4) %>%
  swap(2, 3) %>%
  swap(4, 3)
#raise element 3 to the power of element 6, take last digit
id[3] <- lastd(id[3]^id[6]);
#raise element 5 to the power of element 4, take last digit
id[5] <- lastd(id[5]^id[4]);
#multiply element 1 by 6, add 9, take last digit
id[1] <- lastd(id[1]*6+9);
#multiply element 2 by 2, add 1, take last digit
id[2] <- lastd(id[2]*2+1);
#element 4 is transformed into element 5 multipled by element 6, take last digit
id[4] <- lastd(id[5]*id[6]);
#element 6 equals the sum of elements 1 through to 4, take last digit
id[6] <- lastd(id[1] + id[2] + id[3] + id[4]);
#final element equals 10 - itself, take last digit
id[7] <- lastd(10 - id[7]);
#finally, move last element to the front of the sequence & convert to single sequence
c(id[7], id[1:6]) %>%
  paste(sep = '', collapse = '') %>%
  return()
}
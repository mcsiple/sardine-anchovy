# Mutual information test
devtools::install_github("mdscheuerell/muti")

library(muti)
library(here)
library(tidyverse)
# Calculate mutual information for each sardine-anchvoy pair
# Inputs: two vectors of numeric data
# Requires: getMARSSstates function

# Get RBF dataframe - estimated states from MARSS
load(here("R","DataCleanup","RAM_Barange_FAO_States.RData")) # dataframe RBF

test <- RBF %>% 
  filter(datasource=="Barange" & 
           region=="Benguela" & 
           variable =="ssb")

muti(x = test$Sardine.est,y = test$Anchovy.est)

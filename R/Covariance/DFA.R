# Another way to get at whether the time series are asycnhronous at short time scales is to use a DFA and fit a single trend. If sardine and anchovy are asynchronous, they should have opposite loadings at the time scale you're interested in. This starts the same as the "EstimateCovariance.R" code but uses DFA instead of MARSS. 
library(tidyverse)
source(here::here("R/DataCleanup/getMARSSstates.R"))
source(here::here("R/DataCleanup/Fill_NAs_SA.R"))
load(here::here("R/DataCleanup/allsardineanchovy_3.RData")) # dataframe alldat

region <- c("Benguela","California","Humboldt", "Kuroshio-Oyashio", "NE Atlantic")
variable <- c("rec","ssb","landings")
datasource <- c("Barange")


# Set up y so that observations are lined up after one another:


# Standardize the same way as for the MARSS model

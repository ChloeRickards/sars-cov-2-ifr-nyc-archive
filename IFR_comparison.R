# Preamble -------------------------------------------------------------------------
library(tidyverse)
library(magrittr)

IFR_comparison <- read_csv("data/IFR_comparison.csv") 

IFR_comparison$IFR <- IFR_comparison$IFR*100
IFR_comparison$LowerCI <- IFR_comparison$LowerCI*100
IFR_comparison$UpperCI <- IFR_comparison$UpperCI*100


IFR_stats_adj_comp <- IFR_stats_adj 

IFR_stats_adj_comp$IFR <- gsub("[()]",",",IFR_stats_adj_comp$IFR)
IFR_stats_adj_comp$IFR <- gsub("-",",",IFR_stats_adj_comp$IFR)
IFR_stats_adj_comp$IFR <- sub(" ", "", IFR_stats_adj_comp$IFR)

split <- strsplit(IFR_stats_adj_comp$IFR, ",")
mean_adj <- vector()
lowerCI_adj <- vector()
upperCI_adj <- vector()
for (i in 1:5) {
  vec <- split[[i]]
  mean_adj[i] <- as.numeric(vec[1])
  lowerCI_adj[i] <- as.numeric(vec[2])
  upperCI_adj[i] <- as.numeric(vec[3])
}


IFR_stats_adj_comp <- IFR_stats_adj_comp %>%
  select(-'Age class') %>%
  select(-'Population') %>%
  select(-'Estimated infected') %>%
  select(-'Deaths') %>%
  mutate('IFR' = mean_adj) %>%
  cbind('AgeMidpoint' = c(9, 31, 55, 70, 80)) %>%
  cbind('Dataset' = "Rosenberg, adjusted w/ Pollan") %>%
  cbind('LowerCI' = lowerCI_adj) %>%
  cbind('UpperCI' = upperCI_adj)
  
IFR_stats_adj_comp <- IFR_stats_adj_comp[, c(2, 3, 1, 4, 5)]

# Raw data
# 18-45, 45+
  
IFR_stats_raw_comp <- IFR_stats_raw 

IFR_stats_raw_comp$IFR <- gsub("[()]",",",IFR_stats_raw_comp$IFR)
IFR_stats_raw_comp$IFR <- gsub("-",",",IFR_stats_raw_comp$IFR)
IFR_stats_raw_comp$IFR <- sub(" ", "", IFR_stats_raw_comp$IFR)

split <- strsplit(IFR_stats_raw_comp$IFR, ",")
mean_raw <- vector()
lowerCI_raw <- vector()
upperCI_raw <- vector()
for (i in 1:2) {
  vec <- split[[i]]
  mean_raw[i] <- as.numeric(vec[1])
  lowerCI_raw[i] <- as.numeric(vec[2])
  upperCI_raw[i] <- as.numeric(vec[3])
}

IFR_stats_raw_comp <- IFR_stats_raw_comp %>%
  select(-'Age class') %>%
  select(-'Population') %>%
  select(-'Estimated infected') %>%
  select(-'Deaths') %>%
  mutate('IFR' = mean_raw) %>%
  cbind('AgeMidpoint' = c(31, 50)) %>%
  cbind('Dataset' = "Rosenberg, raw") %>%
  cbind('LowerCI' = lowerCI_raw) %>%
  cbind('UpperCI' = upperCI_raw)

IFR_stats_raw_comp <- IFR_stats_raw_comp[, c(2, 3, 1, 4, 5)]

IFR_comparison <- IFR_comparison %>%
  rbind(IFR_stats_adj_comp) %>%
  rbind(IFR_stats_raw_comp)

IFR_comparison <- arrange(IFR_comparison, IFR_comparison$AgeMidpoint)
  
# Preamble -------------------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(rstan)
library(foreach)

options(mc.cores = 5)
rm(list = ls())

# TO DO
# add toggles for: including probable deaths, whether or not to save as RDS
#

# Here are my steps for generating bootstrapped serodata:
# make rnorm draws from the Rosenberg data. Age classes: 18-35, 35-45, 45-55, 55+
# take the 5000 55+, and scale to 7 age classes 55-60, 60-65, 65-70, 70-75,
#       75-80, 80-85, 85+ using age trend data fit. Each of these 7 age classes has 5000 draws
# create 5000 duplicate draws of the scaled 0-18 serodata, weighted by population
# 

source("utils.R")

redo_fit <- T
by_boro <- F # Toggle all boroughs or by NYC total

suffix <- ("")  # suffix for result file names

n_post <- 100

age_classes <- c("[-Inf, 18)", "[18, 45)", "[45, 65)","[65, 75)","[75, Inf)")

# Epi data ---------------------------------------------------------------------
# Cumulative deaths and cumulative cases in NYC
# March 23 - May 17
# Age groups: 0-18, 18-45, 45-65, 65+
age_epidata <- read_csv("NYCHealth_cumul_data.csv") 

age_epidata <- age_epidata %>% 
  group_by(date, age_class, var) %>% 
  summarise(value = sum(Total)) %>% 
  ungroup() %>% 
  arrange(age_class, date) %>% 
  pivot_wider(id_cols = c("age_class", "date"),
              names_from = "var",values_from = "value")

# Population data --------------------------------------------------------------
# Source: https://www.baruch.cuny.edu/nycdata/population-geography/age_distribution.htm
# Adapted from census data
age_popdata <- read_csv("NYC_pop_age_boro.csv")

age_popdata<- age_popdata %>%
  rbind(data.frame(age_class = "[15, 18)", age_midpoint = 16, 0.6*age_popdata[4,-1:-2])) %>%
  rbind(data.frame(age_class = "[18, 20)", age_midpoint = 19, 0.4*age_popdata[4,-1:-2])) %>%
  filter(age_class != "[15, 20)") %>%
  arrange(age_midpoint)

# grouping age_popdata into serodata age groups
sero_pop <- age_popdata
sero_pop <- filter(sero_pop, sero_pop$age_midpoint > 18)
sero_pop$age_class[sero_pop$age_midpoint > 18 & sero_pop$age_midpoint < 35] <- "[18, 35)"
sero_pop$age_class[sero_pop$age_midpoint > 35 & sero_pop$age_midpoint < 45] <- "[35, 45)" 
sero_pop$age_class[sero_pop$age_midpoint > 45 & sero_pop$age_midpoint < 55] <- "[45, 55)" 
sero_pop$age_class[sero_pop$age_midpoint > 55] <- "[55, Inf)" 

# This gives population by age and by boro, by the serosurvey age categories
sero_pop <- aggregate(cbind(NYC_Total=sero_pop$NYC_Total, Manhattan=sero_pop$Manhattan,
                            Bronx=sero_pop$Bronx, Brooklyn=sero_pop$Brooklyn,
                            Queens=sero_pop$Queens, Staten_Island=sero_pop$Staten_Island),
                      by=list(age_class=sero_pop$age_class), FUN=sum)

#  Seroprevalence estimates ----------------------------------------------------
# Source: Rosenberg et al., specifically Table 2, NYC, age groups
prev_est_age <- read_csv("serosurvey_by_age_draws.csv") %>% 
  select(-X1) %>% 
  select(-age_midpoint) %>%
  mutate(pop = sero_pop$NYC_Total) %>%
  group_by(age_class) 

# to be filled with bootstrapped draws
prev_est_age_draws <- data.frame(age_class = character(),
                                 seropos = double(),
                                 sim = integer(),
                                 pop = double())

for (i in 1:nrow(prev_est_age)){
  draws <- rnorm(n_post,prev_est_age$seropos_mu[i],prev_est_age$seropos_sig[i])
  age_pop <- data.frame(age_class = prev_est_age$age_class,
                        pop = prev_est_age$pop)
  for (j in 1:n_post){
    prev_est_age_draws <- prev_est_age_draws %>%
      rbind(data.frame(age_class = age_pop$age_class[i], seropos = draws[j], sim = j,pop = age_pop$pop[i]))
  }
}

# I'll use this later
over_55_sero <- prev_est_age_draws$seropos[prev_est_age_draws$age_class == "[55, Inf)"]

# Serosurvey + age trend from Pollán -------------------------------------------
age_trend_data <- read.csv("serop_age_spain.csv")

# Creating the (-Inf, 5), [15, 18), [18, 20), and [85, Inf) age groups
# Deleting the (-Inf, 1), [1, 5), [15, 20), [85, 90), and [90, Inf)
# I average 85-90 and 90+ together to get 85+. 
# To do later? find and add Spain population so I can take the weighted average for 85+
# Currently modified to exclude the 0-18 population
age_trend_data<- age_trend_data %>%
  rbind(data.frame(age_class = "[-Inf, 5)", age_midpoint = 3,
                   N_POC = sum(age_trend_data$N_POC[age_trend_data$age_midpoint < 5]),
                   POC = weighted.mean(age_trend_data$POC[age_trend_data$age_midpoint < 5], age_trend_data$N_POC[age_trend_data$age_midpoint < 5]),
                   N_IA = sum(age_trend_data$N_IA[age_trend_data$age_midpoint < 5]),
                   IA = weighted.mean(age_trend_data$IA[age_trend_data$age_midpoint < 5], age_trend_data$N_IA[age_trend_data$age_midpoint < 5])
  )) %>%
  rbind(data.frame(age_class = "[15, 18)", age_midpoint = 16, age_trend_data[5,3:6]*0.6)) %>%
  rbind(data.frame(age_class = "[18, 20)", age_midpoint = 19, age_trend_data[5,3:6]*0.4)) %>%
  rbind(data.frame(age_class = "[85, Inf)", age_midpoint = 90, 
                   N_POC = sum(age_trend_data$N_POC[age_trend_data$age_midpoint > 85]),
                   POC = weighted.mean(age_trend_data$POC[age_trend_data$age_midpoint > 85], age_trend_data$N_POC[age_trend_data$age_midpoint > 85]),
                   N_IA = sum(age_trend_data$N_IA[age_trend_data$age_midpoint > 85]),
                   IA = weighted.mean(age_trend_data$IA[age_trend_data$age_midpoint > 85], age_trend_data$N_IA[age_trend_data$age_midpoint > 85])
  )) %>%
  filter(age_class != "[-Inf, 1)" & age_class != "[1, 5)" & age_class != "[15, 20)" & age_class != "[85, 90)" & age_class != "[90, Inf)") %>%
  transform(age_midpoint = as.double(age_midpoint), N_POC = as.double(N_POC), POC = as.double(POC), N_IA = as.double(N_IA), IA = as.double(IA)) %>%
  arrange(age_midpoint)

# Linear fit for whole dataset
age_sero_fit=lm(IA~age_midpoint+I(age_midpoint^2),data=age_trend_data);summary(age_sero_fit) # fits a line to the data using a linear model; summary of fit
age_trend_data$fitted=fitted(age_sero_fit) # either here or on very end - adds a column to age_trend_data of fitted IA values - basically I'll tack on a fitted column onto the serodata and use that
age_trend_data$scale=age_trend_data$fitted/mean(age_trend_data$IA[age_trend_data$age_midpoint<70]) #either here or on very end -- adds a column for scale. This is the multiplier for the serodata

# Finds ratio from 0-18 to 18-45, to prep for 0-18 conversion
IA_0to18 <- age_trend_data$IA[age_trend_data$age_midpoint < 18]
IA_0to18 <- weighted.mean(IA_0to18, age_popdata$NYC_Total[age_popdata$age_midpoint < 18])

IA_18to45 <- age_trend_data$IA[age_trend_data$age_midpoint > 18 & age_trend_data$age_midpoint < 45]
IA_18to45 <- weighted.mean(IA_18to45, age_popdata$NYC_Total[age_popdata$age_midpoint > 18 & age_trend_data$age_midpoint < 45])

ratio_0to18 <- IA_0to18/IA_18to45

# NEED TO FIND CI FOR SCALE
#rnorm(1, mean = over_55_sero*adj_serodata_55_up$scale[1], sd = sd_55_up)
# f1=lm(serop~age+age^2,data=df)
# 
# 
# predict(f1,se.fit=T

# f1=lm(serop~age+age^2,data=df)
# 
# predict(f1,se.fit=T,newdata=data.frame(age=c(5,10,22,36))
# make sure age --> age midpoint, adjust values, adjust names for data frame

adj_serodata_55_up = age_trend_data[age_trend_data$age_midpoint>55,] %>% #selects 55+
  cbind(pop = age_popdata$NYC_Total[age_popdata$age_midpoint > 55])
for (i in 1:nrow(adj_serodata_55_up)){
  prev_est_age_draws <- prev_est_age_draws %>%
    rbind(data.frame(age_class = adj_serodata_55_up$age_class[i], sim = 1:n_post, seropos = over_55_sero*adj_serodata_55_up$scale[i], pop = adj_serodata_55_up$pop[i]))
}

prev_est_age_draws <- filter(prev_est_age_draws, age_class != "[55, Inf)")

# Adjusting serodata to match desired age classes ------------------------------

# 18-45, 45-65, 65-75, and 75+ classes
draws_18_45 <- prev_est_age_draws %>%
  filter(age_class == "[18, 35)" | age_class == "[35, 45)")
draws_45_65 <- prev_est_age_draws %>%
  filter(age_class == "[45, 55)" | age_class == "[55, 60)" | age_class == "[60, 65)")
draws_65_75 <- prev_est_age_draws %>%
  filter(age_class == "[65, 70)" | age_class == "[70, 75)")
draws_75_up <- prev_est_age_draws %>%
  filter(age_class == "[75, 80)" | age_class == "[80, 85)" | age_class == "[85, Inf)")

all_draws <- list(draws_18_45, draws_45_65, draws_65_75, draws_75_up)
adj_age_classes <- age_classes[age_classes != "[-Inf, 18)" & age_classes != "all"]
dummy_draws <- data.frame(age_class = character(),
                                 seropos = double(),
                                 sim = integer(),
                                 pop = double())

for (i in 1:length(adj_age_classes)){
  new_entry <- getWeightedDraws(all_draws[[i]], n_post, adj_age_classes[i])
  dummy_draws <- dummy_draws %>%
    rbind(new_entry)
  prev_est_age <- prev_est_age %>%
    rbind(data.frame(age_class = adj_age_classes[i], seropos_mu = mean(new_entry$seropos),
          seropos_sig = sd(new_entry$seropos), pop = new_entry$pop[1]))
  }

draws_0_18 <- ratio_0to18*dummy_draws$seropos[dummy_draws$age_class == "[18, 45)"]

prev_est_age_draws <- prev_est_age_draws %>%
  rbind(data.frame(age_class = "[-Inf, 18)",
                   seropos = draws_0_18,
                   sim = 1:n_post,
                   pop = sum(age_popdata$NYC_Total[age_popdata$age_midpoint < 18]))) %>%
  rbind(dummy_draws) %>%
  filter(age_class %in% age_classes) %>%
  ungroup()

prev_est_age <- prev_est_age %>%
  rbind(data.frame(age_class = "[-Inf, 18)", seropos_mu = mean(draws_0_18),
        seropos_sig = sd(draws_0_18), pop = sum(age_popdata$NYC_Total[age_popdata$age_midpoint < 18]))) %>%
  filter(age_class %in% age_classes) %>%
  slice(match(age_classes, age_class))

# Simulation Setup -------------------------------------------------------------
# One more adjustment to age_popdata to make it line up with the desired age classes

# grouping age_popdata into serodata age groups
age_popdata$age_class[age_popdata$age_midpoint < 18] <- "[-Inf, 18)"
age_popdata$age_class[age_popdata$age_midpoint > 18 & age_popdata$age_midpoint < 45] <- "[18, 45)"
age_popdata$age_class[age_popdata$age_midpoint > 45 & age_popdata$age_midpoint < 65] <- "[45, 65)" 
age_popdata$age_class[age_popdata$age_midpoint > 65 & age_popdata$age_midpoint < 75] <- "[65, 75)" 
age_popdata$age_class[age_popdata$age_midpoint > 75] <- "[75, Inf)" 

# This gives population by age and by boro, by the epidata age categories
age_popdata <- aggregate(cbind(NYC_Total=age_popdata$NYC_Total, Manhattan=age_popdata$Manhattan,
                               Bronx=age_popdata$Bronx, Brooklyn=age_popdata$Brooklyn,
                               Queens=age_popdata$Queens, Staten_Island=age_popdata$Staten_Island),
                         by=list(age_class=age_popdata$age_class), FUN=sum)

# I have it set to only NYC_Total right now 
if (!by_boro){
  age_popdata <- data.frame(age_class = age_classes, "pop" = age_popdata$NYC_Total)
}

n_sero_weeks <- 1  # number of serosurvey weeks

# This is only used for serosurvey dates, and then for postprocessing
# Filled in with mean and 95% CI from Rosenburg
serodata <- tribble(
  ~date, ~median, ~q025, ~q975,
  "2020-04-23", 22.7, 21.5, 24.0
) %>% 
  mutate(date = as.Date(date))

# Delay distributions ----------------------------------------------------------

delay_params <- tribble(
  ~delay, ~logmu, ~logmu.sd, ~logmu.low, ~logmu.high, ~logsigma, ~logsigma.sd, ~logsigma.low, ~logsigma.high,
  "inc", 1.57, NA, 1.44, 1.69, 0.65, NA, 0.56, 0.73,     # Bi et al. 2020
  "report", 1.49, 0.065, NA, NA, 0.756, 0.0458, NA, NA,  # Scire et al. 2020
  "symp_sero", 2.34, 0.114, NA, NA, 0.38, 0.26, NA, NA,  # Stringhini et al. 2020
  "report_death", 2.1, 0.055, NA, NA, 0.87,  0.039, NA, NA # DGS
) %>% 
  group_by(delay) %>% 
  mutate(
    logmu.sd = case_when(is.na(logmu.sd) ~ getSD(logmu, logmu.low, logmu.high), T ~ logmu.sd),
    logsigma.sd = case_when(is.na(logsigma.sd) ~ getSD(logsigma, logsigma.low, logsigma.high),  T ~ logsigma.sd)
  ) %>% 
  ungroup() %>% 
  mutate(mean = lnormMean(logmu, logsigma),   # Mean and SD on natural scale
         sd = lnormSD(logmu, logsigma))

# Sample delay PDFs 
n_dist_samples <- n_post # number of PDF samples to draw from
nmax <- 100  # max number of days in PMFs

# Initialize PDF matrices
pdf_inc <- matrix(rep(0, nmax * n_dist_samples), nrow = n_dist_samples)
pdf_report <- pdf_inc
pdf_symp_sero <- pdf_inc
pdf_sero <- pdf_inc
pdf_report_death <- pdf_inc
pdf_symp_death <- pdf_inc
pdf_death <- pdf_inc

# Sample distributions
for (i in 1:n_dist_samples) {
  pdf_inc[i, ] <- setDelayPMF(delay_params, "inc", nmax, rnd = T)
  pdf_report[i, ] <- setDelayPMF(delay_params, "report", nmax, rnd = T)
  pdf_symp_sero[i, ] <- setDelayPMF(delay_params, "symp_sero", nmax, rnd = T)
  pdf_sero[i, ] <- convolve(pdf_inc[i, ], rev(pdf_symp_sero[i, ]), type = "o")[1:nmax]  # convolve with incubation period
  pdf_report_death[i, ] <- setDelayPMF(delay_params, "report_death", nmax, rnd = T)
  pdf_symp_death[i, ] <- convolve(pdf_report[i, ], rev(pdf_report_death[i, ]), type = "o")[1:nmax]  # convolve with delay from symptoms to hospitalization
  pdf_death[i, ] <- convolve(pdf_inc[i, ], rev(pdf_symp_death[i, ]), type = "o")[1:nmax]  # convolve with incubation period
}

# Inference --------------------------------------------------------------------
# Settings for stan
control <- list(adapt_delta = .9, max_treedepth = 12, metric = "dense_e")
nwarmup <- 400
niter <- nwarmup + 100
n_chains <- 5

IFRs <- foreach(agec = age_classes[age_classes != "all"],
                .combine = rbind) %do% 
  { res_file <- paste0("age_stratified_IFRs_distsample_", agec, suffix, ".rds") 
  stanfit_file <- paste0("stanfit_age_stratified_IFRs_distsample_", agec, suffix, ".rds") 
  
  if (!file.exists(res_file) | redo_fit) {
    
    # Compile stan model
    ifr_stan <- stan_model("ifr_betabinomial.stan")
    
    epidata <- filter(age_epidata, age_class == agec)
    pop_age <- age_popdata$pop[age_popdata$age_class == agec]  # population of the age class
    
    # Dates on which serosurveys were done
    ind_date <- which(epidata$date %in% serodata$date)
    
    # Deaths for model
    deaths <- epidata$death_cumul[ind_date]
    n_data <- nrow(epidata)
    
    # Seroprevalence estimates 
    # Posterior draws as matrix with each row corresponding to a given date, 
    # and each column to a sample. 
    thetas <- prev_est_age_draws %>% 
      filter(age_class == agec) %>% 
      select(seropos, sim) %>% 
      spread(sim, seropos) %>% 
      as.matrix()
    
    # Posterior draws of seroprevalent population
    Isero <- pop_age * thetas
    
    # Initialize Istar and phi matrices [n data points x n samples]
    Istar <- matrix(rep(0, n_data * n_dist_samples), ncol = n_dist_samples)
    phi <- Istar
    
    # Loop over samples of delay distributions
    for (i in 1:n_dist_samples) {
      
      # Deconvolve state of cumulative infections up to proportionality (Istar)
      Istar[, i] <- computeIstar(epidata$case_cumul, pdf_inc[i, ], pdf_report[i, ])
      
      # Compute phi(t)
      phi[, i] <- (convolve(Istar[, i], rev(pdf_death[i, ]), type = "o")/
                     convolve(Istar[, i], rev(pdf_sero[i, ]), type = "o"))[1:n_data]
    }
    
    # Get phis for dates of serosurveys
    phis <- phi[ind_date, ]
    print(phis)
    
    # Compute number of Infected at risk of dying
    I <- round(Isero * phis) 
    
    # Set min of I to observed deaths to avoid errors in stan
    # Set max of I to total population in this age class - maybe change to estimated infections (sero * pop)
    # original code does not contain max - probably shouldn't contain max either
    for (i in 1:nrow(I)) {
      I[i, I[i,] < deaths[i]] <- deaths[i] 
      I[i, I[i,] > pop_age] <- pop_age # use this to set the max
      I[i,] <- as.integer(I[i,])
    }
    
    # Data for stan model
    data <- list(N = n_sero_weeks, 
                 M = n_post, 
                 deaths = as.array(deaths), 
                 I = I)
    
    ifr_stanfit <- sampling(ifr_stan,
                            data = data,
                            init = rndInit, 
                            chains = n_chains,
                            warmup = nwarmup, 
                            iter = niter,
                            control = control,
                            save_warmup = FALSE)
    
    saveRDS(ifr_stanfit, file = "stanfit_file")
    
    # Extract and save
    ifrs <- extract(ifr_stanfit, pars = "IFR")$IFR
    saveRDS(ifrs, file = "res_file")

  } else {
    # Load results
    ifrs <- readRDS(res_file)
  }
  tibble(ifr = ifrs, age_class = agec)
  }

save(list = ls(), file = "data_setup_adj.rda")
saveRDS(IFRs, file = paste0("age_stratified_IFRs_adj", suffix, ".rds"))

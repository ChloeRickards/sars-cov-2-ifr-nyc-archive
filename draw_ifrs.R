# Notes ------------------------------------------------------------------------

# This code takes in daily death and case data and seroprevalence data, and
# creates a distribution of potential infection fatality rates based off of age

# Be sure to change this line to set the working directory to the folder 
# that you are working from 
setwd('/Users/chloerickards/Desktop/Work Projects/SARS-CoV-2 IFR/sarscov2-ifr-nyc')

# Preamble ---------------------------------------------------------------------
# These are essential for the function of the code
library(tidyverse)
library(magrittr)
library(rstan)
library(foreach)
library(truncnorm)
options(mc.cores = 5)

# Clear environment
rm(list = ls())

# loads necessary functions
source("utils.R")

# User inputs ------------------------------------------------------------------
# These are inputs that you can adjust to match your own scenario

# NOTE: you are welcome to use your own data to calculate IFR. However, you MUST
# provide the following:
# Age classes
# Age midpoints
# Age-stratified daily deaths and daily cases
# Age-stratified population data
# Age-stratified serodata (if applicable)
# All age classes must be the same

# Are you adjusting your serovalues?
# If so, you MUST provide age-stratified serodata
adjusted <- T
Pollan <- F

# This calibrates the suffixes and loads the adjusting data
# If you are using a different source than Pollan et al. 2020 or Ward et al. 2020
# then you WILL have to adjust the linear model in the main data
if (adjusted & Pollan) {
  age_trend_data <- read.csv("serop_age_spain.csv")
  suffix <- "adj"
} else if (adjusted & !Pollan) {
  age_trend_data <- read.csv("serop_age_ward.csv")
  suffix <- "adj_ward"
} else {
  suffix <- "raw"
}

# List of age classes and approximate age midpoints that you are using
age_classes <- c("0-17", "18-44", "45-64","65-74","75+")
age_midpoints <- c(9, 31, 55, 70, 85)

# Cumulative case and cumulative death data, organized by age and data
# This can be replaced by your own data, but it must contain the following
# columns with the following headings: age_class, date, case_cumul, death_cumul
#
# Default NYC case and death epidata:
# March 23 - May 17
# Age groups: 0-18, 18-45, 45-65, 65+
# Source: NYC Dept. of Health 
age_epidata <- read_csv("NYCHealth_cumul_data.csv") 

# Population by age 
# This can be replaced by your own data, but it must contain the following
# columns with the following headings: age_class, age_midpoint, total
#
# The age classes also MUST match the previously defined age classes
#
# Adapted from NYC population data, which is based on census data
# Source: https://www.baruch.cuny.edu/nycdata/population-geography/age_distribution.htm
age_popdata <- read_csv("NYC_pop_age_total.csv")

# Number of draws
# The higher, the better and more accurate, but the longer it'll take to run
n_post <- 1000

# Redo the entire fit?
redo_fit <- T

# Loading and formatting data --------------------------------------------------

# Default NYC case and death epidata:
# March 23 - May 17
# Age groups: 0-18, 18-45, 45-65, 65+
# Source: NYC Dept. of Health 
# Reformat age epidata
age_epidata <- age_epidata %>% 
  group_by(date, age_class, var) %>% 
  summarise(value = sum(Total)) %>% 
  ungroup() %>% 
  arrange(age_class, date) %>% 
  pivot_wider(id_cols = c("age_class", "date"),
              names_from = "var",values_from = "value") %>%
  mutate(date = as.Date(date))


# Source: Rosenberg et al., specifically Table 2, NYC, age groups
# "pop" column comes from the CUNY resource
# 18-34 and 35-44 were combined in a population-weighted average
prev_est_age <- read_csv("serosurvey_by_age_draws.csv") %>% 
  group_by(age_class) 

# to be filled with bootstrapped draws
prev_est_age_draws <- data.frame(age_class = character(),
                                 age_midpoint = integer(),
                                 seropos = double(),
                                 sim = integer(),
                                 pop = double())

# Creating seroprevalence draws ------------------------------------------------

# This is where I'll adjust the serovalues
# If you just want to use raw values, then set "adjusted" as F above 
# Default: adjusted with values from Pollan et al., I did a bit of Excel manipulation
# before exporting here
# Adjustments to Pollan: split 15-19 to 15-17 and 18-19, deleted infants, made 1-4
# into 0-4, condensed 85-89 and 90+ into 85+ to fit with NYC population data
#
# The adjustment method is specific to the NYC dataset - it might look different
# for other age classes
if (adjusted){
  
  # This first adjustment is for Pollan et al. and uses the seroprev_age_spain dataset
  # It includes a ratio for 0-17:18-44, which allows us to calculate the 0-17 seroprevalence
  # (and later, IFRs). This also creates the linear model used to scale the 55+ ages
  if (suffix == "adj") {
    # 0-18 age class
    # Method: Generate draws to get the 0-17:18-44 ratio
    ratio_draws_0to17 <- rep(0, n_post)
    for (i in 1:n_post){
      # 0-17 draws - draw 1 number from each class mu and sig
      # but what do I do with these 4 numbers? Average them? pick a random one??
      l_0to17 <- length(age_trend_data$age_class[age_trend_data$age_midpoint < 18])
      pollan_draws_0to17 <- rep(0, l_0to17)
      for (j in 1:l_0to17){
        pollan_draws_0to17[j] <- rlnorm(1,meanlog=log(age_trend_data[j,]$MeanLN),sdlog=log(age_trend_data[j,]$SDLN))
      }
      
      # 18-44 draws - draw 1 number from each class mu and sig
      l_18to44 <- length(age_trend_data$age_class[age_trend_data$age_midpoint > 18 & age_trend_data$age_midpoint < 45])
      pollan_draws_18to44 <- rep(0, l_18to44)
      for (j in 1:l_0to17){
        pollan_draws_18to44[j] <- rnorm(1,mean=(age_trend_data[j,]$IA)/100,sd=(age_trend_data[j,]$UpperCI-age_trend_data[j,]$IA)/200)
      }
      
      # creates a single draw for ratio 0-17:18-44
      ratio_draws_0to17[i] <- mean(pollan_draws_0to17)/mean(pollan_draws_18to44)
      
    }
    
    #18-44 class - no adjustment needed, just population-weighted draws from Rosenberg
    draws_18to44 <- data.frame(age_class = "18-44", 
                               age_midpoint = 31,
                               seropos = rtruncnorm(n = n_post, a = 0, b = 1,
                                                    mean = prev_est_age$seropos_mu[prev_est_age$age_midpoint < 45],
                                                    sd = prev_est_age$seropos_sig[prev_est_age$age_midpoint < 45]),
                               sim = 1:n_post,
                               pop = age_popdata$total[age_popdata$age_midpoint > 18 & age_popdata$age_midpoint < 45])
    
    # 0-17 class, based on 18-44 class and ratio from earlier
    draws_0to17 <- data.frame(age_class = "0-17", 
                              age_midpoint = 9,
                              seropos = ratio_draws_0to17*draws_18to44$seropos,
                              sim = 1:n_post,
                              pop = age_popdata$total[age_popdata$age_midpoint < 18])
    
    # Create a fitted line to Pollan data for 45+
    age_trend_data = age_trend_data[age_trend_data$age_midpoint>45,]
    age_sero_fit=lm(IA~age_midpoint+I(age_midpoint^2),data=age_trend_data);summary(age_sero_fit)
    
    # Append fitted values, st dev, scale for values, and st dev for scale
    age_sero_predict=predict(age_sero_fit, se.fit = T)
    age_trend_data$fitted=age_sero_predict$fit
    age_trend_data$sd=age_sero_predict$se.fit
    age_trend_data$scale=age_trend_data$fitted/mean(age_trend_data$IA[age_trend_data$age_midpoint<70])
    age_trend_data$scale_sd=age_trend_data$sd/mean(age_trend_data$UpperCI[age_trend_data$age_midpoint<70] - age_trend_data$IA[age_trend_data$age_midpoint<70])
    
  } else {
    # this block uses data from Ward et al. 2020 to calculate scaling factors for the 55+ class
    # it does not include serodata for the 0-17 class
    
    # this creates dummy draw for the 0-17 class
    draws_0to17 <- data.frame(age_class = "0-17", 
                              age_midpoint =9,
                              seropos = rep(0, n_post),
                              sim = 1:n_post,
                              pop = age_popdata$total[age_popdata$age_midpoint < 18])
    
    #18-44 class - no adjustment needed, just population-weighted draws from Rosenberg
    draws_18to44 <- data.frame(age_class = "18-44", 
                               age_midpoint = 31,
                               seropos = rtruncnorm(n = n_post, a = 0, b = 1,
                                                    mean = prev_est_age$seropos_mu[prev_est_age$age_midpoint < 45],
                                                    sd = prev_est_age$seropos_sig[prev_est_age$age_midpoint < 45]),
                               sim = 1:n_post,
                               pop = age_popdata$total[age_popdata$age_midpoint > 18 & age_popdata$age_midpoint < 45])
    
    # Create a fitted line to Ward data for 45+
    age_trend_data = age_trend_data[age_trend_data$age_midpoint>45,]
    age_sero_fit=lm(mean~age_midpoint+I(age_midpoint^2),data=age_trend_data);summary(age_sero_fit)
    
    # Append fitted values, st dev, scale for values, and st dev for scale
    age_sero_predict=predict(age_sero_fit, se.fit = T)
    age_trend_data$fitted=age_sero_predict$fit
    age_trend_data$sd=age_sero_predict$se.fit
    age_trend_data$scale=age_trend_data$fitted/mean(age_trend_data$mean[age_trend_data$age_midpoint<70])
    age_trend_data$scale_sd=age_trend_data$sd/mean(age_trend_data$UpperCI[age_trend_data$age_midpoint<70] - age_trend_data$mean[age_trend_data$age_midpoint<70])
  }
  
  # Creates raw draws for 45-54
  sero_45to54 <- rnorm(n_post, 
                        mean = prev_est_age$seropos_mu[prev_est_age$age_midpoint > 45 & prev_est_age$age_midpoint < 55],
                        sd = prev_est_age$seropos_sig[prev_est_age$age_midpoint > 45 & prev_est_age$age_midpoint < 55])
  
  # Creates raw draws for 55 and up. This will be used for 55-64, 65-75, 75+
  # Should I do separate draws for each age class? (If so, just make multiple sero_ draws)
  sero_55up <- rnorm(n_post, 
                       mean = prev_est_age$seropos_mu[prev_est_age$age_midpoint > 55],
                       sd = prev_est_age$seropos_sig[prev_est_age$age_midpoint > 55])
  
  # Creates scale draws for 55-64
  scales_55to64 <- rtruncnorm(n_post, a = 0, b = Inf,
                         mean = age_trend_data$scale[age_trend_data$age_midpoint > 55 & age_trend_data$age_midpoint < 65],
                         sd = age_trend_data$scale_sd[age_trend_data$age_midpoint > 55 & age_trend_data$age_midpoint < 65])
  
  # Element-wise multiplication of raw draws with scales for 55-64
  sero_55to64 <- sero_55up*scales_55to64
  
  # Creates population-weighted draw for 45-64
  n_post_45to54 <- as.integer(prev_est_age$pop[prev_est_age$age_midpoint > 45 & prev_est_age$age_midpoint < 55]/
                                (age_popdata$total[age_popdata$age_midpoint > 45 & age_popdata$age_midpoint < 65])*n_post)
  n_post_55to64 <- n_post - n_post_45to54
  pop_weighted_sero_45to65 <- sample(sero_45to54, n_post_45to54) %>%
    append(sample(sero_55to64, n_post_55to64))
  
  # Formatting as draws
  draws_45to64 <- data.frame(age_class = "45-64", 
                             age_midpoint = 55,
                             seropos = pop_weighted_sero_45to65,
                             sim = 1:n_post,
                             pop = age_popdata$total[age_popdata$age_midpoint > 45 & age_popdata$age_midpoint < 65])
  
  # Creates scale draws for 65-74
  scales_65to74 <- rtruncnorm(n = n_post, a = 0, b = Inf,
                         mean = age_trend_data$scale[age_trend_data$age_midpoint > 65 & age_trend_data$age_midpoint < 75],
                         sd = age_trend_data$scale_sd[age_trend_data$age_midpoint > 65 & age_trend_data$age_midpoint < 75])
  
  # Element-wise multiplication of raw draws with scales for 65-74
  draws_65to74 <- data.frame(age_class = "65-74", 
                             age_midpoint = 70,
                             seropos = sero_55up*scales_65to74,
                             sim = 1:n_post,
                             pop = age_popdata$total[age_popdata$age_midpoint > 65 & age_popdata$age_midpoint < 75])
  
  # Creates scale draws for 75+
  scales_75up <- rtruncnorm(n_post, a = 0, b = Inf,
                      mean = age_trend_data$scale[age_trend_data$age_midpoint > 75],
                      sd = age_trend_data$scale_sd[age_trend_data$age_midpoint > 75])
  
  # Element-wise multiplication of raw draws with scales for 75+
  draws_75up <- data.frame(age_class = "75+", 
                             age_midpoint = 85,
                             seropos = sero_55up*scales_75up,
                             sim = 1:n_post,
                             pop = age_popdata$total[age_popdata$age_midpoint > 75])
  
  # Putting all of the age class draws together
  prev_est_age_draws <- prev_est_age_draws %>%
    rbind(draws_0to17) %>%
    rbind(draws_18to44) %>%
    rbind(draws_45to64) %>%
    rbind(draws_65to74) %>%
    rbind(draws_75up)
  
} else {

  # No 0-17 age class - this contains  dummy draws (sero = 0) so that the 
  # age classes line up for the stats later
  draws_0to17 <- data.frame(age_class = "0-17", 
                           age_midpoint =9,
                           seropos = rep(0, n_post),
                           sim = 1:n_post,
                           pop = age_popdata$total[age_popdata$age_midpoint < 18])
  
  
  # raw draws for 18-44 age class
  draws_18to44 <- data.frame(age_class = "18-44", 
                             age_midpoint = 31,
                             seropos = rnorm(n_post, 
                                             mean = prev_est_age$seropos_mu[prev_est_age$age_midpoint < 45],
                                             sd = prev_est_age$seropos_sig[prev_est_age$age_midpoint < 45]),
                             sim = 1:n_post,
                             pop = age_popdata$total[age_popdata$age_midpoint > 18 & age_popdata$age_midpoint < 45])
  
  # raw serodraws for 45-54 class
  sero_45to54 <- rnorm(n_post, 
                       mean = prev_est_age$seropos_mu[prev_est_age$age_midpoint > 45 & prev_est_age$age_midpoint < 55],
                       sd = prev_est_age$seropos_sig[prev_est_age$age_midpoint > 45 & prev_est_age$age_midpoint < 55])
  
  
  # raw draws for 55-64 class (based on 55+)
  sero_55to64 <- rnorm(n_post, 
                     mean = prev_est_age$seropos_mu[prev_est_age$age_midpoint > 55],
                     sd = prev_est_age$seropos_sig[prev_est_age$age_midpoint > 55])
  
  
  # population-weighted 45-64 class
  n_post_45to54 <- as.integer(prev_est_age$pop[prev_est_age$age_midpoint > 45 & prev_est_age$age_midpoint < 55]/
                                (age_popdata$total[age_popdata$age_midpoint > 45 & age_popdata$age_midpoint < 65])*n_post)
  n_post_55to64 <- n_post - n_post_45to54
  pop_weighted_sero_45to65 <- sample(sero_45to54, n_post_45to54) %>%
    append(sample(sero_55to64, n_post_55to64))
  
  # Formatting as draws for 45-64 class
  draws_45to64 <- data.frame(age_class = "45-64", 
                             age_midpoint = 55,
                             seropos = pop_weighted_sero_45to65,
                             sim = 1:n_post,
                             pop = age_popdata$total[age_popdata$age_midpoint > 45 & age_popdata$age_midpoint < 65])
  
  
  # raw draws for 65-74 class (based on 55+)
  draws_65to74 <- data.frame(age_class = "65-74", 
                             age_midpoint = 70,
                             seropos = rnorm(n_post, 
                                             mean = prev_est_age$seropos_mu[prev_est_age$age_midpoint > 55],
                                             sd = prev_est_age$seropos_sig[prev_est_age$age_midpoint > 55]),
                             sim = 1:n_post,
                             pop = age_popdata$total[age_popdata$age_midpoint > 55])
  
  # raw draws for 75+ class (based on 55+)
  draws_75up <- data.frame(age_class = "75+", 
                             age_midpoint = 85,
                             seropos = rnorm(n_post, 
                                             mean = prev_est_age$seropos_mu[prev_est_age$age_midpoint > 55],
                                             sd = prev_est_age$seropos_sig[prev_est_age$age_midpoint > 55]),
                             sim = 1:n_post,
                             pop = age_popdata$total[age_popdata$age_midpoint > 55])
  
  # Puts all of the draws together
  prev_est_age_draws <- prev_est_age_draws %>%
    rbind(draws_0to17) %>%
    rbind(draws_18to44) %>%
    rbind(draws_45to64) %>%
    rbind(draws_65to74) %>%
    rbind(draws_75up)
  
}

# Delay distributions ----------------------------------------------------------

n_sero_weeks <- 1  # number of serosurvey weeks

# Filled in with mean and 95% CI from Rosenburg et al.
serodata <- tribble(
  ~date, ~median, ~q025, ~q975,
  "2020-04-23", 22.7, 21.5, 24.0
) %>% 
  mutate(date = as.Date(date))

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
nwarmup <- 4000
niter <- nwarmup + 1000
n_chains <- 5

# IFR Computation
IFRs <- foreach(agec = age_classes[age_classes != "all"],
                .combine = rbind) %do% 
  { res_file <- paste0("age_stratified_IFRs_distsample_", agec, suffix, ".rds") 
  stanfit_file <- paste0("stanfit_age_stratified_IFRs_distsample_", agec, suffix, ".rds") 
  
  if (!file.exists(res_file) | redo_fit) {
    
    # Compile stan model
    ifr_stan <- stan_model("ifr_betabinomial.stan")
    
    epidata <- filter(age_epidata, age_class == agec)
    pop_age <- age_popdata$total[age_popdata$age_class == agec]  # population of the age class
    
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

# Export -----------------------------------------------------------------------

# Saves outputs to R Data Structure for postprocess
save(list = ls(), file = "data_setup.rda")
saveRDS(IFRs, file = paste0("age_stratified_IFRs_", suffix, ".rds"))

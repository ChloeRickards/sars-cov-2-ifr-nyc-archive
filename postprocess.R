# Notes ------------------------------------------------------------------------
# This code should be run *after* "draw_ifrs.R"
# This produces all of the graphs and figures, plus a few extra

# Preamble ---------------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(gridExtra)
source("utils.R")

# Setup ------------------------------------------------------------------------
# Loading outputs from draw_ifrs
load("data_setup.rda")

# Setting up some conditions
if (adjusted & Pollan) {
  suffix <- "adj"
} else if (adjusted & !Pollan) {
  suffix <- "adj_ward"
} else {
  suffix <- "raw"
}
age_classes <- c("0-17", "18-44", "45-64","65-74","75+", "all")

# Load and format IFRs
IFRs <- readRDS(paste0("age_stratified_IFRs_", suffix, ".rds"))

IFRs <- IFRs %>%
  group_by(age_class) %>% 
  mutate(sim = row_number()) %>% 
  ungroup() %>% 
  inner_join(prev_est_age_draws,
             by = c("sim", "age_class" = "age_class")) %>% 
  group_by(sim) %>% 
  summarise(ifr = weighted.mean(ifr, seropos * pop)) %>% 
  mutate(age_class = "all") %>% 
  select(ifr, age_class) %>% 
  rbind(IFRs) 

# Compute Bayesian p-values ----------------------------------------------------

ref_age <- "45-64"    # reference age class -- this used to be [20, 50)
ref_ifrs <- IFRs$ifr[IFRs$age_class == ref_age] # reference IFR posterior draws

IFR_pvals <- IFRs %>% 
  group_by(age_class) %>% 
  summarise(pval = computePval(ifrs, ref_ifrs)) %>% 
  mutate(pval = case_when(age_class == ref_age ~ as.double(NA), T ~ pval),
         pval_string = case_when(pval == 0 ~ "<0.001",
                                 is.na(pval) ~ "-",
                                 T ~ format(pval, digits = 2)))

IFR_stats <- 
  # Epi data
  age_epidata %>%
  group_by(age_class) %>% 
  arrange(desc(case_cumul)) %>% 
  select(case_cumul, death_cumul) %>% 
  slice(1) %>% ungroup() %>% 
  # Populations by age
  inner_join(age_popdata %>% select(age_class, total)) %>% 
  # Seroprevalence estimates
  inner_join(prev_est_age_draws %>% 
               group_by(age_class) %>% 
               summarise(sero.mean = mean(seropos),
                         sero.025 = quantile(seropos, 0.025),
                         sero.975 = quantile(seropos, 0.975))) %>% 
  mutate(seropop.mean = sero.mean * total,
         seropop.025 = sero.025 * total,
         seropop.975 = sero.975 * total) %>% 
  # IFR estimates
  inner_join(
    IFRs %>% 
      group_by(age_class) %>% 
      mutate(ifr = ifr*1e2) %>% 
      summarise(mean = mean(ifr),
                q025 = quantile(ifr, 0.025),
                q975 = quantile(ifr, 0.975))
  ) %>% 
  inner_join(IFR_pvals) %>% 
  group_by(age_class) %>% 
  mutate_at(vars(mean, q025, q975), 
            function(x) format(x, digits = 2)) %>% 
  mutate_at(vars(contains("seropop")), 
            function(x) 100*round(x/100)) %>% 
  ungroup() %>% 
  mutate(ifr = paste0(mean, " (", q025, "-", q975, ")"),
         seropop = paste0(format(seropop.mean, big.mark = ","),
                          " (", format(seropop.025, big.mark = ","), 
                          "-", 
                          format(seropop.975, big.mark = ","), ")"),
         age_class = factor(age_class, 
                            levels = age_classes[c(1, 2, 3, 4:length(age_classes))])) %>%
  arrange(age_class)

IFR_stats$age_class <- as.character(IFR_stats$age_class)
IFR_stats[1,1] <- "0-17 *"
IFR_stats[2,1] <- "18-44"
IFR_stats[3,1] <- "45-64"
IFR_stats[4,1] <- "65-74"
IFR_stats[5,1] <- "75+"

IFR_pvals <- IFR_stats %>% 
  select(age_class, seropop, death_cumul, ifr, pval_string) %>% 
  rename(`Age class` = age_class,
         `Estimated infected` = seropop,
         `Deaths` = death_cumul,
         IFR = ifr,
         `p-value` = pval_string) 

knitr::kable(IFR_pvals,
             format = "latex",
             booktabs = T) %>% 
  write(file = paste0("IFRs_pvals_", suffix, ".tex"))

# Epidata plot  ----------------------------------------------------------------
Sys.setlocale("LC_TIME", "C")


# - - - - 
# Plot epi data

var_dict <- c(
  "case_cumul" = "Cumulative incidence",
  "death_cumul" = "Cumulative deaths"
)

stratified_data <- read_csv("NYCHealth_cumul_data.csv")

p_data <- ggplot(stratified_data, aes(x = date, y = Total, color = age_class)) +
  geom_point(size = .6) +
  geom_line() +
  theme_bw() +
  facet_wrap(~var, scales = "free") +
  labs(x = "", y = "") +
  scale_color_discrete(labels= c("0-17", "18-44", "45-64", "65-74", "75+"))

ggsave(p_data, filename = "epidata.png", width = 8, height = 7)

# Delay distribution plot ------------------------------------------------------

# Delays for plotting
nmax <- 50  # max number of days in PDFS
dt <- .1
ntot <- nmax/dt
n_dist_samples <- 1000

# Initialize PDF matrices
pdf_inc <- matrix(rep(0, ntot * n_dist_samples), nrow = n_dist_samples)
pdf_report <- pdf_inc
pdf_case <- pdf_inc
pdf_symp_sero <- pdf_inc
pdf_sero <- pdf_inc
pdf_report_death <- pdf_inc
pdf_symp_death <- pdf_inc
pdf_death <- pdf_inc

# Sample distributions
for (i in 1:n_dist_samples) {
  pdf_inc[i, ] <- getDelay(delay_params, "inc", nmax, dt = dt, rnd = T)*dt
  pdf_report[i, ] <- getDelay(delay_params, "report", nmax, dt = dt, rnd = T)*dt
  pdf_case[i, ] <- convolve( pdf_inc[i, ], rev(pdf_report[i, ]), type = "o")[1:ntot]
  pdf_symp_sero[i, ] <- getDelay(delay_params, "symp_sero", nmax, dt = dt, rnd = T)*dt
  pdf_sero[i, ] <- convolve(pdf_inc[i, ], rev(pdf_symp_sero[i, ]), type = "o")[1:ntot]  # convolve with incubation period
  pdf_report_death[i, ] <- getDelay(delay_params, "report_death", nmax, dt = dt, rnd = T)*dt
  pdf_symp_death[i, ] <- convolve(pdf_report[i, ], rev(pdf_report_death[i, ]), type = "o")[1:ntot]  # convolve with delay from symptoms to hospitalization
  pdf_death[i, ] <- convolve(pdf_inc[i, ], rev(pdf_symp_death[i, ]), type = "o")[1:ntot]  # convolve with incubation period
}

delay_dict <- c(
  "i2s" = "Infection to symptom onset",
  "s2c" = "Onset to reporting",
  "s2sero" = "Onset to seroconversion",
  "c2d" = "Reporting to death",
  "i2c" = "Infection to reporting",
  "i2sero" = "Infection to seroconversion",
  "i2d" = "Infection to death"
)

legend_title <- "Delay distribution"

delays_2plot <- map2_dfr(
  list(pdf_inc, pdf_report, pdf_symp_sero, pdf_report_death, pdf_case, pdf_sero, pdf_death),
  c("i2s", "s2c", "s2sero", "c2d", "i2c", "i2sero", "i2d"), 
  ~getQuantiles(.x) %>%
    mutate(days = seq(dt/2, nmax, by = dt),
           var = .y)) %>% 
  mutate(type = case_when(
    var %in% c("i2s", "s2c", "s2sero", "c2d") ~ "literature",
    T ~ "derived"),
    var_name = delay_dict[var]) 

p_lit <- delays_2plot %>% 
  filter(type == "literature") %>% 
  ggplot(aes(x = days, y = cdf.median, ymin = cdf.q025, ymax = cdf.q975, fill = var_name, lty = var_name)) + 
  geom_ribbon(alpha = .25) +
  geom_line(aes(color = var_name)) +
  guides(color = guide_legend(legend_title),
         fill = guide_legend(legend_title),
         lty = guide_legend(legend_title)) +
  theme_bw() +
  labs(x = "time [days]", y = "CDF") +
  ggthemes::scale_color_few() +
  ggthemes::scale_fill_few() +
  theme(legend.position = c(.71, .2))

p_conv <- delays_2plot %>% 
  filter(type == "derived") %>% 
  ggplot(aes(x = days, y = cdf.median, ymin = cdf.q025, ymax = cdf.q975, fill = var_name, lty = var_name)) + 
  geom_ribbon(alpha = .25) +
  geom_line(aes(color = var_name)) +
  guides(color = guide_legend(legend_title),
         fill = guide_legend(legend_title),
         lty = guide_legend(legend_title)) +
  theme_bw() +
  labs(x = "time [days]", y = "CDF") +
  theme(legend.position = c(.71, .2))

p_comb <- arrangeGrob(grobs = list(p_lit, p_conv), nrow = 1)
ggsave(plot = p_comb, filename = "delay_dist.png", width = 9, height = 4.5)


# Timeline for daily cases and death -------------------------------------------
# Based on confirmed daily cases
# Cases and deaths truncated to dates of interest: 2020-03-01 to 2020-06-15
# Source: NYC Dept. of Health

# Formatting daily case data
daily_cases_all_ages <-read_csv("daily_cases_trunc.csv") %>%
  select(DATE_OF_INTEREST, '7-day average')%>%
  cbind(var = "Cases (7-day average)") %>%
  rename('Date' = DATE_OF_INTEREST) %>%
  rename('DailyCount' = '7-day average') %>%
  mutate(Date = as.Date(Date))

# Formatting daily death data
daily_deaths_all_ages <-read_csv("daily_deaths_trunc.csv") %>%
  select(DATE_OF_DEATH, Confirmed) %>%
  cbind(var = "Deaths") %>%
  rename('Date' = DATE_OF_DEATH) %>%
  rename('DailyCount' = 'Confirmed') %>%
  mutate(Date = as.Date(Date))

timeline <- daily_deaths_all_ages %>%
  rbind(daily_cases_all_ages)

# Placing the serosurvey dates on the timeline
serosurvey_start <- as.Date("2020-04-19")
serosurvey <- data.frame(Date = serosurvey_start + 0:9, DailyCount = rep(1500, 10), var = "Serosurvey Dates (04/19/20 - 04/28/20")

# Calculating the range of days when infection likely occurred
inf2sero <- delays_2plot %>%
  filter(var == "i2sero") %>%
  filter(cdf.mean > 0.05) 

# Identifies start, midpoint, and end date for when infection likely occurred
infection_start = round(inf2sero$days[1])
infection_midpoint <- inf2sero %>%
  filter(cdf.mean < 0.51) %>%
  filter(cdf.mean > 0.49)
infection_midpoint = round(mean(infection_midpoint$days))
infection_end = round(inf2sero$days[length(inf2sero$days)])

# Calculating the range of days when deaths (from the likely infections) likely occurred
inf2death <- delays_2plot %>%
  filter(var == "i2d") %>%
  filter(cdf.mean < 0.95)

# Identifies start, midpoint, and end date for when deaths (from the likely infections) likely occurred
death_start = round(inf2death$days[1])
death_midpoint <- inf2death %>%
  filter(cdf.mean < 0.51) %>%
  filter(cdf.mean > 0.49)
death_midpoint = round(mean(death_midpoint$days))
death_end = round(inf2death$days[length(inf2death$days)])

# Reformats the dates
serosurvey_midpoint <- as.Date("2020-04-23")
infection <- data.frame(Date = serosurvey_midpoint + -infection_start:-(infection_end+3), DailyCount = rep(1400, length(infection_start:(infection_end+3))), var = "Estimated Dates of Infection (03/01/2020 - 04/14/20)")
infection_midpoint_day <- data.frame(Date = serosurvey_midpoint + -infection_midpoint, DailyCount = 1400, var = "Estimated Infection Time Period Midpoint")
deaths <- data.frame(Date = infection_midpoint_day$Date + -24:death_end, DailyCount = rep(1600, length(-24:death_end)), var = "Estimated Dates of Death (03/14/20 - 05/26/20)")
deaths_midpoint_day <- data.frame(Date = infection_midpoint_day$Date + death_midpoint, DailyCount = 1600, var = "Probable Death Time Period Midpoint")

# Combines all of the relevant dates
timeline <- daily_deaths_all_ages %>%
  rbind(daily_cases_all_ages) %>%
  rbind(serosurvey) %>%
  rbind(infection) %>%
  #rbind(deaths) %>%
  rename(Legend = var)

#https://stackoverflow.com/questions/17996410/ggplot-specific-thick-line

# Plot
ggplot(timeline, aes(x=Date,y = DailyCount)) +
  geom_line(aes(color = Legend, size = Legend))+
  theme_bw()+
  theme(axis.title=element_text(size=25),legend.position = c(.7, .8),
        axis.text=element_text(size=25),legend.text=element_text(size=18),
        legend.title=element_text(size=20))+
  scale_size_manual(values = c(1, 1, 3, 3, 3))+
  scale_color_manual(values=c("#61ddff", "#2051e3", "#de2509", "#edc01c"))+
  scale_x_date(date_breaks = "1 month", date_labels = "%b")+
  geom_point(aes(x = serosurvey_midpoint, y = 1500), color = "black", size = 5)+
  xlab("Date")+
  ylab("Case/Death Count")

ggsave("nyc_timeline.png", width = 15, height = 12)

# IFR posterior draws plot -----------------------------------------------------

IFRs <- IFRs %>% 
  mutate(age_class = factor(age_class, 
                            levels = age_classes[c(1, 2, 3, 4:length(age_classes))]))

age.labs <- c("0-17", "18-44", "45-64", "65-74", "75+")
names(age.labs) <- c("0-17", "18-44", "45-64", "65-74", "75+")

p_ifrs <- IFRs %>%
  ggplot(aes(x = ifr*100)) +
  geom_histogram() +
  facet_wrap(~age_class, scales = "free", labeller = labeller(age_class = age.labs)) +
  xlab("IFR [%]") +
  theme_bw()

ggsave(p_ifrs, filename = paste0("ifr_posteriors_", suffix, ".png"), width = 8, height = 4)


# IFR & Sero Table Setup -------------------------------------------------------
# Saves R data to a file, to be used later in "additional_postprocess.R"
# These are also the finalized versions of the IFR and the seroprevalence
# If you use different age classes, you DO need to modify this area

# IFR setup
IFR_stats <- IFR_stats %>% 
  select(age_class, total, seropop, death_cumul, ifr) %>% 
  rename(`Age class` = age_class,
         Population = total,
         `Estimated infected` = seropop,
         `Deaths` = death_cumul,
         IFR = ifr)

# In the case of raw numbers, we don't have a 0-17 estimate, so the value is 
# set to 'NA'
if (!adjusted | !Pollan){
  IFR_stats$`Estimated infected`[IFR_stats$`Age class` == '0-17 *'] <- NA
  IFR_stats$IFR[IFR_stats$`Age class` == '0-17 *'] <- NA
}

# saving as a data structure to be used in "additional_postprocess"
saveRDS(IFR_stats, file = paste0("IFR_stats_", suffix, ".rds"))

# Seroprevalence setup
draws_means_and_sds <- data.frame(age_class = "0-17", 
                           age_midpoint = 9,
                           seropos_mu = round(mean(draws_0to17$seropos),3),
                           seropos_sd = round(sd(draws_0to17$seropos),3)) %>%
  rbind(data.frame(age_class = "18-44", 
                   age_midpoint = 31,
                   seropos_mu = round(mean(draws_18to44$seropos),3),
                   seropos_sd = round(sd(draws_18to44$seropos),3))) %>%
  rbind(data.frame(age_class = "45-64", 
                   age_midpoint = 55,
                   seropos_mu = round(mean(draws_45to64$seropos),3),
                   seropos_sd = round(sd(draws_45to64$seropos),3))) %>%
  rbind(data.frame(age_class = "65-74", 
                   age_midpoint = 70,
                   seropos_mu = round(mean(draws_65to74$seropos),3),
                   seropos_sd = round(sd(draws_45to64$seropos),3))) %>%
  rbind(data.frame(age_class = "75+", 
                   age_midpoint = 85,
                   seropos_mu = round(mean(draws_75up$seropos),3),
                   seropos_sd = round(sd(draws_75up$seropos),3))) 

# Creates seroprevalence and 95% CI
seroprevalence <- rep(0,5)
for (i in 1:length(seroprevalence)){
  seroprevalence[i] <- paste0(draws_means_and_sds[i,]$seropos_mu*100, " (", (draws_means_and_sds[i,]$seropos_mu - 2*draws_means_and_sds[i,]$seropos_sd)*100, "-", (draws_means_and_sds[i,]$seropos_mu + 2*draws_means_and_sds[i,]$seropos_sd)*100, ")")
}

# saving as a data structure to be used in "additional_postprocess"
saveRDS(seroprevalence, file = paste0("seroprevalence_", suffix, ".rds"))


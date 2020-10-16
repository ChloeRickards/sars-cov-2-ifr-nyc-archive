# Preamble
library(tidyverse)
library(magrittr)
library(ggplot2)
library(lubridate)

# FINAL IFR/SEROPREV TABLE

IFR_stats_adj <- readRDS("IFR_stats_adj.RDS")
IFR_stats_raw <- readRDS("IFR_stats_raw.RDS")
prev_est_age_adj <- readRDS("prev_est_age_adj.RDS")
prev_est_age_raw <- readRDS("prev_est_age_raw.RDS")

# IFR_stats_adj <- readRDS("IFR_stats_adj_allDeaths.RDS")
# IFR_stats_raw <- readRDS("IFR_stats_raw_allDeaths.RDS")
# prev_est_age_adj <- readRDS("prev_est_age_adj_allDeaths.RDS")
# prev_est_age_raw <- readRDS("prev_est_age_raw_allDeaths.RDS")

resultsTable <- 
  cbind('Age Class' = IFR_stats_adj$`Age class`) %>%
  cbind('Population' = IFR_stats_adj$Population) %>%
  cbind('Deaths' = IFR_stats_adj$Deaths) %>%
  cbind('Seroprevalence (raw)' = c("NA", prev_est_age_raw$`% Seropositive (95% CI)`)) %>%
  cbind('IFR (raw)' = c("NA", IFR_stats_raw$IFR)) %>%
  cbind('Seroprevalence (adjusted with Pollan)' = prev_est_age_adj$`% Seropositive (95% CI)`) %>%
  cbind('IFR (adjusted with Pollan)' = IFR_stats_adj$IFR)

# IFR COMPARISON

#Qs for Marm:
# why remove Italy?

IFR_comparison <- read_csv("data/IFR_comparison.csv") 
IFR_comparison=IFR_comparison[IFR_comparison$Dataset!="Italy",]

require(ggplot2);require(ggthemes)

IFR_comparison[,3:5]=100*IFR_comparison[,3:5]#converting to %

# Midpoints for our age classes: 0-18, 18-45, 45-65, 65-75, 75+
NYC_midpoints = c(9, 32, 55, 70, 85)

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

NYC_midpoints = c(9, 32, 55, 70, 85)

IFR_stats_raw_comp <- IFR_stats_raw

IFR_stats_raw_comp$IFR <- gsub("[()]",",",IFR_stats_raw_comp$IFR)
IFR_stats_raw_comp$IFR <- gsub("-",",",IFR_stats_raw_comp$IFR)
IFR_stats_raw_comp$IFR <- sub(" ", "", IFR_stats_raw_comp$IFR)

split <- strsplit(IFR_stats_raw_comp$IFR, ",")
mean_raw <- vector()
lowerCI_raw <- vector()
upperCI_raw <- vector()
for (i in 1:4) {
  vec <- split[[i]]
  mean_raw[i] <- as.numeric(vec[1])
  lowerCI_raw[i] <- as.numeric(vec[2])
  upperCI_raw[i] <- as.numeric(vec[3])
}

IFR_comparison <- IFR_comparison %>%
  rbind(data.frame(AgeMidpoint = NYC_midpoints[-1], Dataset = 'Rickards (raw)', IFR = mean_raw, LowerCI = lowerCI_raw, UpperCI = upperCI_raw)) %>%
  rbind(data.frame(AgeMidpoint = NYC_midpoints, Dataset = 'Rickards (adjusted with Pollan)', IFR = mean_adj, LowerCI = lowerCI_adj, UpperCI = upperCI_adj))

options(scipen=10000)
IFR_comparison$Dataset <- factor(IFR_comparison$Dataset, levels = c("Rickards (raw)", "Rickards (adjusted with Pollan)", "Pastor*", "Pastor* All cause", "Perez-Saez*", "Russell", "Salje", "Verity", "Ward*", "Ward* All cause"))

ggplot(IFR_comparison, aes(y=IFR, x=AgeMidpoint,col=Dataset)) +
  geom_line(aes(group = Dataset),position=position_dodge(width=2))+theme_bw()+
  geom_point(size=4,position=position_dodge(width=2))+
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI),
                width = 0.25,position=position_dodge(width=2))+
  coord_cartesian(xlim=c(1, 86), ylim = c(.9e-3,1.8e1))+scale_y_log10()+
  theme(axis.title=element_text(size=25),legend.position = c(.2, .80),
        axis.text=element_text(size=25),legend.text=element_text(size=15),
        legend.title=element_text(size=15))+
  scale_color_manual(values=c("#ff1900", "#c92a0a", "#E3D400", "#b1de10", "#24b379", "36b6e0", "#2888bf", "#2d5ba1", "#7679d6", "#ae91eb"))+
  xlab("Age Group midpoint")+
  ylab("IFR log scale (in percent)")

ggsave("results/fig/ifr_comparison.png", width = 12, height = 15)

#https://stackoverflow.com/questions/43577579/multiple-lines-multiple-error-bars-using-ggplot2-in-r

# DAILY CASES AND DEATHS
# Daily deaths are just confirmed here
# May do combined when I do the IFR stats for combined
# I chose to use confirmed deaths here because I'm doing IFR stats based on confirmed

delays_2plot <- readRDS("delays_2plot.RDS")

daily_cases_all_ages <-read_csv("data/daily cases trunc.csv") %>%
  select(DATE_OF_INTEREST, '7-day average')%>%
  cbind(age_class = "all") %>%
  cbind(var = "Cases (7-day average)")
daily_cases_all_ages$DATE_OF_INTEREST <-as.Date(daily_cases_all_ages$DATE_OF_INTEREST, "%m/%d/%y")

daily_deaths_all_ages <-read_csv("data/daily deaths trunc.csv") %>%
  select(DATE_OF_DEATH, Confirmed) %>%
  cbind(age_class = "all") %>%
  cbind(var = "Deaths")
daily_deaths_all_ages$DATE_OF_DEATH <-as.Date(daily_deaths_all_ages$DATE_OF_DEATH, "%m/%d/%y")

daily_cases_by_age <- read_csv("data/NYCHealth_cumul_data.csv") %>%
  filter(var == 'case_cumul') %>%
  arrange(age_class, date)
daily_vec = vector()
# for each age class, take the difference of the two days
for (i in 1:length(unique(daily_cases_by_age$age_class))){
  c <- daily_cases_by_age %>%
    filter(age_class == unique(daily_cases_by_age$age_class)[i])
  for (j in 1:dim(c)[1]){
    if (j == 1){
      daily_vec <- append(daily_vec, 0)
    } else{
      daily_vec <- append(daily_vec,c$Total[j] - c$Total[j-1])
    }
  }
}
daily_cases_by_age <- daily_cases_by_age %>%
  cbind(Daily = daily_vec) %>%
  select(-Total) %>%
  select(-var)
daily_cases_by_age <- filter(daily_cases_by_age, daily_cases_by_age$date != "2020-03-23")

daily_deaths_by_age <- read_csv("data/NYCHealth_cumul_data.csv") %>%
  filter(var == 'death_cumul') %>%
  arrange(age_class, date)
daily_vec = vector()
# for each age class, take the difference of the two days
for (i in 1:length(unique(daily_deaths_by_age$age_class))){
  c <- daily_deaths_by_age %>%
    filter(age_class == unique(daily_deaths_by_age$age_class)[i])
  for (j in 1:dim(c)[1]){
    if (j == 1){
      daily_vec <- append(daily_vec, 0)
    } else{
      daily_vec <- append(daily_vec,c$Total[j] - c$Total[j-1])
    }
  }
}
daily_deaths_by_age <- daily_deaths_by_age %>%
  cbind(Daily = daily_vec) %>%
  select(-Total) %>%
  select(-var)
daily_deaths_by_age <- filter(daily_deaths_by_age, daily_deaths_by_age$date != "2020-03-23")
for (i in 1:dim(daily_deaths_by_age)[1]){
  if (daily_deaths_by_age$Daily[i] <0){
    daily_deaths_by_age$Daily[i] = 0
  }
}

names(daily_deaths_all_ages)[names(daily_deaths_all_ages) == "DATE_OF_DEATH"] <- "Date"
names(daily_deaths_all_ages)[names(daily_deaths_all_ages) == "Confirmed"] <- "DailyCount"
names(daily_cases_all_ages)[names(daily_cases_all_ages) == "DATE_OF_INTEREST"] <- "Date"
names(daily_cases_all_ages)[names(daily_cases_all_ages) == "7-day average"] <- "DailyCount"

timeline <- daily_deaths_all_ages %>%
  rbind(daily_cases_all_ages) %>%
  select(-age_class)

serosurvey_start <- as.Date("2020-04-19")
serosurvey <- data.frame(Date = serosurvey_start + 0:9, DailyCount = rep(1500, 10), var = "Serosurvey Dates (04/19/20 - 04/28/20")

# days from infection to seroconversion
# What I'm going to do: 1) Take middle 95% of cdf, 2) convert to days, 3) convert to dates, using days before 4/23 (midpoint)

inf2sero <- delays_2plot %>%
  filter(var == "i2sero") %>%
  filter(cdf.mean > 0.025) %>%
  filter(cdf.mean < 0.975)

infection_start = round(inf2sero$days[1])

infection_midpoint <- inf2sero %>%
  filter(cdf.mean < 0.51) %>%
  filter(cdf.mean > 0.49)
infection_midpoint = round(mean(infection_midpoint$days))

infection_end = round(inf2sero$days[length(inf2sero$days)])

# taking deaths from midpoint of infections
inf2death <- delays_2plot %>%
  filter(var == "i2d") %>%
  filter(cdf.mean > 0.025) %>%
  filter(cdf.mean < 0.975)

death_start = round(inf2death$days[1])
death_midpoint <- inf2death %>%
  filter(cdf.mean < 0.51) %>%
  filter(cdf.mean > 0.49)
death_midpoint = round(mean(death_midpoint$days))
death_end = round(inf2death$days[length(inf2death$days)])

serosurvey_midpoint <- as.Date("2020-04-23")
infection <- data.frame(Date = serosurvey_midpoint + -infection_start:-infection_end, DailyCount = rep(1400, length(infection_start:infection_end)), var = "Probable Dates of Infection (03/16/2020 - 04/16/20)")
infection_midpoint_day <- data.frame(Date = serosurvey_midpoint + -infection_midpoint, DailyCount = 1400, var = "Probable Infection Time Period Midpoint")
deaths <- data.frame(Date = infection_midpoint_day$Date + death_start:death_end, DailyCount = rep(1600, length(death_start:death_end)), var = "Probable Dates of Death (04/15/20 - 05/27/20)")
deaths_midpoint_day <- data.frame(Date = infection_midpoint_day$Date + death_midpoint, DailyCount = 1600, var = "Probable Death Time Period Midpoint")

timeline <- timeline %>%
  rbind(serosurvey) %>%
  rbind(infection) %>%
  rbind(deaths) %>%
  rename(Variable = var)

#https://stackoverflow.com/questions/17996410/ggplot-specific-thick-line

ggplot(timeline, aes(x=Date,y = DailyCount)) +
  geom_line(aes(color = Variable, size = Variable))+
  theme_bw()+
  theme(axis.title=element_text(size=25),legend.position = c(.7, .8),
        axis.text=element_text(size=25),legend.text=element_text(size=18),
        legend.title=element_text(size=20))+
  scale_size_manual(values = c(1, 1, 3, 3, 3))+
  scale_color_manual(values=c("#61ddff", "#2051e3", "#de2509", "#9718c9", "#edc01c"))+
  scale_x_date(date_breaks = "1 month", date_labels = "%b")+
  xlab("Date")+
  ylab("Case/Death Count")
  
ggsave("results/fig/nyc_timeline.png", width = 15, height = 12)


# SUPPLEMENTAL FIG - AGE TREND DATA

lapply(c("ggplot2","ggthemes","dplyr"),require,character.only=T) #load multiple packages
m=read.csv("Seroprev_by_age.csv")
m=as.data.frame(group_by(m,Dataset) %>%
                  mutate(SeroprevN=Seroprev/mean(Seroprev)))

ggplot(m, aes(y=SeroprevN, x=AgeMidpoint,col=Dataset)) +
  geom_line(aes(group = Dataset))+theme_bw()+
  geom_point(size=4)+
  theme(axis.title=element_text(size=25),legend.position = c(.6, .20),
        axis.text=element_text(size=25),legend.text=element_text(size=15),
        legend.title=element_text(size=15))+
  xlab("Age Group midpoint")
ylab("Relative Seroprevalence")

ggsave(m, filename = "results/fig/ifr_comparison.png", width = 8, height = 4)

# setwd("c:/marm/research/covid-19/IFR from serology/")
# m=read.csv("age-sex-week-est.csv")
# require(ggplot2);require(ggthemes)
# #m[,3:5]=100*m[,3:5]#converting to %
# #options(scipen=10000)
# mg=aggregate(seropos~week+age_cat,data=m,FUN=mean)
# ggplot(mg, aes(y=seropos, x=week,col=age_cat)) + 
#   geom_line(aes(group = age_cat))+#theme_few()+
#   geom_point(size=1)+
#   #  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), 
#   #                width = 0.25,position=position_dodge(width=2))+
#   #  scale_y_continuous(trans='log10',limits = c(.9e-3,1.6e1))+
#   #  coord_cartesian(xlim=c(10, 86), ylim = c(.9e-3,1.6e1))+scale_y_log10()+
#   theme(axis.title=element_text(size=25),legend.position = c(.2, .80),
#         axis.text=element_text(size=25),legend.text=element_text(size=15),
#         legend.title=element_text(size=15))+
#   xlab("Week")+
#   ylab("Seroprev")


# Age fit
age_trend_data =read.csv("data/serop_age_spain.csv") # this reads the file
age_trend_data = age_trend_data[age_trend_data$AgeMidpoint>45,] #this selects for age groups over 45
age_sero_fit=lm(IA~AgeMidpoint+I(AgeMidpoint^2),data=age_trend_data);summary(age_sero_fit) # fits a line to the data using a linear model; summary of fit
age_trend_data$fitted=fitted(age_sero_fit) # either here or on very end - adds a column to age_trend_data of fitted IA values - basically I'll tack on a fitted column onto the serodata and use that
age_trend_data$scale=age_trend_data$fitted/mean(age_trend_data$IA[age_trend_data$AgeMidpoint<70]) #either here or on very end -- adds a column for scale, not sure what scale does tho
plot(IA~AgeMidpoint,data=age_trend_data)
lines(fitted(age_sero_fit)~age_trend_data$AgeMidpoint) # these last two lines to visualize

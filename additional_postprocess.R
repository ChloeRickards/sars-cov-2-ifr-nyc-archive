# Notes ------------------------------------------------------------------------
# This one requires that you have run "draw_ifrs.R" and "postprocess.R" on BOTH
# the "raw" and "adj"settings, since these analyses compare both settings


# Preamble ---------------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(ggplot2)
library(lubridate)

# Final IFR/seroprevalence table -----------------------------------------------
# This compares IFR/seroprevalence between raw and adjusted runs

IFR_stats_adj <- readRDS("IFR_stats_adj.RDS")
IFR_stats_raw <- readRDS("IFR_stats_raw.RDS")
seroprevalence_adj <- readRDS("seroprevalence_adj.RDS")
seroprevalence_raw <- readRDS("seroprevalence_raw.RDS")


resultsTable <- 
  cbind('Age Class' = IFR_stats_adj$`Age class`) %>%
  cbind('Population' = IFR_stats_adj$Population) %>%
  cbind('Deaths' = IFR_stats_adj$Deaths) %>%
  cbind('Seroprevalence (raw)' =  seroprevalence_raw) %>%
  cbind('IFR (raw)' = IFR_stats_raw$IFR) %>%
  cbind('Seroprevalence (adjusted with Pollan)' = seroprevalence_adj) %>%
  cbind('IFR (adjusted with Pollan)' = IFR_stats_adj$IFR)


# IFR comparison ---------------------------------------------------------------
# This compares the raw and adjusted version of our results against multiple
# results from different countries. The CSV contains the results from other countries
IFR_comparison <- read_csv("IFR_comparison.csv") 
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

ggsave("ifr_comparison.png", width = 12, height = 15)

#https://stackoverflow.com/questions/43577579/multiple-lines-multiple-error-bars-using-ggplot2-in-r

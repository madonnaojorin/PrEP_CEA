############ Incidence of HIV 
############ Original: Dec 27 2020, Modified Nov 27 2021

rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
folder="Figure"

### Library 
library(tidyverse)
library(gridExtra)
theme_set(theme_classic(base_size = 20, base_family = "Helvetica"))

source(file = "Function_code/Incidence_function.R")

Initial = 1000
NumYears = 15
data_0 <- Incidence_function(0,Initial,NumYears)
data_30 <- Incidence_function(0.3,Initial,NumYears)
data_50 <- Incidence_function(0.5,Initial,NumYears)
data_70 <- Incidence_function(0.7,Initial,NumYears)
data_100 <- Incidence_function(1,Initial,NumYears)
##### Incidence PrEP 
PrEP_infected_data_0 <- data_0$PrEP_infected_data
PrEP_infected_data_30 <- data_30$PrEP_infected_data
PrEP_infected_data_50 <- data_50$PrEP_infected_data
PrEP_infected_data_70 <- data_70$PrEP_infected_data
PrEP_infected_data_100 <- data_100$PrEP_infected_data

data_PrEP_plot <-
  rbind(PrEP_infected_data_0,PrEP_infected_data_30,PrEP_infected_data_50,
        PrEP_infected_data_70,PrEP_infected_data_100)

##### Incidence No PrEP 
NoPrEP_infected_data_0 <- data_0$NoPrEP_infected_data
NoPrEP_infected_data_30 <- data_30$NoPrEP_infected_data
NoPrEP_infected_data_50 <- data_50$NoPrEP_infected_data
NoPrEP_infected_data_70 <- data_70$NoPrEP_infected_data
NoPrEP_infected_data_100 <- data_100$NoPrEP_infected_data

data_NoPrEP_plot <-
  rbind(NoPrEP_infected_data_0,NoPrEP_infected_data_30,NoPrEP_infected_data_50,
        NoPrEP_infected_data_70,NoPrEP_infected_data_100)

##### Simulated incidence No PrEP 
Sample <- sample(1:1000,100)
PrEP_simulated_data_0 <- data_0$PrEP_simulated_data
PrEP_simulated_data_0 <- PrEP_simulated_data_0[PrEP_simulated_data_0$Trial %in% Sample,]
PrEP_simulated_data_30 <- data_30$PrEP_simulated_data
PrEP_simulated_data_30 <- PrEP_simulated_data_30[PrEP_simulated_data_30$Trial %in% Sample,]
PrEP_simulated_data_50 <- data_50$PrEP_simulated_data
PrEP_simulated_data_50 <- PrEP_simulated_data_50[PrEP_simulated_data_50$Trial %in% Sample,]
PrEP_simulated_data_70 <- data_70$PrEP_simulated_data
PrEP_simulated_data_70 <- PrEP_simulated_data_70[PrEP_simulated_data_70$Trial %in% Sample,]
PrEP_simulated_data_100 <- data_100$PrEP_simulated_data
PrEP_simulated_data_100 <- PrEP_simulated_data_100[PrEP_simulated_data_100$Trial %in% Sample,]


data_PrEP_simulated <-
  rbind(PrEP_simulated_data_0,PrEP_simulated_data_30,PrEP_simulated_data_50,
        PrEP_simulated_data_70,PrEP_simulated_data_100)

##### Simulated incidence No PrEP 
NoPrEP_simulated_data_0 <- data_0$NoPrEP_simulated_data
NoPrEP_simulated_data_0 <- NoPrEP_simulated_data_0[NoPrEP_simulated_data_0$Trial %in% Sample,]
NoPrEP_simulated_data_30 <- data_30$NoPrEP_simulated_data
NoPrEP_simulated_data_30 <- NoPrEP_simulated_data_30[NoPrEP_simulated_data_30$Trial %in% Sample,]
NoPrEP_simulated_data_50 <- data_50$NoPrEP_simulated_data
NoPrEP_simulated_data_50 <- NoPrEP_simulated_data_50[NoPrEP_simulated_data_50$Trial %in% Sample,]
NoPrEP_simulated_data_70 <- data_70$NoPrEP_simulated_data
NoPrEP_simulated_data_70 <- NoPrEP_simulated_data_70[NoPrEP_simulated_data_70$Trial %in% Sample,]
NoPrEP_simulated_data_100 <- data_100$NoPrEP_simulated_data
NoPrEP_simulated_data_100 <- NoPrEP_simulated_data_100[NoPrEP_simulated_data_100$Trial %in% Sample,]

data_NoPrEP_simulated <-
  rbind(NoPrEP_simulated_data_0,NoPrEP_simulated_data_30,NoPrEP_simulated_data_50,
        NoPrEP_simulated_data_70,NoPrEP_simulated_data_100)


####### Plot for whole population #######
Mean <- data_PrEP_plot$Mean + data_NoPrEP_plot$Mean
High <- data_PrEP_plot$High + data_NoPrEP_plot$High
Low <- data_PrEP_plot$Low + data_NoPrEP_plot$Low

data_PrEP_plot$Coverage = c(rep("0% PrEP coverage",(NumYears+1)),rep("30% PrEP coverage",(NumYears+1)),rep("50% PrEP coverage",(NumYears+1)),rep("70% PrEP coverage",(NumYears+1)),rep("100% PrEP coverage",(NumYears+1)))
data_PrEP_plot$Coverage<- factor(data_PrEP_plot$Coverage,levels = c("0% PrEP coverage","30% PrEP coverage","50% PrEP coverage","70% PrEP coverage","100% PrEP coverage"))

data_all_1000 <-
  data.frame(
    Mean = Mean,
    High = High,
    Low = Low,
    Year = data_PrEP_plot$Year,
    Coverage = data_PrEP_plot$Coverage)

data_sim <- data_PrEP_simulated$data + data_NoPrEP_simulated$data

data_all_sim <-
  data.frame(
    data = data_sim,
    Year = data_NoPrEP_simulated$Year,
    Trial = data_NoPrEP_simulated$Trial,
    Coverage = data_NoPrEP_simulated$Coverage)

Whole <- ggplot(data_all_1000) +
  geom_line(
    data = data_all_sim[data_all_sim$Coverage=="00%",],
    aes(x = Year, y = data, group = Trial),
    alpha = 0.2,
    color = "#42B540",
    size = 0.1
  ) +
  geom_line(
    data = data_all_sim[data_all_sim$Coverage=="30%",],
    aes(x = Year, y = data, group = Trial),
    alpha = 0.2,
    color = "#0099B4"
  ) +
  geom_line(
    data = data_all_sim[data_all_sim$Coverage=="50%",],
    aes(x = Year, y = data, group = Trial),
    alpha = 0.2,
    color = "#00468B"
  ) +
  geom_line(
    data = data_all_sim[data_all_sim$Coverage=="70%",],
    aes(x = Year, y = data, group = Trial),
    alpha = 0.2,
    color = "#FDAF91"
  ) +
  geom_line(
    data = data_all_sim[data_all_sim$Coverage=="100%",],
    aes(x = Year, y = data, group = Trial),
    alpha = 0.2,
    color = "#AD172A"
  ) +
  geom_line(aes(x = Year, y = Mean, colour = Coverage), size = 2) +
  scale_colour_manual(values = c("#42B540", "#0099B4","#00468B","#FDAF91","#AD172A"))+
  scale_x_continuous(breaks = 1:15,limits = c(1,15)) +
  ylab("People living with HIV") +
  xlab("Year since the PrEP program initiation")+
  theme(legend.position = c(0.3, 0.8),
        legend.key.size = unit(1.5, "cm"),
        legend.title=element_blank()) +
  labs(fill = "PrEP coverage", color = "PrEP coverage")
Whole
#ggsave(filename = paste0(folder,"/Figure2.png"), Whole, dpi = 300,width = 200,height = 200,units = c("mm"))

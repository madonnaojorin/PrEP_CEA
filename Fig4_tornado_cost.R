############ Tornado plot for cost
############ Original: Dec 27 2020, Modified Nov 27 2021

rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
folder="Figure"
foldername = "data_result"
Nyear=15

### Library 
library(tidyverse)
library(magrittr)
library(reshape2)

source(file = "Function_code/DSA_cost.R")
DSA_cost(Nyear,foldername)
df <- read.csv(paste0(foldername,"/cost_results_",Nyear,".csv"),sep=",")
df <- df[,-1]

source(file ="Function_code/ICER_calculation.R")
ICER_base(Nyear,foldername)
icer <- read.csv(paste0(foldername,"/ICER_results_",Nyear,".csv"))
Base = icer[1,2]*10000

##### Define min and max
# DrugPrice
A <- 2
colnames(df[A])
(DrugPrice_Value_min <- round(min(as.numeric(df[,A]))*10000,digits=1))
(DrugPrice_Value_max <- round(max(as.numeric(df[,A]))*10000,digits=1))
(DrugPrice <- paste0("Cost of PrEP \n","(min:",DrugPrice_Value_min,"; max:",DrugPrice_Value_max,")"))
(DrugPrice_ICER_min <- round(min(as.numeric(df[,(A-1)])),digits=1))
(DrugPrice_ICER_max <- round(max(as.numeric(df[,(A-1)])),digits=1))
(DrugPrice_size <- (Base - DrugPrice_ICER_min) + (DrugPrice_ICER_max - Base))


# CostART
B <- 10
colnames(df[B])
(CostART_Value_min <- round(min(as.numeric(df[,B]))*10000,digits=1))
(CostART_Value_max <- round(max(as.numeric(df[,B]))*10000,digits=1))
(CostART <- paste0("Cost of ART \n","(min:",CostART_Value_min,"; max:",CostART_Value_max,")"))
(CostART_ICER_min <- round(min(as.numeric(df[,(B-1)])),digits=1))
(CostART_ICER_max <- round(max(as.numeric(df[,(B-1)])),digits=1))
(CostART_size <- (Base - CostART_ICER_min) + (CostART_ICER_max - Base))


# Screening Negative
D <- 4
colnames(df[D])
(ScreeningNegative_Value_min <- round(min(as.numeric(df[,D]))*10000,digits=1))
(ScreeningNegative_Value_max <- round(max(as.numeric(df[,D]))*10000,digits=1))
(Negative <- paste0("Cost of Screening Negative \n","(min:",ScreeningNegative_Value_min,"; max:",ScreeningNegative_Value_max,")"))
(ScreeningNegative_ICER_min <- round(min(as.numeric(df[,(D-1)])),digits=1))
(ScreeningNegative_ICER_max <- round(max(as.numeric(df[,(D-1)])),digits=1))
(Negative_size <- (Base - ScreeningNegative_ICER_min) + (ScreeningNegative_ICER_max - Base))


# Screening Positive
E <- 6
colnames(df[E])
(ScreeningPositive_Value_min <- round(min(as.numeric(df[,E]))*10000,digits=1))
(ScreeningPositive_Value_max <- round(max(as.numeric(df[,E]))*10000,digits=1))
(Positive <- paste0("Cost of Screening Positive \n","(min:",ScreeningPositive_Value_min,"; max:",ScreeningPositive_Value_max,")"))
(ScreeningPositive_ICER_min <- round(min(as.numeric(df[,(E-1)])),digits=1))
(ScreeningPositive_ICER_max <- round(max(as.numeric(df[,(E-1)])),digits=1))
(Positive_size <- (Base - ScreeningPositive_ICER_min) + (ScreeningPositive_ICER_max - Base))


# Hospitalization
G <- 8
colnames(df[G])
(Hospitalization_Value_min <- round(min(as.numeric(df[,G]))*10000,digits=1))
(Hospitalization_Value_max <- round(max(as.numeric(df[,G]))*10000,digits=1))
(Hospitalization <- paste0("Cost of Hospitalization \n","(min:",Hospitalization_Value_min,"; max:",Hospitalization_Value_max,")"))
(Hospitalization_ICER_min <- round(min(as.numeric(df[,(G-1)])),digits=1))
(Hospitalization_ICER_max <- round(max(as.numeric(df[,(G-1)])),digits=1))
(Hospitalization_size <- (Base - Hospitalization_ICER_min) + (Hospitalization_ICER_max - Base))

df2 <- data.frame(names = c(DrugPrice, CostART, Negative, Positive, Hospitalization),
                  min = c(DrugPrice_ICER_min,CostART_ICER_min,ScreeningNegative_ICER_min,ScreeningPositive_ICER_min,Hospitalization_ICER_min),
                  max = c(DrugPrice_ICER_max,CostART_ICER_max,ScreeningNegative_ICER_max,ScreeningPositive_ICER_max,Hospitalization_ICER_max))


dat <- melt(df2, id.vars = "names",
            variable.name = "val",
            value.name = "output")

# names val  output
# 1                  Cost of PrEP \n(min:3863.6; max:6000) min  433.10
# 2               Cost of ART \n(min:1801716; max:2573880) min  430.15 <- Flip
# 3      Cost of Screening Negative \n(min:4040; max:8190) min  431.91 
# 4    Cost of Screening Positive \n(min:16430; max:45220) min  387.43 <- Flip
# 5  Cost of Hospitalization \n(min:4339680; max:46113840) min  173.18 <- Flip
# 6                  Cost of PrEP \n(min:3863.6; max:6000) max 1026.59
# 7               Cost of ART \n(min:1801716; max:2573880) max  561.47
# 8      Cost of Screening Negative \n(min:4040; max:8190) max  434.79
# 9    Cost of Screening Positive \n(min:16430; max:45220) max  446.18
# 10 Cost of Hospitalization \n(min:4339680; max:46113840) max  556.13

dat[2,2]<-"max"
dat[7,2]<-"min"
dat[4,2]<-"max"
dat[9,2]<-"min"
dat[5,2]<-"max"
dat[10,2]<-"min"
dat

size <-rep(c(DrugPrice_size,CostART_size,Negative_size,Positive_size,Hospitalization_size),2)
dat <- cbind(dat,size)
dat <- arrange(dat,desc(size))
class(dat) <- c("tornado", class(dat))
attr(dat, "output_name") <- "output"

source(file = "Function_code/tornado_function_code.R")
tornado <- ggplot_tornado(dat, Base)
tornado
#Save
#ggsave(filename = paste0(folder,"/Figure4_cost.png"), tornado, dpi=300, width = 280, height = 150, units = "mm")


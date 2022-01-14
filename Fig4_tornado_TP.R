############ Tornado plot for TP
############ Original: Dec 27 2020, Modified Nov 27 2021
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
folder="Figure"
foldername = "data_result"
Nyear = 15

### Library 
library(tidyverse)
library(magrittr)
library(reshape2)

source(file = "Function_code/DSA_TP.R")
DSA_TP(Nyear,foldername)
df <- read.csv(paste0(foldername,"/TP_results_",Nyear,".csv"))
df <- df[,-1]

source(file ="Function_code/ICER_calculation.R")
ICER_base(Nyear,foldername)
icer <- read.csv(paste0(foldername,"/ICER_results_",Nyear,".csv"))
Base = icer[1,2]*10000

##### Define min and max
# NoPrEP -> HIV
A <- 4
colnames(df[A])
(NoPrEP_HIV_Value_min <- round(min(as.numeric(df[,A])),digits = 2))
(NoPrEP_HIV_Value_max <- round(max(as.numeric(df[,A])),digits = 2))
(NoPrEP_HIV <- paste0("No PrEP to HIV \n","(min:",NoPrEP_HIV_Value_min,"; max:",NoPrEP_HIV_Value_max,")"))
colnames(df[A-1])
(NoPrEP_HIV_ICER_min <- round(min(as.numeric(df[,(A-1)])),digits = 2))
(NoPrEP_HIV_ICER_max <- round(max(as.numeric(df[,(A-1)])),digits = 2))
(NoPrEP_HIV_size <- (Base - NoPrEP_HIV_ICER_min) + (NoPrEP_HIV_ICER_max - Base))


# PrEP -> HIV
B <- 2
colnames(df[B])
(PrEP_HIV_Value_min <- round(min(as.numeric(df[,B])),digits = 2))
(PrEP_HIV_Value_max <- round(max(as.numeric(df[,B])),digits = 2))
(PrEP_HIV <- paste0("PrEP to HIV \n","(min:",PrEP_HIV_Value_min,"; max:",PrEP_HIV_Value_max,")"))
colnames(df[(B-1)])
(PrEP_HIV_ICER_min <- round(min(as.numeric(df[,(B-1)])),digits = 2))
(PrEP_HIV_ICER_max <- round(max(as.numeric(df[,(B-1)])),digits = 2))
(PrEP_HIV_size <- (Base - PrEP_HIV_ICER_min) + (PrEP_HIV_ICER_max - Base))


# AIDS -> death
C <- 10
colnames(df[C])
(AIDS_death_Value_min <- round(min(as.numeric(df[,C])),digits = 2))
(AIDS_death_Value_max <- round(max(as.numeric(df[,C])),digits = 2))
(AIDS_death <- paste0("AIDS to death \n","(min:",AIDS_death_Value_min,"; max:",AIDS_death_Value_max,")"))
colnames(df[(C-1)])
(AIDS_death_ICER_min <- round(min(as.numeric(df[,(C-1)])),digits = 2))
(AIDS_death_ICER_max <- round(max(as.numeric(df[,(C-1)])),digits = 2))
(AIDS_death_size <-  (Base - AIDS_death_ICER_min) + (AIDS_death_ICER_max - Base))


# HIV -> AIDS
D <- 6
colnames(df[D])
(HIV_AIDS_Value_min <- round(min(as.numeric(df[,D])),digits = 2))
(HIV_AIDS_Value_max <- round(max(as.numeric(df[,D])),digits = 2))
(HIV_AIDS <- paste0("HIV to AIDS \n","(min:",HIV_AIDS_Value_min,"; max:",HIV_AIDS_Value_max,")"))
colnames(df[(D-1)])
(HIV_AIDS_ICER_min <- round(min(as.numeric(df[,(D-1)])),digits = 2))
(HIV_AIDS_ICER_max <- round(max(as.numeric(df[,(D-1)])),digits = 2))
(HIV_AIDS_size <-  (Base - HIV_AIDS_ICER_min) + (HIV_AIDS_ICER_max - Base))


# natural death
E <- 12
colnames(df[E])
(natural_death_Value_min <- round(min(as.numeric(df[,E])),digits = 2))
(natural_death_Value_max <- round(max(as.numeric(df[,E])),digits = 2))
(natural_death <- paste0("Natural death \n","(min:",natural_death_Value_min,"; max:",natural_death_Value_max,")"))
colnames(df[(E-1)])
(natural_death_ICER_min <- round(min(as.numeric(df[,(E-1)])),digits = 2))
(natural_death_ICER_max <- round(max(as.numeric(df[,(E-1)])),digits = 2))
(natural_death_size <-  (Base - natural_death_ICER_min) + (natural_death_ICER_max - Base))


# HIV death
G <- 8
colnames(df[G])
(HIV_death_Value_min <- round(min(as.numeric(df[,G])),digits = 2))
(HIV_death_Value_max <- round(max(as.numeric(df[,G])),digits = 2))
(HIV_death <- paste0("HIV to death \n","(min:",HIV_death_Value_min,"; max:",HIV_death_Value_max,")"))
colnames(df[(G-1)])
(HIV_death_ICER_min <- round(min(as.numeric(df[,(G-1)])),digits = 2))
(HIV_death_ICER_max <- round(max(as.numeric(df[,(G-1)])),digits = 2))
(HIV_death_size <-  (Base - HIV_death_ICER_min) + (HIV_death_ICER_max - Base))


df2 <- data.frame(names = c(NoPrEP_HIV, PrEP_HIV, AIDS_death, HIV_AIDS,natural_death,HIV_death),
                  min = c(NoPrEP_HIV_ICER_min,PrEP_HIV_ICER_min,AIDS_death_ICER_min,HIV_AIDS_ICER_min,natural_death_ICER_min,HIV_death_ICER_min),
                  max = c(NoPrEP_HIV_ICER_max,PrEP_HIV_ICER_max,AIDS_death_ICER_max,HIV_AIDS_ICER_max,natural_death_ICER_max,HIV_death_ICER_max))

dat <- melt(df2, id.vars = "names",
            variable.name = "val",
            value.name = "output")

# dat
# names val  output
# 1  No PrEP to HIV \n(min:0.03; max:0.05) min  179.59 <- Flip
# 2        PrEP to HIV \n(min:0; max:0.02) min  424.89
# 3   AIDS to death \n(min:0.03; max:0.05) min  429.36
# 4        HIV to AIDS \n(min:0; max:0.09) min -236.14 <- Flip
# 5      Natural death \n(min:0; max:0.06) min  432.42
# 6       HIV to death \n(min:0; max:0.01) min  432.50
# 7  No PrEP to HIV \n(min:0.03; max:0.05) max  753.79
# 8        PrEP to HIV \n(min:0; max:0.02) max  964.28
# 9   AIDS to death \n(min:0.03; max:0.05) max  441.06
# 10       HIV to AIDS \n(min:0; max:0.09) max  606.19
# 11     Natural death \n(min:0; max:0.06) max  497.06
# 12      HIV to death \n(min:0; max:0.01) max  434.10


dat[1,2]<-"max"
dat[7,2]<-"min"
dat[4,2]<-"max"
dat[10,2]<-"min"
dat


size <-rep(c(NoPrEP_HIV_size, PrEP_HIV_size, AIDS_death_size, HIV_AIDS_size,natural_death_size,HIV_death_size),2)
dat <- cbind(dat,size)
dat <- arrange(dat,desc(size))
class(dat) <- c("tornado", class(dat))
attr(dat, "output_name") <- "output"

source(file = "Function_code/tornado_function_code.R")
ggplot_tornado(dat, Base)

ver1 <- ggplot_tornado(dat, Base)

#Save
#ggsave(filename = paste0(folder,"/Figure4_TP.png"), ver1, dpi=300, width = 280, height = 150, units = "mm")

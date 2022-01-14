############ Tornado plot for qaly
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

source(file = "Function_code/DSA_qaly.R")
DSA_QALY(Nyear,foldername)
df <- read.csv(paste0(foldername,"/QALY_results_",Nyear,".csv"))
df <- df[,-1]

source(file ="Function_code/ICER_calculation.R")
ICER_base(Nyear,foldername)
icer <- read.csv(paste0(foldername,"/ICER_results_",Nyear,".csv"))
Base = icer[1,2]*10000

##### Define min and max
# PrEP 
B <- 2
colnames(df[B])
(PrEP_Value_min <- round(min(as.numeric(df[,B])),digits = 2))
(PrEP_Value_max <- round(max(as.numeric(df[,B])),digits = 2))
(PrEP <- paste0("QALY of PrEP \n","(min:",PrEP_Value_min,"; max:",PrEP_Value_max,")"))
colnames(df[(B-1)])
(PrEP_ICER_min <- round(min(as.numeric(df[,(B-1)])),digits = 2))
(PrEP_ICER_max <- round(max(as.numeric(df[,(B-1)])),digits = 2))
(PrEP_size <- (Base - PrEP_ICER_min) + (PrEP_ICER_max - Base))

(df[df$Target_PrEP == PrEP_Value_min,B])
(df[df$Target_PrEP == PrEP_Value_max,B])

# HIV
C <- 6
colnames(df[C])
(HIV_Value_min <- round(min(as.numeric(df[,C])),digits = 2))
(HIV_Value_max <- round(max(as.numeric(df[,C])),digits = 2))
(HIV <- paste0("QALY of HIV \n","(min:",HIV_Value_min,"; max:",HIV_Value_max,")"))
colnames(df[(C-1)])
(HIV_ICER_min <- round(min(as.numeric(df[,(C-1)])),digits = 2))
(HIV_ICER_max <- round(max(as.numeric(df[,(C-1)])),digits = 2))
(HIV_size <-  (Base - HIV_ICER_min) + (HIV_ICER_max - Base))

(df[df$Target_HIV == HIV_Value_min,C])
(df[df$Target_HIV == HIV_Value_max,C])

# AIDS
D <- 8
colnames(df[D])
(AIDS_Value_min <- round(min(as.numeric(df[,D])),digits = 2))
(AIDS_Value_max <- round(max(as.numeric(df[,D])),digits = 2))
(AIDS <- paste0("QALY of AIDS \n","(min:",AIDS_Value_min,"; max:",AIDS_Value_max,")"))
colnames(df[(D-1)])
(AIDS_ICER_min <- round(min(as.numeric(df[,(D-1)])),digits = 2))
(AIDS_ICER_max <- round(max(as.numeric(df[,(D-1)])),digits = 2))
(AIDS_size <-  (Base - AIDS_ICER_min) + (AIDS_ICER_max - Base))

(df[df$Target_AIDS == AIDS_Value_min,D])  
(df[df$Target_AIDS == AIDS_Value_max,D])  

df2 <- data.frame(names = c( PrEP, HIV, AIDS),
                  min = c(PrEP_ICER_min,HIV_ICER_min,AIDS_ICER_min),
                  max = c(PrEP_ICER_max,HIV_ICER_max,AIDS_ICER_max))

dat <- melt(df2, id.vars = "names",
            variable.name = "val",
            value.name = "output")

# > dat
# names val  output
# 1    QALY of PrEP \n(min:0.93; max:1) min  433.10 <- Flip
# 2  QALY of HIV \n(min:0.58; max:0.99) max  406.91
# 3 QALY of AIDS \n(min:0.29; max:0.62) min  426.84
# 4    QALY of PrEP \n(min:0.93; max:1) max  941.86 <- Flip
# 5  QALY of HIV \n(min:0.58; max:0.99) max 1070.30
# 6 QALY of AIDS \n(min:0.29; max:0.62) max  439.96


dat[1,2]<-"max"
dat[4,2]<-"min"

size <-rep(c(PrEP_size, HIV_size, AIDS_size),2)
dat <- cbind(dat,size)
dat <- arrange(dat,desc(size))
class(dat) <- c("tornado", class(dat))
attr(dat, "output_name") <- "output"

source(file = "Function_code/tornado_function_code.R")
ggplot_tornado(dat, Base)

ver1 <- ggplot_tornado(dat, Base)

#Save
#ggsave(filename = paste0(folder,"/Figure4_qaly.png"), ver1, dpi=300, width = 280, height = 150, units = "mm")


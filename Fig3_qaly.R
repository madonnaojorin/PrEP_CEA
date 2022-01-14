############ QALY for varied PrEP coverage
############ Original: Dec 26 2020, Modified Nov 27 2021

rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
folder="Figure"

### Library 
library(tidyverse)
library(magrittr)
theme_set(theme_classic(base_size = 20,base_family = "Helvetica"))

source(file = "Function_code/Cost_QALY_function.R")

################# Storing Cost ################# 

Discount_rate = 0.02 
Initial = 1000
Numcycles = 15
##### Whole qaly 
Whole_coverage_0 <- Discounted_cost(0,Discount_rate,Initial,Numcycles)$Whole_qaly
Whole_coverage_0 <- Whole_coverage_0[1:((Numcycles+1)*4),]
Whole_coverage_50 <- Discounted_cost(0.5,Discount_rate,Initial,Numcycles)$Whole_qaly
Whole_coverage_50 <- Whole_coverage_50[1:((Numcycles+1)*4),]
Whole_coverage_100 <- Discounted_cost(1,Discount_rate,Initial,Numcycles)$Whole_qaly
Whole_coverage_100 <- Whole_coverage_100[1:((Numcycles+1)*4),]

##### Total qaly and calculating differences between coverages 
qaly_total_0 <- Discounted_cost(0,Discount_rate,Initial,Numcycles)$qaly_total
qaly_total_50 <- Discounted_cost(0.5,Discount_rate,Initial,Numcycles)$qaly_total
qaly_total_100 <- Discounted_cost(1,Discount_rate,Initial,Numcycles)$qaly_total

(difference_0_50<- qaly_total_50 -qaly_total_0)
(difference_0_100<- qaly_total_100 -qaly_total_0)
(Difference_df <- data.frame(Year = rep(qaly_total_0$Year,2),
                             data=c(difference_0_50$Total,difference_0_100$Total),
                             Coverage = c(rep("50",nrow(difference_0_50)),rep("100",nrow(difference_0_100)))))

################# Plotting qaly ################# 
Coverage_0 <- Whole_coverage_0[-(((Numcycles+2)*1):((Numcycles+1)*2)),]
QALY_PrEP_0 <- ggplot(data = Coverage_0)+
  geom_bar(aes(x=Year , y=data, fill=State),
           stat ="identity")+
  theme(axis.text = element_text(face = "bold",colour = "black"),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size=24),
        axis.text.y = element_text(size=24),
        legend.position = c(0.8,0.9))+
  scale_x_continuous(limits=c(0,(Numcycles+1))) +
  scale_y_continuous(breaks=c(0,250,500,750,1000),limits=c(0,1000))+
  ylab("QALY")+
  xlab("Year since the PrEP program initiation")+
  ggtitle("PrEP coverage 0%")+
  scale_fill_manual(values = alpha(c("#AD172A","#42B540","#00468B"),0.7))
QALY_PrEP_0
#ggsave(filename = paste0(folder,"/Figure3_qaly_0.png"), QALY_PrEP_0, dpi = 300,width = 180,height = 200,units = c("mm"))

QALY_PrEP_50 <- ggplot(data = Whole_coverage_50)+
  geom_bar(aes(x=Year , y=data, fill=State),
           stat ="identity")+
  theme(axis.text = element_text(face = "bold",colour = "black"),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size=24),
        axis.text.y = element_text(size=24),
        legend.position = c(0.8,0.9))+
  scale_x_continuous(limits=c(0,(Numcycles+1))) +
  scale_y_continuous(breaks=c(0,250,500,750,1000),limits=c(0,1000))+
  ylab("QALY")+
  xlab("Year since the PrEP program initiation")+
  ggtitle("PrEP coverage 50%")+
  scale_fill_manual(values = alpha(c("#AD172A", "#0099B4","#42B540","#00468B"),0.7))
QALY_PrEP_50
#ggsave(filename = paste0(folder,"/Figure3_qaly_50.png"), QALY_PrEP_50, dpi = 300,width = 180,height = 200,units = c("mm"))

Coverage_100 <- Whole_coverage_100[-(1:((Numcycles+1)*1)),]
QALY_PrEP_100 <- ggplot(data = Coverage_100)+
  geom_bar(aes(x=Year , y=data, fill=State),
           stat ="identity")+
  theme(axis.text = element_text(face = "bold",colour = "black"),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size=24),
        axis.text.y = element_text(size=24),
        legend.position = c(0.8,0.9))+
  scale_x_continuous(limits=c(0,(Numcycles+1))) +
  scale_y_continuous(breaks=c(0,250,500,750,1000),limits=c(0,1000))+
  ylab("QALY")+
  xlab("Year since the PrEP program initiation")+
  ggtitle("PrEP coverage 100%")+
  scale_fill_manual(values = alpha(c("#0099B4","#42B540","#00468B"),0.7))
QALY_PrEP_100
#ggsave(filename = paste0(folder,"/Figure3_qaly_100.png"), QALY_PrEP_100, dpi = 300,width = 180,height = 200,units = c("mm"))

################# Plotting Total Cost ################# 
Coverage_0%<>%mutate(Cumulative = 0)
Coverage_0[Coverage_0$State=="Without PrEP",]$Cumulative =cumsum(Coverage_0[Coverage_0$State=="Without PrEP",]$data)
Coverage_0[Coverage_0$State=="HIV",]$Cumulative =cumsum(Coverage_0[Coverage_0$State=="HIV",]$data)
Coverage_0[Coverage_0$State=="AIDS",]$Cumulative =cumsum(Coverage_0[Coverage_0$State=="AIDS",]$data)

CumQALY_PrEP_0 <- ggplot(data = Coverage_0)+
  geom_bar(aes(x=Year , y=Cumulative, fill=State),
           stat ="identity")+
  theme(axis.text = element_text(face = "bold",colour = "black"),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size=24),
        axis.text.y = element_text(size=24),
        legend.position = c(0.25,0.8))+
  scale_x_continuous(limits=c(0,(Numcycles+1))) +
  scale_y_continuous(breaks=c(0,4000,8000,12000),limits=c(0,14000))+
  ylab("QALY")+
  ggtitle("PrEP coverage 0%")+
  xlab("Year since the PrEP program initiation")+
  scale_fill_manual(values = alpha(c("#AD172A","#42B540","#00468B"),0.7))
CumQALY_PrEP_0
#ggsave(filename = paste0(folder,"/Figure3_cumqaly_0.png"), CumQALY_PrEP_0, dpi = 300,width = 180,height = 200,units = c("mm"))


Whole_coverage_50%<>%mutate(Cumulative = 0)
Whole_coverage_50[Whole_coverage_50$State=="Without PrEP",]$Cumulative =cumsum(Whole_coverage_50[Whole_coverage_50$State=="Without PrEP",]$data)
Whole_coverage_50[Whole_coverage_50$State=="With PrEP",]$Cumulative =cumsum(Whole_coverage_50[Whole_coverage_50$State=="With PrEP",]$data)
Whole_coverage_50[Whole_coverage_50$State=="HIV",]$Cumulative =cumsum(Whole_coverage_50[Whole_coverage_50$State=="HIV",]$data)
Whole_coverage_50[Whole_coverage_50$State=="AIDS",]$Cumulative =cumsum(Whole_coverage_50[Whole_coverage_50$State=="AIDS",]$data)
CumQALY_PrEP_50 <- ggplot(data = Whole_coverage_50)+
  geom_bar(aes(x=Year , y=Cumulative, fill=State),
           stat ="identity")+
  theme(axis.text = element_text(face = "bold",colour = "black"),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size=24),
        axis.text.y = element_text(size=24),
        legend.position = c(0.25,0.8))+
  scale_x_continuous(limits=c(0,(Numcycles+1))) +
  scale_y_continuous(breaks=c(0,4000,8000,12000),limits=c(0,14000))+
  ylab("QALY")+
  xlab("Year since the PrEP program initiation")+
  ggtitle("PrEP coverage 50%")+
  scale_fill_manual(values = alpha(c("#AD172A", "#0099B4","#42B540","#00468B"),0.7))
CumQALY_PrEP_50
#ggsave(filename = paste0(folder,"/Figure3_cumqaly_50.png"), CumQALY_PrEP_50, dpi = 300,width = 180,height = 200,units = c("mm"))

Coverage_100%<>%mutate(Cumulative = 0)
Coverage_100[Coverage_100$State=="With PrEP",]$Cumulative =cumsum(Coverage_100[Coverage_100$State=="With PrEP",]$data)
Coverage_100[Coverage_100$State=="HIV",]$Cumulative =cumsum(Coverage_100[Coverage_100$State=="HIV",]$data)
Coverage_100[Coverage_100$State=="AIDS",]$Cumulative =cumsum(Coverage_100[Coverage_100$State=="AIDS",]$data)

CumQALY_PrEP_100 <- ggplot(data = Coverage_100)+
  geom_bar(aes(x=Year , y=Cumulative, fill=State),
           stat ="identity")+
  theme(axis.text = element_text(face = "bold",colour = "black"),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size=24),
        axis.text.y = element_text(size=24),
        legend.position = c(0.25,0.8))+
  scale_x_continuous(limits=c(0,(Numcycles+1))) +
  scale_y_continuous(breaks=c(0,4000,8000,12000),limits=c(0,14000))+
  ylab("QALY")+
  xlab("Year since the PrEP program initiation")+
  ggtitle("PrEP coverage 100%")+
  scale_fill_manual(values = alpha(c("#0099B4","#42B540","#00468B"),0.7))
CumQALY_PrEP_100
#ggsave(filename = paste0(folder,"/Figure3_cumqaly_100.png"), CumQALY_PrEP_100, dpi = 300,width = 180,height = 200,units = c("mm"))

############ Cost for varied PrEP coverage
############ Original: Dec 27 2020, Modified Nov 27 2021

rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
folder="Figure"

### Library 
library(tidyverse)
library(gridExtra)
library(magrittr)
theme_set(theme_classic(base_size = 20, base_family = "Helvetica"))

source(file = "Function_code/Cost_QALY_function.R")

################# Storing Cost ################# 
Discount_rate = 0.02 
Initial = 1000
Numcycles =15

Rate = 110  # 1 USD = 110 JPY

##### Whole cost 
Whole_coverage_0 <- Discounted_cost(0,Discount_rate,Initial, Numcycles)$Whole_cost
Whole_coverage_0 <- Whole_coverage_0[1:((Numcycles+1)*4),]
Whole_coverage_0$data <- Whole_coverage_0$data/Rate
Whole_coverage_50 <- Discounted_cost(0.5,Discount_rate,Initial, Numcycles)$Whole_cost
Whole_coverage_50 <- Whole_coverage_50[1:((Numcycles+1)*4),]
Whole_coverage_50$data <- Whole_coverage_50$data/Rate
Whole_coverage_100 <- Discounted_cost(1,Discount_rate,Initial, Numcycles)$Whole_cost
Whole_coverage_100 <- Whole_coverage_100[1:((Numcycles+1)*4),]
Whole_coverage_100$data <- Whole_coverage_100$data/Rate

##### Total Cost and calculating differences between coverages 
Cost_total_0 <- Discounted_cost(0,Discount_rate,Initial, Numcycles)$Cost_total
Cost_total_0$Total <- Cost_total_0$Total/Rate
Cost_total_50 <- Discounted_cost(0.5,Discount_rate,Initial, Numcycles)$Cost_total
Cost_total_50$Total <- Cost_total_50$Total/Rate
Cost_total_100 <- Discounted_cost(1,Discount_rate,Initial, Numcycles)$Cost_total
Cost_total_100$Total <- Cost_total_100$Total/Rate

(difference_0_50<- Cost_total_50 -Cost_total_0)
(difference_0_100<- Cost_total_100 -Cost_total_0)
(Difference_df <- data.frame(Year = rep(Cost_total_0$Year,2),
                             data=c(difference_0_50$Total,difference_0_100$Total),
                             Coverage = c(rep("50",nrow(difference_0_50)),rep("100",nrow(difference_0_100)))))

################# Plotting Cost ################# 
Coverage_0 <- Whole_coverage_0[-((Numcycles+1):((Numcycles+1)*2)),]
Cost_PrEP_0  <- ggplot(data = Coverage_0)+
  geom_bar(aes(x=Year , y=data, fill=State),
           stat ="identity")+
  theme(axis.text = element_text(face = "bold",colour = "black"),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size=24),
        axis.text.y = element_text(size=24),
        legend.position = c(0.25,0.9))+
  scale_x_continuous(limits=c(0,(Numcycles+1))) +
  scale_y_continuous(breaks=c(0,500,1000,1500),labels = c(0,5,10,15),limits=c(0,1500))+
  ylab("Costs in million (USD)")+
  ggtitle("PrEP coverage 0%")+
  xlab("Year since the PrEP program initiation")+
  scale_fill_manual(values = alpha(c("#AD172A","#42B540","#00468B"),0.7))
Cost_PrEP_0 
#ggsave(filename = paste0(folder,"/Figure3_cost_0_USD.png"), Cost_PrEP_0, dpi = 300,width = 180,height = 200,units = c("mm"))

Cost_PrEP_50  <- ggplot(data = Whole_coverage_50)+
  geom_bar(aes(x=Year , y=data, fill=State),
           stat ="identity")+
  theme(axis.text = element_text(face = "bold",colour = "black"),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size=24),
        axis.text.y = element_text(size=24),
        legend.position = c(0.25,0.9))+
  scale_x_continuous(limits=c(0,(Numcycles+1))) +
  scale_y_continuous(breaks=c(0,500,1000,1500),labels = c(0,5,10,15),limits=c(0,1500))+
  ylab("Costs in million (USD)")+
  xlab("Year since the PrEP program initiation")+
  ggtitle("PrEP coverage 50%")+
  scale_fill_manual(values = alpha(c("#AD172A", "#0099B4","#42B540","#00468B"),0.7))
Cost_PrEP_50 
#ggsave(filename = paste0(folder,"/Figure3_cost_50_USD.png"), Cost_PrEP_50, dpi = 300,width = 180,height = 200,units = c("mm"))

Coverage_100 <- Whole_coverage_100[-(1:(Numcycles+1)),]
Cost_PrEP_100 <- ggplot(data = Coverage_100)+
  geom_bar(aes(x=Year , y=data, fill=State),
           stat ="identity")+
  theme(axis.text = element_text(face = "bold",colour = "black"),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size=24),
        axis.text.y = element_text(size=24),
        legend.position = c(0.8,0.9))+
  scale_x_continuous(limits=c(0,(Numcycles+1))) +
  scale_y_continuous(breaks=c(0,500,1000,1500),labels = c(0,5,10,15),limits=c(0,1500))+
  ylab("Costs in million (USD)")+
  xlab("Year since the PrEP program initiation")+
  ggtitle("PrEP coverage 100%")+
  scale_fill_manual(values = alpha(c("#0099B4","#42B540","#00468B"),0.7))
Cost_PrEP_100
#ggsave(filename = paste0(folder,"/Figure3_cost_100_USD.png"), Cost_PrEP_100, dpi = 300,width = 180,height = 200,units = c("mm"))


################# Plotting difference ################# 
Difference_cost <- ggplot(data = Difference_df)+
  geom_line(aes(x=Year,y=data,colour=Coverage),size=2)+
  geom_hline(yintercept = 0,linetype=2)+
  theme(axis.text = element_text(face = "bold",colour = "black"),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size=24),
        axis.text.y = element_text(size=24),
        legend.position = c(0.7,0.8))+
  scale_x_continuous(limits=c(0,(Numcycles+1))) +
  scale_y_continuous(breaks=c(-500,0,500,1000,1500),labels = c(-0.5,0,0.5,1.0,1.5),limits=c(-500,1500))+
  ylab("Costs in million (USD)")+
  xlab("Year since the PrEP program initiation")+
  scale_colour_manual(values = c("#ED2200","#00468B"))
Difference_cost

#ggsave(filename = paste0(folder,"/Figure3_cost_diff_USD.png"), Difference_cost, dpi = 300,width = 180,height = 200,units = c("mm"))


################# Plotting Cumulative Cost ################# 
Coverage_0%<>%mutate(Cumulative = 0)
Coverage_0[Coverage_0$State=="Without PrEP",]$Cumulative =cumsum(Coverage_0[Coverage_0$State=="Without PrEP",]$data)
Coverage_0[Coverage_0$State=="HIV",]$Cumulative =cumsum(Coverage_0[Coverage_0$State=="HIV",]$data)
Coverage_0[Coverage_0$State=="AIDS",]$Cumulative =cumsum(Coverage_0[Coverage_0$State=="AIDS",]$data)

CumCost_PrEP_0 <- ggplot(data = Coverage_0)+
  geom_bar(aes(x=Year , y=Cumulative, fill=State),
           stat ="identity")+
  theme(axis.text = element_text(face = "bold",colour = "black"),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size=24),
        axis.text.y = element_text(size=24),
        legend.position = c(0.25,0.8))+
  scale_x_continuous(limits=c(0,(Numcycles+1))) +
  scale_y_continuous(breaks=c(0,5000,10000,15000),labels = c(0,5.0,10.0,15.0),limits=c(0,17000))+
  ylab("Cumulative costs in million (USD)")+
  ggtitle("PrEP coverage 0%")+
  xlab("Year since the PrEP program initiation")+
  scale_fill_manual(values = alpha(c("#AD172A","#42B540","#00468B"),0.7))
CumCost_PrEP_0
#ggsave(filename = paste0(folder,"/Figure3_cumcost_0_USD.png"), CumCost_PrEP_0, dpi = 300,width = 180,height = 200,units = c("mm"))


Whole_coverage_50%<>%mutate(Cumulative = 0)
Whole_coverage_50[Whole_coverage_50$State=="Without PrEP",]$Cumulative =cumsum(Whole_coverage_50[Whole_coverage_50$State=="Without PrEP",]$data)
Whole_coverage_50[Whole_coverage_50$State=="With PrEP",]$Cumulative =cumsum(Whole_coverage_50[Whole_coverage_50$State=="With PrEP",]$data)
Whole_coverage_50[Whole_coverage_50$State=="HIV",]$Cumulative =cumsum(Whole_coverage_50[Whole_coverage_50$State=="HIV",]$data)
Whole_coverage_50[Whole_coverage_50$State=="AIDS",]$Cumulative =cumsum(Whole_coverage_50[Whole_coverage_50$State=="AIDS",]$data)
CumCost_PrEP_50 <- ggplot(data = Whole_coverage_50)+
  geom_bar(aes(x=Year , y=Cumulative, fill=State),
           stat ="identity")+
  theme(axis.text = element_text(face = "bold",colour = "black"),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size=24),
        axis.text.y = element_text(size=24),
        legend.position = c(0.25,0.8))+
  scale_x_continuous(limits=c(0,(Numcycles+1))) +
  scale_y_continuous(breaks=c(0,5000,10000,15000),labels = c(0,5.0,10.0,15.0),limits=c(0,17000))+
  ylab("Cumulative costs in million (USD)")+
  xlab("Year since the PrEP program initiation")+
  ggtitle("PrEP coverage 50%")+
  scale_fill_manual(values = alpha(c("#AD172A", "#0099B4","#42B540","#00468B"),0.7))
CumCost_PrEP_50
#ggsave(filename = paste0(folder,"/Figure3_cumcost_50_USD.png"), CumCost_PrEP_50, dpi = 300,width = 180,height = 200,units = c("mm"))

Coverage_100%<>%mutate(Cumulative = 0)
Coverage_100[Coverage_100$State=="With PrEP",]$Cumulative =cumsum(Coverage_100[Coverage_100$State=="With PrEP",]$data)
Coverage_100[Coverage_100$State=="HIV",]$Cumulative =cumsum(Coverage_100[Coverage_100$State=="HIV",]$data)
Coverage_100[Coverage_100$State=="AIDS",]$Cumulative =cumsum(Coverage_100[Coverage_100$State=="AIDS",]$data)

CumCost_PrEP_100 <- ggplot(data = Coverage_100)+
  geom_bar(aes(x=Year , y=Cumulative, fill=State),
           stat ="identity")+
  theme(axis.text = element_text(face = "bold",colour = "black"),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size=24),
        axis.text.y = element_text(size=24),
        legend.position = c(0.25,0.8))+
  scale_x_continuous(limits=c(0,(Numcycles+1))) +
  scale_y_continuous(breaks=c(0,5000,10000,15000),labels = c(0,5.0,10.0,15.0),limits=c(0,17000))+
  ylab("Cumulative costs in million (USD)")+
  xlab("Year since the PrEP program initiation")+
  ggtitle("PrEP coverage 100%")+
  scale_fill_manual(values = alpha(c("#0099B4","#42B540","#00468B"),0.7))
CumCost_PrEP_100
#ggsave(filename = paste0(folder,"/Figure3_cumcost_100_USD.png"), CumCost_PrEP_100, dpi = 300,width = 180,height = 200,units = c("mm"))

################# Plotting difference ################# 
##### Total Cost and calculating differences between coverages 
Cost_total_0 <- Discounted_cost(0,Discount_rate,Initial,Numcycles)$Cost_total
Cost_total_0$Total<- Cost_total_0$Total/Rate
Cost_total_50 <- Discounted_cost(0.5,Discount_rate,Initial,Numcycles)$Cost_total
Cost_total_50$Total<- Cost_total_50$Total/Rate
Cost_total_100 <- Discounted_cost(1,Discount_rate,Initial,Numcycles)$Cost_total
Cost_total_100$Total<- Cost_total_100$Total/Rate

Cum_Cost_total_0 <-cumsum(Cost_total_0)
Cum_Cost_total_50 <-cumsum(Cost_total_50)
Cum_Cost_total_100 <-cumsum(Cost_total_100)
(difference_0_50_cum <- Cum_Cost_total_50 -Cum_Cost_total_0)
(difference_0_100_cum <- Cum_Cost_total_100 -Cum_Cost_total_0)
(Difference_df <- data.frame(Year = rep(Cost_total_0$Year,2),
                             data=c(difference_0_50_cum$Total,difference_0_100_cum$Total),
                             Coverage = c(rep("50",nrow(difference_0_50_cum)),rep("100",nrow(difference_0_100_cum)))))

Difference_cum <- ggplot(data = Difference_df)+
  geom_line(aes(x=Year,y=data,colour=Coverage),size=2)+
  theme(axis.text = element_text(face = "bold",colour = "black"),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size=24),
        axis.text.y = element_text(size=24),
        legend.position = c(0.2,0.9))+
  scale_x_continuous(limits=c(0,15)) +
  scale_y_continuous(breaks=c(0,3000,6000,9000),labels = c(0,3.0,6.0,9.0),limits=c(0,9000))+
  ylab("Costs in million (USD)")+
  xlab("Year since the PrEP program initiation")+
  scale_colour_manual(values = c("#ED2200","#00468B"))
Difference_cum

#ggsave(filename = paste0(folder,"/Figure3_cumcost_diff_USD.png"), Difference_cum, dpi = 300,width = 180,height = 200,units = c("mm"))


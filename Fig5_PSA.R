############ Probabilistic Sensitivity Analysis
############ Original: Dec 27 2020, Modified Nov 27 2021

rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
folder="Figure"

### Library 
library(tidyverse)
library(fitdistrplus)

theme_set(theme_classic(base_size = 20))
simulation <- 10000
set.seed(100000)
###### Setting parameter for PSA
Rate <- 110
#### Cost of PrEP
PrEP_df <- c(1368031, 1227586, 1179310, 1096551, 1055172, 1062068, 1020689, 986206, 979310, 896551, 852361,1390896)
PrEP_data <- PrEP_df/1000000
PrEP_mean <- mean(PrEP_data)
PrEP_sd <- sd(PrEP_data)

shape = (PrEP_mean)^2/PrEP_sd
scale = PrEP_sd/PrEP_mean
PrEPcost <- rgamma(simulation,shape=shape,scale=scale)/110
PrEP_yearly <- PrEPcost*1000000/10000

ART <- rgamma(simulation,213042/10000*2,2)/110
ART_yearly = ART*12

Hospital <- rgamma(simulation,1480000/10000*0.15,0.15)/110
Hospitalization_yearly = Hospital*12


#### QOL No PrEP
eq5d<-read.csv("analysis_data.csv")
qol_NoPrEP <- eq5d[eq5d$visit == 201,c("visit","EQ_index2005")]
qol_NoPrEP <- na.omit(qol_NoPrEP)
QOL_NoPrEP <- sample(qol_NoPrEP$EQ_index2005,simulation,replace = T)
qol_PrEP <- eq5d[eq5d$visit == 1001,c("visit","EQ_index2005")]
qol_PrEP <- na.omit(qol_NoPrEP)
QOL_PrEP <- QOL_NoPrEP/0.95

source(file = "Function_code/pert.R")
qalys <-read.csv("QALY_value.csv")
row.names(qalys)<-qalys$QALY
QOL_HIV <- rpert(simulation, min=qalys["HIV","MIN"],mod=qalys["HIV","Base"],max=qalys["HIV","MAX"],shape=4)
QOL_AIDS <- rpert(simulation, min=qalys["AIDS","MIN"],mod=qalys["AIDS","Base"],max=qalys["AIDS","MAX"],shape=4)

coverage <- rbeta(simulation,10,10) #0.5
Numcycles <- 15

Simulation_cost_total <- c()
Simulation_qaly_total <- c()
Cost_and_QALY <- data.frame()
Simulation_cost_total_control <- c()
Simulation_qaly_total_control <- c()
Cost_and_QALY_control <- data.frame()

for (p in 1:simulation){
  ###### Setting parameters for COST #####
  cost <- read.csv("Cost_value.csv")
  cost[1:5,2:4] <- cost[1:5,2:4]/Rate 
  rownames(cost) <- cost[1:5,1]
  cost <- cost[,-1]
  
  ### Cost for PrEP ###
  (PrEP_Yearly = PrEP_yearly[p])
  
  ### Cost for ART ###
  (ART_Yearly = ART_yearly[p])
  
  ### Hospitalization ###
  (Hospitalization_Yearly = Hospitalization_yearly[p])
  #(Hospitalization_Monthly = as.numeric(cost["Hospitalization","Base"])/10000)
  #(Hospitalization_Yearly = Hospitalization_Monthly*12)
  
  ### Screening  ###
  Screening_Neg_Monthly = as.numeric(cost["Screening_Negative","Base"])/10000
  Screening_Pos_Monthly = as.numeric(cost["Screening_Positive","Base"])/10000
  
  (Screening_Neg_Yearly = Screening_Neg_Monthly * 4)
  # once per half a year
  (Screening_Pos_Yearly = Screening_Pos_Monthly * 12)
  # once a month
  
  ###### Set parameters for Markov model ######
  cycles <- Numcycles
  Yearth <- cycles + 1
  (QALY_NoPrEP <- QOL_NoPrEP[p])
  (QALY_PrEP <- QOL_PrEP[p])
  (QALY_HIV <- QOL_HIV[p])
  (QALY_AIDS <- QOL_AIDS[p])
  QALY_Death <- 0
  (Cost_NoPrEP <-
      Screening_Neg_Yearly)
  (Cost_PrEP <-
      Screening_Neg_Yearly + PrEP_Yearly)
  (Cost_HIV <-
      Screening_Pos_Yearly + ART_Yearly)
  (Cost_AIDS <-
      Screening_Pos_Yearly + ART_Yearly + Hospitalization_Yearly)
  Cost_Death <- 0
  
  Pop_initial <- 1000 
  
  (Percentage_covered <- coverage[p])  
  dis_rate = 0.02
  
  ###### Set parameters for transition probability
  TP <- read.csv("TP_value.csv")
  rownames(TP) <- TP[1:6,1]
  TP <- TP[,-1]
  
  natural_Death = TP["Natural_death","Base"]
  HIV_AIDS = TP["HIV_AIDS","Base"]
  HIV_Death = TP["HIV_death","Base"]
  AIDS_Death = TP["AIDS_death","Base"]
  NoPrEP_HIV = TP["NoPrEP_HIV","Base"]
  NoPrEP_AIDS = 0
  NoPrEP_Death = natural_Death
  PrEP_HIV = TP["PrEP_HIV","Base"]
  PrEP_AIDS = 0
  PrEP_Death = natural_Death
  
  ###### Calculate initial value
  (NoPrEP_initial <- Pop_initial*(1-Percentage_covered))
  (PrEP_initial <- Pop_initial*Percentage_covered)
  
  ###### Set parameters for WITHOUT PrEP ######
  #1. No PrEP
  #2. HIV positive
  #3. AIDS
  #4. Death
  
  Transition_from1_without <- c((1 - (NoPrEP_HIV + NoPrEP_AIDS + NoPrEP_Death)), NoPrEP_HIV, NoPrEP_AIDS, NoPrEP_Death)
  Transition_from2_without <- c(0, (1-HIV_AIDS-HIV_Death), HIV_AIDS, HIV_Death)
  Transition_from3_without <- c(0, 0, (1-AIDS_Death), AIDS_Death)
  Transition_from4_without <- c(0,0,0,1)
  
  (TP_without <- matrix(data = c(Transition_from1_without,Transition_from2_without,Transition_from3_without,Transition_from4_without), nrow = 4, byrow = TRUE))
  rowSums(TP_without) ### Check if probabilities adding up to 1
  
  trace_NoPrEP <- matrix(data = NA, nrow = cycles+1, ncol = 4)
  trace_NoPrEP_cost <- matrix(data = NA, nrow = cycles+1, ncol = 4)
  trace_NoPrEP_qaly <- matrix(data = NA, nrow = cycles+1, ncol = 4)
  colnames(trace_NoPrEP)<- c("NoPrEP","HIV","AIDS","Death")
  colnames(trace_NoPrEP_cost)<- c("NoPrEP","HIV","AIDS","Death")
  colnames(trace_NoPrEP_qaly)<- c("NoPrEP","HIV","AIDS","Death")
  
  trace_NoPrEP[1,] <- c(NoPrEP_initial,0,0,0)
  trace_NoPrEP_cost[1,] <- c(NoPrEP_initial*Cost_NoPrEP,0,0,0)
  trace_NoPrEP_qaly[1,] <- c(NoPrEP_initial*QALY_NoPrEP,0,0,0)
  
  NoPrEP_cost_initial <- c(Cost_NoPrEP,Cost_HIV,Cost_AIDS,Cost_Death) 
  NoPrEP_qaly_initial <- c(QALY_NoPrEP,QALY_HIV,QALY_AIDS,QALY_Death) 
  
  for (i in 1:cycles){
    (NoPrEP_discounted_cost <- NoPrEP_cost_initial/(1+dis_rate)^(i-1))
    (NoPrEP_discounted_qaly <- NoPrEP_qaly_initial/(1+dis_rate)^(i-1))
    
    (trace_NoPrEP[i+1,] <- trace_NoPrEP[i,] %*% TP_without)
    (trace_NoPrEP_cost[i+1,] <- trace_NoPrEP[i+1,] * (NoPrEP_discounted_cost))
    (trace_NoPrEP_qaly[i+1,] <- trace_NoPrEP[i+1,] * (NoPrEP_discounted_qaly))
    (trace_NoPrEP[i+1,] <- trace_NoPrEP[i+1,])
  }
  
  NoPrEP_cost <- data.frame(data = c(trace_NoPrEP_cost[,"NoPrEP"],trace_NoPrEP_cost[,"HIV"],trace_NoPrEP_cost[,"AIDS"],trace_NoPrEP_cost[,"Death"]),
                            State = c(rep("No PrEP",Yearth),rep("HIV",Yearth),rep("AIDS",Yearth),rep("Death",Yearth)),
                            Year = c(rep(1:Yearth,4)))
  NoPrEP_cost$State <- factor(NoPrEP_cost$State,levels = c("No PrEP","HIV","AIDS","Death"))
  
  NoPrEP_qaly <- data.frame(data = c(trace_NoPrEP_qaly[,"NoPrEP"],trace_NoPrEP_qaly[,"HIV"],trace_NoPrEP_qaly[,"AIDS"],trace_NoPrEP_qaly[,"Death"]),
                            State = c(rep("No PrEP",Yearth),rep("HIV",Yearth),rep("AIDS",Yearth),rep("Death",Yearth)),
                            Year = c(rep(1:Yearth,4)))
  NoPrEP_qaly$State <- factor(NoPrEP_qaly$State,levels = c("No PrEP","HIV","AIDS","Death"))
  
  ###### Set parameters for WITH PrEP starting from WITHOUT PrEP ######
  #1. WITH PrEP
  #2. HIV positive
  #3. AIDS
  #4. Death
  
  Transition_from1_with <- c((1 - (PrEP_HIV + PrEP_AIDS + PrEP_Death)), PrEP_HIV, PrEP_AIDS, PrEP_Death)
  Transition_from2_with <- c(0, (1-HIV_AIDS-HIV_Death), HIV_AIDS, HIV_Death)
  Transition_from3_with <- c(0, 0, (1-AIDS_Death), AIDS_Death)
  Transition_from4_with <- c(0,0,0,1)
  
  (TP_with <- matrix(data = c(Transition_from1_with,Transition_from2_with,Transition_from3_with,Transition_from4_with), nrow = 4, byrow = TRUE))
  rowSums(TP_with) ### Check if probabilities adding up to 1
  
  trace_PrEP <- matrix(data = NA, nrow = cycles+1, ncol = 4)
  trace_PrEP_cost <- matrix(data = NA, nrow = cycles+1, ncol = 4)
  trace_PrEP_qaly <- matrix(data = NA, nrow = cycles+1, ncol = 4)
  colnames(trace_PrEP)<- c("PrEP","HIV","AIDS","Death")
  colnames(trace_PrEP_cost)<- c("PrEP","HIV","AIDS","Death")
  colnames(trace_PrEP_qaly)<- c("PrEP","HIV","AIDS","Death")
  
  (trace_PrEP[1,] <- c(PrEP_initial,0,0,0))
  trace_PrEP_cost[1,] <- c(PrEP_initial*Cost_PrEP,0,0,0)
  trace_PrEP_qaly[1,] <- c(PrEP_initial*QALY_PrEP,0,0,0)
  
  PrEP_cost_initial <- c(Cost_PrEP,Cost_HIV,Cost_AIDS,Cost_Death) 
  PrEP_qaly_initial <- c(QALY_PrEP,QALY_HIV,QALY_AIDS,QALY_Death) 
  
  for (i in 1:cycles){
    (PrEP_discounted_cost <- PrEP_cost_initial/(1+dis_rate)^(i-1))
    (PrEP_discounted_qaly <- PrEP_qaly_initial/(1+dis_rate)^(i-1))
    
    (trace_PrEP[i+1,] <- trace_PrEP[i,] %*% TP_with)
    (trace_PrEP_cost[i+1,] <- trace_PrEP[i+1,] * (PrEP_discounted_cost))
    (trace_PrEP_qaly[i+1,] <- trace_PrEP[i+1,] * (PrEP_discounted_qaly))
    (trace_PrEP[i+1,] <- trace_PrEP[i+1,])
  }
  
  PrEP_cost <- data.frame(data = c(trace_PrEP_cost[,"PrEP"],trace_PrEP_cost[,"HIV"],trace_PrEP_cost[,"AIDS"],trace_PrEP_cost[,"Death"]),
                          State = c(rep("PrEP",Yearth),rep("HIV",Yearth),rep("AIDS",Yearth),rep("Death",Yearth)),
                          Year = c(rep(1:Yearth,4)))
  PrEP_cost$State <- factor(PrEP_cost$State,levels = c("PrEP","HIV","AIDS","Death"))
  
  PrEP_qaly <- data.frame(data = c(trace_PrEP_qaly[,"PrEP"],trace_PrEP_qaly[,"HIV"],trace_PrEP_qaly[,"AIDS"],trace_PrEP_qaly[,"Death"]),
                          State = c(rep("PrEP",Yearth),rep("HIV",Yearth),rep("AIDS",Yearth),rep("Death",Yearth)),
                          Year = c(rep(1:Yearth,4)))
  PrEP_qaly$State <- factor(PrEP_qaly$State,levels = c("PrEP","HIV","AIDS","Death"))
  
  ###### Calculate Cost
  Simulation_cost_NoPrEP <- sum(NoPrEP_cost$data)
  Simulation_cost_PrEP <- sum(PrEP_cost$data)
  Simulation_qaly_NoPrEP <- sum(NoPrEP_qaly$data)
  Simulation_qaly_PrEP <- sum(PrEP_qaly$data)
  (Simulation_cost_total[p] <- Simulation_cost_NoPrEP + Simulation_cost_PrEP)
  (Simulation_qaly_total[p] <- Simulation_qaly_NoPrEP + Simulation_qaly_PrEP)
  
  (Simulation_total <- c(Simulation_cost_total[p],Simulation_qaly_total[p]))
  Cost_and_QALY <- rbind(Cost_and_QALY,Simulation_total)
  
  #############################################################
  ###### Calculate control group  #############################
  #############################################################
  Control = 0
  ###### Calculate initial value and inflow value for each case
  (Control_initial <- Pop_initial*(1-Control))
  
  ###### Set parameters for WITHOUT PrEP ######
  trace_control <- matrix(data = NA, nrow = cycles+1, ncol = 4)
  trace_control_cost <- matrix(data = NA, nrow = cycles+1, ncol = 4)
  trace_control_qaly <- matrix(data = NA, nrow = cycles+1, ncol = 4)
  colnames(trace_control)<- c("Control","HIV","AIDS","Death")
  colnames(trace_control_cost)<- c("Control","HIV","AIDS","Death")
  colnames(trace_control_qaly)<- c("Control","HIV","AIDS","Death")
  
  trace_control[1,] <- c(Control_initial,0,0,0)
  trace_control_cost[1,] <- c(Control_initial*Cost_NoPrEP,0,0,0)
  trace_control_qaly[1,] <- c(Control_initial*QALY_NoPrEP,0,0,0)
  
  Control_cost_initial <- c(Cost_NoPrEP,Cost_HIV,Cost_AIDS,Cost_Death) 
  Control_qaly_initial <- c(QALY_NoPrEP,QALY_HIV,QALY_AIDS,QALY_Death) 
  
  for (i in 1:cycles){
    (Control_discounted_cost <- Control_cost_initial/(1+dis_rate)^(i-1))
    (Control_discounted_qaly <- Control_qaly_initial/(1+dis_rate)^(i-1))
    
    (trace_control[i+1,] <- trace_control[i,] %*% TP_without)
    (trace_control_cost[i+1,] <- trace_control[i+1,] * (Control_discounted_cost))
    (trace_control_qaly[i+1,] <- trace_control[i+1,] * (Control_discounted_qaly))
    (trace_control[i+1,] <- trace_control[i+1,])
  }
  
  Control_cost <- data.frame(data = c(trace_control_cost[,"Control"],trace_control_cost[,"HIV"],trace_control_cost[,"AIDS"],trace_control_cost[,"Death"]),
                             State = c(rep("Control",Yearth),rep("HIV",Yearth),rep("AIDS",Yearth),rep("Death",Yearth)),
                             Year = c(rep(1:Yearth,4)))
  Control_cost$State <- factor(Control_cost$State,levels = c("Control","HIV","AIDS","Death"))
  
  Control_qaly <- data.frame(data = c(trace_control_qaly[,"Control"],trace_control_qaly[,"HIV"],trace_control_qaly[,"AIDS"],trace_control_qaly[,"Death"]),
                             State = c(rep("Control",Yearth),rep("HIV",Yearth),rep("AIDS",Yearth),rep("Death",Yearth)),
                             Year = c(rep(1:Yearth,4)))
  Control_qaly$State <- factor(Control_qaly$State,levels = c("Control","HIV","AIDS","Death"))
  
  ###### Calculate Cost
  Simulation_cost_control <- sum(Control_cost$data)
  Simulation_qaly_control <- sum(Control_qaly$data)
  (Simulation_cost_total_control[p] <- Simulation_cost_control)
  (Simulation_qaly_total_control[p] <- Simulation_qaly_control)
  (Simulation_total_control <- c(Simulation_cost_total_control[p],Simulation_qaly_total_control[p]))
  Cost_and_QALY_control <- rbind(Cost_and_QALY_control,Simulation_total_control)
}

Difference_cost <- Simulation_cost_total - Simulation_cost_total_control
Difference_qaly <- Simulation_qaly_total - Simulation_qaly_total_control

################## PSA #####################################
WTP = 500/110
serial = c(seq(1,simulation))
simdata = as.data.frame(serial)

## calculate difference
simdata$Difference_qaly = Difference_qaly
simdata$Difference_cost = Difference_cost

## calculate ICER
simdata$icer = simdata$Difference_cost/simdata$Difference_qaly

simdata$model = WTP * simdata$Difference_qaly 

simdata$model_true = simdata$model - simdata$Difference_cost 

simdata$CE = ifelse(test = simdata$model_true>0,yes = "1", no = "0" )

simdata$CE_col = ifelse(test =simdata$CE == 0,yes = "2", no = "3" )
table(simdata$CE)

ICER <- ggplot(simdata) +
  geom_point(aes(x = Difference_qaly, y = Difference_cost,color =CE_col,shape=f),size = 1,shape=20)+
  geom_abline(intercept = 0, slope = WTP, size= 2)+
  scale_x_continuous(limits=c(-1000,2000))+
  scale_y_continuous(labels = c(-5000,0,5000,10000,15000),breaks = c(-5000,0,5000,10000,15000),limits = c(-7000,17000))+
  scale_colour_manual(values=c("#ED2200","#00468B"))+
  ylab("Incremental costs (USD)") +
  xlab("Incremental effectiveness (QALY)") +
  theme(legend.position = "none",
        axis.text = element_text(face = "bold",colour = "black"),
        axis.title = element_text(face = "bold",size=20),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))

ICER

Threshold <- list()
for (w in 1:25){
  WTP_sim <- w *100/110
  simdata$model = WTP_sim * simdata$Difference_qaly 
  
  simdata$model_true = simdata$model - simdata$Difference_cost 
  
  simdata$CE = ifelse(test = simdata$model_true>0,yes = "1", no = "0" )
  
  simdata$CE_col = ifelse(test =simdata$CE == 0,yes = "2", no = "3" )
  Probability <- table(simdata$CE)[2]/simulation
  Threshold[[w]]<-c(WTP_sim,Probability)
}

WTP_data <- do.call("rbind", Threshold)
#WTP_data[17:25,2]<-1
colnames(WTP_data)<- c("Threshold","Probability")
WTP_data <- as.data.frame(WTP_data)
WTP_data$Threshold<-WTP_data$Threshold*1000

WTP_plot <- ggplot(WTP_data) +
  geom_area(data=subset(WTP_data, Threshold>=500/110),aes(x=Threshold,y=Probability) ,fill="grey",alpha=0.5)+
  geom_line(aes(x = Threshold, y = Probability),size=1.5)+
  geom_vline(xintercept = 500/110*1000,size=1,linetype="dashed")+
  ylab("Probability of being cost-effective") +
  xlab("Cost-effectiveness threshold (USD)") +
  scale_y_continuous(breaks = seq(0,1,by=0.1),labels = seq("0","1",by=0.1))+
  theme(axis.text = element_text(face = "bold",colour = "black"),
        axis.title = element_text(face = "bold",size=20),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))

WTP_plot
#Save
#ggsave(filename = paste0(folder,"/Figure5_PSA.png"), ICER, dpi=300, width = 200, height = 150, units = "mm")

#ggsave(filename = paste0(folder,"/Figure5_curve2.png"), WTP_plot, dpi=300, width = 200, height = 150, units = "mm")

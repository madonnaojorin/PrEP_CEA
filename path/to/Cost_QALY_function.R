############ Final code after checking ##############
############ Jan 16th: update all the data for cost and trandition probability 
############ Modified: Nov 27 2021 

Discounted_cost <- function (Coverage, Discount_rate, Initial, Numcycles){
  cost <- read.csv("Cost_value.csv")
  rownames(cost) <- cost[1:5,1]
  cost <- cost[,-1]
  ###### Setting parameters for COST #####
  ### Cost for PrEP ###
  (PrEP_pill = as.numeric(cost["PrEP","Base"]))
  PrEP_Monthly = PrEP_pill*30/10000
  (PrEP_Yearly = PrEP_Monthly*12)

  ### Cost for ART ###
  (ART_Monthly = as.numeric(cost["ART","Base"])/10000)
  (ART_Yearly = ART_Monthly*12)
  
  ### Hospitalization ###
  (Hospitalization_Monthly = as.numeric(cost["Hospitalization","Base"])/10000)
  (Hospitalization_Yearly = Hospitalization_Monthly*12)
  
  ### Screening  ###
  Screening_Neg_Monthly = as.numeric(cost["Screening_Negative","Base"])/10000
  Screening_Pos_Monthly = as.numeric(cost["Screening_Positive","Base"])/10000
  
  (Screening_Neg_Yearly = Screening_Neg_Monthly*4)
  # once per half a year
  (Screening_Pos_Yearly = Screening_Pos_Monthly*12)
  # once a month
  
  (Cost_NoPrEP <- Screening_Neg_Yearly); 
  (Cost_PrEP <- Screening_Neg_Yearly+PrEP_Yearly);
  (Cost_HIV <- Screening_Pos_Yearly+ART_Yearly);
  (Cost_AIDS <- Screening_Pos_Yearly+ART_Yearly+Hospitalization_Yearly); 
  Cost_Death <- 0
  
  ###### Set parameters for QALY ######
  QALY_data <- read.csv("QALY_value.csv")
  
  QALY_PrEP <- QALY_data[1,"Base"]
  QALY_NoPrEP <- QALY_data[2,"Base"]
  QALY_HIV <- QALY_data[3,"Base"]
  QALY_AIDS <- QALY_data[4,"Base"]
  QALY_Death <- QALY_data[5,"Base"]

  ###### Set parameters for Markov model ######
  cycles <- Numcycles
  Yearth <- cycles +1  
  Pop_initial <-Initial
  Percentage_covered <- Coverage
  dis_rate <- Discount_rate
  
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

  ###### Calculate initial value and inflow value for each case
  
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
                            State = c(rep("Without PrEP",Yearth),rep("HIV",Yearth),rep("AIDS",Yearth),rep("Death",Yearth)),
                            Year = c(rep(1:Yearth,4)))
  NoPrEP_cost$State <- factor(NoPrEP_cost$State,levels = c("Without PrEP","HIV","AIDS","Death"))
  
  NoPrEP_qaly <- data.frame(data = c(trace_NoPrEP_qaly[,"NoPrEP"],trace_NoPrEP_qaly[,"HIV"],trace_NoPrEP_qaly[,"AIDS"],trace_NoPrEP_qaly[,"Death"]),
                            State = c(rep("Without PrEP",Yearth),rep("HIV",Yearth),rep("AIDS",Yearth),rep("Death",Yearth)),
                            Year = c(rep(1:Yearth,4)))
  NoPrEP_qaly$State <- factor(NoPrEP_qaly$State,levels = c("Without PrEP","HIV","AIDS","Death"))
  
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
  
  ##### Cost #####
  
  PrEP_cost <- data.frame(data = c(trace_PrEP_cost[,"PrEP"],trace_PrEP_cost[,"HIV"],trace_PrEP_cost[,"AIDS"],trace_PrEP_cost[,"Death"]),
                          State = c(rep("With PrEP",Yearth),rep("HIV",Yearth),rep("AIDS",Yearth),rep("Death",Yearth)),
                          Year = c(rep(1:Yearth,4)))
  PrEP_cost$State <- factor(PrEP_cost$State,levels = c("With PrEP","HIV","AIDS","Death"))
  
  cost_NoPrEP <- trace_NoPrEP_cost[,"NoPrEP"]
  cost_PrEP <- trace_PrEP_cost[,"PrEP"]
  cost_HIV <- trace_NoPrEP_cost[,"HIV"] + trace_PrEP_cost[,"HIV"]
  cost_AIDS <- trace_NoPrEP_cost[,"AIDS"] + trace_PrEP_cost[,"AIDS"]
  cost_Death <- trace_NoPrEP_cost[,"Death"] + trace_PrEP_cost[,"Death"]
  
  Whole_cost <- data.frame(data = c(cost_NoPrEP,cost_PrEP,cost_HIV,cost_AIDS,cost_Death),
                           State = c(rep("Without PrEP",Yearth),rep("With PrEP",Yearth),rep("HIV",Yearth),rep("AIDS",Yearth),rep("Death",Yearth)),
                           Year = c(rep(1:Yearth,5)))
  Whole_cost$State <- factor(Whole_cost$State,levels = c("Without PrEP","With PrEP","HIV","AIDS","Death"))
  
  Cost_total <- aggregate(x=Whole_cost[c("data")],by=list(Whole_cost$Year),FUN=sum)
  colnames(Cost_total) <- c("Year","Total")
  
  ##### QALY #####
  PrEP_qaly <- data.frame(data = c(trace_PrEP_qaly[,"PrEP"],trace_PrEP_qaly[,"HIV"],trace_PrEP_qaly[,"AIDS"],trace_PrEP_qaly[,"Death"]),
                          State = c(rep("With PrEP",Yearth),rep("HIV",Yearth),rep("AIDS",Yearth),rep("Death",Yearth)),
                          Year = c(rep(1:Yearth,4)))
  PrEP_qaly$State <- factor(PrEP_qaly$State,levels = c("With PrEP","HIV","AIDS","Death"))
  
  qaly_NoPrEP <- trace_NoPrEP_qaly[,"NoPrEP"]
  qaly_PrEP <- trace_PrEP_qaly[,"PrEP"]
  qaly_HIV <- trace_NoPrEP_qaly[,"HIV"] + trace_PrEP_qaly[,"HIV"]
  qaly_AIDS <- trace_NoPrEP_qaly[,"AIDS"] + trace_PrEP_qaly[,"AIDS"]
  qaly_Death <- trace_NoPrEP_qaly[,"Death"] + trace_PrEP_qaly[,"Death"]
  
  Whole_qaly <- data.frame(data = c(qaly_NoPrEP,qaly_PrEP,qaly_HIV,qaly_AIDS,qaly_Death),
                           State = c(rep("Without PrEP",Yearth),rep("With PrEP",Yearth),rep("HIV",Yearth),rep("AIDS",Yearth),rep("Death",Yearth)),
                           Year = c(rep(1:Yearth,5)))
  Whole_qaly$State <- factor(Whole_qaly$State,levels = c("Without PrEP","With PrEP","HIV","AIDS","Death"))
  
  qaly_total <- aggregate(x=Whole_qaly[c("data")],by=list(Whole_qaly$Year),FUN=sum)
  colnames(qaly_total) <- c("Year","Total")
  
  
  return(list(Cost_total=Cost_total, Whole_cost=Whole_cost,PrEP_cost=PrEP_cost,NoPrEP_cost=NoPrEP_cost,
              qaly_total=qaly_total, Whole_qaly=Whole_qaly,PrEP_qaly=PrEP_qaly,NoPrEP_qaly=NoPrEP_qaly))
  
  
}

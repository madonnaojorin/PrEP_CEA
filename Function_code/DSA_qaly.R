############ Added code  ##############
############ Feb 20 2021
############ Modified Oct 30, 2021, this is for USD

#rm(list = ls())
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

DSA_QALY <- function(numcycles,filename){
  ###### Set parameters for Markov model ######
  dis_rate = 0.02
  Numcycles = numcycles
  cycles <- Numcycles
  Yearth <- cycles +1
  Pop_initial <-1000 
  Percentage_covered <- 0.5 
  Rate <- 110 #
  ###### Calculate initial value
  (NoPrEP_initial <- Pop_initial*(1-Percentage_covered))
  (PrEP_initial <- Pop_initial*Percentage_covered)
  
  ###### Setting parameters for COST #####
  cost <- read.csv("cost_value.csv")
  rownames(cost) <- cost[1:5,1]
  cost[1:5,2:4] <- cost[1:5,2:4]/Rate # This is the only part changed (+ file name)
  cost <- cost[,-1]
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
  
  (Cost_NoPrEP <- Screening_Neg_Yearly)
  (Cost_PrEP <- Screening_Neg_Yearly+PrEP_Yearly)
  (Cost_HIV <- Screening_Pos_Yearly+ART_Yearly)
  (Cost_AIDS <- Screening_Pos_Yearly+ART_Yearly+Hospitalization_Yearly)
  Cost_Death <- 0
  
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
  
  ###### Set parameters for QALY ######
  qaly <- read.csv("QALY_value.csv")
  rownames(qaly) <- qaly[1:5,1]
  qaly <- qaly[,-1]
  
  
  ################################################################
  #####                   No    PrEP                         #####
  ################################################################
  ICER_NoPrEP <- list()
  Target_NoPrEP <- list()
  for(j in 1:100){
    
    ###### Set parameters for QALY ######
    # No PrEP QALY
    QALY_NoPrEP_min = as.numeric(qaly["NoPrEP","MIN"])
    QALY_NoPrEP_max = as.numeric(qaly["NoPrEP","MAX"])
    increment = (QALY_NoPrEP_max - QALY_NoPrEP_min)/99
    QALY_NoPrEP = (QALY_NoPrEP_min + (j-1)*increment)
    
    QALY_PrEP <- as.numeric(qaly["PrEP","Base"])
    QALY_HIV <- as.numeric(qaly["HIV","Base"])
    QALY_AIDS <- as.numeric(qaly["AIDS","Base"])
    QALY_Death <- 0
    
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
    
    ##### ICER calculation ######
    Cost_QALY <- c(sum(PrEP_cost$data),sum(NoPrEP_cost$data),sum(PrEP_qaly$data),sum(NoPrEP_qaly$data))
    
    (ICER_This <- (Cost_QALY[1] - Cost_QALY[2])/(Cost_QALY[3] - Cost_QALY[4]))  
    ICER_NoPrEP[[j]] <- ICER_This
    Target_NoPrEP[[j]] <- QALY_NoPrEP
    
  }
  
  (ICER_NoPrEP_df <- do.call("rbind", ICER_NoPrEP)*10000)
  (Target_NoPrEP_df <- do.call("rbind", Target_NoPrEP))
  (data_NoPrEP <- cbind(ICER_NoPrEP_df,Target_NoPrEP_df))
  colnames(data_NoPrEP) <- c("ICER_NoPrEP","Target_NoPrEP")
  
  ################################################################
  #####                         PrEP                         #####
  ################################################################
  ICER_PrEP<- list()
  Target_PrEP <- list()
  
  for(j in 1:100){
    
    ###### Set parameters for QALY ######
    QALY_NoPrEP <- as.numeric(qaly["NoPrEP","Base"])
    
    # PrEP QALY
    QALY_PrEP_min = as.numeric(qaly["PrEP","MIN"])
    QALY_PrEP_max = as.numeric(qaly["PrEP","MAX"])
    increment = (QALY_PrEP_max - QALY_PrEP_min)/99
    QALY_PrEP = (QALY_PrEP_min + (j-1)*increment)
    
    QALY_HIV <- as.numeric(qaly["HIV","Base"])
    QALY_AIDS <- as.numeric(qaly["AIDS","Base"])
    QALY_Death <- 0
    
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
    
    ##### ICER calculation ######
    Cost_QALY <- c(sum(PrEP_cost$data),sum(NoPrEP_cost$data),sum(PrEP_qaly$data),sum(NoPrEP_qaly$data))
    
    (ICER_This <- (Cost_QALY[1] - Cost_QALY[2])/(Cost_QALY[3] - Cost_QALY[4]))  
    ICER_PrEP[[j]] <- ICER_This
    Target_PrEP[[j]] <- QALY_PrEP
  }
  
  (ICER_PrEP_df <- do.call("rbind", ICER_PrEP)*10000)
  (Target_PrEP_df <- do.call("rbind", Target_PrEP))
  (data_PrEP <- cbind(ICER_PrEP_df,Target_PrEP_df))
  colnames(data_PrEP) <- c("ICER_PrEP","Target_PrEP")
  
  
  ################################################################
  #####                         HIV                          #####
  ################################################################
  ICER_HIV <- list()
  Target_HIV <- list()
  
  for(j in 1:100){
    
    ###### Set parameters for QALY ######
    QALY_NoPrEP <- as.numeric(qaly["NoPrEP","Base"])
    QALY_PrEP <- as.numeric(qaly["PrEP","Base"])
    
    # HIV QALY
    QALY_HIV_min = as.numeric(qaly["HIV","MIN"])
    QALY_HIV_max = as.numeric(qaly["HIV","MAX"])
    increment = (QALY_HIV_max - QALY_HIV_min)/99
    QALY_HIV = (QALY_HIV_min + (j-1)*increment)
    
    QALY_AIDS <- as.numeric(qaly["AIDS","Base"])
    QALY_Death <- 0
    
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
    
    ##### ICER calculation ######
    Cost_QALY <- c(sum(PrEP_cost$data),sum(NoPrEP_cost$data),sum(PrEP_qaly$data),sum(NoPrEP_qaly$data))
    
    (ICER_This <- (Cost_QALY[1] - Cost_QALY[2])/(Cost_QALY[3] - Cost_QALY[4]))  
    ICER_HIV[[j]] <- ICER_This
    Target_HIV[[j]] <- QALY_HIV
    
  }
  
  (ICER_HIV_df <- do.call("rbind", ICER_HIV)*10000)
  (Target_HIV_df <- do.call("rbind", Target_HIV))
  (data_HIV <- cbind(ICER_HIV_df,Target_HIV_df))
  colnames(data_HIV) <- c("ICER_HIV","Target_HIV")
  
  ################################################################
  #####                         AIDS                         #####
  ################################################################
  ICER_AIDS <- list()
  Target_AIDS <- list()
  
  for(j in 1:100){
    
    ###### Set parameters for QALY ######
    QALY_NoPrEP <- as.numeric(qaly["NoPrEP","Base"])
    QALY_PrEP <- as.numeric(qaly["PrEP","Base"])
    QALY_HIV <- as.numeric(qaly["HIV","Base"])
    
    # AIDS QALY
    QALY_AIDS_min = as.numeric(qaly["AIDS","MIN"])
    QALY_AIDS_max = as.numeric(qaly["AIDS","MAX"])
    increment = (QALY_AIDS_max - QALY_AIDS_min)/99
    QALY_AIDS = (QALY_AIDS_min + (j-1)*increment)
    
    QALY_Death <- 0
    
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
    
    ##### ICER calculation ######
    Cost_QALY <- c(sum(PrEP_cost$data),sum(NoPrEP_cost$data),sum(PrEP_qaly$data),sum(NoPrEP_qaly$data))
    
    (ICER_This <- (Cost_QALY[1] - Cost_QALY[2])/(Cost_QALY[3] - Cost_QALY[4]))  
    ICER_AIDS[[j]] <- ICER_This
    Target_AIDS[[j]] <- QALY_AIDS
  }
  
  (ICER_AIDS_df <- do.call("rbind", ICER_AIDS)*10000)
  (Target_AIDS_df <- do.call("rbind", Target_AIDS))
  (data_AIDS <- cbind(ICER_AIDS_df,Target_AIDS_df))
  colnames(data_AIDS) <- c("ICER_AIDS","Target_AIDS")
  
  data = cbind(data_PrEP,data_NoPrEP,data_HIV,data_AIDS)
  dir.create(filename, showWarnings = FALSE)
  write.csv(x = data,paste0(filename,"/QALY_results_",numcycles,".csv"))
}
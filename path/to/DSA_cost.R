############ DSA combined ##############
############ Feb 21 2021, modified Oct 30, 2021
############ Modified Oct 30, 2021, this is for USD

#rm(list = ls())
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

DSA_cost <- function(numcycles,filename){
  ###### Set parameters for Markov model ######
  dis_rate = 0.02
  Numcycles = numcycles
  cycles <- Numcycles
  Yearth <- cycles +1
  Pop_initial <-1000 
  Percentage_covered <- 0.5
  Rate <- 110
  ###### Calculate initial value
  (NoPrEP_initial <- Pop_initial*(1-Percentage_covered))
  (PrEP_initial <- Pop_initial*Percentage_covered)
  
  
  ###### Set parameters for cost ######
  cost <- read.csv("Cost_value.csv")
  cost[1:5,2:4] <- cost[1:5,2:4]/Rate # This is the only few parts changed (+ file name, ICER*10000)
  rownames(cost) <- cost[1:5,1]
  cost <- cost[,-1]
  
  ###### Set parameters for QALY ######
  QALY_data <- read.csv("QALY_value.csv")
  QALY_PrEP <- QALY_data[1,"Base"]
  QALY_NoPrEP <- QALY_data[2,"Base"]
  QALY_HIV <- QALY_data[3,"Base"]
  QALY_AIDS <- QALY_data[4,"Base"]
  QALY_Death <- QALY_data[5,"Base"]
  
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
  
  
  ################################################################
  #####                         PrEP                         #####
  ################################################################
  ICER_PrEP <- list()
  Target_PrEP <- list()
  
  for(j in 1:100){
    
    ###### Setting parameters for COST #####
    ### Cost for PrEP ###
    #(PrEP_pill = as.numeric(cost["PrEP","Base"]))
    PrEP_min = as.numeric(cost["PrEP","MIN"])
    PrEP_max = as.numeric(cost["PrEP","MAX"])
    increment = (PrEP_max - PrEP_min)/99
    PrEP_pill = (PrEP_min + (j-1)*increment)/10000
    PrEP_Monthly = PrEP_pill*30
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
    Target_PrEP[[j]] <- PrEP_pill
    
  }
  
  (ICER_PrEP_df <- do.call("rbind", ICER_PrEP)*10000)
  (Target_PrEP_df <- do.call("rbind", Target_PrEP))
  (data_PrEP <- cbind(ICER_PrEP_df,Target_PrEP_df))
  colnames(data_PrEP) <- c("ICER_PrEP","Target_PrEP")
  
  
  ################################################################
  #####                Pos screening                         #####
  ################################################################
  ICER_pos <- list()
  Target_pos <- list()
  for(j in 1:100){
    
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
    #Screening_Pos_Monthly = as.numeric(cost["Screening_Positive","Base"])/10000
    Positive_min = as.numeric(cost["Screening_Positive","MIN"])
    Positive_max = as.numeric(cost["Screening_Positive","MAX"])
    increment = (Positive_max - Positive_min)/99
    Screening_Pos_Monthly = (Positive_min + (j-1)*increment)/10000
    
    (Screening_Neg_Yearly = Screening_Neg_Monthly*4)
    # once per half a year
    (Screening_Pos_Yearly = Screening_Pos_Monthly*12)
    # once a month
    (Cost_NoPrEP <- Screening_Neg_Yearly); 
    (Cost_PrEP <- Screening_Neg_Yearly+PrEP_Yearly);
    (Cost_HIV <- Screening_Pos_Yearly+ART_Yearly);
    (Cost_AIDS <- Screening_Pos_Yearly+ART_Yearly+Hospitalization_Yearly); 
    Cost_Death <- 0
    
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
    ICER_pos[[j]] <- ICER_This
    Target_pos[[j]] <- Screening_Pos_Monthly
    
  }
  
  (ICER_pos_df <- do.call("rbind", ICER_pos)*10000)
  (Target_pos_df <- do.call("rbind", Target_pos))
  (data_pos <- cbind(ICER_pos_df,Target_pos_df))
  colnames(data_pos) <- c("ICER_pos","Target_pos")
  
  ################################################################
  #####                Neg screening                         #####
  ################################################################
  
  ICER_neg <- list()
  Target_neg <- list()
  for(j in 1:100){
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
    #Screening_Neg_Monthly = as.numeric(cost["Screening_Negative","Base"])/10000
    Negative_min = as.numeric(cost["Screening_Negative","MIN"])
    Negative_max = as.numeric(cost["Screening_Negative","MAX"])
    increment = (Negative_max - Negative_min)/99
    Screening_Neg_Monthly = (Negative_min + (j-1)*increment)/10000
    
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
    ICER_neg[[j]] <- ICER_This
    Target_neg[[j]] <- Screening_Neg_Monthly
    
  }
  
  (ICER_neg_df <- do.call("rbind", ICER_neg)*10000)
  (Target_neg_df <- do.call("rbind", Target_neg))
  (data_neg <- cbind(ICER_neg_df,Target_neg_df))
  colnames(data_neg) <- c("ICER_neg","Target_neg")
  
  ################################################################
  #####                Hospitalization                       #####
  ################################################################
  ICER_hospitalization <- list()
  Target_hospitalization <- list()
  for(j in 1:100){
    
    ###### Setting parameters for COST #####
    ### Cost for PrEP ###
    (PrEP_pill = as.numeric(cost["PrEP","Base"]))
    PrEP_Monthly = PrEP_pill*30/10000
    (PrEP_Yearly = PrEP_Monthly*12)
    
    ### Cost for ART ###
    (ART_Monthly = as.numeric(cost["ART","Base"])/10000)
    (ART_Yearly = ART_Monthly*12)
    
    ### Hospitalization ###
    #(Hospitalization_Monthly = as.numeric(cost["Hospitalization","Base"])/10000)
    Hospitalization_min = as.numeric(cost["Hospitalization","MIN"])
    Hospitalization_max = as.numeric(cost["Hospitalization","MAX"])
    increment = (Hospitalization_max - Hospitalization_min)/99
    Hospitalization_Monthly = (Hospitalization_min + (j-1)*increment)/10000
    
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
    ICER_hospitalization[[j]] <- ICER_This
    Target_hospitalization[[j]] <- Hospitalization_Yearly
    
  }
  
  (ICER_hospitalization_df <- do.call("rbind", ICER_hospitalization)*10000)
  (Target_hospitalization_df <- do.call("rbind", Target_hospitalization))
  (data_hospitalization <- cbind(ICER_hospitalization_df,Target_hospitalization_df))
  colnames(data_hospitalization) <- c("ICER_hospitalization","Target_hospitalization")
  
  ################################################################
  #####                ART                                   #####
  ################################################################
  ICER_ART <- list()
  Target_ART <- list()
  for(j in 1:100){
    ###### Setting parameters for COST #####
    ### Cost for PrEP ###
    (PrEP_pill = as.numeric(cost["PrEP","Base"]))
    PrEP_Monthly = PrEP_pill*30/10000
    (PrEP_Yearly = PrEP_Monthly*12)
    
    ### Cost for ART ###
    #(ART_Monthly = as.numeric(cost["ART","Base"])/10000)
    ### Cost for ART ###
    ART_min = as.numeric(cost["ART","MIN"])
    ART_max = as.numeric(cost["ART","MAX"])
    increment = (ART_max - ART_min)/99
    ART_Monthly = (ART_min + (j-1)*increment)/10000
    
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
    ICER_ART[[j]] <- ICER_This
    Target_ART[[j]] <- ART_Yearly
    
  }
  
  (ICER_ART_df <- do.call("rbind", ICER_ART)*10000)
  (Target_ART_df <- do.call("rbind", Target_ART))
  (data_ART <- cbind(ICER_ART_df,Target_ART_df))
  colnames(data_ART) <- c("ICER_ART","Target_ART")
  
  
  data = cbind(data_PrEP,data_neg,data_pos,data_hospitalization,data_ART)
  dir.create(filename, showWarnings = FALSE)
  write.csv(x = data,paste0(filename,"/cost_results_",numcycles,".csv"))
}

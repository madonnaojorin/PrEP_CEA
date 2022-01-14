############ Final code after checking ##############
############ Jan 16th: updated all the data for cost and trandition probability 
############ Modified: Nov 27 2021 

Incidence_function <- function(Coverage, Initial, NumYears){

  ###### Set parameters for Markov model ######
  cycles <- NumYears
  Yearth <- cycles +1  
  Pop_initial <-Initial
  Percentage_covered <- Coverage
  dis_rate <- 0.2
  
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
  
  simulation = 1000
  set.seed(simulation)
  ###### Set parameters for WITHOUT PrEP ######
  #y1. No PrEP
  #y2. HIV positive
  #y3. AIDS
  #y4. death
  p_NoPrEP_HIV <- NoPrEP_HIV
  p_NoPrEP_AIDS <- NoPrEP_AIDS
  p_NoPrEP_death <- NoPrEP_Death
  p_HIV_AIDS <- HIV_AIDS
  p_HIV_death <- HIV_Death
  p_AIDS_death <- AIDS_Death
  
  mat_1_NoPrEP <- matrix(NA, nrow = Yearth, ncol = simulation)
  mat_2_NoPrEP <- matrix(NA, nrow = Yearth, ncol = simulation)
  mat_3_NoPrEP <- matrix(NA, nrow = Yearth, ncol = simulation)
  mat_4_NoPrEP <- matrix(NA, nrow = Yearth, ncol = simulation)
  
  ###### Set parameters for PrEP ######
  #y1. PrEP
  #y2. HIV positive
  #y3. AIDS
  #y4. death
  p_PrEP_HIV <- PrEP_HIV
  p_PrEP_AIDS <- PrEP_AIDS
  p_PrEP_death <- PrEP_Death
  
  mat_1_PrEP <- matrix(NA, nrow = Yearth, ncol = simulation)
  mat_2_PrEP <- matrix(NA, nrow = Yearth, ncol = simulation)
  mat_3_PrEP <- matrix(NA, nrow = Yearth, ncol = simulation)
  mat_4_PrEP <- matrix(NA, nrow = Yearth, ncol = simulation)
  

  for (i in 1:simulation) {
    #x0 <- 1000
    (NoPrEP_initial <- Pop_initial * (1 - Percentage_covered))

    (PrEP_initial <- Pop_initial * Percentage_covered)

    #y0 <- 0
    y_2_0_NoPrEP <- 0
    y_3_0_NoPrEP <- 0
    y_4_0_NoPrEP <- 0
    
    y_2_0_PrEP <- 0
    y_3_0_PrEP <- 0
    y_4_0_PrEP <- 0
    
    #x <- x0
    y_1_NoPrEP <- NoPrEP_initial
    
    y_1_PrEP <- PrEP_initial
    
    #y <- y0
    y_2_NoPrEP <- y_2_0_NoPrEP
    y_3_NoPrEP <- y_3_0_NoPrEP
    y_4_NoPrEP <- y_4_0_NoPrEP
    
    y_2_PrEP <- y_2_0_PrEP
    y_3_PrEP <- y_3_0_PrEP
    y_4_PrEP <- y_4_0_PrEP
    
    #vec_x <- x
    vec_1_NoPrEP <- y_1_NoPrEP
    vec_2_NoPrEP <- y_2_NoPrEP
    vec_3_NoPrEP <- y_3_NoPrEP
    vec_4_NoPrEP <- y_4_NoPrEP
    
    vec_1_PrEP <- y_1_PrEP
    vec_2_PrEP <- y_2_PrEP
    vec_3_PrEP <- y_3_PrEP
    vec_4_PrEP <- y_4_PrEP
    
    for (j in 1:NumYears) {
      ### From No PrEP compartment to NoPrEP(Stay), HIV, AIDS, death
      random_NoPrEP <- runif(y_1_NoPrEP)
      (num_x2_NoPrEP <-
          length(random_NoPrEP[(0 < random_NoPrEP) &
                                 (random_NoPrEP < p_NoPrEP_HIV)]))
      (num_x3_NoPrEP <-
          length(random_NoPrEP[(p_NoPrEP_HIV < random_NoPrEP) &
                                 (random_NoPrEP < (p_NoPrEP_HIV + p_NoPrEP_AIDS))]))
      (num_x4_NoPrEP <-
          length(random_NoPrEP[((p_NoPrEP_HIV + p_NoPrEP_AIDS) < random_NoPrEP) &
                                 (random_NoPrEP < (p_NoPrEP_HIV + p_NoPrEP_AIDS + p_NoPrEP_death))]))
      
      ### From PrEP compartment to PrEP(Stay), HIV, AIDS, death
      random_PrEP <- runif(y_1_PrEP)
      (num_x2_PrEP <-
          length(random_PrEP[(0 < random_PrEP) &
                               (random_PrEP < p_PrEP_HIV)]))
      (num_x3_PrEP <-
          length(random_PrEP[(p_PrEP_HIV < random_PrEP) &
                               (random_PrEP < (p_PrEP_HIV + p_PrEP_AIDS))]))
      (num_x4_PrEP <-
          length(random_PrEP[((p_PrEP_HIV + p_PrEP_AIDS) < random_PrEP) &
                               (random_PrEP < (p_PrEP_HIV + p_PrEP_AIDS + p_PrEP_death))]))
      
      
      ### From HIV compartment to HIV(Stay), AIDS, death
      random_HIV <- runif(y_2_NoPrEP)
      (num_x3_HIV <-
          length(random_HIV[(0 < random_HIV) &
                              (random_HIV < p_HIV_AIDS)]))
      (num_x4_HIV <-
          length(random_HIV[(p_HIV_AIDS < random_HIV) &
                              (random_HIV < (p_HIV_AIDS + p_HIV_death))]))
      
      ### From HIV compartment to HIV(Stay), AIDS, death
      random_HIV <- runif(y_2_PrEP)
      (num_x3_HIV <-
          length(random_HIV[(0 < random_HIV) &
                              (random_HIV < p_HIV_AIDS)]))
      (num_x4_HIV <-
          length(random_HIV[(p_HIV_AIDS < random_HIV) &
                              (random_HIV < (p_HIV_AIDS + p_HIV_death))]))
      
      
      ### From AIDS compartment to AIDS(Stay), death
      random_AIDS <- runif(y_3_NoPrEP)
      (num_x4_AIDS <-
          length(random_AIDS[p_AIDS_death < random_AIDS]))
      
      ### From AIDS compartment to AIDS(Stay), death
      random_AIDS <- runif(y_3_PrEP)
      (num_x4_AIDS <-
          length(random_AIDS[p_AIDS_death < random_AIDS]))
      
      
      ### Number of people for each compartment
      
      (yy_1_NoPrEP <-
          y_1_NoPrEP - (num_x2_NoPrEP + num_x3_NoPrEP + num_x4_NoPrEP))
      (yy_2_NoPrEP <- y_2_NoPrEP + num_x2_NoPrEP)
      (yy_3_NoPrEP <- y_3_NoPrEP + num_x3_NoPrEP + num_x3_HIV)
      (yy_4_NoPrEP <-
          y_4_NoPrEP + num_x4_NoPrEP + num_x4_HIV + num_x4_AIDS)
      
      
      (vec_1_NoPrEP <- c(vec_1_NoPrEP, yy_1_NoPrEP))
      (vec_2_NoPrEP <- c(vec_2_NoPrEP, yy_2_NoPrEP))
      (vec_3_NoPrEP <- c(vec_3_NoPrEP, yy_3_NoPrEP))
      (vec_4_NoPrEP <- c(vec_4_NoPrEP, yy_4_NoPrEP))
      
      (y_1_NoPrEP <- yy_1_NoPrEP)
      (y_2_NoPrEP <- yy_2_NoPrEP)
      (y_3_NoPrEP <- yy_3_NoPrEP)
      (y_4_NoPrEP <- yy_4_NoPrEP)
      
      ### Number of people for each compartment
      
      (yy_1_PrEP <-
          y_1_PrEP - (num_x2_PrEP + num_x3_PrEP + num_x4_PrEP))
      (yy_2_PrEP <- y_2_PrEP + num_x2_PrEP)
      (yy_3_PrEP <- y_3_PrEP + num_x3_PrEP + num_x3_HIV)
      (yy_4_PrEP <-
          y_4_PrEP + num_x4_PrEP + num_x4_HIV + num_x4_AIDS)
      
      
      (vec_1_PrEP <- c(vec_1_PrEP, yy_1_PrEP))
      (vec_2_PrEP <- c(vec_2_PrEP, yy_2_PrEP))
      (vec_3_PrEP <- c(vec_3_PrEP, yy_3_PrEP))
      (vec_4_PrEP <- c(vec_4_PrEP, yy_4_PrEP))
      
      (y_1_PrEP <- yy_1_PrEP)
      (y_2_PrEP <- yy_2_PrEP)
      (y_3_PrEP <- yy_3_PrEP)
      (y_4_PrEP <- yy_4_PrEP)
    }
    mat_1_NoPrEP[, i] <- vec_1_NoPrEP
    mat_2_NoPrEP[, i] <- vec_2_NoPrEP
    mat_3_NoPrEP[, i] <- vec_3_NoPrEP
    mat_4_NoPrEP[, i] <- vec_4_NoPrEP
    
    mat_1_PrEP[, i] <- vec_1_PrEP
    mat_2_PrEP[, i] <- vec_2_PrEP
    mat_3_PrEP[, i] <- vec_3_PrEP
    mat_4_PrEP[, i] <- vec_4_PrEP
  }
  
  ###### Calculate mean and 95% CI
  res_mean_NoPrEP <-
    list()
  res_low_NoPrEP <-
    list()
  res_high_NoPrEP <- list()
  res_summary_NoPrEP <- list()
  for (k in 1:Yearth) {
      (res_mean_NoPrEP[[k]] <- mean((mat_2_NoPrEP[k,]+mat_3_NoPrEP[k,])))
      (res_low_NoPrEP[[k]] <-  quantile((mat_2_NoPrEP[k,]+mat_3_NoPrEP[k,]), c(0.025)))
      (res_high_NoPrEP[[k]] <- quantile((mat_2_NoPrEP[k,]+mat_3_NoPrEP[k,]), c(0.975)))
      (res_summary_NoPrEP[[k]] <-
          c(res_mean_NoPrEP[[k]], res_low_NoPrEP[[k]], res_high_NoPrEP[[k]]))
  }
  
  NoPrEP_infected <- do.call("rbind", res_summary_NoPrEP)
  colnames(NoPrEP_infected) <- c("Mean", "Low", "High")
  NoPrEP_infected_data <-
    data.frame(
      NoPrEP_infected,
      Year = 1:Yearth,
      State = rep("No PrPEP", Yearth),
      Coverage = rep(paste0(Coverage*10, "0%"), Yearth)
    )
  
  res_mean_PrEP <-
    list()
  res_low_PrEP <-
    list()
  res_high_PrEP <- list()
  res_summary_PrEP <- list()
  for (k in 1:Yearth) {
      (res_mean_PrEP[[k]] <- mean((mat_2_PrEP[k,]+mat_2_PrEP[k,])))
      (res_low_PrEP[[k]] <-  quantile((mat_2_PrEP[k,]+mat_2_PrEP[k,]), c(0.025)))
      (res_high_PrEP[[k]] <- quantile((mat_2_PrEP[k,]+mat_2_PrEP[k,]), c(0.975)))
      (res_summary_PrEP[[k]] <-
          c(res_mean_PrEP[[k]], res_low_PrEP[[k]], res_high_PrEP[[k]]))
  }
  
  PrEP_infected <- do.call("rbind", res_summary_PrEP)
  colnames(PrEP_infected) <- c("Mean", "Low", "High")
  PrEP_infected_data <-
    data.frame(
      PrEP_infected,
      Year = 1:Yearth,
      State = rep("PrPEP", Yearth),
      Coverage = rep(paste0(Coverage*10, "0%"), Yearth)
    )
  
  PrEP_simulated_data <-
    data.frame(
      data = as.vector(mat_2_PrEP+mat_3_PrEP),
      Year = rep(1:Yearth, simulation),
      Trial = as.character(rep(1:simulation, each = Yearth)),
      Coverage = rep(paste0(Coverage*10, "0%"), Yearth)
    )
  
  NoPrEP_simulated_data <-
    data.frame(
      data = as.vector(mat_2_NoPrEP+mat_3_NoPrEP),
      Year = rep(1:Yearth, simulation),
      Trial = as.character(rep(1:simulation, each = Yearth)),
      Coverage = rep(paste0(Coverage*10, "0%"), Yearth)
    )
  return(list(NoPrEP_infected_data=NoPrEP_infected_data, PrEP_infected_data=PrEP_infected_data,NoPrEP_simulated_data=NoPrEP_simulated_data,PrEP_simulated_data=PrEP_simulated_data))
}

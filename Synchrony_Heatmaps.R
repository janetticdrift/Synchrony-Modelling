#####################################
#This script produces figures that show the range in variance ratios of resident communities,
#and their subsequent synchrony dynamics after invasion of a strong and weak invader.
#####################################

setwd("~/Dropbox/Synchrony/Modeling Code")
library(codyn)
runs <- 2000

#load species parameters
source("Final Code/Species_Parameters_Source.R")

# ----------------------------------------------------------------------------------------
###Figure 1: Resident communities' starting variance ratios----

#Set outputs of interest
cv_total_biomass <-matrix(NA,nrow=length(env_condition), ncol=length(beta_range))
VR <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))

for (x in 1:length(env_condition)) {
  for (y in 1:length(beta_range)) {
    
    env <- env_condition[x]
    beta <- beta_range[y]
    
    VR_current <- rep(NA, runs)
    CV_current <- rep(NA, runs)
    
    for (z in 1:runs) {
      
      #State the variables
      time <- 500 
      N1 <- rep(NA, time) 
      N2 <- rep(NA, time)
      
      N1[1] <- 50 
      N2[1] <- 50
      
      #Set parameters; growth rates, carrying capacities, and demographic effects are from the parameter source code
      beta12 <- beta   #effect of species 2 on species 1
      beta21 <- beta   #effect of species 1 on species 2
      sigmaE1 <- -env #env effect on species 1
      sigmaE2 <- -env #env effect on species 2
      miuE <- rnorm(time, mean = 0, sd = 1) #environmental timeseries variation
      miuD1 <- rnorm(time, mean = 0, sd = 1) #dem timeseries for species 1
      miuD2 <- rnorm(time, mean = 0, sd = 1) #dem timeseries for species 2
      
      #Create the model
      for (t in 1:(time-1)) {
        #calculate population sizes for species 1
        N1[t+1] <- N1[t]*exp(r1*(1-(N1[t]/K1) - (beta12*N2[t]/K2))+
                               (sigmaE1*miuE[t])+(sigmaD1*miuD1[t])/sqrt(N1[t]))
        
        #calculate population sizes for species 2
        N2[t+1] <- N2[t]*exp(r2*(1-(N2[t]/K2) - (beta21*N1[t]/K1))+
                               (sigmaE2*miuE[t])+(sigmaD2*miuD2[t])/sqrt(N2[t]))
      }
      counter <- 1 
      if(N1[counter] < 0) {
        N1[counter] <- 0
      }
      if(N2[counter] < 0) {
        N2[counter] <- 0
      }
      
      #Create graph showing biomass through time for each species 
      #and then species through time
      total_biomass <- N1 + N2
      
      #Current coefficient of variation
      CV_current[z] <- sd(total_biomass)/mean(total_biomass)
      
      #calculate variance ratio
      our_data <- matrix(NA, nrow = 2*time, ncol = 3)
      our_data[,1] <- rep(c(1, 2), time)
      our_data[,2] <- rep(seq(1:time), each = 2)
      our_data[,3] <- as.vector(rbind(N1,N2))
      our_data_v2 <- data.frame(our_data)
      colnames(our_data_v2) <- c("species", "time", "abundance")
      
      VR_temp <- variance_ratio(our_data_v2, time.var = "time", species.var = "species", 
                                abundance.var = "abundance", bootnumber = 1)
      VR_current[z] <- VR_temp$VR
    }
    VR[x,y] <- mean(VR_current)
    cv_total_biomass[x,y] <- mean(CV_current)
  }
  print(x/length(env_condition))
}

#Packages for creating plots
library(ggplot2)
library(reshape2)

M1 <- melt(VR) 

ggplot(M1, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#008080", high ="#ca562c", mid = "#f6edbd", midpoint = 1, limit = c(0,2)) +
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), y= expression(paste("Strength of Competititon (", beta, ")")), title = "Resident Communities' Variance Ratios", fill="VR") +
   scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95"))


# ----------------------------------------------------------------------------------------
###Figure 2: Resistance of resident communities to invader growth of a weak and strong invader

#Figure 2.a: Weak Invader-----
#Set outputs of interest
cv_total_biomass <-matrix(NA,nrow=length(env_condition), ncol = length(beta_range))
avg_lambda_poor <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))
success_lambda_poor <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))
burn_in <- 100
time <- 500

for (x in 1:length(env_condition)) {
  for (y in 1:length(beta_range)) {
    
    env <- env_condition[x]
    beta <- beta_range[y]
    
    lambda_poor_a <- rep(NA, runs) 
    lambda_poor <- matrix(NA, nrow = time-burn_in, ncol = runs)
    
    for (z in 1:runs) {
      
      #State the variables
      N1 <- rep(NA, time) 
      N2 <- rep(NA, time)
      N3 <- rep(NA, (time-burn_in))
      
      N1[1] <- 50 
      N2[1] <- 50
      N3[1] <- 1 
      
      invader_abund_poor <- 1 
      
      #Set parameters; growth rates, carrying capacities, and demographic effects are from the parameter source code
      beta12 <- beta  
      beta21 <- beta  
      sigmaE1 <- -env 
      sigmaE2 <- -env
      miuE <- rnorm(time, mean = 0, sd = 1) #environmental timeseries variation
      miuD1 <- rnorm(time, mean = 0, sd = 1) #dem timeseries for species 1
      miuD2 <- rnorm(time, mean = 0, sd = 1) #dem timeseries for species 2
      miuD3 <- rnorm((time-burn_in), mean = 0, sd = 1)
      
      
      #Create the model for resident community
      for (t in 1:(time-1)) {
        
        if(N1[t] < 0) {
          N1[t] <- 0
        }
        if(N2[t] < 0) {
          N2[t] <- 0
        }
        
        #calculate population size for species 1
        N1[t+1] <- N1[t]*exp(r1*(1-(N1[t]/K1) - (beta12*N2[t]/K2))+
                               (sigmaE1*miuE[t])+(sigmaD1*miuD1[t])/sqrt(N1[t]))
        
        #calculate population size for species 2
        N2[t+1] <- N2[t]*exp(r2*(1-(N2[t]/K2) - (beta21*N1[t]/K1))+
                               (sigmaE2*miuE[t])+(sigmaD2*miuD2[t])/sqrt(N2[t]))
      }
      counter <- 1
      
      #Model effects of poor invasion
      for (t in burn_in:(time-1)) {
        
        N3[counter] <- invader_abund_poor*exp(r3*(1-(invader_abund_poor/K3) - (beta31*N1[t]/K1) - (beta32*N2[t]/K2))
                                              + (sigmaE3*miuE[t])+(sigmaD3*miuD3[counter])/sqrt(invader_abund_poor))
        if(N3[counter] < 0) {
          N3[counter] <- 0
        }
        
        #Calculate invader's growth rate
        lambda_poor[counter, z] <- log(N3[counter]/invader_abund_poor)
        
        counter <- counter + 1
      }
      
      lambda_poor_a[z] <- mean(lambda_poor[,z])
    }
    avg_lambda_poor[x,y] <- mean(lambda_poor_a)
    success_lambda_poor[x,y] <- length(which(lambda_poor>0))/((time-burn_in)*runs)
    
  }
  print(x/length(env_condition))
}

#Successful Invasion of Poor Invader Plot
M2 <- melt(success_lambda_poor)

ggplot(M2, aes(x=Var1, y=Var2, fill=value*100)) + 
  geom_tile() + 
  scale_fill_distiller(palette = "RdBu", limits = c(min(success_lambda_poor)*100, 100)) +
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), y= expression(paste("Strength of Competititon (", beta, ")")), title = "Growth Success of Weak Invader", fill="Growth
Success (%)") +
  theme(axis.text = element_text( size = 12)) + 
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95"))

#Figure 2.b: Strong Invader-----
#Set outputs of interest
cv_total_biomass <-matrix(NA,nrow=length(env_condition), ncol = length(beta_range))
avg_lambda_good <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))
success_lambda_good <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))
burn_in <- 100
time <- 500

for (x in 1:length(env_condition)) {
  for (y in 1:length(beta_range)) {
    
    env <- env_condition[x]
    beta <- beta_range[y]
    
    lambda_good_a <- rep(NA, runs) 
    lambda_good <- matrix(NA, nrow = time-burn_in, ncol = runs)
    
    for (z in 1:runs) {
      
      #State the variables
      N1 <- rep(NA, time) 
      N2 <- rep(NA, time)
      N4 <- rep(NA, (time-burn_in))
      
      N1[1] <- 50 
      N2[1] <- 50
      N4[1] <- 1 
      
      invader_abund_good <- 1 
      
      #Set parameters; growth rates, carrying capacities, and demographic effects are from the parameter source code
      beta12 <- beta   
      beta21 <- beta  
      sigmaE1 <- -env 
      sigmaE2 <- -env 
      miuE <- rnorm(time, mean = 0, sd = 1) #environmental timeseries variation
      miuD1 <- rnorm(time, mean = 0, sd = 1) #dem timeseries for species 1
      miuD2 <- rnorm(time, mean = 0, sd = 1) #dem timeseries for species 2
      miuD4 <- rnorm((time-burn_in), mean = 0, sd = 1)
      
      
      #Create the model for resident community
      for (t in 1:(time-1)) {
        
        if(N1[t] < 0) {
          N1[t] <- 0
        }
        if(N2[t] < 0) {
          N2[t] <- 0
        }
        
        #calculate population size for species 1
        if(N1[t] > 0) {
          N1[t+1] <- N1[t]*exp(r1*(1-(N1[t]/K1) - (beta12*N2[t]/K2))+
                                 (sigmaE1*miuE[t])+(sigmaD1*miuD1[t])/sqrt(N1[t]))
        } else {
          N1[t+1] <- 0
        }
        
        #calculate population size for species 2
        if (N2[t] > 0) {
          N2[t+1] <- N2[t]*exp(r2*(1-(N2[t]/K2) - (beta21*N1[t]/K1))+
                                 (sigmaE2*miuE[t])+(sigmaD2*miuD2[t])/sqrt(N2[t]))
        } else {
          N2[t+1] <- 0
        }
        
      }
      counter <- 1
      
      #Model effects of strong invasion
      for (t in burn_in:(time-1)) {
        N4[counter] <- invader_abund_good*exp(r4*(1-(invader_abund_good/K4) - (beta41*N1[t]/K1) - (beta42*N2[t]/K2))
                                              + (sigmaE4*miuE[t])+(sigmaD4*miuD4[counter])/sqrt(invader_abund_good))
        if(N4[counter] < 0) {
          N4[counter] <- 0
        }
        
        #Calculate invader's growth rate
        lambda_good[counter, z] <- log(N4[counter]/invader_abund_good)
        
        counter <- counter + 1
      }
      lambda_good_a[z] <- mean(lambda_good[,z])
    }
    avg_lambda_good[x,y] <- mean(lambda_good_a)
    success_lambda_good[x,y] <- length(which(lambda_good>0))/((time-burn_in)*runs)
    
  }
  print(x/length(env_condition))
}

#Successful Invasion of Good Invader Plot
M3 <- melt(success_lambda_good)

ggplot(M3, aes(x=Var1, y=Var2, fill=value*100)) +
  geom_tile() + 
  scale_fill_distiller(palette = "RdBu", limits = c(min(success_lambda_good)*100, 100)) +
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), y= expression(paste("Strength of Competititon (", beta, ")")), title = "Growth Success of Strong Invader", fill="Growth
Success (%)") +
  theme(axis.text = element_text( size = 12)) + 
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95"))

# ----------------------------------------------------------------------------------------
###Figure 3: Resilience of resident communities to established invaders

#Figure 3.a: Weak Invader in 2 Species Community----
burn_in <- 200
timeseries <- 100
invasion_success <- burn_in + 10

VR_pre_poor_2 <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))
VR_post_poor_2 <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))

for (x in 1:length(env_condition)) {
  for (y in 1:length(beta_range)) {
    
    env <- env_condition[x]
    beta <- beta_range[y]
    
    VR_current_pre_poor_2 <- rep(NA, runs)
    VR_current_post_poor_2 <- rep(NA, runs)
    
    for (z in 1:runs) {
      
      #State the variables
      time <- 300 
      N1 <- rep(NA, time) 
      N2 <- rep(NA, time)
      N3 <- rep(NA, time)
      
      N1[1] <- 50
      N2[1] <- 50
      N3[1] <- 0
      
      #Set parameters
      beta12 <- beta   
      beta21 <- beta   
      sigmaE1 <- -env 
      sigmaE2 <- -env 
      miuE <- rnorm(time, mean = 0, sd = 1) 
      miuD1 <- rnorm(time, mean = 0, sd = 1)
      miuD2 <- rnorm(time, mean = 0, sd = 1)
      miuD3 <- rnorm(time, mean = 0, sd = 1)
      
      #Create the model
      for (t in 1:(time-1)) {
        #Introduce invader at time burn_in at abundance of 1
        if(t == burn_in) {
          N3[t] <- 1
        }
        
        #calculate population sizes for species 1
        N1[t+1] <- N1[t]*exp(r1*(1-(N1[t]/K1) - (beta12*N2[t]/K2)) - (beta13*N3[t]/K3) +
                               (sigmaE1*miuE[t])+(sigmaD1*miuD1[t])/sqrt(N1[t]))
        
        #calculate population sizes for species 2
        N2[t+1] <- N2[t]*exp(r2*(1-(N2[t]/K2) - (beta21*N1[t]/K1)) - (beta23*N3[t]/K3) +
                               (sigmaE2*miuE[t])+(sigmaD2*miuD2[t])/sqrt(N2[t]))
        
        #calculate population sizes for species 3
        N3[t+1] <- N3[t]*exp(r3*(1-(N3[t]/K3) - (beta31*N1[t]/K1)) - (beta32*N2[t]/K2) +
                               (sigmaE3*miuE[t]) + (sigmaD3*miuD3[t])/sqrt(N3[t]))
        
        N1[is.na(N1)] <- 0
        N2[is.na(N2)] <- 0
        N3[is.na(N3)] <- 0
        
        if(N1[t+1] < 1) {
          N1[t+1] <- 0
        }
        if(N2[t+1] < 1) {
          N2[t+1] <- 0
        }
        if(N3[t+1] < 1) {
          N3[t+1] <- 0
        }
      }  
      
      #calculate variance ratio, pre-invasion
      # create an empty matrix
      our_data_pre_poor_2 <- matrix(NA, nrow = 2*(burn_in + timeseries), ncol = 3)
      # fill in the first column with our species ID
      our_data_pre_poor_2[,1] <- rep(c(1, 2), (burn_in + timeseries))
      # fill in the second column with the timepoint
      our_data_pre_poor_2[,2] <- rep(seq(1:(burn_in + timeseries)), each = 2)
      # fill in our third colum with the abundances from our simulation model
      our_data_pre_poor_2[,3] <- as.vector(rbind(N1[(burn_in-timeseries):(burn_in-1)], 
                                               N2[(burn_in-timeseries):(burn_in-1)]))
      
      our_data_pre_poor_v2_2 <- data.frame(our_data_pre_poor_2)
      colnames(our_data_pre_poor_v2_2) <- c("species", "time", "abundance")
      
      # calculate the variance ratio
      VR_temp_pre_poor_2 <- variance_ratio(our_data_pre_poor_v2_2, time.var = "time", species.var = "species", 
                                         abundance.var = "abundance", bootnumber = 1)
      VR_current_pre_poor_2[z] <- VR_temp_pre_poor_2$VR
      
      #calculate variance ratio, post-invasion
      # create an empty matrix
      our_data_post_poor_2 <- matrix(NA, nrow = 2*timeseries, ncol = 3)
      # fill in the first column with our species ID
      our_data_post_poor_2[,1] <- rep(c(1, 2), timeseries)
      # fill in the second column with the timepoint
      our_data_post_poor_2[,2] <- rep(seq(1:timeseries), each = 2)
      # fill in our third colum with the abundances from our simulation model
      our_data_post_poor_2[,3] <- as.vector(rbind(N1[(burn_in+1):time], 
                                                N2[(burn_in+1):time]))

      our_data_post_poor_v2_2 <- data.frame(our_data_post_poor_2)
      colnames(our_data_post_poor_v2_2) <- c("species", "time", "abundance")
      
      if(N3[invasion_success] > 1) {
        
        # calculate the variance ratio
        VR_temp_post_poor_2 <- variance_ratio(our_data_post_poor_v2_2, time.var = "time", 
                                            species.var = "species", abundance.var = "abundance", bootnumber = 1)
        VR_current_post_poor_2[z] <- VR_temp_post_poor_2$VR
      }
    }
    VR_pre_poor_2[x,y] <- mean(VR_current_pre_poor_2)
    VR_post_poor_2[x,y] <- mean(VR_current_post_poor_2, na.rm = TRUE)
  }
  print(x/length(env_condition))
}

#Proportion of VR > 1 in pre invasion communities
mean(VR_pre_poor_2>1, na.rm = TRUE)
#Proportion of VR > 1 in post invasion communities
mean(VR_post_poor_2>1, na.rm = TRUE)

M4 <- melt(VR_post_poor_2)

ggplot(M4, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#008080", high ="#ca562c", mid = "#f6edbd", midpoint = 1, limit = c(0,2)) + 
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), y= expression(paste("Strength of Competititon (", beta, ")")), title = "Two-Species with Weak Invader VRs", fill="VR") +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95"))

#Figure 3.b: Strong Invader in 2 Species Community----
burn_in <- 200
timeseries <- 100
invasion_success <- burn_in + 10

VR_pre_good_2 <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))
VR_post_good_2 <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))

for (x in 1:length(env_condition)) {
  for (y in 1:length(beta_range)) {
    
    env <- env_condition[x]
    beta <- beta_range[y]
    
    VR_current_pre_good_2 <- rep(NA, runs)
    VR_current_post_good_2 <- rep(NA, runs)
    
    for (z in 1:runs) {
      
      #State the variables
      time <- 300 
      N1 <- rep(NA, time) 
      N2 <- rep(NA, time)
      N4 <- rep(NA, time)
      
      N1[1] <- 50
      N2[1] <- 50
      N4[1] <- 0
      
      #Set parameters
      beta12 <- beta   
      beta21 <- beta   
      sigmaE1 <- -env 
      sigmaE2 <- -env 
      miuE <- rnorm(time, mean = 0, sd = 1) #environmental timeseries variation
      miuD1 <- rnorm(time, mean = 0, sd = 1) #dem timeseries for species 1
      miuD2 <- rnorm(time, mean = 0, sd = 1) #dem timeseries for species 2
      miuD4 <- rnorm(time, mean = 0, sd = 1)
      
      #Create the model
      for (t in 1:(time-1)) {
        if(t == burn_in) {
          N4[t] <- 1
        }
        
        #calculate population sizes for species 1
        N1[t+1] <- N1[t]*exp(r1*(1-(N1[t]/K1) - (beta12*N2[t]/K2)) - (beta14*N4[t]/K4) +
                               (sigmaE1*miuE[t])+(sigmaD1*miuD1[t])/sqrt(N1[t]))
        
        #calculate population sizes for species 2
        N2[t+1] <- N2[t]*exp(r2*(1-(N2[t]/K2) - (beta21*N1[t]/K1)) - (beta24*N4[t]/K4) +
                               (sigmaE2*miuE[t])+(sigmaD2*miuD2[t])/sqrt(N2[t]))
        
        #calculate population sizes for species 4
        N4[t+1] <- N4[t]*exp(r4*(1-(N4[t]/K4) - (beta41*N1[t]/K1)) - (beta42*N2[t]/K2) +
                               (sigmaE4*miuE[t])+(sigmaD4*miuD4[t])/sqrt(N4[t]))
        
        N1[is.na(N1)] <- 0
        N2[is.na(N2)] <- 0
        N4[is.na(N4)] <- 0
        
        if(N1[t+1] < 1) {
          N1[t+1] <- 0
        }
        if(N2[t+1] < 1) {
          N2[t+1] <- 0
        }
        if(N4[t+1] < 1) {
          N4[t+1] <- 0
        }
      }  

      #calculate variance ratio, pre-invasion
      # create an empty matrix
      our_data_pre_good_2 <- matrix(NA, nrow = 2*timeseries, ncol = 3)
      # fill in the first column with our species ID
      our_data_pre_good_2[,1] <- rep(c(1, 2), timeseries)
      # fill in the second column with the timepoint
      our_data_pre_good_2[,2] <- rep(seq(1:timeseries), each = 2)
      # fill in our third colum with the abundances from our simulation model
      our_data_pre_good_2[,3] <- as.vector(rbind(N1[(burn_in-timeseries):(burn_in-1)], 
                                               N2[(burn_in-timeseries):(burn_in-1)]))
      
      our_data_pre_good_v2_2 <- data.frame(our_data_pre_good_2)
      colnames(our_data_pre_good_v2_2) <- c("species", "time", "abundance")
      
      # calculate the variance ratio
      VR_temp_pre_good_2 <- variance_ratio(our_data_pre_good_v2_2, time.var = "time", 
                                         species.var = "species", abundance.var = "abundance", bootnumber = 1)
      VR_current_pre_good_2[z] <- VR_temp_pre_good_2$VR
      
      #calculate variance ratio, post-invasion
      # create an empty matrix
      our_data_post_good_2 <- matrix(NA, nrow = 2*timeseries, ncol = 3)
      # fill in the first column with our species ID
      our_data_post_good_2[,1] <- rep(c(1, 2), timeseries)
      # fill in the second column with the timepoint
      our_data_post_good_2[,2] <- rep(seq(1:timeseries), each = 2)
      # fill in our third colum with the abundances from our simulation model
      our_data_post_good_2[,3] <- as.vector(rbind(N1[(burn_in+1):time], 
                                                N2[(burn_in+1):time]))
      
      our_data_post_good_v2_2 <- data.frame(our_data_post_good_2)
      colnames(our_data_post_good_v2_2) <- c("species", "time", "abundance")
      
      if(N4[invasion_success] > 1) {
        # calculate the variance ratio
        VR_temp_post_good_2 <- variance_ratio(our_data_post_good_v2_2, time.var = "time",
                                            species.var = "species", abundance.var = "abundance", bootnumber = 1)
        VR_current_post_good_2[z] <- VR_temp_post_good_2$VR
      }
    }
    VR_pre_good_2[x,y] <- mean(VR_current_pre_good_2)
    VR_post_good_2[x,y] <- mean(VR_current_post_good_2, na.rm = TRUE)
  }
  print(x/length(env_condition))
}

#Proportion of VR > 1 in pre invasion communities
mean(VR_pre_good_2>1, na.rm = TRUE)
#Proportion of VR > 1 in post invasion communities
mean(VR_post_good_2>1, na.rm = TRUE)

M5 <- melt(VR_post_good_2)

ggplot(M5, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#008080", high ="#ca562c", mid = "#f6edbd", midpoint = 1, limit = c(0,2)) + 
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), y= expression(paste("Strength of Competititon (", beta, ")")), title = "Two-Species with Strong Invader VRs", fill="VR") +
  theme(axis.text = element_text( size = 12)) + 
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95"))

#Figure 3.c: Weak Invader in 3 Species Community----
burn_in <- 200
timeseries <- 100
invasion_success <- burn_in + 10 

VR_pre_poor_3 <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))
VR_post_poor_3 <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))

for (x in 1:length(env_condition)) {
  for (y in 1:length(beta_range)) {
    
    env <- env_condition[x]
    beta <- beta_range[y]
    
    VR_current_pre_poor_3 <- rep(NA, runs)
    VR_current_post_poor_3 <- rep(NA, runs)
    
    for (z in 1:runs) {
      
      #State the variables
      time <- 300 
      N1 <- rep(NA, time) 
      N2 <- rep(NA, time)
      N3 <- rep(NA, time)
      
      N1[1] <- 50
      N2[1] <- 50
      N3[1] <- 0
      
      #Set parameters
      beta12 <- beta   
      beta21 <- beta   
      sigmaE1 <- -env 
      sigmaE2 <- -env 
      miuE <- rnorm(time, mean = 0, sd = 1) #environmental timeseries variation
      miuD1 <- rnorm(time, mean = 0, sd = 1) #dem timeseries for species 1
      miuD2 <- rnorm(time, mean = 0, sd = 1) #dem timeseries for species 2
      miuD3 <- rnorm(time, mean = 0, sd = 1)
      
      #Create the model
      for (t in 1:(time-1)) {
        if(t == burn_in) {
          N3[t] <- 1
        }
        
        #calculate population sizes for species 1
        N1[t+1] <- N1[t]*exp(r1*(1-(N1[t]/K1) - (beta12*N2[t]/K2)) - (beta13*N3[t]/K3) +
                               (sigmaE1*miuE[t])+(sigmaD1*miuD1[t])/sqrt(N1[t]))
        
        #calculate population sizes for species 2
        N2[t+1] <- N2[t]*exp(r2*(1-(N2[t]/K2) - (beta21*N1[t]/K1)) - (beta23*N3[t]/K3) +
                               (sigmaE2*miuE[t])+(sigmaD2*miuD2[t])/sqrt(N2[t]))
        
        #calculate population sizes for species 3
        N3[t+1] <- N3[t]*exp(r3*(1-(N3[t]/K3) - (beta31*N1[t]/K1)) - (beta32*N2[t]/K2) +
                               (sigmaE3*miuE[t]) + (sigmaD3*miuD3[t])/sqrt(N3[t]))
        
        N1[is.na(N1)] <- 0
        N2[is.na(N2)] <- 0
        N3[is.na(N3)] <- 0
        
        if(N1[t+1] < 1) {
          N1[t+1] <- 0
        }
        if(N2[t+1] < 1) {
          N2[t+1] <- 0
        }
        if(N3[t+1] < 1) {
          N3[t+1] <- 0
        }
      }  
      
      #calculate variance ratio, pre-invasion
      # create an empty matrix
      our_data_pre_poor_3 <- matrix(NA, nrow = 3*(burn_in + timeseries), ncol = 3)
      # fill in the first column with our species ID
      our_data_pre_poor_3[,1] <- rep(c(1, 2, 3), (burn_in + timeseries))
      # fill in the second column with the timepoint
      our_data_pre_poor_3[,2] <- rep(seq(1:(burn_in + timeseries)), each = 3)
      # fill in our third colum with the abundances from our simulation model
      our_data_pre_poor_3[,3] <- as.vector(rbind(N1[(burn_in-timeseries):(burn_in-1)], 
                                               N2[(burn_in-timeseries):(burn_in-1)],
                                               N3[(burn_in-timeseries):(burn_in-1)]))
      
      our_data_pre_poor_v2_3 <- data.frame(our_data_pre_poor_3)
      colnames(our_data_pre_poor_v2_3) <- c("species", "time", "abundance")
      
      # calculate the variance ratio
      VR_temp_pre_poor_3 <- variance_ratio(our_data_pre_poor_v2_3, time.var = "time", species.var = "species", 
                                         abundance.var = "abundance", bootnumber = 1)
      VR_current_pre_poor_3[z] <- VR_temp_pre_poor_3$VR
      
      #calculate variance ratio, post-invasion
      # create an empty matrix
      our_data_post_poor_3 <- matrix(NA, nrow = 3*timeseries, ncol = 3)
      # fill in the first column with our species ID
      our_data_post_poor_3[,1] <- rep(c(1, 2, 3), timeseries)
      # fill in the second column with the timepoint
      our_data_post_poor_3[,2] <- rep(seq(1:timeseries), each = 3)
      # fill in our third colum with the abundances from our simulation model
      our_data_post_poor_3[,3] <- as.vector(rbind(N1[(burn_in+1):time], 
                                                N2[(burn_in+1):time],
                                                N3[(burn_in+1):time]))
      
      our_data_post_poor_v2_3 <- data.frame(our_data_post_poor_3)
      colnames(our_data_post_poor_v2_3) <- c("species", "time", "abundance")
      
      if(N3[invasion_success] > 1) {
        
        # calculate the variance ratio
        VR_temp_post_poor_3 <- variance_ratio(our_data_post_poor_v2_3, time.var = "time", 
                                            species.var = "species", abundance.var = "abundance", bootnumber = 1)
        VR_current_post_poor_3[z] <- VR_temp_post_poor_3$VR
      }
    }
    VR_pre_poor_3[x,y] <- mean(VR_current_pre_poor_3)
    VR_post_poor_3[x,y] <- mean(VR_current_post_poor_3, na.rm = TRUE)
  }
  print(x/length(env_condition))
}

#3 species community post-invasion heatmap
M6 <- melt(VR_post_poor_3) 

ggplot(M6, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#008080", high ="#ca562c", mid = "#f6edbd", midpoint = 1, limit = c(0,2)) + 
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), y= expression(paste("Strength of Competititon (", beta, ")")), title = "Three-Species with Weak Invader VRs", fill="VR") +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95"))

#Figure 3.d: Strong Invader in 3 Species Community----
burn_in <- 200
timeseries <- 100
invasion_success <- burn_in + 10

VR_pre_good_3 <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))
VR_post_good_3 <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))

for (x in 1:length(env_condition)) {
  for (y in 1:length(beta_range)) {
    
    env <- env_condition[x]
    beta <- beta_range[y]
    
    VR_current_pre_good_3 <- rep(NA, runs)
    VR_current_post_good_3 <- rep(NA, runs)
    
    for (z in 1:runs) {
      
      #State the variables
      time <- 300 
      N1 <- rep(NA, time) 
      N2 <- rep(NA, time)
      N4 <- rep(NA, time)
      
      N1[1] <- 50
      N2[1] <- 50
      N4[1] <- 0
      
      #Set parameters
      beta12 <- beta   
      beta21 <- beta   
      sigmaE1 <- -env 
      sigmaE2 <- -env 
      miuE <- rnorm(time, mean = 0, sd = 1) #environmental timeseries variation
      miuD1 <- rnorm(time, mean = 0, sd = 1) #dem timeseries for species 1
      miuD2 <- rnorm(time, mean = 0, sd = 1) #dem timeseries for species 2
      miuD4 <- rnorm(time, mean = 0, sd = 1)
      
      #Create the model
      for (t in 1:(time-1)) {
        if(t == burn_in) {
          N4[t] <- 1
        }
        
        #calculate population sizes for species 1
        N1[t+1] <- N1[t]*exp(r1*(1-(N1[t]/K1) - (beta12*N2[t]/K2)) - (beta14*N4[t]/K4) +
                               (sigmaE1*miuE[t])+(sigmaD1*miuD1[t])/sqrt(N1[t]))
        
        #calculate population sizes for species 2
        N2[t+1] <- N2[t]*exp(r2*(1-(N2[t]/K2) - (beta21*N1[t]/K1)) - (beta24*N4[t]/K4) +
                               (sigmaE2*miuE[t])+(sigmaD2*miuD2[t])/sqrt(N2[t]))
        
        #calculate population sizes for species 4
        N4[t+1] <- N4[t]*exp(r4*(1-(N4[t]/K4) - (beta41*N1[t]/K1)) - (beta42*N2[t]/K2) +
                               (sigmaE4*miuE[t])+(sigmaD4*miuD4[t])/sqrt(N4[t]))
        
        
        N1[is.na(N1)] <- 0
        N2[is.na(N2)] <- 0
        N4[is.na(N4)] <- 0
        
        if(N1[t+1] < 1) {
          N1[t+1] <- 0
        }
        if(N2[t+1] < 1) {
          N2[t+1] <- 0
        }
        if(N4[t+1] < 1) {
          N4[t+1] <- 0
        }
      }  
      
      #calculate variance ratio
      # create an empty matrix
      our_data_pre_good_3 <- matrix(NA, nrow = 3*timeseries, ncol = 3)
      # fill in the first column with our species ID
      our_data_pre_good_3[,1] <- rep(c(1, 2, 4), timeseries)
      # fill in the second column with the timepoint
      our_data_pre_good_3[,2] <- rep(seq(1:timeseries), each = 3)
      # fill in our third colum with the abundances from our simulation model
      our_data_pre_good_3[,3] <- as.vector(rbind(N1[(burn_in-timeseries):(burn_in-1)], 
                                               N2[(burn_in-timeseries):(burn_in-1)],
                                               N4[(burn_in-timeseries):(burn_in-1)]))
      
      our_data_pre_good_v2_3 <- data.frame(our_data_pre_good_3)
      colnames(our_data_pre_good_v2_3) <- c("species", "time", "abundance")
      
      # calculate the variance ratio
      VR_temp_pre_good_3 <- variance_ratio(our_data_pre_good_v2_3, time.var = "time", 
                                         species.var = "species", abundance.var = "abundance", bootnumber = 1)
      VR_current_pre_good_3[z] <- VR_temp_pre_good_3$VR
      
      #calculate variance ratio, post-invasion
      # create an empty matrix
      our_data_post_good_3 <- matrix(NA, nrow = 3*timeseries, ncol = 3)
      # fill in the first column with our species ID
      our_data_post_good_3[,1] <- rep(c(1, 2, 4), timeseries)
      # fill in the second column with the timepoint
      our_data_post_good_3[,2] <- rep(seq(1:timeseries), each = 3)
      # fill in our third colum with the abundances from our simulation model
      our_data_post_good_3[,3] <- as.vector(rbind(N1[(burn_in+1):time], 
                                                N2[(burn_in+1):time],
                                                N4[(burn_in+1):time]))
      
      our_data_post_good_v2_3 <- data.frame(our_data_post_good_3)
      colnames(our_data_post_good_v2_3) <- c("species", "time", "abundance")
      
      if(N4[invasion_success] > 1) {
        # calculate the variance ratio
        VR_temp_post_good_3 <- variance_ratio(our_data_post_good_v2_3, time.var = "time",
                                            species.var = "species", abundance.var = "abundance", bootnumber = 1)
        VR_current_post_good_3[z] <- VR_temp_post_good_3$VR
      }
    }
    VR_pre_good_3[x,y] <- mean(VR_current_pre_good_3)
    VR_post_good_3[x,y] <- mean(VR_current_post_good_3, na.rm = TRUE)
  }
  print(x/length(env_condition))
}

M7 <- melt(VR_post_good_3)

ggplot(M7, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#008080", high ="#ca562c", mid = "#f6edbd", midpoint = 1, limit = c(0,2)) + 
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), y= expression(paste("Strength of Competititon (", beta, ")")), title = "Three-Species with Strong Invader VRs", fill="VR") +
  theme(axis.text = element_text( size = 12)) + 
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95"))

# ----------------------------------------------------------------------------------------
###Figure 4: Comparing the difference between pre- and post-invasion variance ratios

#Figure 4.a: Weak Invader in 2 Species Community----
M8 <- melt(VR_post_poor_2 - VR_pre_poor_2)

ggplot(M8, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#509b4b", high ="#9669b0", mid = "#f7f7f7", midpoint = 0, limit = c(-1, 1)) + 
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), y= expression(paste("Strength of Competititon (", beta, ")")), title = "Change in Two Species VRs with Weak Invader", fill="VR") +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95"))

#Figure 4.b: Strong Invader in 2 Species Community----
M9 <- melt(VR_post_good_2 - VR_pre_good_2)

ggplot(M9, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#509b4b", high ="#9669b0", mid = "#f7f7f7", midpoint = 0, limit = c(-1, 1)) + 
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), y= expression(paste("Strength of Competititon (", beta, ")")), title = "Change in Two Species VRs with Strong Invader", fill="VR") +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95"))

#Figure 4.c: Weak Invader in 3 Species Community----
M10 <- melt(VR_post_poor_3 - VR_pre_poor_3)

ggplot(M10, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#509b4b", high ="#9669b0", mid = "#f7f7f7", midpoint = 0, limit = c(-1, 1)) + 
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), y= expression(paste("Strength of Competititon (", beta, ")")), title = "Change in Three Species VRs with Weak Invader", fill="VR") +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95"))

#Figure 4.d: Strong Invader in 3 Species Community----
M11 <- melt(VR_post_good_3 - VR_pre_good_3)

ggplot(M11, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#509b4b", high ="#9669b0", mid = "#f7f7f7", midpoint = 0, limit = c(-1.25, 1)) + 
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), y= expression(paste("Strength of Competititon (", beta, ")")), title = "Change in Three Species VRs with Strong Invader", fill="VR") +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95"))



#Figure 5: Growth Rate depending on resident abundance-----





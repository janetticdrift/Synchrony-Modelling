#####################################
# This script produces figures that show the range in 
# variance ratios of resident communities,
# and their subsequent synchrony dynamics after invasion 
# of a strong and weak invader.
#####################################

# Package for calculating the variance ratio
library(codyn)

#Packages for creating plots
library(ggplot2)
library(reshape2)
library(ggpubr)

# Package for setting paths
library(here)

runs <- 2000 # number of runs

#load species parameters
source(here("Species_Parameters_Source.R"))

###Figure 1: Resident communities' starting variance ratios----
time <- 200 
burn_in <- 100

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
      
      #Set state variables
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
        N1[t+1] <- N1[t]*exp(r1*(1-(N1[t]/K1) - (beta12*N2[t]/K2)) +
                               (sigmaE1*miuE[t]) + (sigmaD1*miuD1[t])/sqrt(N1[t]))
        
        if(is.nan(N1[t+1])) {
          N1[t+1] <- 0
        }
        
        if(N1[t+1] < 1) {
          N1[t+1] <- 0
        }
        #calculate population sizes for species 2
        N2[t+1] <- N2[t]*exp(r2*(1-(N2[t]/K2) - (beta21*N1[t]/K1)) +
                               (sigmaE2*miuE[t]) + (sigmaD2*miuD2[t])/sqrt(N2[t]))
        
        if(is.nan(N2[t+1])) {
          N2[t+1] <- 0
        }
        
        if(N2[t+1] < 1) {
          N2[t+1] <- 0
        }
        
      }
      
      #Calculate total biomass through time in the community
      total_biomass <- N1 + N2
      
      #Current coefficient of variation
      CV_current[z] <- sd(total_biomass)/mean(total_biomass)
      
      #calculate variance ratio
      our_data <- matrix(NA, nrow = 2*(time-burn_in), ncol = 3)
      our_data[,1] <- rep(c(1, 2), (time-burn_in))
      our_data[,2] <- rep(seq(1:(time-burn_in)), each = 2)
      our_data[,3] <- as.vector(rbind(N1[(burn_in+1):time],N2[(burn_in+1):time]))
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

# Create heat map
M1 <- melt(VR) 

ggplot(M1, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#008080", high ="#ca562c", mid = "#f6edbd", midpoint = 1, limit = c(0,2)) +
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), y= expression(paste("Strength of Competititon (", beta, ")")), title = "Resident Communities' Variance Ratios", fill="VR") +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95"))

###Figure 2: Resistance of resident communities to invader growth of a weak and strong invader

#Figure 2.a: Weak Invader-----
#Set outputs of interest
avg_lambda_poor <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))
success_lambda_poor <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))
burn_in <- 100
time <- 200 

for (x in 1:length(env_condition)) {
  for (y in 1:length(beta_range)) {
    
    env <- env_condition[x]
    beta <- beta_range[y]
    
    lambda_poor_a <- rep(NA, runs) 
    lambda_poor <- matrix(NA, nrow = time-burn_in, ncol = runs)
    
    for (z in 1:runs) {
      
      #State variables
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
        
        #calculate population size for species 1
        N1[t+1] <- N1[t]*exp(r1*(1-(N1[t]/K1) - (beta12*N2[t]/K2)) +
                               (sigmaE1*miuE[t]) + (sigmaD1*miuD1[t])/sqrt(N1[t]))
        
        #calculate population size for species 2
        N2[t+1] <- N2[t]*exp(r2*(1-(N2[t]/K2) - (beta21*N1[t]/K1)) +
                               (sigmaE2*miuE[t]) + (sigmaD2*miuD2[t])/sqrt(N2[t]))
        
        # Needed for if N1 or N2 are 0 for including demographic stochasticity
        # as sqrt(0) yields NaN
        if(is.nan(N1[t+1])) {
          N1[t+1] <- 0
        }
        if(is.nan(N2[t+1])) {
          N2[t+1] <- 0
        }
        
        if(N1[t+1] < 1) {
          N1[t+1] <- 0
        }
        if(N2[t+1] < 1) {
          N2[t+1] <- 0
        }
        
      }
      
      counter <- 1
      
      #Model effects of poor invasion
      for (t in burn_in:(time-1)) {
        
        N3[counter] <- invader_abund_poor*exp(r3*(1-(invader_abund_poor/K3) - (beta31*N1[t]/K1) - (beta32*N2[t]/K2))
                                              + (sigmaE3*miuE[t]) + (sigmaD3*miuD3[counter])/sqrt(invader_abund_poor))
        # Needed for when N4 is 0 for including demographic stochasticity
        # as sqrt(0) yields NaN

        if(is.nan(N3[t+1])) {
          N3[t+1] <- 0
        }
        if(N3[counter] < 1) {
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


#Figure 2.b: Strong Invader-----
#Set outputs of interest
avg_lambda_good <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))
success_lambda_good <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))

for (x in 1:length(env_condition)) {
  for (y in 1:length(beta_range)) {
    
    env <- env_condition[x]
    beta <- beta_range[y]
    
    lambda_good_a <- rep(NA, runs) 
    lambda_good <- matrix(NA, nrow = time-burn_in, ncol = runs)
    
    for (z in 1:runs) {
      
      #State variables
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
        #calculate population size for species 1
        N1[t+1] <- N1[t]*exp(r1*(1-(N1[t]/K1) - (beta12*N2[t]/K2)) +
                               (sigmaE1*miuE[t]) + (sigmaD1*miuD1[t])/sqrt(N1[t]))
        
        #calculate population size for species 2
        N2[t+1] <- N2[t]*exp(r2*(1-(N2[t]/K2) - (beta21*N1[t]/K1)) +
                               (sigmaE2*miuE[t]) + (sigmaD2*miuD2[t])/sqrt(N2[t]))
        
        # Needed for if N1 or N2 are 0 for including demographic stochasticity
        # as sqrt(0) yields NaN
        if(is.nan(N1[t+1])) {
          N1[t+1] <- 0
        }
        if(is.nan(N2[t+1])) {
          N2[t+1] <- 0
        }

        if(N1[t+1] < 1) {
          N1[t+1] <- 0
        }
        if(N2[t+1] < 1) {
          N2[t+1] <- 0
        }
        
      }

      counter <- 1
      
      #Model effects of strong invasion
      for (t in burn_in:(time-1)) {
        N4[counter] <- invader_abund_good*exp(r4*(1-(invader_abund_good/K4) - (beta41*N1[t]/K1) - (beta42*N2[t]/K2))
                                              + (sigmaE4*miuE[t])+(sigmaD4*miuD4[counter])/sqrt(invader_abund_good))
        
        # Needed for when N4 is 0 for including demographic stochasticity
        # as sqrt(0) yields NaN
        # Needed for when N4 is 0. sqrt(0) yields NaN
        if(is.nan(N4[t+1])) {
          N4[t+1] <- 0
        }
        
        if(N4[counter] < 1) {
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

#Figure 2.a and b: Plotting----
#Successful Invasion of Poor Invader Plot

# set plot limits
min_lim <- min(c(success_lambda_good, success_lambda_poor))
max_lim <- max(c(success_lambda_good, success_lambda_poor))

M2 <- melt(success_lambda_poor)

plot1 <- ggplot(M2, aes(x=Var1, y=Var2, fill=value*100)) + 
  geom_tile() + 
  scale_fill_distiller(palette = "RdBu", limits = c(min_lim*100, max_lim*100)) +
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), y= expression(paste("Strength of Competititon (", beta, ")")), title = "Growth Rates of Weak Invader", fill="Growth
Success (%)") +
  theme(axis.text = element_text( size = 12)) + 
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95"))

#Successful Invasion of Good Invader Plot
M3 <- melt(success_lambda_good)

plot2 <- ggplot(M3, aes(x=Var1, y=Var2, fill=value*100)) +
  geom_tile() + 
  scale_fill_distiller(palette = "RdBu", limits = c(min_lim*100, max_lim*100)) +
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), y= expression(paste("Strength of Competititon (", beta, ")")), title = "Growth Rates of Strong Invader", fill="Growth
Success (%)") +
  theme(axis.text = element_text( size = 12)) + 
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95"))

#Combine plots together, label A and B
quartz(width=11, height=5)

ggarrange(
  plot1, plot2, labels = c("A", "B"),
  common.legend = TRUE, legend = "right"
)

###Figure 3: Resilience of resident communities to established invaders

#Figure 3.a and c: Weak Invader in 2 Species Community----
time <- 300 
burn_in <- 200
timeseries <- 100 # length of timeseries to use for VR calculations
invasion_success <- burn_in + 10

# two species community
VR_pre_poor <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))
VR_post_poor <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))

# three species community
VR_pre_poor_full <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))
VR_post_poor_full <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))

for (x in 1:length(env_condition)) {
  for (y in 1:length(beta_range)) {
    
    env <- env_condition[x]
    beta <- beta_range[y]
    
    VR_current_pre_poor <- rep(NA, runs)
    VR_current_post_poor <- rep(NA, runs)
    
    VR_current_pre_poor_full <- rep(NA, runs)
    VR_current_post_poor_full <- rep(NA, runs)
    
    for (z in 1:runs) {
      
      #State variables
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
        N1[t+1] <- N1[t]*exp(r1*(1-(N1[t]/K1) - (beta12*N2[t]/K2) - (beta13*N3[t]/K3)) +
                               (sigmaE1*miuE[t]) + (sigmaD1*miuD1[t])/sqrt(N1[t]))
        
        #calculate population sizes for species 2
        N2[t+1] <- N2[t]*exp(r2*(1-(N2[t]/K2) - (beta21*N1[t]/K1) - (beta23*N3[t]/K3)) +
                               (sigmaE2*miuE[t]) + (sigmaD2*miuD2[t])/sqrt(N2[t]))
        
        #calculate population sizes for species 3
        N3[t+1] <- N3[t]*exp(r3*(1-(N3[t]/K3) - (beta31*N1[t]/K1) - (beta32*N2[t]/K2)) +
                               (sigmaE3*miuE[t]) + (sigmaD3*miuD3[t])/sqrt(N3[t]))
        
        # Needed for when N is 0. sqrt(0) yields NaN
        if(is.nan(N1[t+1])) {
          N1[t+1] <- 0
        }
        if(is.nan(N2[t+1])) {
          N2[t+1] <- 0
        }
        if(is.nan(N3[t+1])) {
          N3[t+1] <- 0
        }
        
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
      
      #calculate variance ratio, pre-invasion; 2 species community
      # create an empty matrix
      our_data_pre_poor <- matrix(NA, nrow = 2*(burn_in - timeseries), ncol = 3)
      # fill in the first column with our species ID
      our_data_pre_poor[,1] <- rep(c(1, 2), (burn_in - timeseries))
      # fill in the second column with the timepoint
      our_data_pre_poor[,2] <- rep(seq(1:(burn_in - timeseries)), each = 2)
      # fill in our third colum with the abundances from our simulation model
      our_data_pre_poor[,3] <- as.vector(rbind(N1[(burn_in-timeseries):(burn_in-1)], 
                                               N2[(burn_in-timeseries):(burn_in-1)]))
      
      df_our_data_pre_poor <- data.frame(our_data_pre_poor)
      colnames(df_our_data_pre_poor) <- c("species", "time", "abundance")
      
      # calculate the variance ratio
      VR_temp_pre_poor <- variance_ratio(df_our_data_pre_poor, time.var = "time", species.var = "species", 
                                         abundance.var = "abundance", bootnumber = 1)
      VR_current_pre_poor[z] <- VR_temp_pre_poor$VR
      
      #calculate variance ratio, pre-invasion; 3 species community
      our_data_pre_poor_full <- matrix(NA, nrow = 3*(burn_in - timeseries), ncol = 3)
      # fill in the first column with our species ID
      our_data_pre_poor_full[,1] <- rep(c(1, 2, 3), (burn_in - timeseries))
      # fill in the second column with the timepoint
      our_data_pre_poor_full[,2] <- rep(seq(1:(burn_in - timeseries)), each = 3)
      # fill in our third colum with the abundances from our simulation model
      our_data_pre_poor_full[,3] <- as.vector(rbind(N1[(burn_in-timeseries):(burn_in-1)], 
                                                 N2[(burn_in-timeseries):(burn_in-1)],
                                                 N3[(burn_in-timeseries):(burn_in-1)]))
      
      df_our_data_pre_poor_full <- data.frame(our_data_pre_poor_full)
      colnames(df_our_data_pre_poor_full) <- c("species", "time", "abundance")
      
      # calculate the variance ratio
      VR_temp_pre_poor_full <- variance_ratio(df_our_data_pre_poor_full, time.var = "time", species.var = "species", 
                                           abundance.var = "abundance", bootnumber = 1)
      VR_current_pre_poor_full[z] <- VR_temp_pre_poor_full$VR
      
      #calculate variance ratio, post-invasion; 2 species community
      # create an empty matrix
      our_data_post_poor <- matrix(NA, nrow = 2*timeseries, ncol = 3)
      # fill in the first column with our species ID
      our_data_post_poor[,1] <- rep(c(1, 2), timeseries)
      # fill in the second column with the timepoint
      our_data_post_poor[,2] <- rep(seq(1:timeseries), each = 2)
      # fill in our third colum with the abundances from our simulation model
      our_data_post_poor[,3] <- as.vector(rbind(N1[(burn_in+1):time], 
                                                N2[(burn_in+1):time]))

      df_our_data_post_poor <- data.frame(our_data_post_poor)
      colnames(df_our_data_post_poor) <- c("species", "time", "abundance")
      
      # three species community post invasion
      our_data_post_poor_full <- matrix(NA, nrow = 3*timeseries, ncol = 3)
      # fill in the first column with our species ID
      our_data_post_poor_full[,1] <- rep(c(1, 2,3), timeseries)
      # fill in the second column with the timepoint
      our_data_post_poor_full[,2] <- rep(seq(1:timeseries), each = 3)
      # fill in our third colum with the abundances from our simulation model
      our_data_post_poor_full[,3] <- as.vector(rbind(N1[(burn_in+1):time], 
                                                N2[(burn_in+1):time],
                                                N3[(burn_in+1):time]))
      
      df_our_data_post_poor_full <- data.frame(our_data_post_poor_full)
      colnames(df_our_data_post_poor_full) <- c("species", "time", "abundance")
      
      if(N3[invasion_success] > 1) {
        
        # calculate the variance ratio; 2 species
        VR_temp_post_poor <- variance_ratio(df_our_data_post_poor, time.var = "time", 
                                            species.var = "species", abundance.var = "abundance", bootnumber = 1)
        VR_current_post_poor[z] <- VR_temp_post_poor$VR
        
        # calculate the variance ratio; full community
        VR_temp_post_poor_full <- variance_ratio(df_our_data_post_poor_full, time.var = "time", 
                                            species.var = "species", abundance.var = "abundance", bootnumber = 1)
        VR_current_post_poor_full[z] <- VR_temp_post_poor_full$VR
      }
    }
    
    VR_pre_poor[x,y] <- mean(VR_current_pre_poor, na.rm=TRUE)
    VR_post_poor[x,y] <- mean(VR_current_post_poor, na.rm = TRUE)
    
    VR_pre_poor_full[x,y] <- mean(VR_current_pre_poor_full, na.rm=TRUE)
    VR_post_poor_full[x,y] <- mean(VR_current_post_poor_full, na.rm = TRUE)
  }
  print(x/length(env_condition))
}

# Proportion of VR > 1 in pre invasion communities
mean(VR_pre_poor>1, na.rm = TRUE)
mean(VR_pre_poor_full>1, na.rm = TRUE)
#Proportion of VR > 1 in post invasion communities
mean(VR_post_poor>1, na.rm = TRUE)
mean(VR_post_poor_full>1, na.rm = TRUE)

M4 <- melt(VR_post_poor)

ggplot(M4, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#008080", high ="#ca562c", mid = "#f6edbd", midpoint = 1, limit = c(0,2)) + 
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), y= expression(paste("Strength of Competititon (", beta, ")")), title = "Two-Species with Weak Invader VRs", fill="VR") +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95"))

#Figure 3.b and d: Strong Invader in 2 Species Community----
time <- 300 
burn_in <- 200
timeseries <- 100 # length of timeseries to use for VR calculations
invasion_success <- burn_in + 10

# two species community
VR_pre_good <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))
VR_post_good <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))

# three species community
VR_pre_good_full <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))
VR_post_good_full <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))

for (x in 1:length(env_condition)) {
  for (y in 1:length(beta_range)) {
    
    env <- env_condition[x]
    beta <- beta_range[y]
    
    VR_current_pre_good <- rep(NA, runs)
    VR_current_post_good <- rep(NA, runs)
    
    VR_current_pre_good_full <- rep(NA, runs)
    VR_current_post_good_full <- rep(NA, runs)
    
    for (z in 1:runs) {
      
      #State variables
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
      miuE <- rnorm(time, mean = 0, sd = 1) 
      miuD1 <- rnorm(time, mean = 0, sd = 1)
      miuD2 <- rnorm(time, mean = 0, sd = 1)
      miuD4 <- rnorm(time, mean = 0, sd = 1)
      
      #Create the model
      for (t in 1:(time-1)) {
        #Introduce invader at time burn_in at abundance of 1
        if(t == burn_in) {
          N4[t] <- 1
        }
        
        #calculate population sizes for species 1
        N1[t+1] <- N1[t]*exp(r1*(1-(N1[t]/K1) - (beta12*N2[t]/K2) - (beta14*N4[t]/K4)) +
                               (sigmaE1*miuE[t])+(sigmaD1*miuD1[t])/sqrt(N1[t]))
        
        #calculate population sizes for species 2
        N2[t+1] <- N2[t]*exp(r2*(1-(N2[t]/K2) - (beta21*N1[t]/K1) - (beta24*N4[t]/K4)) +
                               (sigmaE2*miuE[t])+(sigmaD2*miuD2[t])/sqrt(N2[t]))
        
        #calculate population sizes for species 4
        N4[t+1] <- N4[t]*exp(r4*(1-(N4[t]/K4) - (beta41*N1[t]/K1) - (beta42*N2[t]/K2)) +
                               (sigmaE4*miuE[t])+(sigmaD4*miuD4[t])/sqrt(N4[t]))
        
        # Needed for when N is 0. sqrt(0) yields NaN
        if(is.nan(N1[t+1])) {
          N1[t+1] <- 0
        }
        if(is.nan(N2[t+1])) {
          N2[t+1] <- 0
        }
        if(is.nan(N4[t+1])) {
          N4[t+1] <- 0
        }
        
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
      
      #calculate variance ratio, pre-invasion; 2 species community
      # create an empty matrix
      our_data_pre_good <- matrix(NA, nrow = 2*(burn_in - timeseries), ncol = 3)
      # fill in the first column with our species ID
      our_data_pre_good[,1] <- rep(c(1, 2), (burn_in - timeseries))
      # fill in the second column with the timepoint
      our_data_pre_good[,2] <- rep(seq(1:(burn_in - timeseries)), each = 2)
      # fill in our third colum with the abundances from our simulation model
      our_data_pre_good[,3] <- as.vector(rbind(N1[(burn_in-timeseries):(burn_in-1)], 
                                               N2[(burn_in-timeseries):(burn_in-1)]))
      
      df_our_data_pre_good <- data.frame(our_data_pre_good)
      colnames(df_our_data_pre_good) <- c("species", "time", "abundance")
      
      # calculate the variance ratio
      VR_temp_pre_good <- variance_ratio(df_our_data_pre_good, time.var = "time", species.var = "species", 
                                         abundance.var = "abundance", bootnumber = 1)
      VR_current_pre_good[z] <- VR_temp_pre_good$VR
      
      #calculate variance ratio, pre-invasion; 3 species community
      our_data_pre_good_full <- matrix(NA, nrow = 3*(burn_in - timeseries), ncol = 3)
      # fill in the first column with our species ID
      our_data_pre_good_full[,1] <- rep(c(1, 2, 3), (burn_in - timeseries))
      # fill in the second column with the timepoint
      our_data_pre_good_full[,2] <- rep(seq(1:(burn_in - timeseries)), each = 3)
      # fill in our third colum with the abundances from our simulation model
      our_data_pre_good_full[,3] <- as.vector(rbind(N1[(burn_in-timeseries):(burn_in-1)], 
                                                    N2[(burn_in-timeseries):(burn_in-1)],
                                                    N4[(burn_in-timeseries):(burn_in-1)]))
      
      df_our_data_pre_good_full <- data.frame(our_data_pre_good_full)
      colnames(df_our_data_pre_good_full) <- c("species", "time", "abundance")
      
      # calculate the variance ratio
      VR_temp_pre_good_full <- variance_ratio(df_our_data_pre_good_full, time.var = "time", species.var = "species", 
                                              abundance.var = "abundance", bootnumber = 1)
      VR_current_pre_good_full[z] <- VR_temp_pre_good_full$VR
      
      #calculate variance ratio, post-invasion; 2 species community
      # create an empty matrix
      our_data_post_good <- matrix(NA, nrow = 2*timeseries, ncol = 3)
      # fill in the first column with our species ID
      our_data_post_good[,1] <- rep(c(1, 2), timeseries)
      # fill in the second column with the timepoint
      our_data_post_good[,2] <- rep(seq(1:timeseries), each = 2)
      # fill in our third colum with the abundances from our simulation model
      our_data_post_good[,3] <- as.vector(rbind(N1[(burn_in+1):time], 
                                                N2[(burn_in+1):time]))
      
      df_our_data_post_good <- data.frame(our_data_post_good)
      colnames(df_our_data_post_good) <- c("species", "time", "abundance")
      
      # three species community post invasion
      our_data_post_good_full <- matrix(NA, nrow = 3*timeseries, ncol = 3)
      # fill in the first column with our species ID
      our_data_post_good_full[,1] <- rep(c(1, 2,3), timeseries)
      # fill in the second column with the timepoint
      our_data_post_good_full[,2] <- rep(seq(1:timeseries), each = 3)
      # fill in our third colum with the abundances from our simulation model
      our_data_post_good_full[,3] <- as.vector(rbind(N1[(burn_in+1):time], 
                                                     N2[(burn_in+1):time],
                                                     N4[(burn_in+1):time]))
      
      df_our_data_post_good_full <- data.frame(our_data_post_good_full)
      colnames(df_our_data_post_good_full) <- c("species", "time", "abundance")
      
      if(N4[invasion_success] > 1) {
        
        # calculate the variance ratio; 2 species
        VR_temp_post_good <- variance_ratio(df_our_data_post_good, time.var = "time", 
                                            species.var = "species", abundance.var = "abundance", bootnumber = 1)
        VR_current_post_good[z] <- VR_temp_post_good$VR
        
        # calculate the variance ratio; full community
        VR_temp_post_good_full <- variance_ratio(df_our_data_post_good_full, time.var = "time", 
                                                 species.var = "species", abundance.var = "abundance", bootnumber = 1)
        VR_current_post_good_full[z] <- VR_temp_post_good_full$VR
      }
    }
    
    VR_pre_good[x,y] <- mean(VR_current_pre_good, na.rm=TRUE)
    VR_post_good[x,y] <- mean(VR_current_post_good, na.rm = TRUE)
    
    VR_pre_good_full[x,y] <- mean(VR_current_pre_good_full, na.rm=TRUE)
    VR_post_good_full[x,y] <- mean(VR_current_post_good_full, na.rm = TRUE)
  }
  print(x/length(env_condition))
}

#Proportion of VR > 1 in pre invasion communities
mean(VR_pre_good>1, na.rm = TRUE)
mean(VR_pre_good_full>1, na.rm = TRUE)
#Proportion of VR > 1 in post invasion communities
mean(VR_post_good>1, na.rm = TRUE)
mean(VR_post_good_full>1, na.rm = TRUE)

M5 <- melt(VR_post_good)

ggplot(M5, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#008080", high ="#ca562c", mid = "#f6edbd", midpoint = 1, limit = c(0,2)) + 
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), y= expression(paste("Strength of Competititon (", beta, ")")), title = "Two-Species with Strong Invader VRs", fill="VR") +
  theme(axis.text = element_text( size = 12)) + 
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95"))

# Plotting Figure 3c d: Strong Invader in 3 Species Community----
# 3 species community post-invasion heatmap
M6 <- melt(VR_post_poor_full) 

ggplot(M6, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#008080", high ="#ca562c", mid = "#f6edbd", midpoint = 1, limit = c(0,2)) + 
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), y= expression(paste("Strength of Competititon (", beta, ")")), title = "Three-Species with Weak Invader VRs", fill="VR") +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95"))

M7 <- melt(VR_post_good_full)

ggplot(M7, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#008080", high ="#ca562c", mid = "#f6edbd", midpoint = 1, limit = c(0,2)) + 
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), y= expression(paste("Strength of Competititon (", beta, ")")), title = "Three-Species with Strong Invader VRs", fill="VR") +
  theme(axis.text = element_text( size = 12)) + 
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95"))

###Figure 4: Comparing the difference between pre- and post-invasion variance ratios

#Figure 4.a: Weak Invader in 2 Species Community----
M8 <- melt(VR_post_poor - VR_pre_poor)

ggplot(M8, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#509b4b", high ="#9669b0", mid = "#f7f7f7", midpoint = 0, limit = c(-1, 1)) + 
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), y= expression(paste("Strength of Competititon (", beta, ")")), title = "Change in Two Species VRs with Weak Invader", fill="VR") +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95"))

#Figure 4.b: Strong Invader in 2 Species Community----
M9 <- melt(VR_post_good - VR_pre_good)

ggplot(M9, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#509b4b", high ="#9669b0", mid = "#f7f7f7", midpoint = 0, limit = c(-1, 1.5)) + 
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), y= expression(paste("Strength of Competititon (", beta, ")")), title = "Change in Two Species VRs with Strong Invader", fill="VR") +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95"))

#Figure 4.c: Weak Invader in 3 Species Community----
M10 <- melt(VR_post_poor_full - VR_pre_poor_full)

ggplot(M10, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#509b4b", high ="#9669b0", mid = "#f7f7f7", midpoint = 0, limit = c(-1, 1)) + 
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), y= expression(paste("Strength of Competititon (", beta, ")")), title = "Change in Three Species VRs with Weak Invader", fill="VR") +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95"))

#Figure 4.d: Strong Invader in 3 Species Community----
M11 <- melt(VR_post_good_full - VR_pre_good_full)

ggplot(M11, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#509b4b", high ="#9669b0", mid = "#f7f7f7", midpoint = 0, limit = c(-1.25, 1)) + 
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), y= expression(paste("Strength of Competititon (", beta, ")")), title = "Change in Three Species VRs with Strong Invader", fill="VR") +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95"))



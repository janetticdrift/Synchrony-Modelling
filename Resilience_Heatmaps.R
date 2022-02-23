#####################################
# This script produces figures that show 
# the subsequent effects of a strong and weak
# invader on the synchrony of the resident
# communities.
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
source(here("Species_Parameters_Source_Fig134.R"))

###Figure 3: Resilience of resident communities to established invaders

#Figure 3.a and c: Weak Invader in 2 and 3 Species Community----

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
#Proportion of VR > 1 in post invasion communities
mean(VR_post_poor>1, na.rm = TRUE)
mean(VR_post_poor_full>1, na.rm = TRUE)

#Figure 3.b and d: Strong Invader in 2 and 3 Species Community----
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
#Proportion of VR > 1 in post invasion communities
mean(VR_post_good>1, na.rm = TRUE)
mean(VR_post_good_full>1, na.rm = TRUE)

# Plotting Figure 3: Invaders in 2 and 3 Species Community----
#A. 2 species community weak invader post-invasion heatmap
M4 <- melt(VR_post_poor)

plot4 <- ggplot(M4, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#008080", high ="#ca562c", mid = "#f6edbd", midpoint = 1, limit = c(0,2)) + 
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), y= expression(paste("Strength of Competititon (", beta, ")")), title = "Two-Species with Weak Invader VRs", fill="VR") +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95"))

#B. 2 species community strong invader post-invasion heatmap
M5 <- melt(VR_post_good)

plot5 <- ggplot(M5, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#008080", high ="#ca562c", mid = "#f6edbd", midpoint = 1, limit = c(0,2)) + 
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), y= expression(paste("Strength of Competititon (", beta, ")")), title = "Two-Species with Strong Invader VRs", fill="VR") +
  theme(axis.text = element_text( size = 12)) + 
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95"))

#C. 3 species community weak invader post-invasion heatmap
M6 <- melt(VR_post_poor_full) 

plot6 <- ggplot(M6, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#008080", high ="#ca562c", mid = "#f6edbd", midpoint = 1, limit = c(0,2)) + 
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), y= expression(paste("Strength of Competititon (", beta, ")")), title = "Three-Species with Weak Invader VRs", fill="VR") +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95"))

#D. 3 species community strong invader post-invasion heatmap
M7 <- melt(VR_post_good_full)

plot7 <- ggplot(M7, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#008080", high ="#ca562c", mid = "#f6edbd", midpoint = 1, limit = c(0,2)) + 
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), y= expression(paste("Strength of Competititon (", beta, ")")), title = "Three-Species with Strong Invader VRs", fill="VR") +
  theme(axis.text = element_text( size = 12)) + 
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95"))

#Combine plots into one figure - TBD
quartz(width=10, height=5)

ggarrange(
  plot4, plot5, plot6, plot7, labels = c("A", "B", "C", "D"),
  common.legend = TRUE, legend = "right", nrow = 2, ncol = 2
)
###Figure 4: Comparing the difference between pre- and post-invasion variance ratios

#Figure 4.a: Weak Invader in 2 Species Community----
M8 <- melt(VR_post_poor - VR_pre_poor)

plot8 <- ggplot(M8, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#509b4b", high ="#9669b0", mid = "#f7f7f7", midpoint = 0, limit = c(-1, 1)) + 
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), y= expression(paste("Strength of Competititon (", beta, ")")), title = "Change in Two Species VRs with Weak Invader", fill=expression(paste("",Delta,"VR"))) +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95"))

#Figure 4.b: Strong Invader in 2 Species Community----
M9 <- melt(VR_post_good - VR_pre_good)

plot9 <- ggplot(M9, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#509b4b", high ="#9669b0", mid = "#f7f7f7", midpoint = 0, limit = c(-1, 1.5)) + 
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), y= expression(paste("Strength of Competititon (", beta, ")")), title = "Change in Two Species VRs with Strong Invader", fill=expression(paste("",Delta,"VR"))) +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95"))

#Figure 4.c: Weak Invader in 3 Species Community----
M10 <- melt(VR_post_poor_full - VR_pre_poor_full)

plot10 <- ggplot(M10, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#509b4b", high ="#9669b0", mid = "#f7f7f7", midpoint = 0, limit = c(-1, 1)) + 
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), y= expression(paste("Strength of Competititon (", beta, ")")), title = "Change in Three Species VRs with Weak Invader", fill=expression(paste("",Delta,"VR"))) +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95"))

#Figure 4.d: Strong Invader in 3 Species Community----
M11 <- melt(VR_post_good_full - VR_pre_good_full)

plot11 <- ggplot(M11, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#509b4b", high ="#9669b0", mid = "#f7f7f7", midpoint = 0, limit = c(-1.25, 1)) + 
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), y= expression(paste("Strength of Competititon (", beta, ")")), title = "Change in Three Species VRs with Strong Invader", fill=expression(paste("",Delta,"VR"))) +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95"))

quartz(width=11, height=5)

ggarrange(
  plot8, plot9, plot10, plot11, labels = c("A", "B", "C", "D"),
  common.legend = TRUE, legend = "right", nrow = 2, ncol = 2
)

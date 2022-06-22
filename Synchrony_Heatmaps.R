#####################################
# This script produces a figure that shows the range in 
# variance ratios of resident communities.
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
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95")) +
  theme(axis.title=element_text(size=14), #change axis title size
        axis.text=element_text(size=12), #change axid text size
        plot.title = element_text(size=18), #change plot title size
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.title = element_text(size=12), #change legend title font size
        legend.text = element_text(size=10)) #change legend text font size)

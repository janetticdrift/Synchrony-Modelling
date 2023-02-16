###
###Resubmission Analyses
###

#Question 5: Speciose Communities#----

# Package for calculating the variance ratio
library(codyn)

#Packages for creating plots
library(ggplot2)
library(reshape2)
library(ggpubr)
library(tidyverse)

#Number of species in the community
species <- 2
a <- species

#Species and Community Parameters
env_condition <- seq(from=0, to=.25, by=.01)
beta_range <- seq(from=0, to=.95, by=.05)
time <- 200 
burn_in <- 100
runs <- 200

# K <- round(runif(species, min=1000, max=1500))
K <- c(1000, 1500)
# r <- round(runif(species, min=0.1, max=0.9), 1)
r <- c(0.5, 0.8)
sigmaD <- 1

#Variance Ratio Empty Output
VR <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))

for (x in 1:length(env_condition)) {
  for (y in 1:length(beta_range)) {
    
    env <- env_condition[x]
    beta <- beta_range[y]
    
    VR_current <- rep(NA, runs)
    
    for (z in 1:runs) {
      
      #Set starting abundances
      N <- matrix(nrow = time, ncol = species)
      colnames(N) <- c(1:species)
      
      N[1,] <- 50 

      #Create competition coefficients
      beta_vector <- rnorm(species*species, mean=beta_range[y], sd=0.0) #Set sd = 0, species to 2
      beta_vector <- ifelse(beta_vector<0, 0, beta_vector)
      beta_matrix <- matrix(data=beta_vector, nrow = species, ncol = species)
      diag(beta_matrix) <- 1
      
      #Create environmental effect
      sigmaE <- -env
      
      #Create environmental variation
      miuE <- rnorm(time, mean = 0, sd = 1)
      
      #Create demographic variation
      miuD <- rnorm(time, mean = 0, sd = 1)
      
    for (t in 1:(time-1)) { # for each species being the focal species
      for (s in 1:species) {
        N[t+1,s] <- N[t,s]*exp(r[s]*(1-sum(beta_matrix[s,]*N[t,]/K[s])) +
                                 (sigmaE*miuE[t]) + (sigmaD*miuD[t])/sqrt(N[t,s]))
      
        if(is.nan(N[t+1,s])) {
          N[t+1,s] <- 0
        }
        if(N[t+1,s] < 1) {
          N[t+1,s] <- 0
        }
       }
     }
      N <- as.data.frame(N)
      N.clean <- N %>%
        slice(-(1:burn_in)) %>%
        gather(key = "species", value = "abundance") %>%
        mutate(time = rep(c(1:(time-burn_in)), times = a))

      VR_temp <- variance_ratio(N.clean, time.var = "time", species.var = "species",
                                abundance.var = "abundance", bootnumber = 1)
      VR_current[z] <- VR_temp$VR
    }
    VR[x,y] <- mean(VR_current)
  }
  print(x/length(env_condition))
}

# N[t+1,s] <- N[t,s]*exp(r[s]*(1-sum(beta_matrix[s,]*N[t,]/K)) +
#                          (sigmaE*miuE[t]) + (sigmaD*miuD[t])/sqrt(N[t,s]))
#With Extra numbers added
# N[t,s]*exp(r[s]*(1-(N[t,s]/K[s]) - sum(beta_matrix[s,]*N[t,]/K[s])) +
#              (sigmaE*miuE[t]) + (sigmaD*miuD[t])/sqrt(N[t,s]))

# Create heat map
M1 <- melt(VR) 

ggplot(M1, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#008080", high ="#ca562c", mid = "#f6edbd", midpoint = 1, limit = c(0,2)) +
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), 
       y= expression(paste("Strength of Competititon (", beta, ")")), 
       title = "Resident Communities' Variance Ratios", fill="VR") +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95")) +
  theme(axis.title=element_text(size=14), #change axis title size
        axis.text=element_text(size=12), #change axis text size
        plot.title = element_text(size=18), #change plot title size
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.title = element_text(size=12), #change legend title font size
        legend.text = element_text(size=10)) #change legend text font size)

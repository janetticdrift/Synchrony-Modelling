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

#Number of species in the community
species <- 15

#Species and Community Parameters
env_condition <- seq(from=0, to=.25, by=.01)
beta_range <- seq(from=0, to=.95, by=.05)
time <- 300
runs <- 50

K <- round(runif(species, min=1000, max=1500))
r <- round(runif(species, min=0.1, max=0.9), 1)
sigmaD <- 1

for (x in 1:length(env_condition)) {
  for (y in 1:length(beta_range)) {
    
    env <- env_condition[x]
    beta <- beta_range[y]
    
    VR_current <- rep(NA, runs)
    
    for (z in 1:runs) {
      
      #Set starting abundances
      N <- matrix(nrow = time, ncol = species)
      
      N[1,] <- 50 

      #Create competition coefficients
      beta_vector <- rnorm(species*species, mean=beta_range[y], sd=0.05)
      ifelse(beta_vector<0, 0, beta_vector)
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
        N[t+1,s] <- N[t,s]*exp(r[s]*(1-sum(beta_matrix[s,]*N[t,]/K)) +
                                 (sigmaE*miuE[t]) + (sigmaD*miuD[t])/sqrt(N[t,s]))
      
        # if(is.nan(N[t+1,s])) {
        #   N[t+1,s] <- 0
        # }
      }
    }  
    }
  }
  print(x/length(env_condition))
}


#Original code
for (s in 1:length(species)) { # for each species being the focal species
  N[t+1,s] <- N[t,s]*exp(r[s]*(1-sum(beta_matrix[s,]*N[t,]/K)) +
                           (sigmaE*miuE[t]) + (sigmaD*miuD[t])/sqrt(N[t,s]))

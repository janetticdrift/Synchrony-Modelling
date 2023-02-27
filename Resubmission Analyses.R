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

#r <- c(0.5, 0.8)
sigmaD <- 1

###Residents Only Code###----

#Variance Ratio Empty Output
VR <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))
number_species <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))

for (x in 1:length(env_condition)) {
  for (y in 1:length(beta_range)) {
    
    env <- env_condition[x]
    beta <- beta_range[y]
    
    VR_current <- rep(NA, runs)
    total_species <- rep(NA, runs)
    
    for (z in 1:runs) {
      
      K <- round(runif(species, min=1000, max=1500))
      r <- round(runif(species, min=0.1, max=0.9), 1)
      
      #Set starting abundances
      N <- matrix(nrow = time, ncol = species)
      colnames(N) <- c(1:species)
      
      N[1,] <- 50 

      #Create competition coefficients
      beta_vector <- rnorm(species*species, mean=beta_range[y], sd=0.5) #Set sd = 0, species to 2
      beta_vector <- ifelse(beta_vector<0, 0, beta_vector)
      beta_matrix <- matrix(data=beta_vector, nrow = species, ncol = species)
      diag(beta_matrix) <- 1
      
      #Create environmental effect
      sigmaE <- -env
      
      #Create environmental variation
      miuE <- rnorm(time, mean = 0, sd = 1)
      
      #Create demographic variation
      miuD <- rnorm(species*time, mean = 0, sd = 1)
      miuD <- matrix(data=miuD, nrow = time, ncol = species)
      
    for (t in 1:(time-1)) { # for each species being the focal species
      for (s in 1:species) {
        N[t+1,s] <- N[t,s]*exp(r[s]*(1-sum(beta_matrix[s,]*N[t,]/K)) +
                                 (sigmaE*miuE[t]) + (sigmaD*miuD[t,s])/sqrt(N[t,s]))
      
        if(is.nan(N[t+1,s])) {
          N[t+1,s] <- 0
        }
        if(N[t+1,s] < 1) {
          N[t+1,s] <- 0
        }
       }
    }
      
      extinct <- which(N[time,]==0) # determine which species go extinct
      number_extinct <- length(extinct)
      
      N <- as.data.frame(N)
      N.clean <- N %>%
        slice(-(1:burn_in)) %>%
        select(-c(extinct)) %>%
        gather(key = "species", value = "abundance") %>%
        mutate(time = rep(c(1:(time-burn_in)), times = (a-number_extinct)))
      
    if (species-number_extinct > 1){
        VR_temp <- variance_ratio(N.clean, time.var = "time", species.var = "species",
                                abundance.var = "abundance", bootnumber = 1)
        total_species[z] <- species-number_extinct
        VR_current[z] <- VR_temp$VR } else {
        VR_current[z] <- NA
        total_species[z] <- NA
        }
      
    }
    
    VR[x,y] <- mean(VR_current, na.rm = TRUE)
    number_species[x,y] <- mean(total_species, na.rm = TRUE)
  }
  print(x/length(env_condition))
}


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

M2 <- melt(number_species) 

ggplot(M2, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#008080", high ="#ca562c", mid = "#f6edbd", midpoint = 2/2, limit = c(0,2)) +
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), 
       y= expression(paste("Strength of Competititon (", beta, ")")), 
       title = "Number of Species", fill="Species") +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95")) +
  theme(axis.title=element_text(size=14), #change axis title size
        axis.text=element_text(size=12), #change axis text size
        plot.title = element_text(size=18), #change plot title size
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.title = element_text(size=12), #change legend title font size
        legend.text = element_text(size=10)) #change legend text font size)

###Add Invader Resistance Code###----
#Growth Rate Empty Output
avg_lambda_good <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))
success_lambda_good <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))

for (x in 1:length(env_condition)) {
  for (y in 1:length(beta_range)) {
    
    env <- env_condition[x]
    beta <- beta_range[y]
    
    lambda_good_a <- rep(NA, runs) 
    lambda_good <- matrix(NA, nrow = time-burn_in, ncol = runs)
    
    for (z in 1:runs) {
      
      K <- round(runif(species, min=1000, max=1500))
      r <- round(runif(species, min=0.1, max=0.9), 1)
      
      #Strong invader parameters
      Ki <- 1000
      ri <- 0.7
      #Weak invader parameters
      # Ki <- 900
      # ri <- 0.4
      
      #Set starting resident abundances
      N <- matrix(nrow = time, ncol = species)
      colnames(N) <- c(1:species)
      
      N[1,] <- 50 
      
      #Set invader vector
      Ni <- rep(NA, (time-burn_in))
      
      Ni[1] <- 1 
      invader_abund_good <- 1
      
      #Create competition coefficients for residents
      beta_vector <- rnorm(species*species, mean=beta_range[y], sd=0.5) #Set sd = 0, species to 2
      beta_vector <- ifelse(beta_vector<0, 0, beta_vector)
      beta_matrix <- matrix(data=beta_vector, nrow = species, ncol = species)
      diag(beta_matrix) <- 1
      
      #Create environmental effect
      sigmaE <- -env
      
      #Create environmental variation
      miuE <- rnorm(time, mean = 0, sd = 1)
      
      #Create demographic variation
      miuD <- rnorm(species*time, mean = 0, sd = 1)
      miuD <- matrix(data=miuD, nrow = time, ncol = species)
      miuDi <- rnorm((time-burn_in), mean = 0, sd = 1)
      
      for (t in 1:(time-1)) { # for each species being the focal species
        for (s in 1:species) {
          N[t+1,s] <- N[t,s]*exp(r[s]*(1-sum(beta_matrix[s,]*N[t,]/K)) +
                                   (sigmaE*miuE[t]) + (sigmaD*miuD[t,s])/sqrt(N[t,s]))
          
          if(is.nan(N[t+1,s])) {
            N[t+1,s] <- 0
          }
          if(N[t+1,s] < 1) {
            N[t+1,s] <- 0
          }
        }
      }
      
      counter <- 1
      
      #Model effects of good invasion
      for (t in burn_in:(time-1)) {
        
        Ni[counter] <- invader_abund_good*exp(ri*(1-(invader_abund_good/Ki) - (beta31*N1[t]/K1) - (beta32*N2[t]/K2))
                                              + (sigmaE3*miuE[t]) + (sigmaD3*miuD3[counter])/sqrt(invader_abund_good))
        # Needed for when N4 is 0 for including demographic stochasticity
        # as sqrt(0) yields NaN
        
        if(is.nan(N3[counter])) {
          N3[counter] <- 0
        }
        if(N3[counter] < 1) {
          N3[counter] <- 0
        }
        
        #Calculate invader's growth rate
        lambda_good[counter, z] <- log(N3[counter]/invader_abund_good)
        
        counter <- counter + 1
      }
      lambda_good_a[z] <- mean(lambda_good[,z])
    }
    
    avg_lambda_good[x,y] <- mean(lambda_good_a)
    success_lambda_good[x,y] <- length(which(lambda_good>0))/((time-burn_in)*runs)
  }
  print(x/length(env_condition))
}


# Create heat map
M1 <- melt(VR) 

ggplot(M1, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#008080", high ="#ca562c", mid = "#f6edbd", midpoint = 1, limit = c(0,5)) +
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

###Add Invader Resilience Code###----

#Number of species in the community
species <- 10
a <- species

#Species and Community Parameters
env_condition <- seq(from=0, to=.25, by=.01)
beta_range <- seq(from=0, to=.95, by=.05)
time <- 200 
burn_in <- 100
runs <- 200

sigmaD <- 1

#Invader Growth Rate Empty Output
VR <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))
number_species <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))

for (x in 1:length(env_condition)) {
  for (y in 1:length(beta_range)) {
    
    env <- env_condition[x]
    beta <- beta_range[y]
    
    VR_current <- rep(NA, runs)
    total_species <- rep(NA, runs)
    
    for (z in 1:runs) {
      
      K <- round(runif(species, min=1000, max=1500))
      #Add invader K
      K <- append(K, 1000)
      r <- round(runif(species, min=0.1, max=0.9), 1)
      #Add invader R
      r <- append(r, 0.7)
      
      #Set starting abundances
      N <- matrix(nrow = time, ncol = species+1)
      colnames(N) <- c(1:(species + 1))
      
      N[1,1:species] <- 50 
      N[1, species+1] <- 0
      
      #Create competition coefficients
      beta_vector <- rnorm(species*species, mean=beta_range[y], sd=0.5)
      beta_matrix <- matrix(data=beta_vector, nrow = species, ncol = species)
      beta_matrix <- ifelse(beta_matrix<0, 0, beta_matrix)
      
      #Add invader competition coefficients
      beta_ri <- sample(seq(from=0.5, to=.95, by=.05), species+1, replace = TRUE) #Effect of invader on residents
      beta_ir <- sample(seq(from=0, to=.5, by=.05), species, replace = TRUE) #Effect of residents on invader
      beta_matrix <- beta_matrix %>%
        rbind(beta_ir) %>%
        cbind(beta_ri)
#        set_colnames(c(1:(species + 1))) 
#        set_rownames(c(1:(species + 1))) maybe unnecessary
      diag(beta_matrix) <- 1
      
      #Create environmental effect
      sigmaE <- -env
      
      #Create environmental variation
      miuE <- rnorm(time, mean = 0, sd = 1)
      
      #Create demographic variation
      miuD <- rnorm((species+1)*time, mean = 0, sd = 1)
      miuD <- matrix(data=miuD, nrow = time, ncol = species+1)
      
      for (t in 1:(time-1)) { # for each species being the focal species
        for (s in 1:(species+1)) {
          if(t == burn_in) {
            N[t,species+1] <- 1
          }
          N[t+1,s] <- N[t,s]*exp(r[s]*(1-sum(beta_matrix[s,]*N[t,]/K)) +
                                   (sigmaE*miuE[t]) + (sigmaD*miuD[t,s])/sqrt(N[t,s]))
          
          if(is.nan(N[t+1,s])) {
            N[t+1,s] <- 0
          }
          if(N[t+1,s] < 1) {
            N[t+1,s] <- 0
          }
        }
      }
      
      extinct <- which(N[time,]==0) # determine which species go extinct
      number_extinct <- length(extinct)
      
      N <- as.data.frame(N)
      N.clean <- N %>%
        slice(-(1:burn_in)) %>%
        select(-c(extinct)) %>%
        gather(key = "species", value = "abundance") %>%
        mutate(time = rep(c(1:(time-burn_in)), times = ((a+1)-number_extinct)))
      
      if ((species+1)-number_extinct > 1){
        VR_temp <- variance_ratio(N.clean, time.var = "time", species.var = "species",
                                  abundance.var = "abundance", bootnumber = 1)
        total_species[z] <- (species+1)-number_extinct
        VR_current[z] <- VR_temp$VR } else {
          VR_current[z] <- NA
          total_species[z] <- NA
        }
      
    }
    
    VR[x,y] <- mean(VR_current, na.rm = TRUE)
    number_species[x,y] <- mean(total_species, na.rm = TRUE)
  }
  print(x/length(env_condition))
}


# Create heat map
M1 <- melt(VR) 

ggplot(M1, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#008080", high ="#ca562c", mid = "#f6edbd", midpoint = 1, limit = c(0,3)) +
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

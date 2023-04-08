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

#Species and Community Parameters
env_condition <- seq(from=0, to=.25, by=.01)
beta_range <- seq(from=0, to=.95, by=.05)
time <- 200 
burn_in <- 100
runs <- 2000

sigmaD <- 1

###Residents Only Code###----

#2 species code
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
      
      K <- c(1000, 1500) #Set to specific values to match the original values used
      r <- c(0.5, 0.8)
      
      #Set starting abundances
      N <- matrix(nrow = time, ncol = species)
      colnames(N) <- c(1:species)
      
      N[1,] <- 50 

      #Create competition coefficients
      beta_vector <- rnorm(species*species, mean=beta_range[y], sd=0) #Set sd = 0, species to 2
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
M_2sp <- melt(VR) 

plot2sp <- ggplot(M_2sp, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#008080", high ="#ca562c", mid = "#f6edbd", midpoint = 1, limit = c(0,2)) +
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), 
       y= expression(paste("Strength of Competititon (", beta, ")")), 
       title = "Synchrony with Two Residents", fill="VR") +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95")) +
  theme(axis.title=element_text(size=10), #change axis title size
        axis.text=element_text(size=12), #change axis tick size
        plot.title = element_text(size=16), #change plot title size
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=10)) #change legend text font size)

M_2num <- melt(number_species) 

plot2num <- ggplot(M_2num, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#008080", high ="#ca562c", mid = "#f6edbd", midpoint = 2/2, limit = c(0,2)) +
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), 
       y= expression(paste("Strength of Competititon (", beta, ")")), 
       title = "Number of Species Remaining", fill="Species") +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95")) +
  theme(axis.title=element_text(size=10), #change axis title size
        axis.text=element_text(size=12), #change axis tick size
        plot.title = element_text(size=16), #change plot title size
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=10)) #change legend text font size)

#5 species code
#Number of species in the community
species <- 5
a <- species
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
M_5sp <- melt(VR) 

plot5sp <- ggplot(M_5sp, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#008080", high ="#ca562c", mid = "#f6edbd", midpoint = 1, limit = c(0,3.021)) +
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), 
       y= expression(paste("Strength of Competititon (", beta, ")")), 
       title = "Synchrony with Five Residents", fill="VR") +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95")) +
  theme(axis.title=element_text(size=10), #change axis title size
        axis.text=element_text(size=12), #change axis tick size
        plot.title = element_text(size=16), #change plot title size
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=10)) #change legend text font size)

M_5num <- melt(number_species) 

plot5num <- ggplot(M_5num, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#008080", high ="#ca562c", mid = "#f6edbd", midpoint = 5/2, limit = c(0,5)) +
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), 
       y= expression(paste("Strength of Competititon (", beta, ")")), 
       title = "Number of Species Remaining", fill="Species") +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95")) +
  theme(axis.title=element_text(size=10), #change axis title size
        axis.text=element_text(size=12), #change axis tick size
        plot.title = element_text(size=16), #change plot title size
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=10)) #change legend text font size)

#5 species code
#Number of species in the community
species <- 10
a <- species
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
M_10sp <- melt(VR) 

plot10sp <- ggplot(M_10sp, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#008080", high ="#ca562c", mid = "#f6edbd", midpoint = 1, limit = c(0,5)) +
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), 
       y= expression(paste("Strength of Competititon (", beta, ")")), 
       title = "Synchrony with Ten Residents", fill="VR") +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95")) +
  theme(axis.title=element_text(size=10), #change axis title size
        axis.text=element_text(size=12), #change axis tick size
        plot.title = element_text(size=16), #change plot title size
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=10)) #change legend text font size)

M_10num <- melt(number_species) 

plot10num <- ggplot(M_10num, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#008080", high ="#ca562c", mid = "#f6edbd", midpoint = 10/2, limit = c(0,10)) +
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), 
       y= expression(paste("Strength of Competititon (", beta, ")")), 
       title = "Number of Species Remaining", fill="Species") +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95")) +
  theme(axis.title=element_text(size=10), #change axis title size
        axis.text=element_text(size=12), #change axis tick size
        plot.title = element_text(size=16), #change plot title size
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=10)) #change legend text font size)

#Combine the 6 plots into one figure
ggarrange(plot2sp, plot2num, plot5sp, plot5num, plot10sp, plot10num,
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 2, nrow = 3)

###Invader Resistance Code###----
sigmaD <- 1 #Residents 
sigmaDi <- 0 #Invaders

#2 species Weak Invader
species <- 2
a <- species
avg_lambda <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))
success_lambda_weak <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))

for (x in 1:length(env_condition)) {
  for (y in 1:length(beta_range)) {
    
    env <- env_condition[x]
    beta <- beta_range[y]
    
    lambda_a <- rep(NA, runs) 
    lambda <- matrix(NA, nrow = time-burn_in, ncol = runs)
    
    for (z in 1:runs) {
      
      #Multi-species resident parameters
      # K <- round(runif(species, min=1000, max=1500))
      # r <- round(runif(species, min=0.1, max=0.9), 1)
      #2 species resident parameters
      K <- c(1000, 1500) 
      r <- c(0.5, 0.8)
      
      #Strong invader parameters
      # Ki <- 1000
      # ri <- 0.7
      #Weak invader parameters
      Ki <- 900
      ri <- 0.4
      
      #Set starting resident abundances
      N <- matrix(nrow = time, ncol = species)
      colnames(N) <- c(1:species)
      
      N[1,] <- 50 
      
      #Set invader vector
      Ni <- rep(NA, (time-burn_in))
      
      Ni[1] <- 1 
      invader_abund <- 1
      
      #Create competition coefficients for residents
      beta_vector <- rnorm(species*species, mean=beta_range[y], sd=0) #Set sd = 0, species to 2
      beta_vector <- ifelse(beta_vector<0, 0, beta_vector)
      beta_matrix <- matrix(data=beta_vector, nrow = species, ncol = species)
      diag(beta_matrix) <- 1
      
      #Create competition coefficients for strong invader
      # beta_res_on_invader <- rep(.5, 2)
      #Create competition coefficients for weak invader
      beta_res_on_invader <- rep(.6, 2)
      
      #Create environmental effect
      sigmaE <- -env
      sigmaEi <- -0.06
      
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
        
        Ni[counter] <- invader_abund*exp(ri*(1-(invader_abund/Ki) - 
                                                    sum(beta_res_on_invader*N[t,]/K)) + 
                                                (sigmaEi*miuE[t]) + 
                                                (sigmaDi*miuDi[counter])/sqrt(invader_abund))
        # Needed for when N4 is 0 for including demographic stochasticity
        # as sqrt(0) yields NaN
        
        if(is.nan(Ni[counter])) {
          Ni[counter] <- 0
        }
        if(Ni[counter] < 1) {
          Ni[counter] <- 0
        }
        
        #Calculate invader's growth rate
        lambda[counter, z] <- log(Ni[counter]/invader_abund)
        
        counter <- counter + 1
      }
      lambda_a[z] <- mean(lambda[,z])
    }
    
    avg_lambda[x,y] <- mean(lambda_a)
    success_lambda_weak[x,y] <- length(which(lambda>0))/((time-burn_in)*runs)
  }
  print(x/length(env_condition))
}

#2 species Strong Invader
#Growth Rate Parameters and Empty Output
avg_lambda <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))
success_lambda_good <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))

for (x in 1:length(env_condition)) {
  for (y in 1:length(beta_range)) {
    
    env <- env_condition[x]
    beta <- beta_range[y]
    
    lambda_a <- rep(NA, runs) 
    lambda <- matrix(NA, nrow = time-burn_in, ncol = runs)
    
    for (z in 1:runs) {
      
      #Multi-species resident parameters
      # K <- round(runif(species, min=1000, max=1500))
      # r <- round(runif(species, min=0.1, max=0.9), 1)
      #2 species resident parameters
      K <- c(1000, 1500) 
      r <- c(0.5, 0.8)
      
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
      invader_abund <- 1
      
      #Create competition coefficients for residents
      beta_vector <- rnorm(species*species, mean=beta_range[y], sd=0) #Set sd = 0, species to 2
      beta_vector <- ifelse(beta_vector<0, 0, beta_vector)
      beta_matrix <- matrix(data=beta_vector, nrow = species, ncol = species)
      diag(beta_matrix) <- 1
      
      #Create competition coefficients for strong invader
      beta_res_on_invader <- rep(.5, 2)
      #Create competition coefficients for weak invader
      # beta_res_on_invader <- rep(.6, 2)
        
      #Create environmental effect
      sigmaE <- -env
      sigmaEi <- -0.1
      
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
        
        Ni[counter] <- invader_abund*exp(ri*(1-(invader_abund/Ki) - 
                                                    sum(beta_res_on_invader*N[t,]/K)) + 
                                                (sigmaEi*miuE[t]) + 
                                                (sigmaDi*miuDi[counter])/sqrt(invader_abund))
        # Needed for when N4 is 0 for including demographic stochasticity
        # as sqrt(0) yields NaN
        
        if(is.nan(Ni[counter])) {
          Ni[counter] <- 0
        }
        if(Ni[counter] < 1) {
          Ni[counter] <- 0
        }
        
        #Calculate invader's growth rate
        lambda[counter, z] <- log(Ni[counter]/invader_abund)
        
        counter <- counter + 1
      }
      lambda_a[z] <- mean(lambda[,z])
    }
    
    avg_lambda[x,y] <- mean(lambda_a)
    success_lambda_good[x,y] <- length(which(lambda>0))/((time-burn_in)*runs)
  }
  print(x/length(env_condition))
}


# Create heat maps for 2species communities
M_2resisstrong <- melt(success_lambda_good) 
M_2resisweak <- melt(success_lambda_weak)

min_lim <- min(success_lambda_good, success_lambda_weak)
max_lim <- max(success_lambda_good, success_lambda_weak)

plot2spstrong <- ggplot(M_2resisstrong, aes(x=Var1, y=Var2, fill=value*100)) + 
  geom_tile() + 
  scale_fill_distiller(palette = "RdBu", limits = c(0, 100)) +
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), 
  y= expression(paste("Strength of Competititon (", beta, ")")), title = "Growth Rates of Strong Invader with 2 Residents", fill="Growth
Success (%)") +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95")) +
  theme(axis.title=element_text(size=10), #change axis title size
        axis.text=element_text(size=12), #change axis tick size
        plot.title = element_text(size=14), #change plot title size
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=10)) #change legend text font size)

plot2spweak <- ggplot(M_2resisweak, aes(x=Var1, y=Var2, fill=value*100)) + 
  geom_tile() + 
  scale_fill_distiller(palette = "RdBu", limits = c(0, 100)) +
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), 
  y= expression(paste("Strength of Competititon (", beta, ")")), title = "Growth Rates of Weak Invader with 2 Residents", fill="Growth
Success (%)") +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95")) +
  theme(axis.title=element_text(size=10), #change axis title size
        axis.text=element_text(size=12), #change axis tick size
        plot.title = element_text(size=14), #change plot title size
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=10)) #change legend text font size)

#5 species Weak Invader
species <- 5
a <- species
avg_lambda <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))
success_lambda_weak <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))

for (x in 1:length(env_condition)) {
  for (y in 1:length(beta_range)) {
    
    env <- env_condition[x]
    beta <- beta_range[y]
    
    lambda_a <- rep(NA, runs) 
    lambda <- matrix(NA, nrow = time-burn_in, ncol = runs)
    
    for (z in 1:runs) {
      
      #Multi-species resident parameters
      K <- round(runif(species, min=1000, max=1500))
      r <- round(runif(species, min=0.1, max=0.9), 1)
      #2 species resident parameters
      # K <- c(1000, 1500) 
      # r <- c(0.5, 0.8)
      
      #Strong invader parameters
      # Ki <- 1000
      # ri <- 0.7
      #Weak invader parameters
      Ki <- 900
      ri <- 0.4
      
      #Set starting resident abundances
      N <- matrix(nrow = time, ncol = species)
      colnames(N) <- c(1:species)
      
      N[1,] <- 50 
      
      #Set invader vector
      Ni <- rep(NA, (time-burn_in))
      
      Ni[1] <- 1 
      invader_abund <- 1
      
      #Create competition coefficients for residents
      beta_vector <- rnorm(species*species, mean=beta_range[y], sd=0.5) #Set sd = 0, species to 2
      beta_vector <- ifelse(beta_vector<0, 0, beta_vector)
      beta_matrix <- matrix(data=beta_vector, nrow = species, ncol = species)
      diag(beta_matrix) <- 1
      
      #Create competition coefficients for strong invader
      # beta_res_on_invader <- rep(.5, 2)
      #Create competition coefficients for weak invader
      beta_res_on_invader <- rep(.6, 2)
      
      #Create environmental effect
      sigmaE <- -env
      sigmaEi <- -0.06
      
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
        
        Ni[counter] <- invader_abund*exp(ri*(1-(invader_abund/Ki) - 
                                               sum(beta_res_on_invader*N[t,]/K)) + 
                                           (sigmaEi*miuE[t]) + 
                                           (sigmaDi*miuDi[counter])/sqrt(invader_abund))
        
        if(is.nan(Ni[counter])) {
          Ni[counter] <- 0
        }
        if(Ni[counter] < 1) {
          Ni[counter] <- 0
        }
        
        #Calculate invader's growth rate
        lambda[counter, z] <- log(Ni[counter]/invader_abund)
        
        counter <- counter + 1
      }
      lambda_a[z] <- mean(lambda[,z])
    }
    
    avg_lambda[x,y] <- mean(lambda_a)
    success_lambda_weak[x,y] <- length(which(lambda>0))/((time-burn_in)*runs)
  }
  print(x/length(env_condition))
}

#5 species Strong Invader
#Growth Rate Parameters and Empty Output
avg_lambda <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))
success_lambda_good <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))

for (x in 1:length(env_condition)) {
  for (y in 1:length(beta_range)) {
    
    env <- env_condition[x]
    beta <- beta_range[y]
    
    lambda_a <- rep(NA, runs) 
    lambda <- matrix(NA, nrow = time-burn_in, ncol = runs)
    
    for (z in 1:runs) {
      
      #Multi-species resident parameters
      K <- round(runif(species, min=1000, max=1500))
      r <- round(runif(species, min=0.1, max=0.9), 1)
      #2 species resident parameters
      # K <- c(1000, 1500) 
      # r <- c(0.5, 0.8)
      
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
      invader_abund <- 1
      
      #Create competition coefficients for residents
      beta_vector <- rnorm(species*species, mean=beta_range[y], sd=0.5) #Set sd = 0, species to 2
      beta_vector <- ifelse(beta_vector<0, 0, beta_vector)
      beta_matrix <- matrix(data=beta_vector, nrow = species, ncol = species)
      diag(beta_matrix) <- 1
      
      #Create competition coefficients for strong invader
      beta_res_on_invader <- rep(.5, 2)
      #Create competition coefficients for weak invader
      # beta_res_on_invader <- rep(.6, 2)
      
      #Create environmental effect
      sigmaE <- -env
      sigmaEi <- -0.1
      
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
        
        Ni[counter] <- invader_abund*exp(ri*(1-(invader_abund/Ki) - 
                                               sum(beta_res_on_invader*N[t,]/K)) + 
                                           (sigmaEi*miuE[t]) + 
                                           (sigmaDi*miuDi[counter])/sqrt(invader_abund))
        # Needed for when N4 is 0 for including demographic stochasticity
        # as sqrt(0) yields NaN
        
        if(is.nan(Ni[counter])) {
          Ni[counter] <- 0
        }
        if(Ni[counter] < 1) {
          Ni[counter] <- 0
        }
        
        #Calculate invader's growth rate
        lambda[counter, z] <- log(Ni[counter]/invader_abund)
        
        counter <- counter + 1
      }
      lambda_a[z] <- mean(lambda[,z])
    }
    
    avg_lambda[x,y] <- mean(lambda_a)
    success_lambda_good[x,y] <- length(which(lambda>0))/((time-burn_in)*runs)
  }
  print(x/length(env_condition))
}


# Create heat maps for 2species communities
M_5resisstrong <- melt(success_lambda_good) 
M_5resisweak <- melt(success_lambda_weak)

min_lim <- min(success_lambda_good, success_lambda_weak)
max_lim <- max(success_lambda_good, success_lambda_weak)

plot5spstrong <- ggplot(M_5resisstrong, aes(x=Var1, y=Var2, fill=value*100)) + 
  geom_tile() + 
  scale_fill_distiller(palette = "RdBu", limits = c(0, 100)) +
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), 
       y= expression(paste("Strength of Competititon (", beta, ")")), title = "Growth Rates of Strong Invader with 5 Residents", fill="Growth
Success (%)") +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95")) +
  theme(axis.title=element_text(size=10), #change axis title size
        axis.text=element_text(size=12), #change axis tick size
        plot.title = element_text(size=14), #change plot title size
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=10)) #change legend text font size)

plot5spweak <- ggplot(M_5resisweak, aes(x=Var1, y=Var2, fill=value*100)) + 
  geom_tile() + 
  scale_fill_distiller(palette = "RdBu", limits = c(0, 100)) +
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), 
       y= expression(paste("Strength of Competititon (", beta, ")")), title = "Growth Rates of Weak Invader with 5 Residents", fill="Growth
Success (%)") +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95")) +
  theme(axis.title=element_text(size=10), #change axis title size
        axis.text=element_text(size=12), #change axis tick size
        plot.title = element_text(size=14), #change plot title size
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=10)) #change legend text font size)

#10 species Weak Invader
species <- 10
a <- species
avg_lambda <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))
success_lambda_weak <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))

for (x in 1:length(env_condition)) {
  for (y in 1:length(beta_range)) {
    
    env <- env_condition[x]
    beta <- beta_range[y]
    
    lambda_a <- rep(NA, runs) 
    lambda <- matrix(NA, nrow = time-burn_in, ncol = runs)
    
    for (z in 1:runs) {
      
      #Multi-species resident parameters
      K <- round(runif(species, min=1000, max=1500))
      r <- round(runif(species, min=0.1, max=0.9), 1)
      #2 species resident parameters
      # K <- c(1000, 1500) 
      # r <- c(0.5, 0.8)
      
      #Strong invader parameters
      # Ki <- 1000
      # ri <- 0.7
      #Weak invader parameters
      Ki <- 900
      ri <- 0.4
      
      #Set starting resident abundances
      N <- matrix(nrow = time, ncol = species)
      colnames(N) <- c(1:species)
      
      N[1,] <- 50 
      
      #Set invader vector
      Ni <- rep(NA, (time-burn_in))
      
      Ni[1] <- 1 
      invader_abund <- 1
      
      #Create competition coefficients for residents
      beta_vector <- rnorm(species*species, mean=beta_range[y], sd=0.5) #Set sd = 0, species to 2
      beta_vector <- ifelse(beta_vector<0, 0, beta_vector)
      beta_matrix <- matrix(data=beta_vector, nrow = species, ncol = species)
      diag(beta_matrix) <- 1
      
      #Create competition coefficients for strong invader
      # beta_res_on_invader <- rep(.5, 2)
      #Create competition coefficients for weak invader
      beta_res_on_invader <- rep(.6, 2)
      
      #Create environmental effect
      sigmaE <- -env
      sigmaEi <- -0.06
      
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
        
        Ni[counter] <- invader_abund*exp(ri*(1-(invader_abund/Ki) - 
                                               sum(beta_res_on_invader*N[t,]/K)) + 
                                           (sigmaEi*miuE[t]) + 
                                           (sigmaDi*miuDi[counter])/sqrt(invader_abund))
        
        if(is.nan(Ni[counter])) {
          Ni[counter] <- 0
        }
        if(Ni[counter] < 1) {
          Ni[counter] <- 0
        }
        
        #Calculate invader's growth rate
        lambda[counter, z] <- log(Ni[counter]/invader_abund)
        
        counter <- counter + 1
      }
      lambda_a[z] <- mean(lambda[,z])
    }
    
    avg_lambda[x,y] <- mean(lambda_a)
    success_lambda_weak[x,y] <- length(which(lambda>0))/((time-burn_in)*runs)
  }
  print(x/length(env_condition))
}

#10 species Strong Invader
#Growth Rate Parameters and Empty Output
avg_lambda <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))
success_lambda_good <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))

for (x in 1:length(env_condition)) {
  for (y in 1:length(beta_range)) {
    
    env <- env_condition[x]
    beta <- beta_range[y]
    
    lambda_a <- rep(NA, runs) 
    lambda <- matrix(NA, nrow = time-burn_in, ncol = runs)
    
    for (z in 1:runs) {
      
      #Multi-species resident parameters
      K <- round(runif(species, min=1000, max=1500))
      r <- round(runif(species, min=0.1, max=0.9), 1)
      #2 species resident parameters
      # K <- c(1000, 1500) 
      # r <- c(0.5, 0.8)
      
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
      invader_abund <- 1
      
      #Create competition coefficients for residents
      beta_vector <- rnorm(species*species, mean=beta_range[y], sd=0.5) #Set sd = 0, species to 2
      beta_vector <- ifelse(beta_vector<0, 0, beta_vector)
      beta_matrix <- matrix(data=beta_vector, nrow = species, ncol = species)
      diag(beta_matrix) <- 1
      
      #Create competition coefficients for strong invader
      beta_res_on_invader <- rep(.5, 2)
      #Create competition coefficients for weak invader
      # beta_res_on_invader <- rep(.6, 2)
      
      #Create environmental effect
      sigmaE <- -env
      sigmaEi <- -0.1
      
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
        
        Ni[counter] <- invader_abund*exp(ri*(1-(invader_abund/Ki) - 
                                               sum(beta_res_on_invader*N[t,]/K)) + 
                                           (sigmaEi*miuE[t]) + 
                                           (sigmaDi*miuDi[counter])/sqrt(invader_abund))
        # Needed for when N4 is 0 for including demographic stochasticity
        # as sqrt(0) yields NaN
        
        if(is.nan(Ni[counter])) {
          Ni[counter] <- 0
        }
        if(Ni[counter] < 1) {
          Ni[counter] <- 0
        }
        
        #Calculate invader's growth rate
        lambda[counter, z] <- log(Ni[counter]/invader_abund)
        
        counter <- counter + 1
      }
      lambda_a[z] <- mean(lambda[,z])
    }
    
    avg_lambda[x,y] <- mean(lambda_a)
    success_lambda_good[x,y] <- length(which(lambda>0))/((time-burn_in)*runs)
  }
  print(x/length(env_condition))
}


# Create heat maps for 2species communities
M_10resisstrong <- melt(success_lambda_good) 
M_10resisweak <- melt(success_lambda_weak)

min_lim <- min(success_lambda_good, success_lambda_weak)
max_lim <- max(success_lambda_good, success_lambda_weak)

plot10spstrong <- ggplot(M_10resisstrong, aes(x=Var1, y=Var2, fill=value*100)) + 
  geom_tile() + 
  scale_fill_distiller(palette = "RdBu", limits = c(0, 100)) +
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), 
       y= expression(paste("Strength of Competititon (", beta, ")")), title = "Growth Rates of Strong Invader with 10 Residents", fill="Growth
Success (%)") +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95")) +
  theme(axis.title=element_text(size=10), #change axis title size
        axis.text=element_text(size=12), #change axis tick size
        plot.title = element_text(size=14), #change plot title size
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=10)) #change legend text font size)

plot10spweak <- ggplot(M_10resisweak, aes(x=Var1, y=Var2, fill=value*100)) + 
  geom_tile() + 
  scale_fill_distiller(palette = "RdBu", limits = c(0, 100)) +
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), 
       y= expression(paste("Strength of Competititon (", beta, ")")), title = "Growth Rates of Weak Invader with 10 Residents", fill="Growth
Success (%)") +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95")) +
  theme(axis.title=element_text(size=10), #change axis title size
        axis.text=element_text(size=12), #change axis tick size
        plot.title = element_text(size=14), #change plot title size
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=10)) #change legend text font size)

#Create 6-panel figure
ggarrange(plot2spweak, plot2spstrong, plot5spweak, plot5spstrong, plot10spweak, plot10spstrong,
          labels = c("A", "B", "C", "D", "E", "F"), ncol = 2, nrow = 3)

###Invader Resilience Code###----

#Species and Community Parameters
env_condition <- seq(from=0, to=.25, by=.01)
beta_range <- seq(from=0, to=.95, by=.05)
time <- 200 
burn_in <- 100
runs <- 2000
sigmaD <- 1

#Only calculate VR if growth persists after 10 timesteps past invasion
invasion_success <- burn_in + 10

#Number of species in the community
species <- 2
a <- species

#Variance Ratio Empty Output
VR <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))
number_species <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))

#2 Species Strong
for (x in 1:length(env_condition)) {
  for (y in 1:length(beta_range)) {
    
    env <- env_condition[x]
    beta <- beta_range[y]
    
    VR_current <- rep(NA, runs)
    total_species <- rep(NA, runs)
    
    for (z in 1:runs) {
      
      K <- c(1000, 1500) #Set to specific values to match the original values used
        #Add strong invader K
        K <- append(K, 1000)
      r <- c(0.5, 0.8)
        #Add invader R
       r <- append(r, 0.7)
      
      #Set starting abundances
      N <- matrix(nrow = time, ncol = species+1)
      colnames(N) <- c(1:(species + 1))
      
      N[1,1:species] <- 50 
      N[1, species+1] <- 0
      
      #Create competition coefficients
      beta_vector <- rnorm(species*species, mean=beta_range[y], sd=0)
      beta_matrix <- matrix(data=beta_vector, nrow = species, ncol = species)
      beta_matrix <- ifelse(beta_matrix<0, 0, beta_matrix)
      
      #Add strong invader competition coefficients
      beta_ri <- rep(0.5, species + 1) #Effect of invader on residents
      beta_ir <- rep(0.5, species) #Effect of residents on invader
      beta_matrix <- beta_matrix %>%
        rbind(beta_ir) %>%
        cbind(beta_ri)
      diag(beta_matrix) <- 1
      
      #Create environmental effect
      sigmaE <- rep(-env, species)
      sigmaE <- append(sigmaE, -0.1)
      
      #Create environmental variation
      miuE <- rnorm(time, mean = 0, sd = 1)
      
      #Create demographic variation
      miuD <- rnorm((species+1)*time, mean = 0, sd = 1)
      miuD <- matrix(data=miuD, nrow = time, ncol = species+1)
      
      for (t in 1:(time-1)) { #for each species being the focal species
        for (s in 1:(species+1)) {
          if(t == burn_in) {
            N[t,species+1] <- 1
          }
          N[t+1,s] <- N[t,s]*exp(r[s]*(1-sum(beta_matrix[s,]*N[t,]/K)) +
                                   (sigmaE[s]*miuE[t]) + (sigmaD*miuD[t,s])/sqrt(N[t,s]))
          
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
      
      if ((species+1)-number_extinct > 1 && N$`3`[invasion_success] > 1){
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
M2resilstrong <- melt(VR) 

#Variance Ratio Empty Output
VR <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))
number_species <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))

#2 Species Weak
for (x in 1:length(env_condition)) {
  for (y in 1:length(beta_range)) {
    
    env <- env_condition[x]
    beta <- beta_range[y]
    
    VR_current <- rep(NA, runs)
    total_species <- rep(NA, runs)
    
    for (z in 1:runs) {
      
      K <- c(1000, 1500) #Set to specific values to match the original values used
      #Add weak invader K
      K <- append(K, 900)
      r <- c(0.5, 0.8)
      #Add invader R
      r <- append(r, 0.4)
      
      #Set starting abundances
      N <- matrix(nrow = time, ncol = species+1)
      colnames(N) <- c(1:(species + 1))
      
      N[1,1:species] <- 50 
      N[1, species+1] <- 0
      
      #Create competition coefficients
      beta_vector <- rnorm(species*species, mean=beta_range[y], sd=0)
      beta_matrix <- matrix(data=beta_vector, nrow = species, ncol = species)
      beta_matrix <- ifelse(beta_matrix<0, 0, beta_matrix)
      
      #Add weak invader competition coefficients
      beta_ri <- rep(0.4, species + 1) #Effect of invader on residents
      beta_ir <- rep(0.6, species) #Effect of residents on invader
      beta_matrix <- beta_matrix %>%
        rbind(beta_ir) %>%
        cbind(beta_ri)
      diag(beta_matrix) <- 1
      
      #Create environmental effect
      sigmaE <- rep(-env, species)
      sigmaE <- append(sigmaE, -0.06)
      
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
                                   (sigmaE[s]*miuE[t]) + (sigmaD*miuD[t,s])/sqrt(N[t,s]))
          
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
      
      if ((species+1)-number_extinct > 1 && N$`3`[invasion_success] > 1){
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
M2resilweak <- melt(VR)

#2 species graphs
plot2resilstrong <- ggplot(M2resilstrong, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#008080", high ="#ca562c", mid = "#f6edbd", midpoint = 1, limit = c(0,2)) +
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), 
       y= expression(paste("Strength of Competititon (", beta, ")")), 
       title = "Post Strong Invasion of 2 Species Communities' VRs", fill="VR") +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95")) +
  theme(axis.title=element_text(size=10), #change axis title size
        axis.text=element_text(size=12), #change axis tick size
        plot.title = element_text(size=14), #change plot title size
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=10)) #change legend text font size)

plot2resilweak <- ggplot(M2resilweak, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#008080", high ="#ca562c", mid = "#f6edbd", midpoint = 1, limit = c(0,2)) +
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), 
       y= expression(paste("Strength of Competititon (", beta, ")")), 
       title = "Post Weak Invasion of 2 Species Communities' VRs", fill="VR") +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95")) +
  theme(axis.title=element_text(size=10), #change axis title size
        axis.text=element_text(size=12), #change axis tick size
        plot.title = element_text(size=14), #change plot title size
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=10)) #change legend text font size)

#Number of species in the community
species <- 5
a <- species

#Variance Ratio Empty Output
VR <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))
number_species <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))

#5 Species Strong
for (x in 1:length(env_condition)) {
  for (y in 1:length(beta_range)) {
    
    env <- env_condition[x]
    beta <- beta_range[y]
    
    VR_current <- rep(NA, runs)
    total_species <- rep(NA, runs)
    
    for (z in 1:runs) {
      
      K <- round(runif(species, min=1000, max=1500))
      #Add strong invader K
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
      
      #Add strong invader competition coefficients
      beta_ri <- rep(0.5, species + 1) #Effect of invader on residents
      beta_ir <- rep(0.5, species) #Effect of residents on invader
      beta_matrix <- beta_matrix %>%
        rbind(beta_ir) %>%
        cbind(beta_ri)
      diag(beta_matrix) <- 1
      
      #Create environmental effect
      sigmaE <- rep(-env, species)
      sigmaE <- append(sigmaE, -0.1)
      
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
                                   (sigmaE[s]*miuE[t]) + (sigmaD*miuD[t,s])/sqrt(N[t,s]))
          
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
      
      if ((species+1)-number_extinct > 1 && N$`6`[invasion_success] > 1){
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
M5resilstrong <- melt(VR) 

#Variance Ratio Empty Output
VR <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))
number_species <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))

#5 Species Weak
for (x in 1:length(env_condition)) {
  for (y in 1:length(beta_range)) {
    
    env <- env_condition[x]
    beta <- beta_range[y]
    
    VR_current <- rep(NA, runs)
    total_species <- rep(NA, runs)
    
    for (z in 1:runs) {
      
      K <- round(runif(species, min=1000, max=1500))
      #Add weak invader K
      K <- append(K, 900)
      r <- round(runif(species, min=0.1, max=0.9), 1)
      #Add invader R
      r <- append(r, 0.4)
      
      #Set starting abundances
      N <- matrix(nrow = time, ncol = species+1)
      colnames(N) <- c(1:(species + 1))
      
      N[1,1:species] <- 50 
      N[1, species+1] <- 0
      
      #Create competition coefficients
      beta_vector <- rnorm(species*species, mean=beta_range[y], sd=0.5)
      beta_matrix <- matrix(data=beta_vector, nrow = species, ncol = species)
      beta_matrix <- ifelse(beta_matrix<0, 0, beta_matrix)
      
      #Add weak invader competition coefficients
      beta_ri <- rep(0.4, species + 1) #Effect of invader on residents
      beta_ir <- rep(0.6, species) #Effect of residents on invader
      beta_matrix <- beta_matrix %>%
        rbind(beta_ir) %>%
        cbind(beta_ri)
      diag(beta_matrix) <- 1
      
      #Create environmental effect
      sigmaE <- rep(-env, species)
      sigmaE <- append(sigmaE, -0.06)
      
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
                                   (sigmaE[s]*miuE[t]) + (sigmaD*miuD[t,s])/sqrt(N[t,s]))
          
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
      
      if ((species+1)-number_extinct > 1 && N$`6`[invasion_success] > 1){
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
M5resilweak <- melt(VR)

#5 species graphs
plot5resilstrong <- ggplot(M5resilstrong, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#008080", high ="#ca562c", mid = "#f6edbd", midpoint = 1, limit = c(0,3)) +
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), 
       y= expression(paste("Strength of Competititon (", beta, ")")), 
       title = "Post Strong Invasion of 5 Species Communities' VRs", fill="VR") +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95")) +
  theme(axis.title=element_text(size=10), #change axis title size
        axis.text=element_text(size=12), #change axis tick size
        plot.title = element_text(size=14), #change plot title size
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=10)) #change legend text font size)

plot5resilweak <- ggplot(M5resilweak, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#008080", high ="#ca562c", mid = "#f6edbd", midpoint = 1, limit = c(0,3)) +
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), 
       y= expression(paste("Strength of Competititon (", beta, ")")), 
       title = "Post Weak Invasion of 5 Species Communities' VRs", fill="VR") +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95")) +
  theme(axis.title=element_text(size=10), #change axis title size
        axis.text=element_text(size=12), #change axis tick size
        plot.title = element_text(size=14), #change plot title size
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=10)) #change legend text font size)

#Number of species in the community
species <- 10
a <- species

#Variance Ratio Empty Output
VR <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))
number_species <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))

#10 Species Strong
for (x in 1:length(env_condition)) {
  for (y in 1:length(beta_range)) {
    
    env <- env_condition[x]
    beta <- beta_range[y]
    
    VR_current <- rep(NA, runs)
    total_species <- rep(NA, runs)
    
    for (z in 1:runs) {
      
      K <- round(runif(species, min=1000, max=1500))
      #Add strong invader K
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
      
      #Add strong invader competition coefficients
      beta_ri <- rep(0.5, species + 1) #Effect of invader on residents
      beta_ir <- rep(0.5, species) #Effect of residents on invader
      beta_matrix <- beta_matrix %>%
        rbind(beta_ir) %>%
        cbind(beta_ri)
      diag(beta_matrix) <- 1
      
      #Create environmental effect
      sigmaE <- rep(-env, species)
      sigmaE <- append(sigmaE, -0.1)
      
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
                                   (sigmaE[s]*miuE[t]) + (sigmaD*miuD[t,s])/sqrt(N[t,s]))
          
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
      
      if ((species+1)-number_extinct > 1 && N$`11`[invasion_success] > 1){
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
M10resilstrong <- melt(VR) 

#Variance Ratio Empty Output
VR <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))
number_species <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))

#10 Species Weak
for (x in 1:length(env_condition)) {
  for (y in 1:length(beta_range)) {
    
    env <- env_condition[x]
    beta <- beta_range[y]
    
    VR_current <- rep(NA, runs)
    total_species <- rep(NA, runs)
    
    for (z in 1:runs) {
      
      K <- round(runif(species, min=1000, max=1500))
      #Add weak invader K
      K <- append(K, 900)
      r <- round(runif(species, min=0.1, max=0.9), 1)
      #Add invader R
      r <- append(r, 0.4)
      
      #Set starting abundances
      N <- matrix(nrow = time, ncol = species+1)
      colnames(N) <- c(1:(species + 1))
      
      N[1,1:species] <- 50 
      N[1, species+1] <- 0
      
      #Create competition coefficients
      beta_vector <- rnorm(species*species, mean=beta_range[y], sd=0.5)
      beta_matrix <- matrix(data=beta_vector, nrow = species, ncol = species)
      beta_matrix <- ifelse(beta_matrix<0, 0, beta_matrix)
      
      #Add weak invader competition coefficients
      beta_ri <- rep(0.4, species + 1) #Effect of invader on residents
      beta_ir <- rep(0.6, species) #Effect of residents on invader
      beta_matrix <- beta_matrix %>%
        rbind(beta_ir) %>%
        cbind(beta_ri)
      diag(beta_matrix) <- 1
      
      #Create environmental effect
      sigmaE <- rep(-env, species)
      sigmaE <- append(sigmaE, -0.06)
      
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
                                   (sigmaE[s]*miuE[t]) + (sigmaD*miuD[t,s])/sqrt(N[t,s]))
          
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
      
      if ((species+1)-number_extinct > 1 && N$`11`[invasion_success] > 1){
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
M10resilweak <- melt(VR)

#10 species graphs
plot10resilstrong <- ggplot(M10resilstrong, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#008080", high ="#ca562c", mid = "#f6edbd", midpoint = 1, limit = c(0,5)) +
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), 
       y= expression(paste("Strength of Competititon (", beta, ")")), 
       title = "Post Strong Invasion of 10 Species Communities' VRs", fill="VR") +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95")) +
  theme(axis.title=element_text(size=10), #change axis title size
        axis.text=element_text(size=12), #change axis tick size
        plot.title = element_text(size=13.5), #change plot title size
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=10)) #change legend text font size)

plot10resilweak <- ggplot(M10resilweak, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low="#008080", high ="#ca562c", mid = "#f6edbd", midpoint = 1, limit = c(0,5)) +
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), 
       y= expression(paste("Strength of Competititon (", beta, ")")), 
       title = "Post Weak Invasion of 10 Species Communities' VRs", fill="VR") +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95")) +
  theme(axis.title=element_text(size=10), #change axis title size
        axis.text=element_text(size=12), #change axis tick size
        plot.title = element_text(size=13.5), #change plot title size
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=10)) #change legend text font size)


#Create 6-panel figure
ggarrange(plot2resilweak, plot2resilstrong, plot5resilweak, plot5resilstrong, plot10resilweak, 
          plot10resilstrong,
          labels = c("A", "B", "C", "D", "E", "F"), ncol = 2, nrow = 3)

#Question 7: Geometric Mean----

env_condition <- seq(from=0, to=.25, by=.01)
beta_range <- seq(from=0, to=.95, by=.05)
time <- 110 
burn_in <- 100
runs <- 2000

#Number of species in the community
species <- 2
a <- species

#Geometric mean Empty Output
avg_geom <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))

#2 Species Strong
for (x in 1:length(env_condition)) {
  for (y in 1:length(beta_range)) {
    
    env <- env_condition[x]
    beta <- beta_range[y]
    
    geom <- matrix(NA, nrow = 1, ncol = runs)
    
    for (z in 1:runs) {
      
      K <- c(1000, 1500) #Set to specific values to match the original values used
      #Add strong invader K
      K <- append(K, 1000)
      r <- c(0.5, 0.8)
      #Add invader R
      r <- append(r, 0.7)
      
      #Set starting abundances
      N <- matrix(nrow = time, ncol = species+1)
      colnames(N) <- c(1:(species + 1))
      
      N[1,1:species] <- 50 
      N[1, species+1] <- 0
      
      #Create competition coefficients
      beta_vector <- rnorm(species*species, mean=beta_range[y], sd=0)
      beta_matrix <- matrix(data=beta_vector, nrow = species, ncol = species)
      beta_matrix <- ifelse(beta_matrix<0, 0, beta_matrix)
      
      #Add strong invader competition coefficients
      beta_ri <- rep(0.5, species + 1) #Effect of invader on residents
      beta_ir <- rep(0.5, species) #Effect of residents on invader
      beta_matrix <- beta_matrix %>%
        rbind(beta_ir) %>%
        cbind(beta_ri)
      diag(beta_matrix) <- 1
      
      #Create environmental effect
      sigmaE <- rep(-env, species)
      sigmaE <- append(sigmaE, -0.1)
      
      #Create environmental variation
      miuE <- rnorm(time, mean = 0, sd = 1)
      
      #Create demographic variation
      miuD <- rnorm((species+1)*time, mean = 0, sd = 1)
      miuD <- matrix(data=miuD, nrow = time, ncol = species+1)
      sigmaD <- c(1, 1, 0)
      
      for (t in 1:(time-1)) { #for each species being the focal species
        for (s in 1:(species+1)) {
          if(t == burn_in) {
            N[t,species+1] <- 1
          }
          N[t+1,s] <- N[t,s]*exp(r[s]*(1-sum(beta_matrix[s,]*N[t,]/K)) +
                                   (sigmaE[s]*miuE[t]) + (sigmaD[s]*miuD[t,s])/sqrt(N[t,s]))
          
          if(is.nan(N[t+1,s])) {
            N[t+1,s] <- 0
          }
          if(N[t+1,s] < 1) {
            N[t+1,s] <- 0
          }
        }
      }
      geom[,z] <- exp(mean(log(N[burn_in:(time-1),3])))
    }
    avg_geom[x,y] <- mean(geom)
  }
  print(x/length(env_condition))
}

#Melt together data
geom2strong <- melt(avg_geom) 

#Variance Ratio Empty Output
avg_geom <- matrix(NA,nrow=length(env_condition), ncol=length(beta_range))

#2 Species Weak
for (x in 1:length(env_condition)) {
  for (y in 1:length(beta_range)) {
    
    env <- env_condition[x]
    beta <- beta_range[y]
    
    geom <- matrix(NA, nrow = 1, ncol = runs)
    
    for (z in 1:runs) {
      
      K <- c(1000, 1500) #Set to specific values to match the original values used
      #Add weak invader K
      K <- append(K, 900)
      r <- c(0.5, 0.8)
      #Add invader R
      r <- append(r, 0.4)
      
      #Set starting abundances
      N <- matrix(nrow = time, ncol = species+1)
      colnames(N) <- c(1:(species + 1))
      
      N[1,1:species] <- 50 
      N[1, species+1] <- 0
      
      #Create competition coefficients
      beta_vector <- rnorm(species*species, mean=beta_range[y], sd=0)
      beta_matrix <- matrix(data=beta_vector, nrow = species, ncol = species)
      beta_matrix <- ifelse(beta_matrix<0, 0, beta_matrix)
      
      #Add weak invader competition coefficients
      beta_ri <- rep(0.4, species + 1) #Effect of invader on residents
      beta_ir <- rep(0.6, species) #Effect of residents on invader
      beta_matrix <- beta_matrix %>%
        rbind(beta_ir) %>%
        cbind(beta_ri)
      diag(beta_matrix) <- 1
      
      #Create environmental effect
      sigmaE <- rep(-env, species)
      sigmaE <- append(sigmaE, -0.06)
      
      #Create environmental variation
      miuE <- rnorm(time, mean = 0, sd = 1)
      
      #Create demographic variation
      miuD <- rnorm((species+1)*time, mean = 0, sd = 1)
      miuD <- matrix(data=miuD, nrow = time, ncol = species+1)
      sigmaD <- c(1, 1, 0)
      
      for (t in 1:(time-1)) { # for each species being the focal species
        for (s in 1:(species+1)) {
          if(t == burn_in) {
            N[t,species+1] <- 1
          }
          N[t+1,s] <- N[t,s]*exp(r[s]*(1-sum(beta_matrix[s,]*N[t,]/K)) +
                                   (sigmaE[s]*miuE[t]) + (sigmaD[s]*miuD[t,s])/sqrt(N[t,s]))
          
          if(is.nan(N[t+1,s])) {
            N[t+1,s] <- 0
          }
          if(N[t+1,s] < 1) {
            N[t+1,s] <- 0
          }
        }
      }
      
      geom[,z] <- exp(mean(log(N[burn_in:(time-1),3])))
    }
    avg_geom[x,y] <- mean(geom)
  }
  print(x/length(env_condition))
}

# Create heat map
geom2weak <- melt(avg_geom) 

max_lim <- max(geom2strong$value, geom2weak$value)

plotgeomweak <- ggplot(geom2weak, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_distiller(palette = "RdBu", limits = c(0, 5)) +
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), 
       y= expression(paste("Strength of Competititon (", beta, ")")), 
       title = "Growth Rate of Weak Invader with 2 Residents", fill="Geometric Mean") +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95")) +
  theme(axis.title=element_text(size=10), #change axis title size
        axis.text=element_text(size=12), #change axis tick size
        plot.title = element_text(size=14), #change plot title size
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=10)) #change legend text font size)

plotgeomstrong <- ggplot(geom2strong, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_distiller(palette = "RdBu", limits = c(0, 5)) +
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), 
       y= expression(paste("Strength of Competititon (", beta, ")")), 
       title = "Growth Rate of Strong Invader with 2 Residents", fill="Geometric Mean") +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95")) +
  theme(axis.title=element_text(size=10), #change axis title size
        axis.text=element_text(size=12), #change axis tick size
        plot.title = element_text(size=14), #change plot title size
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=10)) #change legend text font size)

ggarrange(plotgeomweak, plotgeomstrong,
          common.legend = TRUE, legend = "right", labels = c("A", "B"), ncol = 2, nrow = 1)

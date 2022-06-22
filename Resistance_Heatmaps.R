#####################################
# This script produces a figure that shows the
# growth success of a strong and weak invader
# in the resident communities directly after invasion.
#####################################

# Package for calculating the variance ratio
library(codyn)

# Packages for creating plots
library(ggplot2)
library(reshape2)
library(ggpubr)

# Package for setting paths
library(here)

runs <- 2000 # number of runs

# Load species parameters
source(here("Species_Parameters_Source_Fig2.R"))

###Figure 2: Resistance of resident communities to invader growth of a weak and strong invader----

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
        
        if(is.nan(N3[counter])) {
          N3[counter] <- 0
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
        if(is.nan(N4[counter])) {
          N4[counter] <- 0
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

#Set plot limits
min_lim <- min(success_lambda_good, success_lambda_poor)
max_lim <- max(success_lambda_good, success_lambda_poor)

M2 <- melt(success_lambda_poor)

plot1 <- ggplot(M2, aes(x=Var1, y=Var2, fill=value*100)) + 
  geom_tile() + 
  scale_fill_distiller(palette = "RdBu", limits = c(min_lim*100, max_lim*100)) +
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), y= expression(paste("Strength of Competititon (", beta, ")")), title = "Growth Rates of Weak Invader", fill="Growth
Success (%)") +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95")) +
  theme(axis.title=element_text(size=16), #change axis title size
        axis.text=element_text(size=14), #change axis text size
        plot.title = element_text(size=18), #change plot title size
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.title = element_text(size=16), #change legend title font size
        legend.text = element_text(size=14)) #change legend text font size)

#Successful Invasion of Good Invader Plot
M3 <- melt(success_lambda_good)

plot2 <- ggplot(M3, aes(x=Var1, y=Var2, fill=value*100)) +
  geom_tile() + 
  #geom_contour(aes(z=value*100), stat="contour") +
  scale_fill_distiller(palette = "RdBu", limits = c(min_lim*100, max_lim*100)) +
  labs(x= expression(paste("Effect of Environmental Variability (", sigma[E],")")), y= expression(paste("Strength of Competititon (", beta, ")")), title = "Growth Rates of Strong Invader", fill="Growth
Success (%)") +
  scale_x_continuous(breaks = seq(1, 26, 5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25")) +
  scale_y_continuous(breaks = seq(0, 20, 5), labels = c("0", "0.2", "0.45", "0.7", "0.95")) +
  theme(axis.title=element_text(size=16), #change axis title size
        axis.text=element_text(size=14), #change axis text size
        plot.title = element_text(size=18), #change plot title size
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.title = element_text(size=16), #change legend title font size
        legend.text = element_text(size=14)) #change legend text font size)

#Combine plots together, label A and B
quartz(width=11, height=5)

ggarrange(
  plot1, plot2, labels = c("A", "B"),
  common.legend = TRUE, legend = "right"
)

#Make common axes titles
# Figure.2 <- ggarrange(
#   plot1, plot2, labels = c("A", "B"),
#   common.legend = TRUE, legend = "right"
# )
# annotate_figure(Figure.2,
#                 left = text_grob("Strengh of Competition", rot = 90),
#                 bottom = text_grob("Strength of Environmnental Variability", size = 10)
# )
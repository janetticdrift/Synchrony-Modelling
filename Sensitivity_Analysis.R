#####################################
#This script produces figures that show the range in variance ratios of resident communities,
#and their subsequent synchrony dynamics after invasion of a strong and weak invader.
#####################################
library(here)
library(codyn)
set.seed(13)

#Read in species parameters
source(here("Species_Parameters_Source_Fig134.R"))

###Figure 5: Effects of parameters on invader growth rate

#Create sequences for sigma and beta
beta_min <-  0
beta_max <- 0.95
r_invader_min <- 0
r_invader_max <- 1
k_invader_min <- 500
k_invader_max <- 2000
sigma_e_min <- 0
sigma_e_max <- 0.25
sigma_d_min <- 0
sigma_d_max <- 1

runs <- 2000

#Set outputs of interest
lambda_invader <- rep(NA, runs)         #Growth rate of the invader

#Record parameter set for each run
beta_residents_final <- rep(NA, runs)
beta_invader_final <- rep(NA, runs)
r_invader_final <- rep(NA, runs)
k_invader_final <- rep(NA, runs)
sigma_residents_final <- rep(NA, runs)
sigma_invader_final <- rep(NA, runs)
sigmad_residents_final <- rep(NA, runs)
sigmad_invader_final <- rep(NA, runs)

burn_in <- 100
time <- 200

#Run model#------
    
for (z in 1:runs) {
  env <- runif(1, sigma_e_min, sigma_e_max)
  beta <- runif(1, beta_min, beta_max)
  sigmaE3_range <- runif(1, sigma_e_min, sigma_e_max)
  sigmaD12 <- runif(1, sigma_d_min, sigma_d_max)
  sigmaD3 <- runif(1, sigma_d_min, sigma_d_max)
  beta_invader <- runif(1, beta_min, beta_max)
  r_invader <- runif(1, r_invader_min, r_invader_max)
  k_invader <-runif(1, k_invader_min, k_invader_max)
  
  
  #Set starting conditions
  N1 <- rep(NA, time) 
  N2 <- rep(NA, time)
  N3 <- rep(NA, (time-burn_in))
  
  N1[1] <- 50 
  N2[1] <- 50
  N3[1] <- 1 
  
  invader_abund <- 1
  
  #Set parameters, resident parameters remained constant with Species_Parameters_Source, 
  #while invader parameters were allowed to vary
  r1 <- 0.5
  r2 <- 0.8
  K1 <- 1000
  K2 <- 1500
  beta12 <- beta   
  beta21 <- beta  
  sigmaE1 <- -env 
  sigmaE2 <- -env
  sigmaE3 <- -sigmaE3_range
  r3 <- r_invader       #intrinsic rate of growth of species 4
  K3 <- k_invader     #carrying capacity for species 4
  beta31 <- beta_invader   #effect of species 1 (resident) on species 4
  beta32 <- beta_invader   #effect of species 2 on species 4
  sigmaD1 <- sigmaD12
  sigmaD2 <- sigmaD12
  
  #create environmental and demographic timeseries
  miuE <- rnorm(time, mean = 0, sd = 1) #environmental timeseries variation
  miuD1 <- rnorm(time, mean = 0, sd = 1) #dem timeseries for species 1
  miuD2 <- rnorm(time, mean = 0, sd = 1) #dem timeseries for species 2
  miuD3 <- rnorm((time-burn_in), mean = 0, sd = 1)
  
  #Run the model for resident community
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
  
  counter <- 1
  
  #Model effects of invasion
  for (t in burn_in:(time-1)) {
    N3[counter] <- invader_abund*exp(r3*(1-(invader_abund/K3) - (beta31*N1[t]/K1) - (beta32*N2[t]/K2))
                                          + (sigmaE3*miuE[t]) + (sigmaD3*miuD3[counter])/sqrt(invader_abund))
    if(is.nan(N3[counter])) {
      N3[counter] <- 0
    }
    
    if(N3[counter] < 0) {
      N3[counter] <- 0
    }
    
    counter <- counter + 1
  }
  lambda_per_time <- log(N3/invader_abund)
  lambda_invader[z] <-mean(lambda_per_time) #record results
  
  # record parameters for that run
  beta_residents_final[z] <- beta
  beta_invader_final[z] <- beta_invader
  r_invader_final[z] <- r_invader
  k_invader_final[z] <- k_invader
  sigma_residents_final[z] <- env
  sigma_invader_final[z] <- sigmaE3
  sigmad_residents_final[z] <- sigmaD12
  sigmad_invader_final[z] <- sigmaD3
  
}

#Analyze Results:
beta_residents_final_std <- scale(beta_residents_final)
beta_invader_final_std <- scale(beta_invader_final)
r_invader_final_std <- scale(r_invader_final)
k_invader_final_std <- scale(k_invader_final)
sigma_residents_final_std <- scale(sigma_residents_final)
sigma_invader_final_std <- scale(sigma_invader_final)
sigmad_residents_final_std <- scale(sigmad_residents_final)
sigmad_invader_final_std <- scale(sigmad_invader_final)

#Model of effect sizes
param_effects <- lm(lambda_invader ~ beta_residents_final_std + sigma_residents_final_std + 
                      sigmad_residents_final_std +beta_invader_final_std + r_invader_final_std + 
                      k_invader_final_std +sigma_invader_final_std  + sigmad_invader_final_std)
summary(param_effects)

##########
#Bar Plot#
##########
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

par(mar=c(4.1,3.5,0.7,0.7), tcl=-0.4, mgp=c(1,0.5,0))

#Create sensitivity plot
barplot(param_effects$coefficients[-1], 
        ylim=c(-0.21, 0.2), names.arg=c(expression(paste("Resident ", beta)),
                                         expression(paste("Resident ", sigma[E])),
                                         expression(paste("Resident ", sigma[D])),
                                         expression(paste("Invader ", beta)), 
                                         "Invader R",  
                                         "Invader K",
                                         expression(paste("Invader ", sigma[E])),
                                         expression(paste("Invader ", sigma[D]))), 
        cex.names = .8, las=2, cex.axis = .8, beside=T,
        col=rep(c("grey50"),6))
abline(h=0)
abline(v=3.7, lty = "dashed")
error.bar(seq(from=.7,by=1.2,length.out = 8),param_effects$coefficients[-1],
          summary(param_effects)$coefficients[-1,2],length=0.03)

mtext("Effect size",2,line=2.2)


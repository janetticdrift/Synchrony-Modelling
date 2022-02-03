#####################################
# This source code lists the community parameters 
# of the two resident and two invader species.
#####################################

#######Parameter Abbreviations#######
#r - intrinsic growth rate
#K - carrying capacity
#sigmaD - demographic variance of species 
#sigmaE - Environmental effect on species 
#beta - competitive effect between species
#####################################

#--Resident Community Parameters-----
env_condition <- seq(from=0, to=.25, by=.01)
beta_range <- seq(from=0, to=.95, by=.05)
time <- 300

r1 <- 0.5
r2 <- 0.8
K1 <- 1000      
K2 <- 1500      
sigmaD1 <- 1    
sigmaD2 <- 1    

#--Weak Invader Parameters-----
r3 <- 0.4
K3 <- 900
beta31 <- 0.6  
beta32 <- 0.6  
beta13 <- 0.4  
beta23 <- 0.4  
sigmaE3 <- -0.06 
sigmaD3 <- 1  

#--Strong Invader Parameters-----
r4 <- 0.7    
K4 <- 1000      
beta41 <- 0.5  
beta42 <- 0.5  
beta14 <- 0.5  
beta24 <- 0.5 
sigmaE4 <- -0.1
sigmaD4 <- 1    





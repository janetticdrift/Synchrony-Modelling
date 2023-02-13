###
###Resubmission Analyses
###

#Question 5: Speciose Communities#----

species <- 15
K <- round(runif(species, min=1000, max=1500))
for (x in 1:length(env_condition)) {
  for (y in 1:length(beta_range)) {
    #beta_vector <- runif(species*species, min=0, max=1)
    beta_vector <- rnorm(species*species, mean=beta_range[y], sd=0.05)
    ifelse(beta_vector<0, 0, beta_vector)
    beta_matrix <- matrix(data=beta_vector, nrow = species, ncol = species)
    diag(beta_matrix) <- 1
    for (s in 1:species) { # for each species being the focal species
      N[t+1,s] <- N[t,s]*exp(r[s]*(1-sum(beta[s,]*N[t,]/K)-) +
                               (sigmaE1*miuE[t]) + (sigmaD1*miuD1[t])/sqrt(N[t,s]))
    }
  }
}


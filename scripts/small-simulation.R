##############################################
# This code runs the simulation study for 
# "Informative Bayesian Survival Analysis 
# to Handle Heavy Censoring in Lifetime Data".
# The simulation generates 100 data sets, 
# imposing 3 difference censoring mechanisms 
# on the data for each iteration and fits 3
# different Bayesian Weibull models to each 
# of the 3 censored data sets and stores the 
# posterior of the shape and scale in a list 
# object called posts.
##############################################

set.seed(64467364)

library(dplyr)
library(rstan)
library(stringr)
source("scripts/dataSim_Simple.R")

# simulate data sets

windowLength = 15    # length of observation period
n = 100              # number of idlers

simulate_set <- function(n, windowLength) {
   # sample lifetimes
   samples <- rweibull(n = (n * 99), 
                       shape = 1.15, 
                       scale = 5.253)
   # impose censoring mechanism where t_start = 0
   dfT0 <- dataSim_Simple(t_end = windowLength, 
                          nFrames = n, 
                          samples = samples)
   
   # impose censoring mechanism where t_start = 1 mean life
   dfT5 <- dataSim_Simple(t_start = 5, 
                          t_end = (windowLength + 5), 
                          nFrames = n, 
                          samples = samples)
   
   # impose censoring mechanism where t_start = 4 mean life
   dfT20 <- dataSim_Simple(t_start = 20, 
                           t_end = (windowLength + 20), 
                           nFrames = n, 
                           samples = samples)
   
   return(list(T0 = dfT0, T5 = dfT5, T20 = dfT20))
   
}

# generate 100 * 3 censored data sets and store in list.
sims <- lapply(1:100, function(i) simulate_set(n, windowLength))

# load Stan models
## Objective prior model
m.Improper_R <- stan_model("stan/improper_righ-censored.stan")
m.Improper_Int <- stan_model("stan/improper_interval-censored.stan")

## Independent marginal prior model
m.Independent_R <- stan_model("stan/independent_right-censored.stan")
m.Independent_Int <- stan_model("stan/independent_interval-censored.stan")

## Joint prior model
m.Joint_R <- stan_model("stan/joint_right-cenosored.stan")
m.Joint_Int <- stan_model("stan/joint_interval-cenosored.stan")

# Fit the models
## Functions to prep data for Stan
stanPrepR <-  function(df){
   obs <- df$cens == 0
   rCens <- df$cens == 1
   
   d <- list(N_obs = sum(obs),
             N_cens = sum(rCens),
             lifetime_obs = df$lifetime[obs],
             lifetime_cens = df$lifetime[rCens])
   
   return(d)
}

# function to fit models to generated data sets
stanPrepInt <- function(df, t_start){
   obs <- df$cens == 0
   rCens <- df$endCens == 1
   intCens <- (df$startCens == 1) & (df$endCens != 1)
   
   d <- list(N_obs = sum(obs),
             N_Rcens = sum(rCens),
             N_Icens = sum(intCens),
             lifetime_obs = df$lifetime[obs],
             lifetime_Rcens = df$lifetime[rCens],
             lifetime_Icens_Upper = (df$lifetime[intCens] + t_start),
             lifetime_Icens_Lower = df$lifetime[intCens])
   
   return(d)
}

## Fit each model to each data set
posts <- lapply(1:100, 
                function(i) 
                   lapply(c("T0", "T5", "T20"), 
                          function(d) {
                             modelList <- c("m.Improper",
                                            "m.Independent",
                                            "m.Joint")
               
                             if (d == "T0") {
                                modelList <- str_c(modelList, "_R")
                                data <-stanPrepR(sims[[i]][[d]])
                             } else if (d == "T5") {
                                modelList <- str_c(modelList, "_Int")
                                data <-stanPrepInt(sims[[i]][[d]], 5)
                             } else {
                                modelList <- str_c(modelList, "_Int")
                                data <-stanPrepInt(sims[[i]][[d]], 20)
                             }
                             
                             lapply(modelList, function(model) {
                                
                                if(str_detect(pattern = "m.Joint", 
                                              string = model)){
                                   data$t1 <- 3.82
                                   data$t2 <- 15
                                }
                                
                                model <- sampling(get(model), 
                                                  data = data, 
                                                  chains = 4, 
                                                  cores = 4)
                                list(shape = extract(model)$beta, 
                                     scale = extract(model)$eta)
                                
                             })
                                
                          })
                )

# save list object that contains the posteriors of eta (scale)
# and beta (shape) parameters.
save(posts, file = "results/simulation-posteriors.RData")


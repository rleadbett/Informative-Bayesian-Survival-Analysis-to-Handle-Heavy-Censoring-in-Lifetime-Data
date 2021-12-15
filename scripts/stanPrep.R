#' A function to automate the maniputaltion of a dataset to
#' feed into stan.
#' 
#' @param df A dataframe with one column `lifetime` containing 
#' the lifetimes and another column `cens` containing the censoring
#' indicator for each of the lifetimes.
#' 
#' @return A list element which is ready to feed into rstan

stanPrep <- function(df){
  
  censored <- df$cens
  y <- df$Time
  
  # Stan Bayesian estimation
  N_obs <- sum(censored == 0)
  N_cens <- sum(censored == 1)
  lifetime_obs <- y[censored == 0]
  lifetime_cens <- y[censored == 1]
  
  inputData <- list(N_obs = N_obs,
                    N_cens = N_cens,
                    lifetime_obs = lifetime_obs,
                    lifetime_cens = lifetime_cens)
  
  return(inputData)
}
// Independant-informative-prior 
//
// Stan model for censoring paper.
//
// This is a Bayesian model for right-censored lifetime data 
// which followes a single two-paramater Weibull distribution 
// with informative indapendent marginal priors placed on 
// both the shape and scale paramaters.
//
// The informative priors are constructed from the domain 
// knowledge that
// 1) the manufacturer claims that the mean lifetime is 5 
//    years with sd of 2 years.
// 2) the beta value of steel rolling element bearings is 
//    between 1.1 and 1.5
//

data {
	int N_obs;
	int N_cens;
	real lifetime_obs[N_obs];
	real lifetime_cens[N_cens];
}

parameters {
	real<lower = 0> beta;
	real<lower = 0> lambda;
}

transformed parameters {
  real log_lambda;
  real<lower = 0> eta;
  
  log_lambda = log(lambda);
  eta = exp((-1 / beta) * log_lambda);
}

model{

// Likelihood
// non-censored portion
for(i in 1:N_obs){
	target += weibull_lpdf(lifetime_obs[i]|beta, eta);
}
// censored portion
for(j in 1:N_cens){
	target += weibull_lccdf(lifetime_cens[j]|beta, eta);
}
  
// Prior models
// These two priors place 95% or the probability
// mass between lambda = (, ) and beta = (, )
lambda ~ gamma(6.25, 42.59);
beta ~ gamma(156.25, 120.19);

}


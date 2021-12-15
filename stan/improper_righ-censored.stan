// Independant-uninformative-prior 
//
// Stan model for censoring paper.
//
// This is a Bayesian model for right-censored lifetime data 
// which followes a single two-paramater Weibull distribution 
// with an uninformative indapendent marginal priors placed on 
// both the shape and scale paramaters.
//
// A gamma(0.001, 0.001) prior is placed on both of the 
// paramaters. This prior is dominated by the likelihood as
// long as n >= 1.
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
lambda ~ gamma(0.0001, 0.0001);
beta ~ gamma(0.0001, 0.0001);
}




# Informative Bayesian Survival Analysis to Handle Heavy Censoring in Lifetime Data
This repository contains the `R` code required to reproduce the small scale simulation study from the paper with the same name [leadbetter et al. 2021](). This includes all the Stan files which implement the HMC sampling routines for an objective, and two different informative prior Weibull models.

## Reproducing work
### Run order
To reproduce the results of the simulation study one just needs to run the `RunAll.sh` file in the main directory. Otherwise, if this does not work then the file `small-simulation.R` can be run from within the `R` project by running the code `source('scripts/small-simulation.R', echo = TRUE)` in the console.

### Reproducability 
#### RStan
This code requires the use of `Stan` which is a probabalistic programing language writen in `C++`. `R` comunicates to `Stan` through the package `RStan`. Therefor, to run the code in this repository the `RStan` must be installed. Please follow the instruction provided by the stan-dev team [here](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).
#### renv
For reproducability this repository uses `renv` which is automaticaly activated when you open the main `R` project. For more information about `renv` please see the [vignette](https://rstudio.github.io/renv/articles/collaborating.html).

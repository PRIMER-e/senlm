## --- R code to test senlm2 package

## --- Library
library("devtools")
library ("senlm")

## --- Load all
load_all()

## --- Set models
Models <- set_models (mean_class="all", err_class="all", method="crossed", binomial_n=50)

## --- Set default parameters list
Pars <- create_default_par_list (Models)

## --- Simulate data
Data <- create_simulated_datasets (Pars=Pars, x=runif(500,0,100), seed=12345)

## --- Fit multiple models
Fits   <- senlm.multi (Models, Data)

## --- Save workspace
save.image (file="Run.RData")
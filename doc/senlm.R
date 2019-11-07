## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  ##comment = "#>"  
  comment = ""  
)

## ----set-options, echo=FALSE, cache=FALSE---------------------------------------------------------
options(width = 100)

## ----installation, eval=FALSE---------------------------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("PRIMER-e/senlm")

## ----setup----------------------------------------------------------------------------------------
library(senlm)

## ---- echo=FALSE----------------------------------------------------------------------------------
knitr::kable(mean_functions())

## -------------------------------------------------------------------------------------------------
get_mean_fun_parnames("constant")

## -------------------------------------------------------------------------------------------------
get_mean_fun_parnames("uniform")

## -------------------------------------------------------------------------------------------------
get_mean_fun_parnames("beta")

## -------------------------------------------------------------------------------------------------
get_mean_fun_parnames("sech")

## -------------------------------------------------------------------------------------------------
get_mean_fun_parnames("sech_p1")

## -------------------------------------------------------------------------------------------------
get_mean_fun_parnames("sech_r0p1")

## -------------------------------------------------------------------------------------------------
get_mean_fun_parnames("gaussian")

## -------------------------------------------------------------------------------------------------
get_mean_fun_parnames("mixgaussian")

## -------------------------------------------------------------------------------------------------
get_mean_fun_parnames("mixgaussian_equal")

## -------------------------------------------------------------------------------------------------
get_mean_fun_parnames("hofV")

## -------------------------------------------------------------------------------------------------
get_mean_fun_parnames("hofII")

## -------------------------------------------------------------------------------------------------
get_mean_fun_parnames("hofIV")

## -------------------------------------------------------------------------------------------------
get_mean_fun_parnames("hofIVb")

## -------------------------------------------------------------------------------------------------
get_mean_fun_parnames("hofVb")

## ---- echo=FALSE----------------------------------------------------------------------------------
knitr::kable(error_distributions())

## -------------------------------------------------------------------------------------------------
get_constant_parnames("binomial_count")

## -------------------------------------------------------------------------------------------------
get_constant_parnames("binomial_prop")

## -------------------------------------------------------------------------------------------------
get_err_dist_parnames("negbin")

## -------------------------------------------------------------------------------------------------
get_err_dist_parnames("zip")

## -------------------------------------------------------------------------------------------------
get_err_dist_parnames("zinb")

## -------------------------------------------------------------------------------------------------
get_err_dist_parnames("zipl")

## -------------------------------------------------------------------------------------------------
get_err_dist_parnames("zinbl")

## -------------------------------------------------------------------------------------------------
get_err_dist_parnames("zipl.mu")

## -------------------------------------------------------------------------------------------------
get_err_dist_parnames("zinbl.mu")

## -------------------------------------------------------------------------------------------------
get_err_dist_parnames("gaussian")

## -------------------------------------------------------------------------------------------------
get_err_dist_parnames("tweedie")

## -------------------------------------------------------------------------------------------------
get_err_dist_parnames("zig")

## -------------------------------------------------------------------------------------------------
get_err_dist_parnames("zigl")

## -------------------------------------------------------------------------------------------------
get_err_dist_parnames("zigl.mu")

## -------------------------------------------------------------------------------------------------
get_err_dist_parnames("tab")
get_constant_parnames("tab")

## -------------------------------------------------------------------------------------------------
get_err_dist_parnames("zitab")
get_constant_parnames("zitab")

## -------------------------------------------------------------------------------------------------
## All possible combinations of mean functions and error distrbiutions
set_models (mean_fun=c("gaussian", "beta"), err_dist=c("zinb", "zinbl"),
            method="crossed")
## Fit specific models
set_models (mean_fun=c("gaussian", "beta"), err_dist=c("zinb", "zinbl"),
            method="paired")
## Fit binomial model - with n=40
set_models (mean_fun=c("gaussian"), 
            err_dist=c("poisson", "binomial_prop", "binomial_count"),
            binomial_n=40, method="paired")
## Can also specify via mean and error class
set_models (mean_class="main", err_class="percentage")

## ----eval=FALSE-----------------------------------------------------------------------------------
#  ## --- Set models
#  Models <- set_models (mean_fun=c("gaussian"), err_dist=c("zip", "zinb"))
#  
#  ## --- Set default parameters list
#  Pars <- create_default_par_list (Models)
#  
#  ## --- Simulate data
#  Data <- create_simulated_datasets (Pars, N=200, xmin=0, xmax=100, seed=12345)

## -------------------------------------------------------------------------------------------------
data(haul, package="senlm")
head (haul)

## -------------------------------------------------------------------------------------------------
Fit <- senlm (data=haul, xvar="x_depth", yvar="y_Sebastolobus.altivelis", 
              mean_fun="gaussian", err_dist="zip")
summary (Fit)

## ----eval=FALSE-----------------------------------------------------------------------------------
#  Models <- set_models (mean_class="test", err_class="binary")
#  Pars   <- create_default_par_list (models=Models)
#  Data   <- create_simulated_datasets (Pars, seed=12345)
#  Fits   <- senlm.multi (models=Models, data=Data, xvar="x_depth", yvar="y_Sebastolobus.altivelis")
#  print (Models)
#  print (Fits)

## ----eval=FALSE-----------------------------------------------------------------------------------
#  Models <- set_models (mean_fun="gaussian", err_dist=c("zip", "zinb"))
#  Data <- vector (length=2, mode="list")
#  Data[[1]]$x <- haul$x_depth
#  Data[[1]]$y <- haul$y_Albatrossia.pectoralis
#  Data[[2]]$x <- haul$x_depth
#  Data[[2]]$y <- haul$y_Sebastolobus.altivelis
#  Fits   <- senlm.multi (Models, Data)


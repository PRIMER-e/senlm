## --- Model fitting : Maximum likelihood estimation


## ---- ---- ----

#' Fit species-environment non-linear model via maximum likelihood
#' 
#' 'senlm' fits a species-environment non-linear model via maximum likelihood.
#'
#' @param model Model to fit.
#' @param data Data frame of x and y.
#' @param xvar Name of x (domain) variable (or column number).
#' @param yvar Name of y (response) variable (or column number).
#' @param mean_fun Mean function (if model not supplied).
#' @param err_dist Error distribution (if model not supplied).
#' @param binomial_n The binomial n parameter, if error distribution
#' is "binomial_count" or "binomial_percent" (if model not supplied).
#' @param x Vector of x (domain) values (if data not supplied).
#' @param y Repsonse variable vector (if data not supplied).
#' 
#' @return Object containg model fit to data y and x.
#' 
#' @keywords fit senlm model, mle
#'
#' @examples
#'
#' \dontrun{
#' ## Simulate data
#' Pars <- create_default_par_list (mean_fun="gaussian", err_dist="poisson")
#' Data <- create_simulated_datasets (Pars, N=100, xmin=0, xmax=100, seed=12345)
#' ## Fit model
#' Fit <- senlm (mean_fun="gaussian", err_dist="poisson", x=Data[[1]]$x, y=Data[[1]]$y)
#'
#' ## Real data
#' Model <- set_models (mean_fun="gaussian", err_dist=c("zip"))
#' Fit <- senlm (model=Model, data=haul, xvar="x_depth", yvar="y_Sebastolobus.altivelis")
#' }
#' 
#' @export
#' 
senlm <- function (model=NULL, data=NULL, xvar=NULL, yvar=NULL, 
                   mean_fun=NULL, err_dist=NULL, binomial_n=NULL, x=NULL, y=NULL) {
  ## --- Fit model using maximum likelihood

  ## --- Model not specified
  if (is.null(model)) {
    ## Check if mean function and error distribution supplied
    if (is.null(mean_fun) | is.null(err_dist)) {
      stop ("Must specify either model argument or mean function and error distribution!")
    }
    ## Check length of mean function and error distribution
    if ( (length(mean_fun)>1) | (length(err_dist)>1) | (length(binomial_n)>1) ) {
      stop ("Length of mean function, error distribution, and binomial_n must be 1!")
    }
    ## --- Set model
    model <- set_models(mean_fun, err_dist, binomial_n)
  } else {
    ## Check if only one model given
    if (nrow(model)!=1) { stop ("Only one model should be specified!") }
  }
  
  ## --- Check data
  if (is.null(data)) {
    ## --- Data not specified
    ## Check x and y specified
    if (is.null(x) | is.null(y)) { stop ("Must specify x and y if data not given!") }
    xname <- "x"
    yname <- "y"
  } else {
    ## --- Data specified
    ## Check x and y specified
    if (is.null(xvar) | is.null(yvar)) { stop ("Must specify xvar and yvar if data given!") }
    
    if (!is.character(xvar)) { stop ("xvar must be character!") }
    if (length(xvar)>1) { stop ("Too many x variable names specified!") }
    if (any(is.na(match (xvar, names(data))))) { stop ("xvar not valid!") }
    
    if (!is.character(yvar)) { stop ("yvar must be character!") }
    if (length(yvar)>1) { stop ("Too many y variable names specified!") }
    if (any(is.na(match (yvar, names(data))))) { stop ("Some yvar not valid!") }
    
    ## Explanatory variable
    x <- data[,xvar]
    xname <- xvar
    
    ## Response variables
    y <- data[,yvar]
    yname <- yvar
  }
  
  ## --- Create data set
  Dat <- data.frame (x=x, y=y)
  
  ## --- Create model info object
  ModelInfo <- set_model_info (model=model, data=Dat)
  
  ## --- Set negative log-likelihood function
  NLL <- set_nll (ModelInfo=ModelInfo)
  
  ## --- Set mle fit function
  if (ModelInfo$model == "constant-poisson") {
    ## Poisson with constant mean function
    estimate_mle <- mle_constant_poisson
  } else {
    ## Default method : SANN/Nelder-Mead
    estimate_mle <- mle_default
  }
  
  ## --- Fit mle
  Fit.theta <- try( estimate_mle (ModelInfo, Dat) )
  
  ## --- Store model fit and fail flag
  if (class(Fit.theta) == "try-error") { Fail <- TRUE } else { Fail <- FALSE }
  
  
  ## --- Store fit
  Fit <- list()
  
  ## --- Model name
  Fit$model <- ModelInfo$model
  
  ## --- Model info
  Fit$model_info <- ModelInfo
  
  ## --- Data
  Fit$y <- y
  Fit$x <- x
  
  ## --- Names
  Fit$xname <- xname
  Fit$yname <- yname
  
  ## --- Was fit successful?
  if (Fail == FALSE) {
    
    ## --- Fit successful
    Fit$Fail <- FALSE
    
    ## --- Fitted parameters
    Fit$theta <- Fit.theta
    
    ## --- Goodness of fit
    Fit$IC <- fit_information_criteria (ModelInfo, Dat, Fit.theta)

    ## --- Fitted values
    Fit$fitted <- senlm::mu_meanfunction (ModelInfo=Fit$model_info, theta=Fit$theta, x=Fit$x)
    ## --- Residuals
    Fit$residuals <- Fit$y - Fit$fitted    
    
  } else {
    
    ## --- Fit unsuccessful
    Fit$Fail <- TRUE
  }
  
  ## --- Set class
  class(Fit) <- "senlm"
  
  ## --- Return fit
  return (Fit)
}


#' Predict from senlm model fit
#' 
#' Predict from senlm model fit using x values in newdata, or actual x values.
#'
#' @param object Model fit.
#' @param newdata x values to predict y values from. Optional, data values used by default.
#' @param ... additional optional arguments.
#' 
#' @return Vector of predicted y values.
#' 
#' @keywords predict senlm model fit
#'
#' @examples
#'
#' \dontrun{
#'
#' ## Real data
#' Model <- set_models (mean_fun="gaussian", err_dist=c("zip"))
#' Fit <- senlm (model=Model, data=haul, xvar="x_depth", yvar="y_Sebastolobus.altivelis")
#' predict (Fit)
#' }
#' 
#' @export
#' 
predict.senlm <- function (object, newdata, ...) {
  ## --- Predict senlm model
  
  ## --- Check if object is a senlm object
  if (class(object) != "senlm") {
    stop ("object not a semlm object!")
  }
  
  ## --- Set x if missing to x variable from data
  if (missing(newdata) || is.null(newdata)) {
    newdata <- object$x
  }

  ## --- Calculate predication
  pred <- senlm::mu_meanfunction (ModelInfo=object$model_info, theta=object$theta, x=newdata)
  
  ## --- Return predication
  return (pred)
}


#' Simulate data from senlm model fit
#' 
#' Simulate from senlm model fit, nsim times, using x values in newdata, or actual x values.
#'
#' @param object Model fit.
#' @param nsim Number of simulations to perform. Default value of one.
#' @param seed Random number seed.
#' @param newdata x values to predict y values from. Optional, data values used by default.
#' @param ... additional optional arguments.
#' 
#' @return Data frame of predicted y values.
#' 
#' @keywords simulate senlm model fit
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#'
#' ## Simulate data
#' Model <- set_models (mean_fun="gaussian", err_dist=c("zip"))
#' Fit <- senlm (model=Model, data=haul, xvar="x_depth", yvar="y_Sebastolobus.altivelis")
#' sim(Fit, nsim=2, newdata=c(400,600,800))
#' }
#' 
#' 
sim <- function (object, nsim=1, seed=NULL, newdata=NULL, ...) {
  ## --- Simulate data nsim times from senlm model
  
  ## --- Check if object is a senlm object
  if (class(object) != "senlm") {
    stop ("object not a semlm object!")
  }

  ## --- Set x if missing to x variable from data
  if (missing(newdata) || is.null(newdata)) {
    newdata <- object$x
  }
  
  ## --- Create parameter Fit
  Par <- list ()
  ## Model names
  Par$model_name <- object$model
  Par$mean_fun   <- object$model_info$mean_fun
  Par$err_dist   <- object$model_info$err_dist
  ## Split paramters by type
  if (!is.null(object$model_info$thetaM)) {
    Par$thetaM <- object$theta[object$model_info$thetaM]
  } else {
    Par$thetaM <- NULL
  }
  if (!is.null(object$model_info$thetaE)) {
    Par$thetaE <- object$theta[object$model_info$thetaE]
  } else {
    Par$thetaE <- NULL
  }
  if (!is.null(object$model_info$thetaC)) {
    Par$thetaC <- object$theta[object$model_info$thetaC]
  } else {
    Par$thetaC <- NULL
  }

  ## --- Create empty data object
  Dat <- as.data.frame(matrix(NA, ncol=nsim, nrow=length(newdata)))
  names(Dat) <- paste("sim_", 1:nsim, sep="")
  
  ## --- Simulate data
  for (i in 1:nsim) {
    SimData <- simulate_data (x=newdata, Par=Par, seed=seed)
    Dat[,i] <- SimData$y
  }
  
  ## --- Return data
  return (Dat)
}

#' Print summary of senlm model fit
#' 
#' Print summary information for senlm model fit.
#'
#' @param object Model fit.
#' @param ... additional optional arguments.
#' 
#' @return Summary of senlm model fit.
#' 
#' @keywords summary senlm model fit
#'
#' @examples
#'
#' \dontrun{
#' 
#' ## Summarise data
#' Model <- set_models (mean_fun="gaussian", err_dist=c("zip"))
#' Fit <- senlm (model=Model, data=haul, xvar="x_depth", yvar="y_Sebastolobus.altivelis")
#' summary(Fit)
#' }
#' 
#' @export
#' 
summary.senlm <- function (object, ...) {
  ## --- Display summary of fit of senlm model

  ## --- Create summary object
  Summary       <- list()
  Summary$model <- object$model
  Summary$theta <- object$theta
  Summary$IC    <- object$IC
  
  ## --- Print summary
  print (Summary)
}

#' Plot senlm model fit over data
#' 
#' Produces a scatterplot of y vs x with fitted mean curve overlaid.
#'
#' @param x Model fit.
#' @param ... additional optional arguments.
#' 
#' @return Plot of senlm model fit over data
#' 
#' @keywords plot senlm model fit
#'
#' @examples
#'
#' \dontrun{
#' 
#' ## Plot fitted model
#' Model <- set_models (mean_fun="gaussian", err_dist=c("zip"))
#' Fit <- senlm (model=Model, data=haul, xvar="x_depth", yvar="y_Sebastolobus.altivelis")
#' plot(Fit)
#' }
#' 
#' @export
#' 
plot.senlm <- function (x, ...) {
  ## --- Plot fit of senlm model

  ## --- Grab model fit
  fit <- x
  
  ## --- Sort x variable and predict
  x <- sort(fit$x)
  p <- predict.senlm(fit, newdata=sort(fit$x))
  
  ## --- Plot data and add fitted mean function
  graphics::plot(fit$y ~ fit$x, xlab=fit$xname, ylab=fit$yname, main=fit$model)
  graphics::points(p ~ x, type="l", lwd=2, col="red4")
}


fit_information_criteria <- function (ModelInfo, Dat, fit.theta) {
  ## --- Calculate information criteria from model fit
  
  ## --- Set negative log-likelihood function
  NLL <- set_nll (ModelInfo=ModelInfo)
  
  ## --- Transform parameter to unbounded
  u.theta <- make_unbounded (ModelInfo, fit.theta)
  
  ## --- Value of negative log-likelikehood
  fit.value <- NLL (u.theta, ModelInfo=ModelInfo, Dat=Dat)
  names(fit.value) <- "nll"
  
  ## --- Negative log-likelihood
  nll <- fit.value
 
  ## --- Sample size
  n <- length (Dat$y)

  ## --- Number of parameters
  npar <- length(u.theta)
  
  ## --- AIC
  AIC <- 2*nll + 2*npar
  
  ## --- AICc
  AICc <- AIC + (2*npar^2 + 2*npar)/(n - npar - 1)
  
  ## --- BIC
  BIC <- 2*nll + log(n)*npar

  ## --- Store results
  IC <- c(npar=npar, nll=nll, AIC=AIC, AICc=AICc, BIC=BIC)
  names(IC) <- c("npar", "nll", "AIC", "AICc", "BIC")
  
  ## --- Return results
  return (IC)
}


## ---- ---- ----


mle_default <- function (ModelInfo, Dat, theta0=NULL) {
  ## --- Default maximum likelihood fitting method

  ## --- Estimate paramters using MOM from splines
  if (is.null(theta0)) {
    theta0 <- init_mle (ModelInfo, Dat)
  }
  ## --- Transform estimated parameters to unbounded
  u.theta0 <- make_unbounded (ModelInfo, theta0)
  
  ## --- Set negative log-likelihood (with parameter names)
  NLL <- nll_wrapper (ModelInfo=ModelInfo, Dat=Dat)
  bbmle::parnames (NLL) <- names(u.theta0)
  
  ## --- SANN : Find mle using simulated annealing
  Fit.sann <- try ( suppressWarnings (
    bbmle::mle2 (minuslogl=NLL, optimizer="optim", method="SANN", 
                 vecpar=TRUE, start=u.theta0) ))
  
  ## --- Test if sann model fit
  if (class(Fit.sann) != "try-error") {
    
    ## --- Fit using nlminb from sann estimate
    Fit.nlsann <- try ( suppressWarnings (
      bbmle::mle2 (minuslogl=NLL, optimizer="nlminb",
                   vecpar=TRUE, start=bbmle::coef(Fit.sann)) ))
    
    ## --- Test if nlminb from sann fits
    if (class(Fit.nlsann) != "try-error") {
      ## --- nlminb based on sann succeeded
      FitType   <- "nl-sann"
      Fit.theta <- make_bounded (ModelInfo, bbmle::coef(Fit.nlsann))
    } else {
      ## --- nlminb failed - use sann
      FitType   <- "sann"
      Fit.theta <- make_bounded (ModelInfo, bbmle::coef(Fit.sann))
    }
    
  } else {
    
    ## --- Fit using nlminb from init estimate
    Fit.nl <- try ( suppressWarnings (
      bbmle::mle2 (minuslogl=NLL, optimizer="nlminb",
                   vecpar=TRUE, start=u.theta0) ))
    
    ## --- Test if nlminb from init fits
    if (class(Fit.nl) != "try-error") {
      ## --- nlminb  - use init
      FitType   <- "nl"
      Fit.theta <- make_bounded (ModelInfo, bbmle::coef(Fit.nl))
    } else {
      ## --- Fail
      FitType <- "fail"
      Fit.theta <- NA * theta0
    }
  }
  
  ## Return fit
  return (Fit.theta)
}

mle_constant_bernoulli <- function (ModelInfo, Dat) {
  ## --- MLE for constant-bernoulli mean function

  ## MLE
  Fit.theta        <- mean (Dat$y)
  names(Fit.theta) <- c("H")
  
  ## Return MLE fit
  return (Fit.theta)
}

mle_constant_poisson <- function (ModelInfo, Dat) {
  ## --- MLE for constant-poisson mean function

  ## Fit constant mean function
  Fit.theta <- mle_constant_bernoulli (ModelInfo, Dat)

  ## Return MLE fit
  return (Fit.theta)
}

mle_uniform_bernoulli <- function (ModelInfo, Dat) {
  ## --- MLE for uniform-bernoulli mean function

  ## --- Grab count data
  y <- Dat$y
  x <- Dat$x
  
  ## Sort y & x by x
  Y <- y[order(x)]
  X <- x[order(x)]
  ## Remove data where y is 0
  xlim <- range(X[Y>0])
  ## Find bounds of uniform
  xminpos <- which(X==xlim[1])
  xmaxpos <- which(X==xlim[2])
  if (xminpos > 1)         { xminpos <- c(xminpos-1, xminpos)   }
  if (xmaxpos < length(x)) { xmaxpos <- c(xmaxpos,   xmaxpos+1) }
  
  ## --- MLE
  
  ## Mean parameters
  c <- mean(X[xminpos])
  d <- mean(X[xmaxpos])
  H <- mean(Y[(X>=c) & (X<=d)])
  
  ## Store parameters
  Fit.theta        <- c(H, c, d)
  names(Fit.theta) <- c("H", "c", "d")

  ## Return MLE fit
  return (Fit.theta)
}


## ---- ---- ----

  
#' Fit multiple species-environment non-linear models via maximum likelihood
#' 
#' 'senlm.multi' fits species-environment non-linear models via maximum likelihood.
#'
#' @param models Object listing models to fit (from set_models function).
#' @param data A data frame containing 'x' (explanatory) and 'y' (response) variables.
#' @param xvar Name of explanatory variable (must be univariate).
#' @param yvar Names of response variables.
#' 
#' @return Object containg model fits to data y and x.
#' 
#' @keywords fit senlm model, mle
#'
#' @examples
#'
#' \dontrun{
#'
#' models <- set_models (mean_class="test", err_dist=c("zip","zinb"))
#' Fits <- senlm.multi (models=models, data=haul, xvar="x_depth",
#'                      yvar=c("y_Albatrossia.pectoralis", "y_Sebastolobus.altivelis"))
#' }
#' @export
#' 
senlm.multi <- function (models=NULL, data=NULL, xvar=NULL, yvar=NULL) {
  ## --- Fit models using maximum likelihood
  
  ## --- Check inputs


  ## Check x and y specified
  if (is.null(xvar) | is.null(yvar)) { stop ("Must specify xvar and yvar!") }

  if (!is.character(xvar)) { stop ("xvar must be character!") }
  if (length(xvar)>1) { stop ("Too many x variable names specified!") }
  if (any(is.na(match (xvar, names(data))))) { stop ("xvar not valid!") }
  
  if (!is.character(yvar)) { stop ("yvar must be character!") }
  if (any(is.na(match (yvar, names(data))))) { stop ("Some yvar not valid!") }
  
  ## Explanatory variable
  x <- data[,xvar]
  xname <- xvar
  
  ## Response variables
  y <- data[,yvar]
  yname <- yvar
  
  ## --- Fit models

  ## --- Create fit object
  Fits <- vector (mode="list", length=length(yvar))
  names(Fits) <- yvar

  ## --- Loop through response variables
  for (i in 1:length(Fits)) {

    ## Display iteration
    print (i)
    
    ## Create object to store model fits to data with ith response variable
    ModelFits <-  vector (mode="list", length=nrow(models))
    ModelNames <- rep (NA, length=nrow(models))

    ## --- Loop throough models
    for (j in 1:length(ModelFits)) {

      ## --- Fit model
      Fit <- senlm (model=models[j,], data=data, xvar=xvar, yvar=yvar[i])

      ## --- Store model
      ModelFits[[j]] <- Fit
      ModelNames[j] <- Fit$model
    }

    ## --- Add model names
    names(ModelFits) <- ModelNames

    ## --- Store model fits
    Fits[[i]] <- ModelFits
  }

  ## --- Set class
  class (Fits) <- "senlm.multi"

  ## --- Return fits
  return (Fits)
}


summary.senlm.multi <- function (Fits) {
  ## --- Print summary of model fits 

  ## Loop through response variables
  for (i in 1:length(Fits)) {
    ## Loop through models
    for (j in 1:length(Fits[[i]])) {
      ## Print summary
      summary(Fits[[i]][[j]])
    }
  }
}


plot.senlm.multi <- function (Fits) {
  ## --- Print summary of model fits 

  ## Loop through response variables
  for (i in 1:length(Fits)) {
    ## Loop through models
    for (j in 1:length(Fits[[i]])) {
      ## Plot model fit
      plot.senlm (Fits[[i]][[j]])
    }
  }
}
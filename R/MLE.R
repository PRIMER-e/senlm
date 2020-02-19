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
#' is "binomial.count" or "binomial.percent" (if model not supplied).
#' @param y Repsonse variable vector (if data not supplied).
#' @param x Vector of x (domain) values (if data not supplied).
#'
#' @return Object containg model fit to data y and x.
#'
#' @keywords fit senlm model, mle
#'
#' @examples
#'
#' \dontrun{
#' ## Simulate data
#' 
#' dat  <- create_simulated_datasets(pars, N=100, xmin=0, xmax=100, seed=12345)
#' ## Fit model
#' fit <- senlm(mean_fun="gaussian", err_dist="poisson", data=dat,
#'              xvar="x", yvar="gaussian_poisson")
#'
#' ## Real data
#' model <- set_models(mean_fun="gaussian", err_dist="zip")
#' fit <- senlm(model=model, data=haul, xvar="depth", yvar="Sebastolobus.altivelis")
#' }
#'
#' @export
#'
senlm <- function (model=NULL, data=NULL, xvar=NULL, yvar=NULL,
                   mean_fun=NULL, err_dist=NULL, binomial_n=NULL, y=NULL, x=NULL) {
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

  ## --- Store fit
  Fit <- list()

  ## --- Model name
  Fit$model <- gsub ("-", "_", ModelInfo$model) ## *** This line should not be needed!
  
  ## --- Model info
  Fit$model_info <- ModelInfo

  ## --- Data
  Fit$y <- y
  Fit$x <- x

  ## --- Names
  Fit$yname <- yname
  Fit$xname <- xname
   
  ## --- Fit mle
  FitMLE <- try( mle_default (ModelInfo, Dat) )

  ## --- Store fitted values  
  Fit$convergence <- FitMLE$convergence
  Fit$theta <- FitMLE$theta
  Fit$lb    <- FitMLE$lb
  Fit$ub    <- FitMLE$ub
  Fit$u.hessian <- FitMLE$u.hessian
  
  ## --- Was fit successful?
  if (Fit$convergence  == 0) {
    
    ## --- Goodness of fit
    Fit$IC <- fit_information_criteria (ModelInfo, Dat, Fit$theta)
     
    ## --- Fitted values
    Fit$fitted <- mu_meanfunction (ModelInfo=Fit$model_info, theta=Fit$theta, x=Fit$x)
    ## --- Residuals
    Fit$residuals <- Fit$y - Fit$fitted
    ## --- Quantile residuals
    Fit$qresiduals <- qres (Fit)

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
#' model <- set_models(mean_fun="gaussian", err_dist=c("zip"))
#' fit <- senlm(model=model, data=haul, xvar="depth", yvar="Sebastolobus.altivelis")
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
#' @importFrom stats simulate
#' @export
#'
#' @examples
#'
#' \dontrun{
#'
#' ## Simulate data
#' model <- set_models(mean_fun="gaussian", err_dist="zip")
#' fit <- senlm(model=model, data=haul, xvar="depth", yvar="Sebastolobus.altivelis")
#' simulate(fit)
#'
#' ## Simulate two data sets with x=400,600,and 800
#' simulate(fit, nsim=2, newdata=c(400,600,800))
#' }
#'
simulate.senlm <- function (object, nsim=1, seed=NULL, newdata=NULL, ...) {
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
#' model <- set_models(mean_fun="gaussian", err_dist="zip")
#' fit <- senlm(model=model, data=haul, xvar="depth", yvar="Sebastolobus.altivelis")
#' summary(fit)
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
#' model <- set_models(mean_fun="gaussian", err_dist="zip")
#' fit <- senlm(model=model, data=haul, xvar="depth", yvar="Sebastolobus.altivelis")
#' plot(fit)
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
  FitMLE <- try ( suppressWarnings (
    bbmle::mle2 (minuslogl=NLL, optimizer="optim", method="SANN", vecpar=TRUE, start=u.theta0) ))
  
  ## --- Test if sann model fit
  if (class(FitMLE) != "try-error") {
    
    ## --- Fit using nlminb from sann estimate
    FitMLE <- try ( suppressWarnings (
      bbmle::mle2 (minuslogl=NLL, optimizer="nlminb", vecpar=TRUE, start=bbmle::coef(FitMLE)) ))
    
    ## --- Test if nlminb from sann fits
    if (class(FitMLE) != "try-error") {
      ## --- nlminb based on sann succeeded
      FitType   <- "nl-sann"
    } else {
      ## --- nlminb failed - use sann
      FitType   <- "sann"
    }

  } else {

    ## --- Fit using nlminb from init estimate
    FitMLE <- try ( suppressWarnings (
      bbmle::mle2 (minuslogl=NLL, optimizer="nlminb", vecpar=TRUE, start=u.theta0) ))
    
    ## --- Test if nlminb from init fits
    if (class(FitMLE) != "try-error") {
      ## --- nlminb  - use init
      FitType   <- "nl"
    } else {
      ## --- Fail
      FitType <- "fail"
    }
  }

  ## --- Initialise fit object
  Fit <- list()

  ## --- Set fitted parameters and confidence intervals on bounded parameter space
  if ( FitType =="fail" ) {

    ## --- Fit failed
    Fit$convergence <- 1
    Fit$theta <- NA * theta0
    Fit$lb    <- NA * theta0
    Fit$ub    <- NA * theta0
    Fit$u.hessian <- NA

  } else {
    
    ## --- Fit successful
    Fit$convergence <- 0
    
    ## --- Calculate confidence interval on unbounded parameter space
    u.theta   <- bbmle::coef(FitMLE)
    u.hessian <- attributes(FitMLE)$details$hessian
    colnames(u.hessian) <- names(u.theta)
    rownames(u.hessian) <- names(u.theta)

    ## --- Is standard error ok?
    StdErrOK <- TRUE
    
    ## --- Check if hessian is composed of finite numbers
    if (any(is.nan(u.hessian)))     { StdErrOK <- FALSE }
    if (any(is.na(u.hessian)))      { StdErrOK <- FALSE }
    if (any(!is.finite(u.hessian))) { StdErrOK <- FALSE }
    
    ## --- Check invertibility of hessian
    if (StdErrOK) {
      if (det(u.hessian) <= 0 ) { StdErrOK <- FALSE }
    }
    
    ## --- Check if standard errors are non-negative
    if (StdErrOK) { 
      u.stderr  <- sqrt(diag(solve(u.hessian)))
      if (any(u.stderr < 0)) { StdErrOK  <- FALSE }
      }
      
    ## --- Calculate approximate confidence intervals
    if (StdErrOK) {
      u.lb <- u.theta - 2*u.stderr
      u.ub <- u.theta + 2*u.stderr
    } else {
      u.lb <- rep (NA, length(u.theta))
      u.ub <- rep (NA, length(u.theta))
    }

    ## --- Calculate bounded fit / confidence interval
    Fit$theta <- make_bounded (ModelInfo, u.theta)
    Fit$lb    <- make_bounded (ModelInfo, u.lb)
    Fit$ub    <- make_bounded (ModelInfo, u.ub)
    ## --- Store unbounded hessian
    Fit$u.hessian <- u.hessian

  }

  ## --- Return fit
  return (Fit)
}

mle_constant_bernoulli <- function (ModelInfo, Dat) {
  ## --- MLE for constant-bernoulli mean function
  ## *** THIS FUNCTION IS OUTDATED / NOT CALLED ***
  
  ## MLE
  Fit.theta        <- mean (Dat$y)
  names(Fit.theta) <- c("H")

  ## Return MLE fit
  return (Fit.theta)
}

mle_constant_poisson <- function (ModelInfo, Dat) {
  ## --- MLE for constant-poisson mean function
  ## *** THIS FUNCTION IS OUTDATED / NOT CALLED ***

  ## Fit constant mean function
  Fit.theta <- mle_constant_bernoulli (ModelInfo, Dat)

  ## Return MLE fit
  return (Fit.theta)
}

mle_uniform_bernoulli <- function (ModelInfo, Dat) {
  ## --- MLE for uniform-bernoulli mean function
  ## *** THIS FUNCTION IS OUTDATED / NOT CALLED ***

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
#' 'msenlm' fits multiple species-environment non-linear models via maximum likelihood.
#'
#' @param models Object listing models to fit (from set_models function).
#' @param data A data frame containing 'x' (explanatory) and 'y' (response) variables.
#' @param xvar Name of explanatory variable (must be univariate).
#' @param yvar Names of response variables.
#' @param method If "crossed", fit all models to all response variables. If "paired",
#' fit first model to first response variables, etc.
#' 
#' @return Object containg model fits to data y and x.
#'
#' @keywords fit senlm model, mle
#'
#' @examples
#'
#' \dontrun{
#'
#' models <- set_models(mean_fun=c("gaussian","beta"), err_dist=c("zip","zinb"))
#' fits <- msenlm(models=models, data=haul, xvar="depth",
#'                yvar=c("Albatrossia.pectoralis", "Sebastolobus.altivelis"))
#' }
#' @export
#'
msenlm <- function (models=NULL, data=NULL, xvar=NULL, yvar=NULL, method="crossed") {
  ## --- Fit multiple senlm models using maximum likelihood to multiple response variables

  ## --- Check inputs
  
  ## Check x and y specified
  if (is.null(xvar) | is.null(yvar)) { stop ("Must specify xvar and yvar!") }

  if (!is.character(xvar)) { stop ("xvar must be character!") }
  if (length(xvar)>1) { stop ("Too many x variable names specified!") }
  if (any(is.na(match (xvar, names(data))))) { stop ("xvar not valid!") }

  if (!is.character(yvar)) { stop ("yvar must be character!") }
  if (any(is.na(match (yvar, names(data))))) { stop ("Some yvar not valid!") }

  if ( (method=="paired") & (nrow(models)!=length(yvar)) ) {
    stop ("Length of yvar and number of models must match if method='paired'!")
  }
  
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

  
  ## --- Initialise model names if method is paired
  if (method == "paired") { ModelNames <- rep (NA, length=nrow(models)) }

  ## --- Loop through response variables
  for (i in 1:length(Fits)) {
    
    ## --- Models and y variables paired
    if (method == "paired") {
   
      ## Create object to store model fits to data with ith response variable
      ModelFits <-  vector (mode="list", length=1)
 
      ## --- Fit model
      ModelFits[[1]]   <- senlm (model=models[i,], data=data, xvar=xvar, yvar=yvar[i])
      names(ModelFits) <- ModelFits[[1]]$model
    }
    
    ## --- Model and y variables crossed
    if (method == "crossed") {

      ## --- Models and y variables crossed
      
      ## Create object to store model fits to data with ith response variable
      ModelFits <-  vector (mode="list", length=nrow(models))
      ModelNames <- rep (NA, length=nrow(models))
      
      ## --- Loop throough models
      for (j in 1:length(ModelFits)) {
        
        ## --- Fit model
        Fit <- senlm (model=models[j,], data=data, xvar=xvar, yvar=yvar[i])
        
        ## --- Store model
        ModelFits[[j]] <- Fit
        ModelNames[j]  <- Fit$model
      }
    
      ## --- Add model names
      names(ModelFits) <- ModelNames
    }
    
    ## --- Store model fits
    Fits[[i]] <- ModelFits
  }

  ## --- Set class
  class (Fits) <- "msenlm"

  ## --- Return fits
  return (Fits)
}




#' Print summary of senlm multiple model fit
#'
#' Print summary information for mulitple senlm model fits. Select the
#' best model fit for each response variable using the given criteria.
#'
#' @param object Model fit.
#' @param best Select best model for each response variable according to
#' the following criteria: "nll", "AIC", "AICc", "BIC".
#' @param ... Other arguments that will be ignored.
#'
#' @return Data frame summarising senlm model fits.
#'
#' @keywords summary multiple senlm model fit
#'
#' @examples
#'
#' \dontrun{
#'
#' ## Summarise data
#' models <- set_models(mean_fun=c("gaussian", "beta"),
#'                      err_dist=c("zinb", "zip"), method="crossed")
#' fits <- msenlm(models=models, data=haul, xvar="depth",
#'                yvar=c("Albatrossia.pectoralis", "Sebastolobus.altivelis"))
#' summary(fits)
#' summary(fits, best="AICc")
#' }
#'
#' @export
#'
summary.msenlm <- function (object, best=NULL, ...) {
  ## --- Print summary of model fits

  ## --- Check best value
  if (!is.null(best)) {
    ## Possible values of best
    GOF <- c("nll", "AIC", "AICc", "BIC")
    ## Stop if best value is illegal
    if (all(best!=GOF)) { stop ('best option must be equal to "nll", "AIC", "AICc", "BIC"!') }
  }
  
  ## --- Grab parameter names of all models
  
  ## Initialise parameter names
  parnames <- c()
  ## Initalise model counts
  NModels <- 0
  ## Loop through response variables
  for (i in 1:length(object)) {
    ## Loop through models
    for (j in 1:length(object[[i]])) {
      ## Grab parameter names
      parnames <- c(parnames, names(object[[i]][[j]]$theta))
      ## Increment model counts
      NModels <- NModels +  1
    }
  } 
  ## Extract unique parameter names
  parnames <- unique(parnames)
  
  ## --- Set variable names for summary object
  varnames <- c("convergence", "y", "x", "model", "mean_fun", "err_dist",
                parnames, "npar", "nll", "AIC", "AICc", "BIC")
  
  ## --- Initialise summary object
  SDat <- as.data.frame (matrix(NA, ncol=length(varnames), nrow=NModels))
  names(SDat) <- varnames

  ## --- Loop through response variables

  ## Increment row counter
  Row <- 1

  ## Loop through response variables
  for (i in 1:length(object)) {
    ## Loop through models
    for (j in 1:length(object[[i]])) {

      ## Did model fit fail
      Convergence <- object[[i]][[j]]$convergence
      SDat[Row,]$convergence <- Convergence
      
      ## Grab x and y variable names
      SDat[Row,]$y <- object[[i]][[j]]$yname
      SDat[Row,]$x <- object[[i]][[j]]$xname

      ## Grab mean function and error distribution
      MeanErr  <- unlist(strsplit(object[[i]][[j]]$model, "_"))
      SDat[Row,]$model    <- object[[i]][[j]]$model
      SDat[Row,]$mean_fun <- MeanErr[1]
      SDat[Row,]$err_dist <- MeanErr[2]
      
      ## --- Was fit successful?
      if (Convergence == 0) {

        ## Grab fitted parameter values
        SDat[Row,names(object[[i]][[j]]$theta)] <- object[[i]][[j]]$theta

        ## Grab goodness-of-fit variables
        SDat[Row,c("npar", "nll", "AIC", "AICc", "BIC")] <-
          object[[i]][[j]]$IC[c("npar", "nll", "AIC", "AICc", "BIC")]
      }

      ## Increment row number
      Row <- Row + 1
    }
  }
  
  ## --- Find best models
  if (!is.null(best)) {
    
    ## Initialise row counter
    Row <- 0
    ## Initialise best model index for each response variable
    BestModel <- rep(NA,length(object))

    ## Loop through response variables
    for (i in 1:length(object)) {

      ## --- Put goodness-of-fit values in matrix
      NModels <- length (object[[i]])
      ICMat   <- as.data.frame(matrix(NA, ncol=4, nrow=NModels))
      names(ICMat) <- c("nll", "AIC", "AICc", "BIC")
      ## Loop through models
      for (j in 1:length(object[[i]])) {
        ## Grab goodness-of-fit values if available
        if (object[[i]][[j]]$convergence == 0) {
          ICMat[j,] <- object[[i]][[j]]$IC[grep("npar", names(object[[i]][[j]]$IC), invert=TRUE)]
        }
      }
      
      ## --- Grab select goodness-of-fit metric
      GOF <- ICMat[,best]
      
      ## --- Find row of summary object of best methods
      BestModel[i] <- Row + which(GOF==min(GOF, na.rm=T))
      
      ## Increment row counter
      Row <- Row + length(object[[i]])
    }
    
    ## --- Select best models
    SDat <- SDat[BestModel,]
  }
  
  ## Display summary object
  print (SDat)
}


plot.msenlm <- function (Fits) {
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


#' Find the best model for each response variable from an msenlm object
#'
#' Find the best model for each response variable from an msenlm object
#' using the given criteria.
#'
#' @param object msenlm multiple model fit.
#' @param best Select best model for each response variable according to
#' the following criteria: "nll", "AIC", "AICc", "BIC".
#'
#' @return msenlm object only containing best fitted models.
#' 
#' @keywords best multiple senlm model fit
#'
#' @examples
#'
#' \dontrun{
#'
#' ## Summarise data
#' models <- set_models(mean_fun=c("gaussian", "beta"),
#'                      err_dist=c("zinb", "zip"), method="crossed")
#' fits <- msenlm(models=models, data=haul, xvar="depth",
#'                yvar=c("Albatrossia.pectoralis", "Sebastolobus.altivelis"))
#' best <- msenlm.best(fits, best="AICc")
#' summary(best)
#' }
#'
#' @export
#'
msenlm.best <- function (object, best="AICc" ) {
  ## --- Extract best model for each response variable
   
  ## --- Check best value
  if (!is.null(best)) {
    ## Possible values of best
    GOF <- c("nll", "AIC", "AICc", "BIC")
    ## Stop if best value is illegal
    if (all(best!=GOF)) { stop ('best option must be equal to "nll", "AIC", "AICc", "BIC"!') }
  }
  
  ## --- Create best fits object
  BFits <- vector (mode="list", length=length(object))
  names(BFits) <- names(object)
  for (i in 1:length(BFits)) {
    BFits[[i]] <- vector (mode="list", length=1)
  }
  
  ## Loop through response variables
  for (i in 1:length(object)) {

    ## --- Put goodness-of-fit values in matrix
    NModels <- length (object[[i]])
    ICMat   <- as.data.frame(matrix(NA, ncol=4, nrow=NModels))
    names(ICMat) <- c("nll", "AIC", "AICc", "BIC")
    ## Loop through models
    for (j in 1:length(object[[i]])) {
      ## Grab goodness-of-fit values if available
      if (object[[i]][[j]]$convergence == 0) {
        ICMat[j,] <- object[[i]][[j]]$IC[grep("npar", names(object[[i]][[j]]$IC), invert=TRUE)]
      }
    }
    
    ## --- Grab select goodness-of-fit metric
    GOF <- ICMat[,best]
    
    ## --- Find row of summary object of best methods
    BFits[[i]][[1]] <- object[[i]][[which(GOF==min(GOF, na.rm=T))]]
    names(BFits[[i]]) <- BFits[[i]][[1]]$model
  }
  
  ## --- Set class
  class (BFits) <- "msenlm"

  ## --- Return best fits object
  return (BFits)
}

## --- SIMULATE DATA SET ---

## --- Set parameters for one model

#' Set simulated model parameters
#' 
#' Set parameters for one model.
#'
#' @param mean_fun Mean function of model.
#' @param err_dist Error distribution of model.
#' @param Dat Data frame containg x and y to estimate tail-adjusted beta "delta" parameter.
#' @param thetaM Vector of parameters values for mean function.
#' @param thetaE Vector of parameters values for error distribution.
#' @param thetaC Vector of parameters values for constants.
#' 
#' @return List containing parameters to simulated data set.
#' 
#' @keywords simulation parameters
#'
#' @examples
#'
#' ## Sim simulation parameters for gaussian-negbin model
#' set_sim_par(mean_fun="gaussian", err_dist="negbin",
#'             thetaM=c(H=60,m=50,s=10), thetaE=c(phi=1.5))
#' ## Sim simulation parameters for binomial model
#' set_sim_par(mean_fun="gaussian", err_dist="binomial_count",
#'             thetaM=c(H=60,m=50,s=10), thetaC=c(binomial_n=100))
#' 
#' @export
#'
set_sim_par <- function (mean_fun=NULL, err_dist=NULL, Dat=NULL,
                         thetaM=NULL, thetaE=NULL, thetaC=NULL) {
  ## --- Set simulated model parameters
  
  ## --- Stop if error distribution and error class are not specified
  if ( is.null(err_dist) | is.null(mean_fun) ) { stop ("Must specify err_dist and mean_fun!")}

  ## --- Check mean functions and error distributions
  check_mean_functions_and_error_distributions (mean_fun, err_dist)
  
  ## --- Only one mean function, and error distribution allowed
  if (length(mean_fun) > 1) { stop ("Only one mean function allowed!") }
  if (length(err_dist) > 1) { stop ("Only one error distribution allowed!") }
  
  ## --- Set error class
  err_class <- error_distributions (err_dist=err_dist)
  
  ## --- Set model information
  ModelInfo <- set_model_info (mean_fun=mean_fun, err_dist=err_dist, thetaC=thetaC)
  
  ## --- Initialise object
  Par <- list()

  ## --- Model name
  Par$model_name <- ModelInfo$model_name 
  
  ## --- Constants defined via thetaC or err_dist
  thetaC <- estimate_constant (ModelInfo, data=NULL)
  
  ## --- Get default parameters
  ParDefault <- get_default_par (ModelInfo)
  
  ## --- Mean function
  Par$mean_fun <- mean_fun

  ## --- Data model 
  Par$err_dist <- err_dist
  
  ## --- Error class
  Par$err_class <- err_class
  
  ## --- Use default parameters if none others supplied
  if (is.null(thetaM)) { thetaM <- ParDefault$thetaM }
  if (is.null(thetaE)) { thetaE <- ParDefault$thetaE }
  if (is.null(thetaC)) { thetaC <- ParDefault$thetaC }
  
  ## --- Store parameters
  Par$thetaM <- thetaM       ## Mean parameters
  Par$thetaE <- thetaE       ## Error parameters
  Par$thetaC <- thetaC       ## "Constants" parameters
  
  ## --- Combine parameters
  Par$theta <- c(thetaM, thetaE, thetaC)
  
  ## --- Check if parameter values are within bounds
  parslegal <- check_par_values (ModelInfo=ModelInfo, theta=Par$theta)
  if (!parslegal) { stop ("Parameter values are out of bounds!") }
  
  ## --- Set class
  class (Par) <- "senlm_par"
  
  ## --- Return parameter object
  return (Par)
}


get_default_par <- function (ModelInfo) {
  ## --- Get default parameter values for simulated data
  
  ## --- Grab model name, error distribution, and mean function
  model_name <- ModelInfo$model_name
  err_dist <- ModelInfo$err_dist
  mean_fun <- ModelInfo$mean_fun
  
  ## --- Grab constants
  thetaC <- estimate_constant (ModelInfo, data=NULL)

  ## --- Error distributions
  if (err_dist == "bernoulli")         { thetaE <- NULL }
  if (err_dist == "binomial_count")    { thetaE <- NULL }
  if (err_dist == "binomial_prop")     { thetaE <- NULL }
  if (err_dist == "poisson")           { thetaE <- NULL }
  if (err_dist == "negbin")            { thetaE <- c(phi=1.5) }
  if (err_dist == "zip")               { thetaE <- c(pi=0.3) }
  if (err_dist == "zinb")              { thetaE <- c(pi=0.3, phi=1.5) }
  if (err_dist == "zipl")              { thetaE <- c(g0=-0.9, g1=-0.5) }
  if (err_dist == "zipl.mu")           { thetaE <- c(g0=-0.9, g1=-0.05) }
  if (err_dist == "zinbl")             { thetaE <- c(g0=-0.9, g1=-0.5,  phi=1.5) }
  if (err_dist == "zinbl.mu")          { thetaE <- c(g0=-0.9, g1=-0.05, phi=1.5) }
  if (err_dist == "gaussian")          { thetaE <- c(sigma=5) }
  if (err_dist == "tweedie")           { thetaE <- c(phi=1.5, rho=1.1) }
  if (err_dist == "zig")               { thetaE <- c(pi=0.3, phi=1.5) }
  if (err_dist == "zigl")              { thetaE <- c(g0=-0.9, g1=-0.5, phi=0.5) }
  if (err_dist == "zigl.mu")           { thetaE <- c(g0=-0.9, g1=-0.5, phi=0.5) }
  if (err_dist == "tab")               { thetaE <- c(sigma=0.1) }
  if (err_dist == "zitab")             { thetaE <- c(pi=0.3, sigma=0.1) }
  
  ## --- Mean functions
  if (mean_fun == "constant")          { thetaM <- c(H=60) }
  if (mean_fun == "uniform")           { thetaM <- c(H=60, c=30, d=70) }
  if (mean_fun == "gaussian")          { thetaM <- c(H=60, m=50, s=10) }
  if (mean_fun == "mixgaussian")       { thetaM <- c(H=60, a=0.7, m1=35, m2=65, s1=10, s2=5) }
  if (mean_fun == "mixgaussian_equal") { thetaM <- c(H=60, a=0.7, m1=35, m2=65, s=5) }
  if (mean_fun == "beta")              { thetaM <- c(H=60, c=30, d=70, u=0.4, v=0.2) }
  if (mean_fun == "sech")              { thetaM <- c(H=60, m=50, s=5, r=0.75, p=1.5) }
  if (mean_fun == "sech_p1")           { thetaM <- c(H=60, m=50, s=5, r=0.75) }
  if (mean_fun == "sech_r0p1")         { thetaM <- c(H=60, m=50, s=5) }
  if (mean_fun == "hofII")             { thetaM <- c(H=60, m=50, w0=-0.3) }
  if (mean_fun == "hofIV")             { thetaM <- c(H=60, m=50, w=0.3, k=10) }
  if (mean_fun == "hofIVb")            { thetaM <- c(H=60, m=50, w=0.3) }
  if (mean_fun == "hofV")              { thetaM <- c(H=60, m=50, w1=0.2, w2=0.4, k=10) }
  if (mean_fun == "hofVb")             { thetaM <- c(H=60, m=50, w1=0.2, w2=0.4, k=10) }
  
  ## --- Constant
  if (is.null(thetaC)) {
    ## binomial_n
    if ( (err_dist == "binomial_count") | (err_dist == "binomial_prop") ) {
      thetaC <- c(binomial_n=40)
    }
    ## delta
    if ( (err_dist == "tab") | (err_dist == "zitab") ) {
      thetaC <- c(delta=0.01)
    }
  }

  ## --- Rescale H for bernoulli/percentage models
  if ( (err_dist == "bernoulli") | (err_dist == "binomial_prop") |
       (err_dist == "tab")       | (err_dist == "zitab")         ) {
    if (any(names(thetaM) == "H"))  { thetaM["H"]  <- 0.9 }
  }
  
  ## --- Binomial count
  if (err_dist == "binomial_count") {
    ## Make H 0.9 times binomial n parameter
    if (any(names(thetaM) == "H"))  { thetaM["H"] <- round(0.9*thetaC["binomial_n"]) }
  }
  
  ## --- Set error class
  err_class <- error_distributions (err_dist=err_dist)
  
  ## --- Store parameter names
  Par <- list ()
  Par$model_name <- model_name
  Par$err_class <- err_class
  Par$thetaM <- thetaM
  Par$thetaE <- thetaE
  Par$thetaC <- thetaC
  
  ## --- Return parameter names
  return (Par)
}


#' Create default parameter list
#' 
#' Create parameter object using default values.
#' 
#' @param models Models defined by set_models() object.
#' @param mean_fun Mean function.
#' @param mean_class String that indicates the class of mean functions to list.
#' @param err_dist Error distribution.
#' @param err_class String that indicates the class of error distributions to list.
#' @param binomial_n Value of binomial n parameter for binomial models.
#' @param method How mean function and error distributions should be combined;
#' "paired" or "crossed".
#'
#' @return List of parameter simulation objects.
#' 
#' @keywords Simulation parameters
#'
#' @examples
#'
#' ## Create default parameters from set_models() object
#' Models <- set_models (mean_fun=c("gaussian", "beta", "sech"),
#'                       err_dist=c("poisson","zip"), method="crossed")
#' Pars <- create_default_par_list (Models)
#'
#' ## Create parameter object by pairing mean function and error distribution
#' Pars <- create_default_par_list (method="paired",
#'                                  mean_fun=c("gaussian", "beta"), 
#'                                  err_dist=c("zip", "zinb"))
#'
#' ## Create parameter object by crossing mean function and error distribution
#' Pars <- create_default_par_list (method="crossed",
#'                                  mean_fun=c("gaussian", "beta"), 
#'                                  err_dist=c("poisson", "zip", "zinb"))
#' 
#' @export
#' 
create_default_par_list <- function (models=NULL,
                                     mean_fun=NULL, mean_class=NULL,
                                     err_dist=NULL, err_class=NULL,
                                     binomial_n=NA,
                                     method="paired") {
  ## --- Create default model parameter object
  
  if (!is.null(models)) {
    if (any(class(models)=="model_list")) {
      ## --- Parameters are specified by a set_models() object
      
      ## --- Grab mean functions, error distributions, and binomial_n
      mf <- as.character(models$mean_fun)          
      ed <- as.character(models$err_dist)
      bn <- models$binomial_n
    } else {
      stop ("Models argument should be specified with set_models()!" )
    }
  } else {
    ## --- Mean functions, error distributions, and binomial_n are
    ##     supplied separately in vectors
    
    ## --- Stop if mean function and mean class are both specified
    if (!is.null(mean_fun) & !is.null(mean_class)) {
      stop ("Do not specify mean_fun **and** mean_class!")
    }
    ## --- Stop if error distribution and error class both specified
    if (!is.null(err_dist) & !is.null(err_class)) {
      stop ("Do not specify err_dist **and** err_class!")
    }
    
    ## --- If error distribution and error class not specifed set error class to "all"
    if ( is.null(err_dist) & is.null(err_class)) { err_class <- "all"; method <- "crossed" }
    ## --- Set error distribution if error class supplied
    if (!is.null(err_class)) { 
      err_dist <- error_distributions (err_class=err_class)
    }
    
    ## --- If mean function and mean class not specifed set class to "all"
    if ( is.null(mean_fun) & is.null(mean_class) ) { mean_class <- "all";  method <- "crossed" }
    ## --- Set mean function if class supplied
    if (!is.null(mean_class)) { 
      mean_fun <- mean_functions (mean_class=mean_class)
    }
    
    ## --- Set mean functions to all if not set  
    if (is.null(mean_fun)) { mean_fun <- mean_functions(mean_class="all") }
    
    ## --- Check mean functions and error distributions
    check_mean_functions_and_error_distributions (mean_fun, err_dist)
    
    ## --- Remove NA from binomial_n, unless all are NA then set to singular NA 
    ##     (binomial_n must be NA, an integer, or a vector of integers each each binomial model)
    if (all(is.na(binomial_n))) {
      binomial_n <- NA
    } else {
      binomial_n <- binomial_n[!is.na(binomial_n)]
    }
    
    ## --- Create data frame of models
    if (is.na(binomial_n)) {
      ## Reset binomial n to NA (it will be set by set_sim_par)
      ## Set temporarily to 1 to avoid error from input check
      DF <- set_models (mean_fun=mean_fun, err_dist=err_dist, binomial_n=1, method=method)
      DF$binomial_n <- NA
    } else {
      ## Use supplied binomial_n
      DF <- set_models (mean_fun=mean_fun, err_dist=err_dist, binomial_n=binomial_n, method=method)
    }
    
    ## --- Grab mean functions and error distributions
    mf <- as.character(DF$mean_fun)          
    ed <- as.character(DF$err_dist)
    bn <- DF$binomial_n
  }
  
  ## --- Initalise parameter list
  Pars <- vector (mode="list", length=length(mf))

  ## --- Create default parameter objects
  for (i in 1:length(Pars)) {

    ## Set binomial n parameter
    if (is.na(bn[i])) {
      ## Use default n
      thetaC <- NULL
    } else {
      ## Use supplied n
      thetaC <- c(binomial_n=bn[i])
    }
    
    ## Create default parameter object
    Pars[[i]] <- set_sim_par (mean_fun=mf[i], err_dist=ed[i], thetaC=thetaC)
  }
  
  ## --- Return parameter objects
  return (Pars)
}


simulate_data <- function (x, Par, seed=NULL) {
  ## --- Simulate data
  
  ## --- Grab error distribution
  err_dist <- Par$err_dist
  
  ## Create data object
  Dat <- list()

  ## Random number seed
  if ( is.null (seed) ) {
    ## Define seed if not given
    seed <- as.numeric(paste(c(sample(1:9, 1), sample(0:9,(5-1))), collapse=""))
  }
  ## Set seed
  set.seed (seed)
  ## Store seed
  Dat$seed <- seed
  
  ## Calculate mean function
  mu <- do.call (paste("mu_", Par$mean_fun, sep=""), list(Par$thetaM, x))
  
  ## --- SIMULATE DATA

  ## --- Bernoulli
  if (err_dist == "bernoulli" ) {
    ## --- Bernoulli (0/1 Prescence/Abscence data).
    y <- as.numeric(stats::runif (n=length(x)) <= mu)
  }
  
  ## --- Binomial - count
  if (err_dist == "binomial_count" ) {
    n <- as.list(Par$thetaC)$binomial_n
    p <- mu/n
    y <- stats::rbinom (n=length(x), prob=p, size=n)
  }
  
  ## --- Binomial - prop
  if (err_dist == "binomial_prop" ) {
    n <- as.list(Par$thetaC)$binomial_n
    p <- mu
    y <- (stats::rbinom (n=length(x), prob=p, size=n)) / n
  }
  
  ## --- Poisson
  if (err_dist == "poisson" ) {
    y <- stats::rpois (n=length(x), lambda=mu)
  }
  
  ## --- Negative binomial (Gamma-Poisson)
  if (err_dist == "negbin" ) {
    ## --- Generate a Gamma (shape=alpha, rate=beta) [mean=alpha/beta]
    ## --- then a Poisson with the mean given by the Gamma rv.
    ## --- Mean of the NegBin is mu = alpha/beta,
    ## --- and variance is sigam = mu + phi*mu^2.
    ## --- mu is given by the mean function eta.

    ## If phi=0 then data=poisson
    if (as.list(Par$thetaE)$phi == 0) {
      poisson.means <- mu
    } else {
      ## Convert mean and dispersion parameter to shape parameters of Gamma distribution
      alpha <- 1/as.list(Par$thetaE)$phi
      beta  <- alpha/mu
      
      ## Simulate data
      poisson.means <- stats::rgamma(n=length(x), shape=alpha, rate=beta)
    }
    
    ## Sample Poisson given lambda
    y <- rep(0, length(x))
    for (i in 1:length(x)) {
      y[i] <- stats::rpois(1, lambda=poisson.means[i])
    }
  }

  ## --- Zero-inflated Poisson
  if (err_dist == "zip" ) {

    ## Generate Poisson data
    y <- stats::rpois (n=length(x), lambda=mu)
    
    ## Randomly change pi proportion to zero (extra zeros)
    y[stats::runif(length(x)) < as.list(Par$thetaE)$pi] <- 0
  }
    
  ## --- Zero-inflated negative binomial
  if (err_dist == "zinb" ) {
    
    ## Convert mean and dispersion parameter to shape parameters of Gamma distribution
    alpha <- 1/as.list(Par$thetaE)$phi
    beta  <- alpha/mu
    
    ## Generate poisson means
    poisson.means <- stats::rgamma(n=length(x), shape=alpha, rate=beta)
    
    ## Sample Poisson given lambda
    y <- rep(0, length(x))
    for (i in 1:length(x)) {
      y[i] <- stats::rpois(1, lambda=poisson.means[i])
    }
 
    ## Randomly change pi proportion to zero (extra zeros)
    y[stats::runif(length(x)) < as.list(Par$thetaE)$pi] <- 0
  }
  
  ## --- Zero-linked Poisson
  if ( (err_dist == "zipl") | (err_dist == "zipl.mu") ) {

    ## Grab parameters
    g0 <- as.list(Par$theta)$g0
    g1 <- as.list(Par$theta)$g1
    NonZero <- (mu > 0)

    ## Create y vector
    y <- rep (NA, length(mu))
    ## y is zero if mean is zero
    y[!NonZero] <- 0

    ## Grab non-zero means
    mu1 <- mu[NonZero]

    ## Create spike parameter
    if (err_dist == "zipl") {
      logitpi <- g0 + g1*log(mu1)
    }
    if (err_dist == "zipl.mu") {
      logitpi <- g0 + g1*mu1
    }
    pi <- exp(logitpi) / ( 1 + exp(logitpi) )
 
    ## Generate Poisson data
    y[NonZero] <- stats::rpois (n=sum(NonZero), lambda=mu1)
    
    ## Randomly change pi proportion to zero (extra zeros)
    y[stats::runif(sum(NonZero)) < pi] <- 0
  }

  ## --- Zero-linked negative binomial
  if ( (err_dist == "zinbl") | (err_dist == "zinbl.mu") ) {
    
    ## Grab parameters
    g0  <- as.list(Par$theta)$g0
    g1  <- as.list(Par$theta)$g1
    phi <- as.list(Par$theta)$phi
    NonZero <- (mu > 0)

    ## Convert mean and dispersion parameter to shape parameters of Gamma distribution
    alpha <- 1/phi
    beta  <- alpha/mu

    ## Generate poisson means
    poisson.means <- stats::rgamma(n=length(x), shape=alpha, rate=beta)

    ## Create y vector
    y <- rep (NA, length(poisson.means))
    ## y is zero if mean is zero
    y[!NonZero] <- 0

    ## Grab non-zero means
    mu1 <- poisson.means[NonZero]

    ## Create spike parameter
    if (err_dist == "zinbl") {
      logitpi <- g0 + g1*log(mu1)
    }
    if (err_dist == "zinbl.mu") {
      logitpi <- g0 + g1*mu1
    }
    pi <- exp(logitpi) / ( 1 + exp(logitpi) )
    
    ## Generate Poisson data
    y[NonZero] <- stats::rpois (n=sum(NonZero), lambda=mu1)
    
    ## Randomly change pi proportion to zero (extra zeros)
    y[stats::runif(sum(NonZero)) < pi] <- 0
  }

  ## --- Tail-adjusted beta error and Zero-inflated tail-adjusted beta
  if ( (err_dist == "tab" ) | (err_dist == "zitab" ) ) {

    ## Grab parameters
    delta <- as.list(Par$thetaC)$delta
    sigma <- as.list(Par$thetaE)$sigma

    ## Set beta distribution parameters
    Alpha <- mu / sigma
    Beta  <- (1 - mu) / sigma

    ## Generate beta data
    y <- stats::rbeta (n=length(x), shape1=Alpha, shape2=Beta)

    ## Remove values greater than 1 or less than 0
    y[y<0] <- 0;  y[y>1] <- 1

    ## If mean is 0 or 1, set y to 0 or 1 respectively
    y[mu==0] <- 0 ; y[mu==1] <- 1

    ## Set values < delta to 0 and > (1-delta) to 1
    y[y < delta] <- 0  ; y[y > (1-delta)] <- 1
  }
  
  ## --- Zero-inflated tail-adjusted beta
  if (err_dist == "zitab") {
    
    ## Grab parameters
    pi <- as.list(Par$thetaE)$pi
    
    ## Randomly change pi proportion to zero (extra zeros)
    y[stats::runif(length(x)) < pi] <- 0
  } 
  
  ## --- Gaussian error (data non-negative, zero when mean zero)
  if (err_dist == "gaussian" ) {
    y <- stats::rnorm (n=length(x), mean=mu, sd=as.list(Par$thetaE)$sigma)
    y[y<0] <- 0
    y[mu==0] <- 0
  }
  
  ## --- Tweedie error (data non-negative, zero when mean zero)
  if (err_dist == "tweedie" ) {
    ## --- Find where mean is zero
    Index <- (mu > 0)
    ## Ignore observations where mean is zero
    Mu    <- mu[Index]
    y     <- rep (0, length(x))
    ## Simulate data
    y[Index] <- tweedie::rtweedie (n=sum(Index), power=as.list(Par$thetaE)$rho,
                                   mu=Mu, phi=as.list(Par$thetaE)$phi)
    ## Make sure observation are zero when mean is zero
    y[y < 0] <- 0
    y[mu==0] <- 0
  }
  
  ## --- Zero-inflated Gamma error (data non-negative, zero when mean zero)
  if (err_dist == "zig" ) {

    ## Grab parameters
    pi  <- as.list(Par$theta)$pi
    phi <- as.list(Par$theta)$phi

    ## --- Find where mean is positive
    NonZero <- (mu > 0)

    ## Ignore observations where mean is zero
    Mu    <- mu[NonZero]
    y     <- rep (0, length(x))
    
    ## If phi=0 then data=mean
    if (phi == 0) {
      y[NonZero] <- Mu
    } else {
      ## Convert mean and dispersion parameter to shape parameters of Gamma distribution
      alpha <- 1/phi
      beta  <- alpha/Mu
      
      ## Simulate data
      y[NonZero] <- stats::rgamma(n=sum(NonZero), shape=alpha, rate=beta)
    }
    
    ## Make sure observation are zero when mean is zero
    y[y < 0] <- 0
    y[mu==0] <- 0

    ## Randomly change pi proportion to zero (extra zeros)
    y[stats::runif(length(x)) < pi] <- 0
  }
  
  
  ## --- Zero-linked gamma
  if ( (err_dist == "zigl") | (err_dist == "zigl.mu") ) {

    ## Grab parameters
     g0  <- as.list(Par$theta)$g0
    g1  <- as.list(Par$theta)$g1
    phi <- as.list(Par$theta)$phi

    ## --- Find where mean is positive
    NonZero <- (mu > 0)
    
    ## Ignore observations where mean is zero
    Mu    <- mu[NonZero]
    y     <- rep (0, length(x))

    ## If phi=0 then data=mean
    if (phi == 0) {
      y[NonZero] <- Mu
    } else {
      ## Convert mean and dispersion parameter to shape parameters of Gamma distribution
      alpha <- 1/phi
      beta  <- alpha/Mu
      
      ## Simulate data
      y[NonZero] <- stats::rgamma(n=sum(NonZero), shape=alpha, rate=beta)
    }
    
    ## Make sure observation are zero when mean is zero
    y[y < 0] <- 0
    y[!NonZero] <- 0

    ## Create spike parameter
    if (err_dist == "zigl")    { logitpi <- g0 + g1*log(Mu) }
    if (err_dist == "zigl.mu") { logitpi <- g0 + g1*Mu      }
    pi <- exp(logitpi) / ( 1 + exp(logitpi) )
    
    ## Randomly change pi proportion to zero (extra zeros)
    y[stats::runif(sum(NonZero)) < pi] <- 0
  }
  
  ## Store data
  Dat$data_name  <- "Simulation"
  Dat$model_name <- Par$model_name
  Dat$err_dist   <- err_dist
  Dat$mean_fun   <- Par$mean_fun
  Dat$transform  <- "raw"
  Dat$N          <- length(x)
  Dat$y          <- y
  Dat$x          <- x
  Dat$mu         <- mu
  Dat$thetaE     <- Par$thetaE
  Dat$thetaM     <- Par$thetaM
  
  ## Return data object
  return (Dat)
}


#' Create list of simulated data sets
#' 
#' Simulate data sets defined by the given parameter list.
#'
#' @param Pars List containg parameters for the simulated data sets.
#' @param x Vector of x values to simulated data. If not given, N, xmin, and xmax must be).
#' @param N (optional) Sample size of data set to simulate. Default value 500.
#' @param xmin (optional) Minimum value of x variable. Default value 0.
#' @param xmax (optional) Maximum value of x variable. Default value 100.
#' @param seed Random number seed.
#' @param echo If TRUE display list describing data sets simulated
#' 
#' @return Data frame containing one x variable and simulated data sets from Pars object.
#' 
#' @keywords simulated data sets
#' 
#' @examples
#'
#'
#' ## Set model
#' Models <- set_models (mean_fun="gaussian", err_dist="negbin")
#'
#' ## Set list of default parameters for each model
#' Pars <- create_default_par_list (Models)
#'
#' ## Simulate data
#' Data <- create_simulated_datasets (Pars, seed=12345)
#' 
#' @export
#' 
create_simulated_datasets <- function (Pars, x=NULL, N=500, xmin=0, xmax=100,
                                       seed=NULL, echo=FALSE) {
  ## --- Create simulated data set object
  
  ## --- Define random number seed
  if ( is.null (seed) ) {
    seed <- as.numeric(paste(c(sample(1:9, 1), sample(0:9,(5-1))), collapse=""))
  }
  ## Set random seed
  set.seed (seed)
 
  ## Store given x values or simulate new ones
  if (is.null(x)) {
    ## --- Simulate x values
    
    ## --- Check if N, xmin, and xmax are univariate, numeric
    if (!((N[1]==round(N[1])) & (length(N)==1))) { stop ("N must be univariate integer/numeric!") }
    if (!(is.numeric(xmin) & (length(xmin)==1))) { stop ("xmin must be univariate numeric!") }
    if (!(is.numeric(xmax) & (length(xmax)==1))) { stop ("xmax must be univariate numeric!") }
    if (xmin > xmax) { stop ("xmin must be less than xmax!") }
    
    ## Create random x between xmin and xmax
    x <- sort(stats::runif(N, xmin, xmax))
    
  } else {
    ## --- A data frame of x values is given
    if (!is.vector(x))  { stop ("X must be a vector!") }
    N    <- length(x)
    xmin <- NULL
    xmax <- NULL
  }
  

  ## --- Create data sets

  ## --- Set data object
  Data <- vector (length=length(Pars), mode="list")
  
  ## Loop through parameter sets
  for (i in 1:length(Pars)) {

    ## Generate data
    Data[[i]] <- simulate_data (x=x, Par=Pars[[i]])

    ## Data name
    Data[[i]]$data_name <- paste("Simulation: (N=",Data[[i]]$N,
                                 ", ",Data[[i]]$model_name,")", sep="")
  } 
  
  ## Print data models
  if (echo) {
    Mat <- as.matrix(sapply (Data, '[[', "data_name"))
    row.names(Mat) <- paste (1:length(Data), ":", sep="")
    colnames(Mat)  <- "Data sets"
    print (Mat)
  }

  ## --- Convert data to data frame

  ## Initialise data frame
  df <- as.data.frame (matrix(NA, nrow=length(x), ncol=length(Pars)+1))

  ## Store x
  df[,1] <- x

  ## Create variable names
  Names <- rep ("", ncol(df))
  Names[1] <- "x"

  ## Store data and grab y variable names
  for (i in 1:length(Pars)) {
    ## Response variable
    df[,i+1] <- Data[[i]]$y
    ## Variabe name - convert "-" to "_"
    Names[i+1] <- gsub ("-", "_", Data[[i]]$model_name)
  }
  ## Add names to data frame
  names(df) <- Names
  
  ## Return data object
  return (df)
}


## --- SUMMARISE DATA SET ---

summarise_datasets <- function (Data) {
  ## --- Summarise data sets
  
  ## Initialise summary statistics
  Summary <- matrix (NA, nrow=length(Data), ncol=4)
  DataNames <- rep ("", length(Data))
  
  ## Calculate summary statistics
  for (i in 1:length(Data)) {
    ## Grab data name
    DataNames[i] <- Data[[i]]$data_name
    ## Sample size
    Summary[i,1] <- length (Data[[i]]$y)
    ## Mean
    Summary[i,2] <- mean (Data[[i]]$y, na.rm=TRUE)
    ## Standard deviaton
    Summary[i,3] <- stats::sd (Data[[i]]$y, na.rm=TRUE)
    ## Variance
    Summary[i,4] <- stats::var (Data[[i]]$y, na.rm=TRUE)
  }
  
  ## Add row and column names
  colnames (Summary) <- c("N", "mean", "sd", "var")
  rownames (Summary) <- DataNames
  
  ## Convert to data frame
  Summary <- as.data.frame (Summary)
  
  ## Return summary object
  return (Summary)
}


classify_dataset <- function (Dat) {
  ## --- Classify data set as;
  ## binary     : Prescence/absence {0,1}
  ## percentage : Proportions [0,1]  (includes binomial proportions
  ## count      : Discrete counts; 0,1,2,...
  ## abundance  : Non-negative continuous; >= 0
  ## real       : Unbounded continuous (no models for this at this time)
  ##
  ## Note: Data cannot be classified as binomial count data based on data alone.
  
  ## Grab y data
  y <- Dat$y
  
  ## Classify y variable
  if ( all ( (y==0) | (y==1) ) ) {
    ## Prescence/Abscence (Binary data - 0,1)
    err_class <- "binary"
  } else if ( (min(y) >= 0) & (max(y) <= 1) ) {
    ## Percentage
    err_class <- "percentage"
  } else if ( all (y==round(y,0)) & (min(y)>=0) ) {
    ## Count data (non-negative integer)
    err_class <- "count"
  } else if ( min(y)>=0 ) {
    ## Abundance data (non-negative real)
    err_class <- "abundance"
  } else {
    ## Real data (unbounded)
    err_class <- "real"
  }

  ## Return data classification
  return (err_class)
}


## --- SELECT DATA SETS ---

select_datasets <- function (Data, Index) {
  ## --- Set selected simulated data sets given by Index

  ## Store original data sets
  Data0 <- Data

  ## Create new list to store subsetted data sets
  Data <- vector (length=length(Index), mode="list")

  ## Loop through selected data sets
  for (k in 1:length(Index)) {
    ## Select kth data set
    Data[[k]] <- Data0[[Index[k]]]
  }
  
  ## Return selected data sets
  return (Data)
  
}

select_dataset <- function (Data, k) {
  ## --- Set one simulated data set

  ## Select kth data set
  Dat <- Data[[k]]
  
  ## Return data set
  return (Dat)
}


## --- TRANSFORM DATA SETS ---

transform_datasets <- function (Data, log10trans=NULL, logtrans=NULL, sqrttrans=NULL) {
  ## --- Transform selected repsonse variables

  ## --- Log10 function (modified base 10 log function)
  ## y' = y            , y <= 1 
  ## y' = log10(y) + 1 , y  > 1
  if (!is.null(log10trans)) {
    for (i in 1:length(log10trans)) {
      k <- log10trans[i]
      Data[[k]]$transform <- "log10(y)+1"
      Data[[k]]$data_name <- paste (Data[[k]]$data_name, "-", Data[[k]]$transform)
      Data[[k]]$y         <- Log10 (Data[[k]]$y)
    }
  }
  
  ## --- log function
  ## y' = log(y + 1) 
  if (!is.null(logtrans)) {
    for (i in 1:length(logtrans)) {
      k <- logtrans[i]
      Data[[k]]$transform <- "log (y + 1)"
      Data[[k]]$data_name <- paste (Data[[k]]$data_name, "-", Data[[k]]$transform)
      Data[[k]]$y         <- log (Data[[k]]$y + 1)
    }
  }

  ## --- Sqrt function
  ## y' = sqrt(y)
  if (!is.null(sqrttrans)) {
    for (i in 1:length(sqrttrans)) {
      k <- sqrttrans[i]
      Data[[k]]$transform <- "sqrt"
      Data[[k]]$data_name <- paste (Data[[k]]$data_name, "-", Data[[k]]$transform)
      Data[[k]]$y         <- sqrt (Data[[k]]$y)
    }
  }
  
  ## --- Return data
  return (Data)
}

Log10 <- function (y) {
  ## --- Modified base 10 log function
  ## y' = y            , y <= 1 
  ## y' = log10(y) + 1 , y  > 1

  ## Transform data
  y[y>1] <- log10 (y[y>1]) + 1

  ## Return data
  return (y)
}
  

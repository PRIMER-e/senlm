## ---
## --- Functions for model / parameter definitions
## ---

## --- --- --- ---

#' Mean functions and error distributions
#'
#' Display all mean functions and error distributions of given classes.
#'
#' @param mean_class String that indicates the class of mean functions to list.
#' Possible values are:
#' \describe{
#'   \item{"test"}{Constant and uniform mean functions.}
#'   \item{"main"}{Main mean functions.}
#'   \item{"sub"}{Special cases of main mean functions.}
#'   \item{"all"}{All mean functions.}
#' }
#' 
#' @param err_class String that indicates the class of error distributions to list.
#' \describe{
#'   \item{"binary"}{0/1 data; ("bernoulli"). }
#'   \item{"binomial"}{Binomial data represented as counts or proportions;
#'       ("binomial_count", "binomial_prop"). }
#'   \item{"percentage"}{[0-1] data;
#'       ("tab", "zitab").}
#'   \item{"count"}{Count data;
#'       ("poisson", "negbin", "zip", "zinb", "zipl", "zinbl", "zipl.mu", "zinbl.mu").}
#'   \item{"abundance"}{Non-negative continous data;
#'       ("gaussian", "tweedie", "zig", "zigl", "zigl.mu").}
#' }
#' 
#' @return List containing vectors of mean functions and error distributions of given classes.
#' 
#' @keywords mean function, error distributions
#'
#' @examples
#'
#' ## Display all mean functions and error distributions
#' print_models()
#' 
#' ## Display main mean functions and count data error distributions
#' print_models(mean_class="main", err_class="count")
#'
#' @export
#' 
print_models <- function (mean_class="all", err_class="all") {
  ## --- Definitions of all possible error distributions and mean functions

  ## --- Define mean functions
  MF <- mean_functions (mean_class=mean_class)
  
  ## --- Define error distributions
  ED <- error_distributions (err_class=err_class)
  
  ## --- Store possible models
  Models    <- list()
  Models$mean_functions      <- MF
  Models$error_distributions <- ED
  
  ## --- Return possible models
  return (Models)
}


#' Mean functions
#' 
#' Lists all possible mean functions of given class.
#'
#' @param mean_fun Vector of mean function(s).
#' @param mean_class String that indicates the class of mean functions to list.
#' Possible values are:
#' \describe{
#'   \item{"test"}{Constant and uniform mean functions.}
#'   \item{"main"}{Main mean functions.}
#'   \item{"sub"}{Special cases of main mean functions.}
#'   \item{"all"}{All mean functions.}
#' }
#'
#' @return Vector of mean functions of given class.
#' @return If no arguments are supplied, return a data frame listing all possible mean functions
#' with corresponding classes; if a vector of classes is supplied, return corresponding mean
#' functions; if a vector of mean functions is supplied, return corresponding classes.
#' 
#' @keywords mean function
#'
#' @examples
#'
#' ## Print data frame of all mean functions with corresponding classes
#' mean_functions()
#' 
#' ## Constant and uniform mean functions to test code
#' mean_functions(mean_class="test")
#'
#' ## Main mean functions
#' mean_functions(mean_class="main")
#'
#' ## Classes of given mean functions
#' mean_functions(mean_fun=c("beta", "hofIV"))
#'
#' @export
#' 
mean_functions <- function (mean_fun=NULL, mean_class=NULL) {
  ## Define mean functions

  ## --- Define mean functions
  MF <- c("constant", "uniform",
          "beta", "sech", "gaussian", "mixgaussian", "hofV",
          "sech_p1", "sech_r0p1", "mixgaussian_equal", 
          "hofII", "hofIV", "hofIVb", "hofVb")
  
  ## --- Define corresponding mean function classes
  MC <- c("test", "test",
          "main", "main", "main", "main", "main",
          "sub", "sub", "sub", "sub", "sub", "sub", "sub")
  
  ## --- Create data frame to store mean functions with corresponding mean class
  DF <- data.frame(mean_fun=MF, mean_class=MC)

  ## --- Grab appropriate output

  ## If mean function(s) is not null find matching mean class(es)
  if (!is.null(mean_fun)) {
    Match <- as.character(DF[which(!is.na(match(DF$mean_fun, mean_fun))),]$mean_class)
  }
  
  ## If mean class is not null find matching mean functions
  if (!is.null(mean_class)) {
    if (any(mean_class == "all")) {
      Match <- as.character(DF$mean_fun)
    } else {
      Match <- as.character(DF[which(!is.na(match(DF$mean_class, mean_class))),]$mean_fun)
    }
  }
  
  ## If mean function and mean class are both null return data frame
  if ( is.null(mean_fun) & is.null(mean_class) ) { Match <- DF }
  
  ## --- Return match
  return (Match)
}


#' Error distributions
#'
#' Display error distributions of given class.
#'
#' @param err_dist Vector of error distribution(s).
#' @param err_class String that indicates the which error distributions to list.
#' \describe{
#'   \item{"binary"}{0/1 data; ("bernoulli"). }
#'   \item{"binomial"}{Binomial data represented as counts or proportions;
#'       ("binomial_count", "binomial_prop"). }
#'   \item{"percentage"}{[0-1] data;
#'       ("tab", "zitab").}
#'   \item{"count"}{Count data;
#'       ("poisson", "negbin", "zip", "zinb", "zipl", "zinbl", "zipl.mu", "zinbl.mu").}
#'   \item{"abundance"}{Non-negative continous data;
#'       ("gaussian", "tweedie", "zig", "zigl", "zigl.mu").}
#' }
#' 
#' @return If no arguments supplied, return a data frame listing all possible error distributions
#' with corresponding error classes; if a vector of error classes is supplied, return
#' corresponding error distributions; if a vector of error distributions is supplied, return
#' corresponding error classes.
#' 
#' @keywords error distribution
#'
#' @examples
#'
#' ## Print data frame of all error distributions with corresponding error classes
#' error_distributions()
#'
#' ## Print all error distributions
#' error_distributions(err_class="all")
#'
#' ## Print error distributions for count data
#' error_distributions(err_class="count")
#'
#' ## Print error classes for poisson and gaussian data
#' error_distributions(err_dist=c("poisson", "gaussian"))
#' 
#' @export
#' 
error_distributions <- function (err_dist=NULL, err_class=NULL) {
  ## --- Match error distribution with given error class, or vice versa.
  ##     If no input is supplied, return error class, error distribution list
  
  ## --- Define error distributions
  ED <- c("bernoulli",  "binomial_count", "binomial_prop",
          "poisson",    "negbin",         "zip",           "zinb",
          "zipl",       "zinbl",          "zipl.mu",       "zinbl.mu",
          "gaussian",   "tweedie",        "zig",           "zigl",        "zigl.mu",
          "tab",        "zitab")

  ## --- Define corresponding error classes
  DC <- c("binary",     "binomial",        "binomial",
          "count",      "count",           "count",         "count",
          "count",      "count",           "count",         "count",
          "abundance",  "abundance",       "abundance",     "abundance",  "abundance",
          "percentage", "percentage")
  
  ## --- Create data frame to store error distributions with corresponding error class
  DF <- data.frame(err_dist=ED, err_class=DC)
  
  ## --- Grab appropriate output

  ## If error distribution(s) is not null find matching error class(es)
  if (!is.null(err_dist)) {
    Match <- as.character(DF[which(!is.na(match(DF$err_dist, err_dist))),]$err_class)
  }
  
  ## If error class is not null find matching error distribution
  if (!is.null(err_class)) {
    if (any(err_class == "all")) {
      Match <- as.character(DF$err_dist)
    } else {
      Match <- as.character(DF[which(!is.na(match(DF$err_class, err_class))),]$err_dist)
    }
  }
  
  ## If error distribution and error class are both null return data frame
  if ( is.null(err_dist) & is.null(err_class) ) { Match <- DF }
  
  ## --- Return match
  return (Match)
}


## --- --- --- ---


#' Set models
#'
#' Create data frame of models to fit.
#' Models will be formed by either pairing or crossing the elements of
#' the mean function and error distribution parameters vector.
#'
#' @param mean_fun List of mean functions to create data frame of models from.
#' @param err_dist List of error distributions to create data frame of models from.
#' @param mean_class Class of mean functions to create data frame of models from (optional).
#' @param err_class Class of error distributions to create data frame of models from (optional).
#' @param binomial_n Value of binomial n parameter for binomial models (optional).
#' @param method Method of combining the mean functions and error distributions;
#' "crossed" results in all possible mean function / error distribution combinations; "paired"
#' matches the ith mean functions with the ith error distribution.
#'
#' @return Data frame contaning the mean functions, error distribution, and binomial n parameter
#' of all created models.
#' 
#' @keywords model information
#'
#' @examples
#' 
#' ## Create all possible models with supplied mean functions and error distributions
#' set_models(mean_fun=c("beta","sech"), err_dist=c("zip","zinb"), method="crossed")
#'
#' ## Combine mean functions and error distributions directly
#' set_models(mean_fun=c("beta","sech"), err_dist=c("zip","zinb"), method="paired")
#'
#' ## Create a model with a binomial n parameter
#' set_models(mean_fun=c("beta","sech"), err_dist=c("zip","binomial_count"), binomial_n=40)
#'
#' @export
set_models   <- function (mean_fun=NULL, err_dist=NULL, mean_class=NULL, err_class=NULL,
                        binomial_n=NA, method="crossed") {
  ## --- Create data frame of models to fit.
  ## --- Models will be formed by either pairing or crossing the elements of
  ## --- the mean function and error distribution parameters vector.
  
  ## --- Set default values
  if (is.null(mean_fun) & is.null(mean_class)) { mean_class <- "main"  }
  if (is.null(err_dist) & is.null(err_class))  { err_class  <- "count" }
  
  ## --- Check if mean function and error distribution are defined correctly
  Check <- check_mean_error_inputs (mean_fun=mean_fun, err_dist=err_dist,
                                    mean_class=mean_class, err_class=err_class, method=method)
  
  ## --- Grab check model inputs
  method <- Check$method; mean_fun <- Check$mean_fun; err_dist <- Check$err_dist
  
  ## --- Select method to create data frame
  if (method == "paired") {
    ## --- Paired mean functions and error distribution
    Models <- set_models_paired (mean_fun=mean_fun, err_dist=err_dist, binomial_n=binomial_n)
  } else if (method == "crossed") {
    ## --- Crossed mean functions and error distribution
    Models <- set_models_crossed (mean_fun=mean_fun, err_dist=err_dist, binomial_n=binomial_n)
  } else {
    ## --- Error
    stop ("Method must be 'paired' or 'crossed'!")
  }
  
  ## --- Set class
  class(Models) <- append(class(Models), "model_list")
  
  ## --- Return models data frame
  return (Models)  
}


set_models_paired <- function (mean_fun=NULL, err_dist=NULL, binomial_n=NA) {
  ## --- Create data frame of models to fit.
  ## --- Models will be formed by pairing the ith elements of
  ## --- the mean function and error distribution parameters.
  
  ## --- Check binomial distribution
  BIndex <- check_binomial_n (binomial_n, err_dist) 
  
  ## --- Create models data frame
  DF <- data.frame (mean_fun=mean_fun, err_dist=err_dist, binomial_n=NA, stringsAsFactors=FALSE)

  ## --- Add binomial size parameter if given
  if (any(BIndex)) { DF[BIndex,]$binomial_n <- binomial_n }

  ## --- Sort data frame by mean vector
  DF <- DF[order(DF$mean_fun),]
  rownames(DF) <- 1:nrow(DF)
  
  ## Return models names
  return (DF)
}


set_models_crossed <- function (mean_fun=NULL, err_dist=NULL, binomial_n=NA) {
  ## --- Create data frame of models to fit.
  ## --- Models will be formed by crossing the elements of
  ## --- the mean function and error distribution parameters vector.
  
  ## --- Remove duplicates
  err_dist <- unique (err_dist)
  mean_fun <- unique (mean_fun)

  ## --- Check binomial distribution
  BIndex <- check_binomial_n (binomial_n, err_dist) 

  ## --- Combine mean functions and error distributions
  DF1 <- expand.grid (mean_fun=mean_fun, err_dist=err_dist, stringsAsFactors=FALSE)
  
  ## --- Create vector of length(mean_fun) with binomial n values
  ##     for binomial models and NA else where
  if (length(binomial_n) == 1) { binomial_n <- rep(binomial_n, sum(BIndex)) }
  if (length(binomial_n) == sum(BIndex)) {
    n <- binomial_n
    binomial_n <- rep(NA, length(err_dist))
    binomial_n[BIndex] <- n
    DF2 <- expand.grid (mean_fun=mean_fun, binomial_n=binomial_n)
  } else {
    stop ("** Length of binomial_n does not match number of binomial distributions!")
  }
  
  ## --- Combine mean functions, error distributions, and binomial n parameters
  DF <- cbind (DF1, binomial_n=DF2$binomial_n)

  ## --- Sort data frame by mean vector
  DF <- DF[order(DF$mean_fun),]
  rownames(DF) <- 1:nrow(DF)

  ## --- Return models names
  return (DF)
}


## --- --- --- ---


#' Set model info
#' 
#' Create model information data frame (one row). 
#'
#' @param mean_fun Mean function.
#' @param err_dist Error distribution.
#' @param binomial_n Binomial size parameter, n.
#' @param delta Tail-adjusted beta delta parameter.
#' @param thetaC Alternate specification of constant parameters, either
#' the binomial parameter, thetaC=c(n=40);
#' or the tail-adjusted beta delta parameter, thetaC=c(delta=0.01).
#' @param model A row from a set_models() data frame.  
#' @param data Data (x,y) set used to estimate tail-adjusted beta delta parameter.
#'
#' @return Data frame (one row) containing model information.
#' 
#' @keywords Model information
#'
#' @examples
#'
#' Model <- set_models (mean_fun=c("gaussian"), err_dist=c("zitab"), method="crossed")
#' set_model_info (model=Model)
#' set_model_info (mean_fun="gaussian", err_dist="negbin")
#' set_model_info (mean_fun="gaussian", err_dist="binomial_count", binomial_n=50)
#' set_model_info (mean_fun="gaussian", err_dist="tab", delta=0.01)
#'
#' @export
#' 
set_model_info <- function (model=NULL,  mean_fun=NULL, err_dist=NULL, 
                            binomial_n=NA, delta=NA, thetaC=NULL, data=NULL) {
  ## --- Set model info
  
  ## --- Use parameters specified in set_models() argument
  if (!is.null(model)) {
    ## Stop if model more than one model specified
    if (nrow(model) != 1) { stop ("Specify only one model!") }
    ## Stop if set_models() was not used to specify model
    if (!(any(class(model)=="model_list"))) { stop ("Use set_modes() to specify model!") }

    ## Grab parameters
    mean_fun   <- model$mean_fun
    err_dist   <- model$err_dist
    binomial_n <- model$binomial_n
  }
  
  ## --- Check mean functions and error distributions
  check_mean_functions_and_error_distributions (mean_fun, err_dist)
  
  ## --- Create data frame of model information
  ModelInfo <- list (model_name=NA, mean_fun=mean_fun, err_dist=err_dist, binomial_n=NA, delta=NA,
                     theta=NA, u.theta=NA, trans=NA, theta.lb=NA, theta.ub=NA)
  
  ## --- Add model constants: Tail-adjusted beta delta
  if ( (err_dist == "tab") | (err_dist == "zitab") ) {

    ## delta parameter specification priority: data > thetaC > delta
    ## Replace delta if supplied via delta
    if (!is.na(delta)) {
      ModelInfo$delta <- delta
    }
    
    ## Replace delta if supplied via thetaC
    if (!is.null(thetaC)) {
      if (!is.na(thetaC["delta"])) { ModelInfo$delta <- thetaC["delta"] }
    }
    
    ## If data is supplied, estimate delta
    if (!is.null(data)) {
      ModelInfo$delta <- estimate_delta (y=data$y)
    }
  }
  
  ## --- Add model constants: Binomial n
  ## binomial_n parameter specification priority: thetaC > binomial_n
  if ((err_dist=="binomial_count") | (err_dist=="binomial_prop")) {
    ## Use supplied argument
    ModelInfo$binomial_n <- binomial_n
    ## Replace binomial_n if supplied via thetaC
    if (!is.null(thetaC)) {
      if (!is.na(thetaC["binomial_n"])) { ModelInfo$binomial_n <- thetaC["binomial_n"] }
    }
  }
  
  ## --- Create model name
  ModelInfo$model_name <- paste (ModelInfo$mean_fun, ModelInfo$err_dist, sep="-")
  
  ## --- Grab parameter information
  parinfo  <- get_parinfo (ModelInfo)  
  ModelInfo$theta    <- parinfo$theta
  ModelInfo$u.theta  <- parinfo$u.theta
  ModelInfo$trans    <- parinfo$trans
  ModelInfo$theta.lb <- parinfo$lb
  ModelInfo$theta.ub <- parinfo$ub

  ## --- Grab parameter names 
  parnames <- get_parnames (ModelInfo$model_name)
  ModelInfo$thetaM   <- parnames$thetaM
  ModelInfo$thetaE   <- parnames$thetaE
  ModelInfo$thetaC   <- parnames$thetaC
  
  ## --- Return object
  return (ModelInfo)
}


estimate_constant <- function (ModelInfo, data=NULL) {
  ## --- Estimate either the binomial n parameter or the delta parameter
  
  ## --- Grab model info
  err_dist   <- ModelInfo$err_dist
  mean_fun   <- ModelInfo$mean_fun
  binomial_n <- ModelInfo$binomial_n
  delta      <- ModelInfo$delta

  ## --- Estimate delta if data given
  if (!is.null(data) & ( (err_dist=="tab") | (err_dist=="zitab") ) ) {
    ## --- Estimate delta
    delta <- estimate_delta (y=data$y)
  }

  ## --- Create constant vector
  names(binomial_n) <- NULL; names(delta) <- NULL
  thetaC <- c(binomial_n=binomial_n, delta=delta)
  
  ## --- Remove missing variables 
  if (all(is.na(thetaC))) {
    thetaC <- NULL
  } else {
    thetaC <- thetaC[!is.na(thetaC)]
  }
  
  ## --- Return constant
  return (thetaC)
}


estimate_delta <- function (y) {
  ## --- Estimate delta

  ## --- Set delta
  ## --- Smallest positive distances to zero
  delta0 <- min (y[y>0])
  ## --- Smallest positive distances to one
  delta1 <- min ((1-y)[y<1])
  ## --- Set delta (0.5 if data all 0,1)
  delta  <- min ( c(delta0, delta1, 0.5))

  ## --- Return delta
  return (delta)
}


#' Model parameters
#'
#' Get parameter names for given model.
#'
#' @param model_name Model_name in the form of "mean_fun-err_dist".
#' 
#' @return List containing model names, character vector of mean function parameters,
#' error distribution parameters, and constant parameters.
#' 
#' @keywords model parameters
#'
#' @examples
#'
#' ## Get parameters for a "gaussian-zinb" model
#' get_parnames ("gaussian-zinb")
#'
#' @export
#'
get_parnames <- function (model_name) {
  ## --- Define parameters for each model

  ## --- Grab error distribution, mean function, and model name
  ## --- model_name the "mean_fun-err_dist" model name
  Names      <- unlist (strsplit(model_name,"-"))
  mean_fun   <- Names[1]
  err_dist   <- Names[2]
  
  ## --- Get parameter names
  thetaM <- get_mean_fun_parnames (mean_fun)
  thetaE <- get_err_dist_parnames (err_dist)
  thetaC <- get_constant_parnames (err_dist)

  ## --- Store parameter names
  ParNames <- list ()
  ParNames$model_name <- model_name
  ParNames$thetaM <- thetaM
  ParNames$thetaE <- thetaE
  ParNames$thetaC <- thetaC

  ## MLE parameters
  ParNames$theta <- c(thetaM, thetaE)

  ## --- Return parameter names
  return (ParNames)
}


#' Mean function parameters names
#'
#' Get parameter names for given mean function.
#'
#' @param mean_fun Mean function name.
#' 
#' @return Vector containing parameter names of mean function.
#' 
#' @keywords mean function parameters
#'
#' @examples
#'
#' ## Get parameters for a "gaussian" mean function
#' get_mean_fun_parnames ("gaussian")
#'
#' @export
#'
get_mean_fun_parnames <- function (mean_fun) {
  ## --- Get parameter names for mean function

  ## --- Check to make sure mean function is a character vector of length one
  if ( !((class(mean_fun) == "character") & (length(mean_fun) == 1)) ) {
    stop ("mean function must be a character vector of length one!")
  }

  ## --- Set default mean function parameter names to null
  thetaM <- NULL
  
  ## --- Define parameters: mean functions
  if (mean_fun == "constant")          { thetaM <- c("H") }
  if (mean_fun == "uniform")           { thetaM <- c("H","c","d") }
  if (mean_fun == "gaussian")          { thetaM <- c("H","m","s") }
  if (mean_fun == "mixgaussian")       { thetaM <- c("H","a","m1","m2","s1","s2") }
  if (mean_fun == "mixgaussian_equal") { thetaM <- c("H","a","m1","m2","s") }
  if (mean_fun == "beta")              { thetaM <- c("H","c","d","u","v") }
  if (mean_fun == "sech")              { thetaM <- c("H","m","s","r","p") }
  if (mean_fun == "sech_p1")           { thetaM <- c("H","m","s","r") }
  if (mean_fun == "sech_r0p1")         { thetaM <- c("H","m","s") }
  if (mean_fun == "hofII")             { thetaM <- c("H","m","w0") }
  if (mean_fun == "hofIV")             { thetaM <- c("H","m","w","k") }
  if (mean_fun == "hofIVb")            { thetaM <- c("H","m","w") }
  if (mean_fun == "hofV")              { thetaM <- c("H","m","w1","w2","k") }
  if (mean_fun == "hofVb")             { thetaM <- c("H","m","w1","w2") }

  ## --- Return mean function parameter names
  return (thetaM)
}


#' Error distribution parameters names
#'
#' Get parameter names for given error distribution.
#'
#' @param err_dist Error distribution name.
#' 
#' @return Vector containing parameter names of error distribution.
#' 
#' @keywords error distribution parameters
#'
#' @examples
#'
#' ## Get parameters for a "zinb" er
#' get_err_dist_parnames ("zinb")
#'
#' @export
#'
get_err_dist_parnames <- function (err_dist) {
  ## --- Get parameter names for mean function

  ## --- Check to make sure error distribution is a character vector of length one
  if ( !((class(err_dist) == "character") & (length(err_dist) == 1)) ) {
    stop ("err_dist must be a character vector of length one!")
  }

  ## --- Set default error distribution parameter names to null
  thetaE <- NULL

  ## --- Define parameters: error distributions
  if (err_dist == "bernoulli")      { thetaE <- NULL }
  if (err_dist == "binomial_count") { thetaE <- NULL }
  if (err_dist == "binomial_prop")  { thetaE <- NULL }
  if (err_dist == "poisson")        { thetaE <- NULL }
  if (err_dist == "negbin")         { thetaE <- c("phi") }
  if (err_dist == "zip")            { thetaE <- c("pi") }
  if (err_dist == "zinb")           { thetaE <- c("pi","phi") }
  if (err_dist == "zipl")           { thetaE <- c("g0","g1") }
  if (err_dist == "zipl.mu")        { thetaE <- c("g0","g1") }
  if (err_dist == "zinbl")          { thetaE <- c("g0","g1","phi") }
  if (err_dist == "zinbl.mu")       { thetaE <- c("g0","g1","phi") }
  if (err_dist == "tweedie")        { thetaE <- c("phi","rho") }
  if (err_dist == "zig")            { thetaE <- c("pi","phi") }
  if (err_dist == "zigl")           { thetaE <- c("g0", "g1","phi") }
  if (err_dist == "zigl.mu")        { thetaE <- c("g0", "g1","phi") }
  if (err_dist == "gaussian")       { thetaE <- c("sigma") }
  if (err_dist == "tab")            { thetaE <- c("sigma") }
  if (err_dist == "zitab")          { thetaE <- c("pi","sigma") }
 
  ## --- Return error distribution parameter names
  return (thetaE)
}


#' Constant parameters names
#'
#' Get constant parameter names for given error distribution.
#'
#' @param err_dist Error distribution name.
#' 
#' @return Vector containing names of constant parameters.
#' 
#' @keywords constant parameters
#'
#' @examples
#'
#' ## Get parameters for a "zinb" er
#' get_constant_parnames ("zinb")
#'
#' @export
#'
get_constant_parnames <- function (err_dist) {
  ## --- Get parameter names for constants
  
  ## --- Check to make sure error distribution is a character vector of length one
  if ( !((class(err_dist) == "character") & (length(err_dist) == 1)) ) {
    stop ("err_dist must be a character vector of length one!")
  }
  
  ## --- Set default error distribution parameter names to null
  thetaC <- NULL
  
  ## --- Define constant parameters: error distributions
  if (err_dist == "binomial_count") { thetaC <- "binomial_n" }
  if (err_dist == "binomial_prop")  { thetaC <- "binomial_n" }
  if (err_dist == "tab")            { thetaC <- "delta" }
  if (err_dist == "zitab")          { thetaC <- "delta" }
  
  ## --- Return constant parameter names
  return (thetaC)
}


#' List all model parameters
#' 
#' List the names of all possible model parameters.
#'
#' @return List of all parameter names.
#' 
#' @keywords parameter names
#'
#' @examples
#' 
#' list_all_parnames()
#' 
#' @export
#' 
list_all_parnames <- function () {
  ## --- Produce a character vector listing all possible parameters
  
  ## --- Grab all error distributions and mean functions
  M <- print_models (mean_class="all", err_class="all")
  
  ## --- Create list of all possible models
  Models <- set_models (err_dist=M$error_distributions, mean_fun=M$mean_functions,
                        binomial_n=1, method="crossed")
  
  ## --- Initialise vectors containing parameter names
  ##     for error distributions, mean functions, and constants
  thetaE <- NULL; thetaM <- NULL; thetaC <- NULL
  
  ## --- Loop through models
  for (i in 1:nrow(Models)) {
    ## --- Set model name
    model_name <- paste (Models$mean_fun, Models$err_dist, sep="-")
    ## --- Get parameter names for model
    Par <- get_parnames (model_name)
    ## --- Store models parameters
    thetaE <- c(thetaE, Par$thetaE)
    thetaM <- c(thetaM, Par$thetaM)
    thetaC <- c(thetaC, Par$thetaC)
  }
  
  ## --- Remove duplicates
  thetaE <- unique(sort(thetaE))
  thetaM <- unique(sort(thetaM))
  thetaC <- unique(sort(thetaC))
  
  ## --- Combine parameters in one list
  theta <- c(thetaM, thetaE, thetaC)
  
  ## --- Return parameter list
  return (theta)
}


check_mean_error_inputs <- function (mean_fun=NULL,   err_dist=NULL,
                                     mean_class=NULL, err_class=NULL, method=NULL) {
  ## --- Check if mean functions and error distribution
  ## --- specifications for set_models() are legal.

  ## --- Stop if method not paired or crossed
  if (!((method == "paired") | (method=="crossed"))) { 
    stop ("Method must be 'paired' or 'crossed'!")
  }
  
  ## --- If mean function or error distribution is a singleton, use crossed method
  if ( (length(mean_fun)==1) | (length(err_dist)==1) ) { method <- "crossed" }
  
  ## --- Paired
  if (method == "paired") {
    ## --- Stop if error distribution or mean function are nulll
    if (is.null(err_dist) | is.null(mean_fun)) {
      stop ("Must specify mean_fun **and** err_dist if method=='paired'!")
    }
    ## --- Stop if supplied mean functions and error distributions are of different lengths
    if (length(mean_fun) != length(err_dist)) {
      stop ("Length of mean_fun and err_dist must match if method=='paired'!")
    }
  }
  
  ## --- Crossed
  if (method == "crossed") {
    ## --- Stop if error distribution and error class are both specified
    if (!is.null(err_dist) & !is.null(err_class)) {
      stop ("Do not specify err_dist **and** err_class!")
    }
    ## --- Stop if mean function and mean class are both specified
    if (!is.null(mean_fun) & !is.null(mean_class)) {
      stop ("Do not specify mean_fun **and** mean_class!")
    }
    
    ## --- If error distribution and error class not specifed set error class to "all"
    if ( is.null(err_dist) & is.null(err_class) ) { err_class <- "all" }
    ## --- Set distribution if error class supplied
    if (!is.null(err_class)) { 
      err_dist <- error_distributions (err_class=err_class)
    }
    
    ## --- If mean function and mean class not specifed set mean class to "all"
    if ( is.null(mean_fun) & is.null(mean_class) ) { mean_class <- "all" }
    ## --- Set mean function if mean_class supplied
    if (!is.null(mean_class)) { 
      mean_fun <- mean_functions (mean_class=mean_class)
    }    
    ## --- Set mean functions to all if not set  
    if (is.null(mean_fun)) { mean_fun <- mean_functions(mean_class="all") }
  }

  ## --- Check if mean functions and error distributions are legal names
  check_mean_functions_and_error_distributions (mean_fun, err_dist)

  ## --- Store method, mean functions, and error distributions
  Check <- list()
  Check$method <- method
  Check$mean_fun <- mean_fun
  Check$err_dist <- err_dist

  ## --- Return method
  return (Check)
}


check_binomial_n <- function (binomial_n, err_dist) {
  ## --- Check if the binomial n parameter is specied for binomial models

  ## --- Binomial must be a univariate integer/numeric
  if (length(binomial_n) != 1) {
    stop ("binomial_n must be **univariate** integer/numeric!")
  }
  if (!is.na(binomial_n)) {
    if (binomial_n != round(binomial_n)) {
      stop ("binomial_n must be univariate **integer/numeric**!")
    }
  }
  
  ## --- Binomial distribution
  BIndex <- (err_dist=="binomial_count") | (err_dist=="binomial_prop")
  ## --- Stop if n parameter is not given for binomial distributions
  if (any(BIndex) & is.na(binomial_n)) {
    stop ("binomial_n is not specified for binomial distributions!")
  }

  ## --- Return binomial index variable
  return (BIndex)
}


check_mean_functions_and_error_distributions <- function (mean_fun, err_dist) {
  ## --- Check mean functions and error distributions

  ## --- Check mean functions 
  Fail <- check_mean_functions (mean_fun)
  if (Fail) { stop ("Illegal mean function supplied!") }

  ## --- Check error distributions
  Fail <- check_error_distributions (err_dist)
  if (Fail) { stop ("Illegal error distribution supplied!") }
}


check_mean_functions <- function (mean_fun) {
  ## --- Check if mean functions are legal
  
  ## --- List all possible mean functions
  MF <- mean_functions (mean_class="all")
  
  ## --- Check if any mean functions are illegal
  Fail <- any (is.na(match(mean_fun, MF)))

  ## --- Retun fail flag
  return (Fail)
}


check_error_distributions <- function (err_dist) {
  ## --- Check if error distributions are legal
  
  ## --- List all possible error distributions
  ED <- error_distributions(err_class="all")
  
  ## --- Check if any mean functions are illegal
  Fail <- any (is.na(match(err_dist, ED)))

  ## --- Retun fail flag
  return (Fail)
}

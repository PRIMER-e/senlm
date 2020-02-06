## --- Mean functions


#' Evaluate mean function
#'
#' Evalautes the mean function defined in ModelInfo, with the
#' parameters given in theta, at the points given by the vector x.
#'
#' @param ModelInfo A row from the output of set_models(), containing the
#' mean_function and error distribution of the model in question.
#' @param theta Vector containing parameter values of model.
#' @param x Values at which to calculate the mean function.
#'
#' @return Calculates the constant mean functions defined by thetaM at points x.
#' 
#' @keywords mean function, constant
#'
#' @examples
#'
#' Models <- set_models (mean_fun=c("gaussian"), err_dist=c("negbin"))
#' mu_meanfunction (ModelInfo=Models[1,], theta=c(H=60, m=50, s=10, phi=1.5), x=0:100)
#'
#' @export
mu_meanfunction <- function (ModelInfo, theta, x) {
  ## --- Calculate mean function at values x
  
  ## --- Create ModelInfo object if character string "err_dist-mean_fun" is supplied
  if (class(ModelInfo)=="character") { ModelInfo <- set_model_info (ModelInfo) }
  
  ## --- Set model name
  model_name <- paste (ModelInfo$mean_fun, ModelInfo$err_dist, sep="-")
  
  ## --- Get par names
  ParNames <- get_parnames (model_name)

  ## --- Grab mean function parameters
  thetaM <- theta[ParNames$thetaM]
  
  ## --- Get mean function name
  MuName <- ModelInfo$mean_fun

  ## --- Calculate mean function
  Mu <- do.call (paste("mu_",MuName, sep=""), list(thetaM, x))

  ## --- Return mean function
  return (Mu)
}


#' Constant mean function.
#'
#' @param thetaM Vector containing parameter values of mean function.
#' @param x Values at which to calculate the mean function.
#'
#' @return Calculates the constant mean functions defined by thetaM at points x.
#' 
#' @keywords mean function, constant
#'
#' @examples
#' 
#' ## Constant and uniform mean functions to test code
#' thetaM <- c(H=40)
#' mu_constant(thetaM, x=0:100)
#'
#' @export
mu_constant <- function (thetaM, x) {
  ## --- Constant mean function

  ## Grab parameters
  H <- as.list(thetaM)$H
  
  ## Set mean function
  mu <- rep (H, length(x))

  ## Return mu
  return (mu)
}


## --- Uniform (constant between c and d)

mu_uniform <- function (thetaM, x) {
  ## --- Uniform mean function

  ## Grab parameters
  H <- as.list(thetaM)$H
  c <- as.list(thetaM)$c
  d <- as.list(thetaM)$d

  ## Set mean function
  mu        <- rep (H, length(x))
  mu[x < c] <- 0
  mu[x > d] <- 0
  
  ## Return mu
  return (mu)
}


## --- Gaussian

mu_gaussian <- function (thetaM, x) {
  ## --- Gaussian mean function

  ## Grab parameters
  H <- as.list(thetaM)$H
  m <- as.list(thetaM)$m
  s <- as.list(thetaM)$s

  ## Standardised x
  z <- (x - m)/s
  
  ## Set mu
  mu <- H * exp (-0.5*z^2)

  ## Return mu
  return (mu)
}


## --- Mixture of two Gaussians: Common standard deviation

mu_mixgaussian_equal <- function (thetaM, x) {
  ## --- Gaussian mean function

  ## Grab parameters
  H  <- thetaM['H'];  a  <- thetaM['a'];
  m1 <- thetaM['m1']; m2 <- thetaM['m2']
  s  <- thetaM['s'];

  ## Standardised x
  z1 <- (x - m1)/s
  z2 <- (x - m2)/s

  ## Find max of unscaled function (plus 
  z <- seq(m1,m2, length=10000)
  Tolerance <- 0.00001
  H_mode <- max (a * exp(-0.5*((z-m1)/s)^2) + (1-a)*exp(-0.5*((z-m2)/s)^2)) + Tolerance
  
  ## Set mu
  mu <- (H/H_mode) * ( a * exp (-0.5*z1^2) + (1-a) * exp (-0.5*z2^2) )
  
  ## Return mu
  return (mu)
}


## --- Mixture of two Gaussians : Different standard deviations

mu_mixgaussian <- function (thetaM, x) {
  ## --- Gaussian mean function

  ## Grab parameters
  H  <- thetaM['H'];  a  <- thetaM['a'];
  m1 <- thetaM['m1']; m2 <- thetaM['m2']
  s1 <- thetaM['s1']; s2 <- thetaM['s2']
  
  ## Standardised x
  z1 <- (x - m1)/s1
  z2 <- (x - m2)/s2

  ## Find max of unscaled function (plus 
  z <- seq(m1,m2, length=10000)
  Tolerance <- 0.00001
  H_mode <- max (a * exp(-0.5*((z-m1)/s1)^2) + (1-a)*exp(-0.5*((z-m2)/s2)^2)) + Tolerance
  
  ## Set mu
  mu <- (H/H_mode) * ( a * exp (-0.5*z1^2) + (1-a) * exp (-0.5*z2^2) )
  
  ## Return mu
  return (mu)
}


## --- Modified Beta

mu_beta <- function (thetaM, x) {
  ## --- Beta mean function
  
  ## Grab parameters
  H <- as.list(thetaM)$H
  c <- as.list(thetaM)$c
  d <- as.list(thetaM)$d
  u <- as.list(thetaM)$u
  v <- as.list(thetaM)$v
  
  ## Convert beta mean and shape parameters
  a <- u / v ;     names(a) <- "a"
  b <- (1 - u)/v ; names(b) <- "b"
  
  ## Find region where mean is positive
  I <- (x>c) & (x<d)
  
  ## Set x values standardised by endpoints
  z0 <- (x[I] - c)/(d-c)
  z1 <- (d - x[I])/(d-c)
  
  ## Set maximum H
  if ( ((a==1) & (b==1)) | ((a==1) & (b>1)) | ((a>1) & (b==1)) ) {
    ## Uniform (a=b=1), or triangular (a=1,b>1 or a>1,b=1)
    log_H_mode <- 0
  } else if ( (a>1) & (b>1) ) {
    ## Unimodal
    x_mode  <- (d-c)*(a-1)/(a+b-2) + c
    log_H_mode  <- (a-1)*log((x_mode - c)/(d-c)) + (b-1)*log((d - x_mode)/(d-c))
  } else {
    ## Not uniform, triangular, or unimodal
    log_H_mode <- lbeta (a,b)
  }
  
  ## Log mean
  log.z0  <- log(z0); log.z0[z0==0] <- 0
  log.z1  <- log(z1); log.z1[z1==0] <- 0
  logmean <- log (H) - log_H_mode + (a-1)*log.z0 + (b-1)*log.z1
  
  ## Set mean function
  mu    <- rep (0, length(x))
  mu[I] <- exp(logmean)
  
  ## Return mu
  return (mu)
}


## --- Sech function

mu_sech <- function (thetaM, x) {
  ## --- Sech mean function

  ## --- Set sech function
  logsech <- function (x, s, p) {
    ## Standardise x
    X <- x/s
    ## Calculate log ( sech(x/s)^p )
    Y <- p*( log( 2 ) - ( abs(X) + log( 1 + exp(-2*abs(X))) ) )
    ## Return log sech
    return (Y)
  }
  
  ## Grab parameters
  H <- thetaM['H']; m <- thetaM['m']; s <- thetaM['s']; r <- thetaM['r']; p <- thetaM['p']
  
  ## Set x_mode
  x_mode <- s * 0.5 * log ( (1+r)/(1-r) )

  ## Set H_mode (log scale)
  logH_mode <- x_mode*r*p/s  + 0.5*p*log(1-r^2)
  
  ## Set x value
  X <- (x  - (m - x_mode))
  
  ## Log mean
  logmean <- log(H) - logH_mode + ((r*p/s)*X) + logsech(X,s,p)
  ## Set mean function
  mu    <- exp(logmean)
  
  ## Return mu
  return (mu)
}


## --- Sech (p=1) function

mu_sech_p1 <- function (thetaM, x) {
  ## --- Sech mean function
  ## Peakedness parameter set to 1 (default)

  ## Add p=1 to thetaM
  Theta   <- as.list (thetaM)
  Theta$p <- 1
  thetaM  <- unlist (Theta)
  
  ## Return mu
  mu_sech (thetaM, x)
}


## --- Sech (r=0, p=1) function

mu_sech_r0p1 <- function (thetaM, x) {
  ## --- Sech mean function
  ## Peakedness parameter set to 1 (default)
  ## Skewness parameter set to 0 (symmetric)

  ## Add r=0, p=1 to thetaM
  Theta   <- as.list (thetaM)
  Theta$r <- 0
  Theta$p <- 1
  thetaM  <- unlist (Theta)
  
  ## Return mu
  mu_sech (thetaM, x)
}


## --- HOF II function
 
mu_hofII <- function (thetaM, x) {
  ## --- HOF II function
  ## mu = H * (1/( 1 + exp(-w0*(x-m)) ))
  
  ## Grab parameters
  H  <- as.list(thetaM)$H
  m  <- as.list(thetaM)$m
  w0 <- as.list(thetaM)$w0
  
  ## Standardised variables
  z1 <- w0*(x - m)
  
  ## Calculate HOF
  Hof <- 1/(1 + exp(-z1))

  ## H mode
  H_mode <- 1

  ## Calculate mean of HOF function
  mu <- (H/H_mode) * Hof
  
  ## Return mu
  return (mu)
}


## --- HOF IV function

mu_hofIV  <- function (thetaM, x) {
  ## --- HOF IV function
  ## mu = (H/H_mode) * (1/( 1 + exp(-w(x-(m-k))) )) * (1/( 1 + exp(w(x-(m+k)))))
  
  ## Grab parameters
  H <- as.list(thetaM)$H
  m <- as.list(thetaM)$m
  w <- as.list(thetaM)$w
  k <- as.list(thetaM)$k
  
  ## Standardised variables
  z1 <- w*(x - (m - k))
  z2 <- w*(x - (m + k))
  
  ## Calculate HOF
  Hof <- ( 1/(1 + exp(-z1)) ) * ( 1/(1 + exp(z2)) )
  
  ## H_mode
  H_mode <- (1 + exp(-w*k))^(-2)

  ## Calculate mean of HOF function
  mu <- (H/H_mode) * Hof
  
  ## Return mu
  return (mu)
}


## --- HOF IVb function

mu_hofIVb  <- function (thetaM, x) {
  ## --- HOF IV function
  ## mu = (H/H_mode) * (1/( 1 + exp(w(x-m)) )) * (1/( 1 + exp(-w(x-m))))
  
  ## Grab parameters
  H <- as.list(thetaM)$H
  m <- as.list(thetaM)$m
  w <- as.list(thetaM)$w
  
  ## Standardised variables
  z1 <- w*(x - m)
  
  ## Calculate HOF
  Hof <- ( 1/(1 + exp(-z1)) ) * ( 1/(1 + exp(z1)) )

  ## H mode
  H_mode <- 0.25
  
  ## Calculate mean of HOF function
  mu <- (H/H_mode) * Hof
  
  ## Return mu
  return (mu)
}


## --- HOF V function

mu_hofV  <- function (thetaM, x) {
  ## --- HOF V function
  ## mu = (H/H_mode) * (1/( 1 + exp(-w1(x-(m1-k))) )) * (1/( 1 + exp(w2(x-(m+k)))))
  
  ## Grab parameters
  H  <- as.list(thetaM)$H
  m  <- as.list(thetaM)$m
  w1 <- as.list(thetaM)$w1
  w2 <- as.list(thetaM)$w2
  k  <- as.list(thetaM)$k
  
  ## Standardised variables
  z1 <- w1*(x - (m - k))
  z2 <- w2*(x - (m + k))
  
  ## Calculate HOF
  Hof <- ( 1/(1 + exp(-z1)) ) * ( 1/(1 + exp(z2)) )

  ## --- Find maximum H
  ## Lower limit of mode location
  x0 <- (m-k) - 2/(w1 + 0.01)
  ## Upper limit of mode location
  x1 <- (m+k) + 2/(w2 + 0.01)
  ## Range of x values
  xx <- c(seq(x0,x1, length.out=length(x)), m)
  ## Create component curves
  zz1 <- w1*(xx - (m - k))
  zz2 <- w2*(xx - (m + k))
  ## Calculate HOF
  Hof2 <- ( 1/(1 + exp(-zz1)) ) * ( 1/(1 + exp(zz2)) )
  ## H mode
  H_mode=max(Hof2)
  
  ## Calculate mean of HOF function
  mu <- (H/H_mode) * Hof
  
  ## Return mu
  return (mu)
}
 

## --- HOF Vb function

mu_hofVb  <- function (thetaM, x) {
  ## --- HOF V function
  ## mu = (H/H_mode) * (1/( 1 + exp(-w1(x-m)) )) * (1/( 1 + exp(w2(x-m))))
  
  ## Grab parameters
  H  <- as.list(thetaM)$H
  m  <- as.list(thetaM)$m
  w1 <- as.list(thetaM)$w1
  w2 <- as.list(thetaM)$w2
  
  ## Standardised variables
  z1 <- w1*(x - m)
  z2 <- w2*(x - m)
  
  ## Calculate HOF
  Hof <- ( 1/(1 + exp(-z1)) ) * ( 1/(1 + exp(z2)) )

  ## --- Find maximum H
  ## Lower limit of mode location
  x0 <- m - 2/(w1 + 0.01)
  ## Upper limit of mode location
  x1 <- m + 2/(w2 + 0.01)
  ## Range of x values
  xx <- c(seq(x0,x1, length.out=length(x)), m)
  ## Create component curves
  zz1 <- w1*(xx - m)
  zz2 <- w2*(xx - m)
  ## Calculate HOF
  Hof2 <- ( 1/(1 + exp(-zz1)) ) * ( 1/(1 + exp(zz2)) )
  ## H mode
  H_mode=max(Hof2)
  
  ## Calculate mean of HOF function
  mu <- (H/H_mode) * Hof

  ## Return mu
  return (mu)
}

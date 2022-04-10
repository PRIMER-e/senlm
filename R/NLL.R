## --- Negative log-likelihoods


## --- --- ---

set_nll <- function (ModelInfo) {
  ## --- Set negative log-likelihood
  
  ## --- Get error class
  err_dist  <- ModelInfo$err_dist
  err_class <- error_distributions (err_dist=err_dist)
  
  ## --- Continuous negative log-likelihood
  if ( (err_class == "abundance") | (err_class == "percentage") ) {
    NLL <- NLL.continuous
  }
  
  ## --- Discrete negative log-likelihood
  if ( (err_class == "binary") | (err_class == "binomial") | (err_class == "count") ) {
    NLL <- NLL.discrete
  }
  
  ## Return function
  return (NLL)
}

nll_wrapper <- function (ModelInfo, Dat) {
  ## --- Negative log-likelihood wrapper function for models that are fit by bblme::mle2

  ## --- Get error class
  err_dist  <- ModelInfo$err_dist
  err_class <- error_distributions (err_dist=err_dist)
  
  ## --- Continuous negative log-likelihood
  if ( (err_class == "abundance") | (err_class == "percentage") ) {
    NLL <- function (u.theta) { return (NLL.continuous (u.theta, ModelInfo, Dat)) }
  }
  
  ## --- Discrete negative log-likelihood
  if ( (err_class == "binary") | (err_class == "binomial") | (err_class == "count") ) {
    NLL <- function (u.theta) { return (NLL.discrete (u.theta, ModelInfo, Dat)) }
  }
  
  ## --- Return negative log-likelihood function
  return (NLL)
}


get_parinfo  <- function (ModelInfo) {
  ## --- Get list of transformations for model parameters
  
  ## --- Get error distribution
  err_dist   <- ModelInfo$err_dist
  
  ## --- Get parameter names
  theta <- get_parnames (ModelInfo$model_name)$theta
  
  ## --- Set parameter transformations
  trans <- rep("raw", length(theta))
  
  for (i in 1:length(theta)) {

    ## --- Group parameters by transformation
    LogPars    <- c("phi","k","s","s1","s2","w","w1","w2","v","p","sigma","b") ##    > 0
    LogitPars  <- c("pi","u","a", "q")                                         ##  0 - 1 **q=r modskurt 
    cLogitPars <- c("r")                                                       ## -1 - 1
    rLogitPars <- c("rho")                                                     ##  1 - 2
    nLogitPars <- c()                                                          ##  0 - n
    
    ## --- Transformation for H depends on model
    if ( (err_dist == "bernoulli") | (err_dist == "binomial.prop") |
         (err_dist == "tab")       | (err_dist == "zitab")         ) {
      ## 0-1 : logit
      LogitPars <- c("H", LogitPars)
    } else if (err_dist == "binomial.count") {
      ## 0-n : n logit
      nLogitPars <- c("H")
    } else {
      ## log
      LogPars <- c("H", LogPars)
    }
    
    ## --- Determine transformation for each parameter
    if (!is.na(match(theta[i],  LogPars)))   { trans[i] <- "log"    }
    if (!is.na(match(theta[i],  LogitPars))) { trans[i] <- "logit"  }
    if (!is.na(match(theta[i], cLogitPars))) { trans[i] <- "clogit" }
    if (!is.na(match(theta[i], rLogitPars))) { trans[i] <- "rlogit" }
    if (!is.na(match(theta[i], nLogitPars))) { trans[i] <- "nlogit" }
  }

  ## --- Set transformed parameter name
  u.theta <- paste(trans, theta, sep=".")

  ## --- Set lower and upper bounds
  lb <- rep (-Inf, length(theta))  
  ub <- rep ( Inf, length(theta))
  lb[match(LogPars, theta)]    <-  0; ub[match(LogPars, theta)]    <- Inf
  lb[match(LogitPars, theta)]  <-  0; ub[match(LogitPars, theta)]  <- 1
  lb[match(cLogitPars, theta)] <- -1; ub[match(cLogitPars, theta)] <- 1
  lb[match(rLogitPars, theta)] <-  1; ub[match(rLogitPars, theta)] <- 2
  lb[match(nLogitPars, theta)] <-  0; ub[match(nLogitPars, theta)] <- ModelInfo$binomial_n  

  ## --- Store theta, transformations, and bounds into data frame
  pars <- data.frame (theta=theta, u.theta=u.theta, trans=trans,
                      lb=lb, ub=ub, stringsAsFactors=FALSE)
  
  ## --- Return
  return (pars)
}


make_unbounded <- function (ModelInfo, theta) {
  ## --- Convert parameters so they are unbounded

  ## Transform parameters
  ## Zero :  pi
  ## Dispersion : phi
  ## Tweedie : rho
  ## Height : H
  ## Mixture : a
  ## Kurtosis : k
  ## Standard deviation, Spread : s, s1, s2
  ## Precision (or root precision) : w, w1, w2
  ## Beta : u, v
  ## Sech / Sech-p1 - Skewness : r
  ## ModSkurt / Sech - Peakedness : p
  ## ModSkurt - Flatness : b
  ## ModSkurt - Skewness (Symmetry q=0.5): q (r) 
  ## Gaussian : sigma
  
  ## --- Split model name into error distribution and mean function
  err_dist   <- ModelInfo$err_dist
  mean_fun   <- ModelInfo$mean_fun
  binomial_n <- ModelInfo$binomial_n
  delta      <- ModelInfo$delta
  
  ## --- Make sure theta has parameter names
  names(theta) <- ModelInfo$theta
  
  ## --- Get transformations
  trans <- ModelInfo$trans
  
  ## --- Create unbounded parameters
  u.theta <- rep (NA, length(theta))
  names(u.theta) <- ModelInfo$u.theta

  ## --- Grab lower and upper bounds
  lb <- ModelInfo$theta.lb
  ub <- ModelInfo$theta.ub
  
  ## --- Transform parameters to unbounded
  for (i in 1:length(u.theta)) {
    if (trans[i] == "raw") { u.theta[i] <- theta[i]      }
    if (trans[i] == "log") { u.theta[i] <- log(theta[i]) }
    if ((trans[i] == "logit")  | (trans[i] == "clogit") |
        (trans[i] == "rlogit") | (trans[i] == "nlogit") ) {
      u.theta[i] <- log((theta[i] - lb[i])/(ub[i] - theta[i]))
    }
  }
  
  ## --- Return unbounded parameters
  return (u.theta)
}


make_bounded <- function (ModelInfo, u.theta) {
  ## --- Convert parameters so they are unbounded
  
  ## --- Split model name into error distribution and mean function
  err_dist   <- ModelInfo$err_dist
  mean_fun   <- ModelInfo$mean_fun
  binomial_n <- ModelInfo$binomial_n
  delta      <- ModelInfo$delta
  
  ## --- Make sure theta has parameter names
  names(u.theta) <- ModelInfo$u.theta
  
  ## --- Get transformations
  trans <- ModelInfo$trans
  
  ## --- Create unbounded parameters
  theta <- rep (NA, length(u.theta))
  names(theta) <- ModelInfo$theta
  
  ## --- Grab lower and upper bounds
  lb <- ModelInfo$theta.lb
  ub <- ModelInfo$theta.ub
  
  ## --- Transform parameters to unbounded
  for (i in 1:length(u.theta)) {
    if (trans[i] == "raw") { theta[i] <- u.theta[i] }
    if (trans[i] == "log") { theta[i] <- exp(u.theta[i]) }
    if ((trans[i] == "logit")  | (trans[i] == "clogit") |
        (trans[i] == "rlogit") | (trans[i] == "nlogit") ) {
      theta[i] <- (lb[i] + ub[i]*exp(u.theta[i]))/(1 + exp(u.theta[i]))
    }
  }

  ## --- Return unbounded parameters
  return (theta)
}


check_par_values <- function (ModelInfo, theta) {
  ## --- Check if parameter values are within boundaries
  
  ## --- Get parameter lower and upper bounds
  ParInfo <- get_parinfo (ModelInfo)
  
  ## --- Make sure parameters are in correct order
  theta <- theta[ParInfo$theta]
  
  ## --- Parameters legal if all between lower and upper bounds
  Inbounds <- (theta > ParInfo$lb) & (theta < ParInfo$ub)
  legal <- all(Inbounds)
  
  ## --- Print out-of-bounds parameters
  if (!legal) { print (theta[!Inbounds]) }    
  
  ## --- Return parameters legal flag
  return (legal)  
}


check_par_valid <- function (u.theta, ModelInfo, Dat) {
  ## --- Check if model parameters are valid
  
  ## --- Return invalid if any paramters are NaN, or NA
  if ( any(is.nan(u.theta)) | any(is.na(u.theta)) ) {
    valid <- FALSE
    return (valid)
  } 
  
  ## --- Make sure u.theta has model names
  names (u.theta) <- ModelInfo$u.theta
  
  ## --- Grab model info
  err_dist   <- ModelInfo$err_dist
  mean_fun   <- ModelInfo$mean_fun
  binomial_n <- ModelInfo$binomial_n
  delta      <- ModelInfo$delta
  
  ## --- Grab data
  x <- Dat$x
  y <- Dat$y  
  
  ## Parameter is valid by default
  valid <- TRUE

  ## Transform parameter to be bounded
  theta <- make_bounded (ModelInfo, u.theta)
  
  ## Invalid if any parameters missing/NaN/or infinite
  if (any(is.nan(theta)))     { valid <- FALSE }
  if (any(is.na(theta)))      { valid <- FALSE }
  if (any(!is.finite(theta))) { valid <- FALSE }
  
  ## --- Check mean function parameters
  if (valid) { 
  
    ## --- Constant
    if (mean_fun == "constant") {
      ## Nothing to check
    }
    
    ## --- Uniform
    if (mean_fun == "uniform") {
      ## Grab parameters
      c <- theta['c']; d <- theta['d']
      
      ## Check endpoints
      if (c >= d) { valid <- FALSE }
      
      ## Check if any non-zero counts outside of uniform range
      if (valid) { if ( any( y[ (x < c) | (x > d)] != 0) ) { valid <- FALSE } }
    }
    
    ## --- Gaussian
    if (mean_fun == "gaussian") {
      ## Nothing to check
    }
    
    ## --- Mix gaussian
    if ( (mean_fun == "mixgaussian.equal") | (mean_fun == "mixgaussian") ) {
      ## Make sure m1 is less than m2 to prevent label switching
      
      ## Grab parameters
      m1 <- theta['m1']; m2 <- theta['m2']
      
      ## Check that m1 < m2
      if (m1 > m2) { valid <- FALSE }
    }
    
    ## --- Beta
    if (mean_fun == "beta") {
      
      ## Grab parameters
      c <- theta['c']; d <- theta['d']; u <- theta['u']; v <- theta['v']
      ## Convert beta mean and shape parameters
      a <- u / v ;     names(a) <- "a"
      b <- (1 - u)/v ; names(b) <- "b"

      
      ## Check endpoints
      if (c >= d) { valid <- FALSE }
      
      ## Check if any non-zero counts outside of uniform range
      if (valid) { if ( any( y[ (x < c) | (x > d)] != 0) ) { valid <- FALSE } }
      
      ## No u-shaped or j-shaped parameterisation allowed
      if (valid) { if ( (a < 1) | (b < 1)) { valid <- FALSE } }
    }
    
    ## --- Sech
    if (mean_fun == "sech") {
      ## Grab parameters
      r <- theta['r']; s <- theta['s']; p <- theta['p']

      ## Check endpoints 
      if (abs(r) >= 1) { valid <- FALSE }
      if (s <= 0)      { valid <- FALSE }
      if (p <= 0)      { valid <- FALSE }
    }

    ## --- ModSkurt ****
    if (mean_fun == "modskurt") {
      ## Grab parameters
      q <- theta['q']; b <- theta['b']; s <- theta['s']; p <- theta['p']
      
      ## Check endpoints 
      if (abs(q) >= 1) { valid <- FALSE }
      if (b <= 0)      { valid <- FALSE }
      if (s <= 0)      { valid <- FALSE }
      if (p <= 0)      { valid <- FALSE }
    }
    
    ## --- HOF II
    if (mean_fun == "hofII") {
      ## Nothing to check
    }
    
    ## --- HOF IV / HOF IVb
    if ( (mean_fun == "hofIV") | (mean_fun == "hofIVb") ) {
      ## Grab parameters
      w <- theta['w']

      ## w finite & non-zero
      if (w <= 0) { valid <- FALSE }
    }
    
    ## --- HOF V / HOF Vb
    if ( (mean_fun == "hofV") | (mean_fun == "hofVb") ) {
      ## Grab parameters
      w1 <- theta['w1']
      w2 <- theta['w2']
      ## w1 and w2 finite & non-zero
      if (w1 <= 0) { valid <- FALSE }
      if (w2 <= 0) { valid <- FALSE }
    }
  }
  
  ## Return if parameter is valid
  return (valid)
}


NLL.bernoulli <- function (mu.mean, y) {
  ## --- Caclculate negative log-likelihood for Bernoulli data

  ## Grab probabilities
  p <- mu.mean

  ## Calculate negative log-likelihood
  NegLogLike <- -sum(log(1 - p[y==0])) - sum(log(p[y==1]))
  
  ## Return negative log-likelihood
  return (NegLogLike)
}


NLL.binomial.count <- function (n, mu.mean, y) {
  ## --- Caclculate negative log-likelihood for binomial data with n given

  ## Grab probabilities
  p <- mu.mean/n

  
  ## Calculate negative log-likelihood
  NegLogLike <- -sum (stats::dbinom (x=y, size=n, prob=p, log=TRUE))
  
  ## Return negative log-likelihood
  return (NegLogLike)
}


NLL.binomial.prop <- function (n, mu.mean, y) {
  ## --- Caclculate negative log-likelihood for binomial data with n given

  ## Grab probabilities
  p <- mu.mean

  ## Make n an integer
  n <- round(n,0); if (n < 1) { n <- 1 }
  
  ## Convert data to counts [0-n]
  yc <- trunc(y * n); yc[yc < 0] <- 0; yc[yc > n] <- n

  ## Calculate negative log-likelihood
  NegLogLike <- -sum (stats::dbinom (x=yc, size=n, prob=p, log=TRUE))
  
  ## Return negative log-likelihood
  return (NegLogLike)
}


NLL.negbin <- function (phi, mu.mean, y) {
  ## --- Caclculate negative log-likelihood for negative binomial / Poisson data
  
  if (phi == 0) {
    ## Flag to inidicate if mean is zero
    MZ <- mu.mean == 0
    ## --- Poissons
    NegLogLike <- sum(mu.mean) -sum(y[!MZ]*log(mu.mean[!MZ])) -
      sum(log(mu.mean[MZ]^y[MZ])) + sum(lfactorial(y))
  } else {
    ## Alpha
    alpha <- 1/phi

    ## Split data by whether count is zero
    IsZero <- (y == 0)
    Y1 <- y[!IsZero]; E1 <- mu.mean[!IsZero]; MZ1 <- (mu.mean == 0)[!IsZero]
    Y0 <- y[ IsZero]; E0 <- mu.mean[ IsZero]; MZ0 <- (mu.mean == 0)[ IsZero]
    
    ## Caclculate negative log-likelihood - Positive data
    if (sum(!IsZero) > 0) {
      ## Positive data / mean zero not allowed
      if (any(MZ1)) {
        NegLogLikePos <- Inf
      } else {
        ## Reparameterise mean
        p <- alpha/(E1 + alpha)
        ## NLL for positive data / positive mean 
        NegLogLikePos <- -sum(lgamma(alpha+Y1)) + length(Y1)*lgamma(alpha) +
          sum(lgamma(Y1 + 1)) - sum(alpha*log(p)) - sum(Y1*log(1-p))
      }
    } else {
      NegLogLikePos <- 0
    }
    
    ## Caclculate negative log-likelihood - Zero data
    if (sum(IsZero) > 0) {
      ## Reparameterise mean
      p     <- alpha/(E0 + alpha)
      ## NLL for zero data
      NegLogLikeZero <- -sum(log(p^alpha))
    } else {
      NegLogLikeZero <- 0
    }
 
    ## Combine negative log-likelihood for zero and non-zero data
    NegLogLike <- NegLogLikePos + NegLogLikeZero
  }
  
  ## Return negative log-likelihood
  return (NegLogLike)
}


NLL.discrete <- function (u.theta, ModelInfo, Dat) {
  ## --- Negative log-likelihood for discrete data
  
  ## --- Grab data
  x <- Dat$x
  y <- Dat$y  
  
  ## --- Grab model info
  err_dist   <- ModelInfo$err_dist
  mean_fun   <- ModelInfo$mean_fun
  binomial_n <- ModelInfo$binomial_n
  delta      <- ModelInfo$delta

  ## --- Make sure u.theta has model names
  names (u.theta) <- ModelInfo$u.theta

  ## Calculate if parameter is legal
  valid <- check_par_valid (u.theta=u.theta, ModelInfo=ModelInfo, Dat=Dat)
  
  ## Calculate negative log-likelihood
  if (!valid) {
    
    ## --- Par invalid : NLL infinite
    NegLogLike <- Inf

  } else { 

    ## --- Data OK
    
    ## --- MODEL PARAMETERS

    ## --- Transform to bounded
    theta <- make_bounded (ModelInfo, u.theta)
    
    ## Split parameters into error distribution and mean function parameters
    thetaE   <- theta[ModelInfo$thetaE]
    thetaM   <- theta[ModelInfo$thetaM]
    u.thetaM <- u.theta[match (ModelInfo$thetaM , names(theta))]
    
    ## --- Negative binomial
    if ( err_dist == "negbin" ) {
      phi <- theta['phi']
    }
    ## --- Zero-inflated Poisson
    if ( err_dist == "zip" ) {
      pi  <- theta['pi']
    }
    ## --- Zero-inflated negative binomial
    if ( err_dist == "zinb" ) {
      pi  <- theta['pi']
      phi <- theta['phi']
    }
    ## --- Zero-linked Poisson
    if ( (err_dist == "zipl") | (err_dist == "zipl.mu") ) {
      g0  <- theta['g0']
      g1  <- theta['g1']
    }
    ## --- Zero-linked negative binomial
    if ( (err_dist == "zinbl") | (err_dist == "zinbl.mu") ) {
      g0  <- theta['g0']
      g1  <- theta['g1']
      phi <- theta['phi']
    }
    
    ## --- MEAN FUNCTION
    
    ## --- Uniform & modified beta
    ##     Do nothing for constant and sech/modskurt mean functions
    if ( (mean_fun == "uniform") | (mean_fun == "beta") ) {
      
      ## Remove observations outside of range of uniform
      y <- y[ (x >= theta['c']) & (x <= theta['d']) ]
      x <- x[ (x >= theta['c']) & (x <= theta['d']) ]
    }
    
    ## Calculate mean function
    mu <- do.call (paste("mu_", mean_fun, sep=""), list(thetaM, x))
    
    ## --- NEGATIVE LOG-LIKELIHOOD
    
    ## --- Bernoulli
    if ( err_dist == "bernoulli" ) {
      ## Caclculate negative log-likelihood
      NegLogLike <- NLL.bernoulli (mu.mean=mu, y=y)
    }
    
    ## --- Binomial - count
    if ( err_dist == "binomial.count" ) {
      ## Caclculate negative log-likelihood
      NegLogLike <- NLL.binomial.count (n=binomial_n, mu.mean=mu, y=y)
    }
    
    ## --- Binomial - prop
    if ( err_dist == "binomial.prop" ) {
      ## Caclculate negative log-likelihood
      NegLogLike <- NLL.binomial.prop (n=binomial_n, mu.mean=mu, y=y)
    }
    
    ## --- Poisson
    if ( err_dist == "poisson" ) {
      ## Caclculate negative log-likelihood
      NegLogLike <- NLL.negbin (phi=0, mu.mean=mu, y=y)
    }
    
    ## --- Negative binomial
    if ( err_dist == "negbin" ) {
      ## Caclculate negative log-likelihood
      NegLogLike <- NLL.negbin (phi=phi, mu.mean=mu, y=y)
    }
    
    ## --- Zero-inflated Poisson / Zero-inflated negative binomial
    if ( (err_dist == "zip") | (err_dist == "zinb") ) {
      
      ## Split data by whether count is zero
      IsZero <- (y == 0)
      Y1 <- y[!IsZero]; X1 <- x[!IsZero]; E1 <- mu[!IsZero]
      Y0 <- y[ IsZero]; X0 <- x[ IsZero]; E0 <- mu[ IsZero]
      
      ## Caclculate negative log-likelihood - Positive data
      if (sum(!IsZero) > 0) {
        if ( err_dist == "zip" ) {
          NegLogLikePos <- NLL.negbin (phi=0,   mu.mean=E1, y=Y1) - length(Y1)*log(1-pi)
        }
        if ( err_dist == "zinb" ) {
          NegLogLikePos <- NLL.negbin (phi=phi, mu.mean=E1, y=Y1) - length(Y1)*log(1-pi)
        }
      } else {
        NegLogLikePos <- 0
      }
      
      ## Caclculate negative log-likelihood - Zero data
      if (sum(IsZero) > 0) {
        if ( err_dist == "zip" ) {
          NegLogLikeZero <- -sum(log(pi + exp(-E0)*(1-pi)))
        }
        if ( err_dist == "zinb" ) {
          alpha <- 1/phi
          p     <- alpha/(E0 + alpha)
          NegLogLikeZero <- -sum(log(pi + (p^alpha)*(1-pi)))
        }
      } else {
        NegLogLikeZero <- 0
      }
      
      ## Combine negative log-likelihood for zero and non-zero data
      NegLogLike <- NegLogLikePos + NegLogLikeZero
    }

    ## --- Zero-linked Poisson / Zero-linked negative binomial
    if ( (err_dist == "zipl")  | (err_dist == "zipl.mu")  |
         (err_dist == "zinbl") | (err_dist == "zinbl.mu") ) {

      ## Check if mu is zero anywhere
      MeanNotZero <- (mu > 0)

      ## Check if all y=0 when mu=0
      if (any (MeanNotZero == FALSE) ) {
        ym0 <- y[which(!MeanNotZero)]
        ## Is y always zero when mean is zero?
        AllZero <- all (ym0==0)
      } else {
        AllZero <- TRUE  
      }
      
      ## --- Calculate NLL
      if (AllZero == FALSE) {

        ## --- Model failed - y>0 when mu=0
        NegLogLike <- Inf

      } else {

        ## Remove values where mu and y are both zero
        y <- y[MeanNotZero]; x <- x[MeanNotZero]; E <- mu[MeanNotZero] 
        
        ## Split data by whether count is zero
        IsZero <- (y == 0)
        Y1 <- y[!IsZero]; X1 <- x[!IsZero]; E1 <- E[!IsZero]
        Y0 <- y[ IsZero]; X0 <- x[ IsZero]; E0 <- E[ IsZero]

        ## Logit zero spike
        if ( (err_dist == "zipl")  | (err_dist == "zinbl")  ) { 
          logitpi1 <- g0 + g1*log(E1)
          logitpi0 <- g0 + g1*log(E0)
        }
        if ( (err_dist == "zipl.mu") | (err_dist == "zinbl.mu") ) { 
          logitpi1 <- g0 + g1*E1
          logitpi0 <- g0 + g1*E0
        }
        ## Zero spike probability
        pi1 <- exp(logitpi1) / ( 1 + exp(logitpi1) )
        pi0 <- exp(logitpi0) / ( 1 + exp(logitpi0) )
        
        ## Caclculate negative log-likelihood - Positive data
        if (sum(!IsZero) > 0) {
          if ( (err_dist == "zipl")  | (err_dist == "zipl.mu")  ) {
            NegLogLikePos <- NLL.negbin (phi=0, mu.mean=E1, y=Y1) - sum(log(1-pi1))
          }
          if ( (err_dist == "zinbl") | (err_dist == "zinbl.mu") ) {
            NegLogLikePos <- NLL.negbin (phi=phi, mu.mean=E1, y=Y1) - sum(log(1-pi1))
          }
        } else {
          NegLogLikePos <- 0
        }
        
        ## Caclculate negative log-likelihood - Zero data
        if (sum(IsZero) > 0) {
          if ( (err_dist == "zipl")  | (err_dist == "zipl.mu")  ) {
            NegLogLikeZero <- -sum(log(pi0 + exp(-E0)*(1-pi0)))
          }
          if ( (err_dist == "zinbl") | (err_dist == "zinbl.mu") ) {
            alpha <- 1/phi
            p     <- alpha/(E0 + alpha)
            NegLogLikeZero <- -sum(log(pi0 + (p^alpha)*(1-pi0)))
          }
        } else {
          NegLogLikeZero <- 0
        }
        
        ## Combine negative log-likelihood for zero and non-zero data
        NegLogLike <- NegLogLikePos + NegLogLikeZero
      }
    }
    
  }
  
  ## Return negative log-likelihood
  return (NegLogLike)  
}


NLL.continuous <- function (u.theta, ModelInfo, Dat) {
  ## --- Negative log-likelihood for continuous distribution
  ## --- Gaussian, Tweedie, Zero-inflated gamma
  
  ## --- Grab data
  x <- Dat$x
  y <- Dat$y
  
  ## --- Grab model info
  err_dist <- ModelInfo$err_dist
  mean_fun <- ModelInfo$mean_fun
  delta    <- ModelInfo$delta

  ## --- Make sure u.theta has model names
  names (u.theta) <- ModelInfo$u.theta
  
  ## Calculate if parameter is legal
  valid <- check_par_valid (u.theta=u.theta, ModelInfo=ModelInfo, Dat=Dat)
  
  ## Calculate negative log-likelihood
  if (!valid) {
    
    ## --- Par invalid : NLL infinite
    NegLogLike <- Inf

  } else { 
    
    ## --- Data OK
    
    ## --- MODEL PARAMETERS

    ## --- Transform to bounded
    theta <- make_bounded (ModelInfo, u.theta)
    
    ## Split parameters into error distribution and mean function parameters
    thetaE   <- theta[ModelInfo$thetaE]
    thetaM   <- theta[ModelInfo$thetaM]
    u.thetaM <- u.theta[match (ModelInfo$thetaM , names(theta))]
    
    ## --- Separate model and error distribution parameters

    ## --- Tail-adjusted beta
    if ( err_dist == "tab" ) {
      sigma <- theta['sigma']
    }
    ## --- Zero-inflated beta
    if ( err_dist == "zitab" ) {
      pi    <- theta['pi']
      sigma <- theta['sigma']
    }
    ## --- Gaussian
    if ( err_dist == "gaussian" ) {
      sigma <- theta['sigma']
    }
    ## --- Tweedie
    if ( err_dist == "tweedie" ) {
      phi   <- theta['phi']
      rho   <- theta['rho']
    }
    ## --- Zero-inflated gamma
    if ( err_dist == "zig" ) {
      pi    <- theta['pi']
      phi   <- theta['phi']
    }
    ## --- Zero-linked gamma
    if ( (err_dist == "zigl") | (err_dist == "zigl.mu") ) {
      g0    <- theta['g0']
      g1    <- theta['g1']
      phi   <- theta['phi']
    }
    ## --- Zero-inflated inverse gaussian
    if ( err_dist == "ziig" ) {
      pi    <- theta['pi']
      phi   <- theta['phi']
    }
    ## --- Zero-linked inverse gaussian
    if ( (err_dist == "ziigl") | (err_dist == "ziigl.mu") ) {
      g0    <- theta['g0']
      g1    <- theta['g1']
      phi   <- theta['phi']
    }
    
    
    ## --- MEAN FUNCTION
    
    ## --- Uniform & modified beta
    ##     Do nothing for constant and sech/modskurt mean functions
    if ( (mean_fun == "uniform") | (mean_fun == "beta") ) {
      
      ## Remove observations outside of range of uniform
      y <- y[ (x >= theta['c']) & (x <= theta['d'])]
      x <- x[ (x >= theta['c']) & (x <= theta['d'])]
    }
    
    ## --- Calculate mean function
    mu <- do.call (paste("mu_", mean_fun, sep=""), list(thetaM, x))

    ## Make sure means not negative or NaN's
    mu[mu <= 0]    <- 0
    mu[is.nan(mu)] <- 0 
    
    ## --- NEGATIVE LOG-LIKELIHOOD

    ## Tailed-adjusted beta/ Zero-inflated tail-adjusted beta
    if ( (err_dist == "tab") | (err_dist == "zitab") ) {
      
      ## --- Create index variables : 0, < delta, between delta and 1-delta, > 1-delta 

      ## Lower section (< delta)
      ID0 <-  y < delta                            ## Less than delta
      ## Upper section (> 1-delta)
      ID1 <-  y > (1-delta)                        ## Greater than 1-delta
      ## Middle section (>= delta and =< 1-delta)
      ID2 <- (y >= delta) & (y <= (1-delta))       ## Between delta and 1-delta
      
      ## Split data and means:  Lower section
      n <- length(y)
      YD0 <- y[ID0] ; XD0 <- x[ID0] ; ED0 <- mu[ID0] ; nD0 <- length(ED0)
      ## Split data and means: Upper section
      YD1 <- y[ID1] ; XD1 <- x[ID1] ; ED1 <- mu[ID1] ; nD1 <- length(ED1)
      ## Split data and means: Middle section
      YD2 <- y[ID2] ; XD2 <- x[ID2] ; ED2 <- mu[ID2] ; nD2 <- length(ED2)
      
      ## Grab beta parameters
      Alpha <- mu / sigma
      Beta  <- (1 - mu) / sigma
      
      ## Calculate probability < delta
      log.p0 <- stats::pbeta (delta, shape1=Alpha, shape2=Beta, log.p=TRUE)
      log.p1 <- stats::pbeta (1-delta, shape1=Alpha, shape2=Beta, lower.tail=FALSE, log.p=TRUE)
      p0 <- exp (log.p0) ; p1 <- exp (log.p1) ; p2 <- 1 - p1 - p0
      log.p2 <- log (p2)
    }
    
    ## --- NLL : Tail adjusted beta
    if ( err_dist == "tab" ) {
      NegLogLike <- -sum (log.p0[ID0]) - sum (log.p1[ID1]) -
        sum(stats::dbeta(x=YD2, shape1=Alpha[ID2], shape2=Beta[ID2], log=TRUE))
    }
    
    ## --- Zero-inflated tail-adjusted beta
    if ( err_dist == "zitab" ) {
      
      ## --- Negative log-likelihood components

      ## Zero data
      NegLogLike.0  <- -sum ( log(pi + (1 - pi)*p0[ID0]) )
      ## One data
      NegLogLike.1  <- -sum ( log((1 - pi)*p1[ID1]))
      ## Beta section
      NegLogLike.2 <- -sum(stats::dbeta (x=YD2, shape1=Alpha[ID2], shape2=Beta[ID2], log=TRUE))
      ## Prob of beta
      NegLogLike.3 <- - nD2 * log (1-pi)
      if ( (nD2==0) & (pi=1) ) { NegLogLike.3 <- 0 }
      
      ## --- NLL : zitab
      NegLogLike <- NegLogLike.0 + NegLogLike.1 + NegLogLike.2 + NegLogLike.3
    }
    
    ## --- Gaussian
    if ( err_dist == "gaussian" ) {
      NegLogLike <- -sum(stats::dnorm(y, mean=mu, sd=sigma, log=TRUE))
    }

    ## --- Tweedie
    if ( err_dist == "tweedie" ) {

      ## If any positive y's with zero mu -> fail
      if (any(y[mu == 0] > 0)) {

        ## Negative log-likelihood illegal
        NegLogLike <- Inf

      } else {

        ## Remove observation where y's and mu's both zero
        Mu <- mu[mu > 0]
        Y  <- y[mu > 0]
        NegLogLike <- -sum(log(tweedie::dtweedie(y=Y, xi=rho, mu=Mu, phi=phi)))
      }
    }
    
    ## --- Zero-inflated gamma
    if ( (err_dist == "zig") | (err_dist == "zigl") | (err_dist == "zigl.mu") ) { 
      
      ## Split data by whether count is zero
      IsZero <- (y == 0)
      Y1 <- y[!IsZero]; X1 <- x[!IsZero]; E1 <- mu[!IsZero]
      Y0 <- y[ IsZero]; X0 <- x[ IsZero]; E0 <- mu[ IsZero]

      ## Parameters
      alpha <- 1/phi
      beta0 <- alpha/E0
      beta1 <- alpha/E1
      
      ## Set spike parameters on logit scale
      if ( err_dist == "zigl" ) { 
        logitpi1 <- g0 + g1*log(E1)
        logitpi0 <- g0 + g1*log(E0)
      }
      if ( err_dist == "zigl.mu" ) { 
        logitpi1 <- g0 + g1*E1
        logitpi0 <- g0 + g1*E0
      }
      
      ## Untransform spike parameters
      if ( (err_dist == "zigl") | (err_dist == "zigl.mu") ) { 
        pi1 <- exp(logitpi1) / ( 1 + exp(logitpi1) )
        pi0 <- exp(logitpi0) / ( 1 + exp(logitpi0) )
      }
      ## Convert pi to vector for zig
      if ( err_dist == "zig" ) {
        pi1 <- rep(pi, length(Y1))
        pi0 <- rep(pi, length(Y0))
      }
      
      ## Caclculate negative log-likelihood - Positive data
      if (sum(!IsZero) > 0) {
        alpha <- 1/phi
        beta  <- alpha/E1
        NegLogLikePos <- -sum(stats::dgamma(x=Y1, shape=alpha, rate=beta1, log=TRUE) + log(1-pi1))
      } else {
        NegLogLikePos <- 0
      }
      
      ## Caclculate negative log-likelihood - Zero data
      if (sum(IsZero) > 0) {
        alpha <- 1/phi
        beta  <- alpha/E0
        NegLogLikeZero <- -sum(log(pi0))
      } else {
        NegLogLikeZero <- 0
      }
      
      ## Combine negative log-likelihood for zero and non-zero data
      NegLogLike <- NegLogLikePos + NegLogLikeZero
    }

    
    ## --- Zero-inflated inverse gaussian
    if ( (err_dist == "ziig") | (err_dist == "ziigl") | (err_dist == "ziigl.mu") ) { 
      
      ## Split data by whether count is zero
      IsZero <- (y == 0)
      Y1 <- y[!IsZero]; X1 <- x[!IsZero]; E1 <- mu[!IsZero]
      Y0 <- y[ IsZero]; X0 <- x[ IsZero]; E0 <- mu[ IsZero]

      ## Set spike parameters on logit scale
      if ( err_dist == "ziigl" ) { 
        logitpi1 <- g0 + g1*log(E1)
        logitpi0 <- g0 + g1*log(E0)
      }
      if ( err_dist == "ziigl.mu" ) { 
        logitpi1 <- g0 + g1*E1
        logitpi0 <- g0 + g1*E0
      }
      
      ## Untransform spike parameters
      if ( (err_dist == "ziigl") | (err_dist == "ziigl.mu") ) { 
        pi1 <- exp(logitpi1) / ( 1 + exp(logitpi1) )
        pi0 <- exp(logitpi0) / ( 1 + exp(logitpi0) )
      }
      ## Convert pi to vector for ziig
      if ( err_dist == "ziig" ) {
        pi1 <- rep(pi, length(Y1))
        pi0 <- rep(pi, length(Y0))
      }
      
      ## Caclculate negative log-likelihood - Positive data
      if (sum(!IsZero) > 0) {
        NegLogLikePos <- -sum(dinversegaussian(y=Y1, mu=E1, phi=phi, log=TRUE) + log(1-pi1))
      } else {
        NegLogLikePos <- 0
      }
      
      ## Caclculate negative log-likelihood - Zero data
      if (sum(IsZero) > 0) {
        NegLogLikeZero <- -sum(log(pi0))
      } else {
        NegLogLikeZero <- 0
      }
      
      ## Combine negative log-likelihood for zero and non-zero data
      NegLogLike <- NegLogLikePos + NegLogLikeZero
    }
  }

  
  ## Return negative log-likelihood
  return (NegLogLike)  
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MLE.R
\name{msenlm}
\alias{msenlm}
\title{Fit multiple species-environment non-linear models via maximum likelihood}
\usage{
msenlm(
  models = NULL,
  data = NULL,
  xvar = NULL,
  yvar = NULL,
  method = "crossed",
  conf.level = 0.95
)
}
\arguments{
\item{models}{Object listing models to fit (from set_models function).}

\item{data}{A data frame containing 'x' (explanatory) and 'y' (response) variables.}

\item{xvar}{Name of explanatory variable (must be univariate).}

\item{yvar}{Names of response variables.}

\item{method}{If "crossed", fit all models to all response variables. If "paired",
fit first model to first response variables, etc.}

\item{conf.level}{Confidence level for parameter confidence intervals. Default is 0.95.}
}
\value{
Object containg model fits to data y and x.
}
\description{
'msenlm' fits multiple species-environment non-linear models via maximum likelihood.
}
\examples{

\dontrun{

models <- set_models(mean_fun=c("gaussian","beta"), err_dist=c("zip","zinb"))
fits <- msenlm(models=models, data=haul, xvar="depth",
               yvar=c("Albatrossia.pectoralis", "Sebastolobus.altivelis"))
}
}
\keyword{fit}
\keyword{mle}
\keyword{model,}
\keyword{senlm}

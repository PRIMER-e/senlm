% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Models.R
\name{get_constant_parnames}
\alias{get_constant_parnames}
\title{Constant parameters names}
\usage{
get_constant_parnames(err_dist)
}
\arguments{
\item{err_dist}{Error distribution name.}
}
\value{
Vector containing names of constant parameters.
}
\description{
Get constant parameter names for given error distribution.
}
\examples{

## Get parameters for a "zinb" er
get_constant_parnames ("zinb")

}
\keyword{constant}
\keyword{parameters}

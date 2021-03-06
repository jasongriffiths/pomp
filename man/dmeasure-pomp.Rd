\name{dmeasure-pomp}
\docType{methods}
\alias{dmeasure}
\alias{dmeasure,pomp-method}
\alias{dmeasure-pomp}
\keyword{internal}
\title{Evaluate the probability density of observations given underlying states in a partially-observed Markov process}
\description{
  The method \code{dmeasure} evaluates the probability density of a set of measurements given the state of the system.
  This function is part of the low-level interface to \code{pomp} objects.
  This help page does not give instructions on the implementation of models: see \code{\link{pomp}} for instructions.
}
\usage{
\S4method{dmeasure}{pomp}(object, y, x, times, params, log = FALSE, \dots)
}
\arguments{
  \item{object}{an object of class \code{pomp}.}
  \item{y}{
    a rank-2 array containing observations.
    The dimensions of \code{y} are \code{nobs} x \code{ntimes}, where \code{nobs} is the number of observables and \code{ntimes} is the length of \code{times}.
  }
  \item{x}{
    a rank-3 array containing the states of the unobserved process.
    The dimensions of \code{x} are \code{nvars} x \code{nreps} x \code{ntimes}, where \code{nvars} is the number of state variables, \code{nreps} is the number of replicates, and \code{ntimes} is the length of \code{times}.
  }
  \item{times}{
    a numeric vector containing the times at which the observations were made.
  }
  \item{params}{
    a rank-2 array of parameters with columns corresponding to the columns of \code{x}.
    Note that the \code{x} and \code{params} must agree in the number of their columns.
  }
  \item{log}{if TRUE, log probabilities are returned.}
  \item{\dots}{at present, these are ignored.}
}
\value{
  Returns a matrix of dimensions \code{nreps} x \code{ntimes}.
  If \code{d} is the returned matrix, \code{d[j,k]} is the likelihood of the observation \code{y[,k]} at time \code{times[k]} given the state \code{x[,j,k]}.
}
\details{
  This function is essentially a wrapper around the user-supplied \code{dmeasure} slot of the \code{pomp} object.
  For specifications on writing such a function, see \code{\link{pomp}}.
}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\seealso{\code{\link{pomp-class}}, \code{\link{pomp}}}
\keyword{models}
\keyword{ts}

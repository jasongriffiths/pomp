\name{rw2}
\alias{rw2}
\docType{data}
\title{Two-dimensional random-walk process}
\description{
  \code{rw2} is a \code{pomp} object encoding a 2-D normal random walk.
}
\usage{data(rw2)}
\details{
  The random-walk process is fully but noisily observed.
}
\examples{
data(rw2)
plot(rw2)
x <- simulate(rw2,nsim=10,seed=20348585L,params=c(x1.0=0,x2.0=0,s1=1,s2=3,tau=1))
plot(x[[1]])
}
\seealso{\code{\link{pomp-class}} and the vignettes}
\keyword{datasets}
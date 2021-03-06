\name{eulermultinom}
\alias{eulermultinom}
\alias{reulermultinom}
\alias{deulermultinom}
\alias{rgammawn}
\title{Euler-multinomial death process}
\description{
  Density and random-deviate generation for the Euler-multinomial death process with parameters \code{size}, \code{rate}, and \code{dt}.
}
\usage{
reulermultinom(n = 1, size, rate, dt)
deulermultinom(x, size, rate, dt, log = FALSE)
rgammawn(n = 1, sigma, dt)
}
\arguments{
  \item{n}{integer; number of random variates to generate.}
  \item{size}{scalar integer; number of individuals at risk.}
  \item{rate}{numeric vector of hazard rates.}
  \item{sigma}{numeric scalar; intensity of the Gamma white noise process.}
  \item{dt}{numeric scalar; duration of Euler step.}
  \item{x}{matrix or vector containing number of individuals that have succumbed to each death process.}
  \item{log}{logical; if TRUE, return logarithm(s) of probabilities.}
}
\details{
  If \eqn{N} individuals face constant hazards of death in \eqn{k} ways at rates \eqn{r_1, r_2, \dots, r_k}{r1,r2,\dots,rk}, then in an interval of duration \eqn{\Delta t}{dt}, the number of individuals remaining alive and dying in each way is multinomially distributed:
  \deqn{(N-\sum_{i=1}^k \Delta n_i, \Delta n_1, \dots, \Delta n_k) \sim \mathrm{multinomial}(N;p_0,p_1,\dots,p_k),}{(N-\sum(dni), dn1, \dots, dnk) ~ multinomial(N;p0,p1,\dots,pk),}
  where \eqn{\Delta n_i}{dni} is the number of individuals dying in way \eqn{i} over the interval, the probability of remaining alive is \eqn{p_0=\exp(-\sum_i r_i \Delta t)}{p0=exp(-\sum(ri dt))}, and the probability of dying in way \eqn{j} is
  \deqn{p_j=\frac{r_j}{\sum_i r_i} (1-\exp(-\sum_i r_i \Delta t)).}{pj=(1-exp(-sum(ri dt))) rj/(\sum(ri)).}
  In this case, we can say that
  \deqn{(\Delta n_1, \dots, \Delta n_k) \sim \mathrm{eulermultinom}(N,r,\Delta t),}{(dn1,\dots,dnk)~eulermultinom(N,r,dt),}
  where \eqn{r=(r_1,\dots,r_k)}{r=(r1,\dots,rk)}.
  Draw \eqn{m} random samples from this distribution by doing
  
  \code{dn <- reulermultinom(n=m,size=N,rate=r,dt=dt)},

  where \code{r} is the vector of rates.
  Evaluate the probability that \eqn{x=(x_1,\dots,x_k)}{x=(x1,\dots,xk)} are the numbers of individuals who have died in each of the \eqn{k} ways over the interval \eqn{\Delta t=}{}\code{dt}, by doing

  \code{deulermultinom(x=x,size=N,rate=r,dt=dt)}.
  
  Bret\'o & Ionides discuss how an infinitesimally overdispersed death process can be constructed by compounding a binomial process with a Gamma white noise process.
  The Euler approximation of the resulting process can be obtained as follows.
  Let the increments of the equidispersed process be given by

  \code{reulermultinom(size=N,rate=r,dt=dt)}.
  
  In this expression, replace the rate \eqn{r} with \eqn{r {\Delta W}/{\Delta t}},
  where \eqn{\Delta W ~ Gamma(dt/\sigma^2,\sigma^2)} is the increment of an integrated Gamma white noise process with intensity \eqn{\sigma}.
  That is, \eqn{\Delta W} has mean \eqn{\Delta t} and variance \eqn{\sigma^2 \Delta t}.
  The resulting process is overdispersed and converges (as \eqn{\Delta t} goes to zero) to a well-defined process.
  The following lines of \R code accomplish this:
  \preformatted{
    dW <- rgammawn(sigma=sigma,dt=dt)
    dn <- reulermultinom(size=N,rate=r,dt=dW)
  }
  or
  \preformatted{
    dn <- reulermultinom(size=N,rate=r*dW/dt,dt=dt).
  }
  He et al. use such overdispersed death processes in modeling measles.
  
  For all of the functions described here, direct access to the underlying C routines is available:
  see the header file \dQuote{pomp.h}, included with the package.
}
\value{
  \item{reulermultinom}{
    Returns a \code{length(rate)} by \code{n} matrix.
    Each column is a different random draw.
    Each row contains the numbers of individuals succumbed to the corresponding process.
  }
  \item{deulermultinom}{
    Returns a vector (of length equal to the number of columns of \code{x}) containing the probabilities of observing each column of \code{x} given the specified parameters (\code{size}, \code{rate}, \code{dt}).
  }
  \item{rgammawn}{
    Returns a vector of length \code{n} containing random increments of the integrated Gamma white noise process with intensity \code{sigma}.
  }
}
\author{Aaron A. King \email{kingaa at umich dot edu}}
\examples{
print(dn <- reulermultinom(5,size=100,rate=c(a=1,b=2,c=3),dt=0.1))
deulermultinom(x=dn,size=100,rate=c(1,2,3),dt=0.1)
## an Euler-multinomial with overdispersed transitions:
dt <- 0.1
dW <- rgammawn(sigma=0.1,dt=dt)
print(dn <- reulermultinom(5,size=100,rate=c(a=1,b=2,c=3),dt=dW))
}
\references{
  C. Bret\'o & E. L. Ionides,
  Compound Markov counting processes and their applications to modeling infinitesimally over-dispersed systems.
  Stoch. Proc. Appl., 121:2571--2591, 2011.

  D. He, E. L. Ionides, & A. A. King,
  Plug-and-play inference for disease dynamics: measles in large and small populations as a case study.
  J. R. Soc. Interface, 7:271--283, 2010.
  }
\keyword{distribution}

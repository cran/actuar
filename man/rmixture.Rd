\name{rmixture}
\alias{rmixture}
\title{Simulation from Discrete Mixtures}
\description{
  Generate random variates from a discrete mixture of distributions.
}
\usage{
rmixture(n, probs, models, shuffle = TRUE)
}
\arguments{
  \item{n}{number of random variates to generate. If \code{length(n) >
    1}, the length is taken to be the number required.}
  \item{probs}{numeric non-negative vector specifying the probability
    for each model; is internally normalized to sum 1. Infinite
    and missing values are not allowed. Values are recycled as necessary
    to match the length of \code{models}.}
  \item{models}{vector of expressions specifying the simulation models
    with the number of variates omitted; see Details. Models are
    recycled as necessary to match the length of \code{probs}.}
  \item{shuffle}{logical; should the random variates from the
    distributions be shuffled?}
}
\details{
  \code{rmixture} generates variates from a discrete mixture, that is
  random variable with a probability density function of the form
  \deqn{f(x) = p_1 f_1(x) + ... + p_n f_n(x),}
  where \eqn{f_1, \dots, f_n} are densities and \eqn{\sum_{i = 1}^n p_i
    = 1}{p_1 + \dots + p_n = 1}.

  The values in \code{probs} will be internally normalized to be
  used as probabilities \eqn{p_1 + \dots + p_n}.

  The specification of simulation models uses the syntax of
  \code{\link{rcomphierarc}}. Models \eqn{f_1, \dots, f_n} are expressed in a
  semi-symbolic fashion using an object of mode
  \code{\link[base]{expression}} where each element is a complete call
  to a random number generation function, with the number of variates
  omitted.

  The argument of the random number generation functions for the number
  of variates to simulate \strong{must} be named \code{n}.

  If \code{shuffle} is \code{FALSE}, the output vector contains all the
  random variates from the first model, then all the random variates
  from the second model, and so on. If the order of the variates is
  irrelevant, this cuts the time to generate the variates roughly in
  half.
}
\value{
  A vector of random variates from the mixture with density \eqn{f(x)}.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca}
}
\seealso{
  \code{\link{rcompound}} to simulate from compound models.

  \code{\link{rcomphierarc}} to simulate from compound hierarchical models.
}
\examples{
## Mixture of two exponentials (with means 1/3 and 1/7) with equal
## probabilities.
rmixture(10, 0.5, expression(rexp(3), rexp(7)))
rmixture(10, 42, expression(rexp(3), rexp(7))) # same

## Mixture of two lognormals with different probabilities.
rmixture(10, probs = c(0.55, 0.45),
         models = expression(rlnorm(3.6, 0.6),
                             rlnorm(4.6, 0.3)))
}
\keyword{datagen}

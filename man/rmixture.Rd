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
  the random variable with a probability density function of the form
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
\note{
  Building the expressions in \code{models} from the arguments of
  another function is delicate. The expressions must be such that
  evaluation is possible in the frame of \code{rmixture} or its parent.
  See the examples.  
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

## Building the model expressions in the following example
## works as 'rate' is defined in the parent frame of
## 'rmixture'.
probs <- c(2, 5)
g <- function(n, p, rate)
    rmixture(n, p, expression(rexp(rate[1]), rexp(rate[2])))
g(10, probs, c(3, 7))

## The following example does not work: 'rate' does not exist
## in the evaluation frame of 'rmixture'.
f <- function(n, p, model) rmixture(n, p, model)
h <- function(n, p, rate)
    f(n, p, expression(rexp(rate[1]), rexp(rate[2])))
\dontrun{h(10, probs, c(3, 7))}

## Fix: substitute the values in the model expressions.
h <- function(n, p, rate)
{
    models <- eval(substitute(expression(rexp(a[1]), rexp(a[2])),
                              list(a = rate)))
    f(n, p, models)
}
h(10, probs, c(3, 7))
}
\keyword{datagen}

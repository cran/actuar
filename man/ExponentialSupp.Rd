\name{ExponentialSupp}
\alias{ExponentialSupp}
\alias{mexp}
\alias{levexp}
\alias{mgfexp}
\title{Moments and Moment Generating Function of the Exponential Distribution}
\description{
  Raw moments, limited moments and moment generating function for the
  exponential distribution with rate \code{rate} (i.e., mean
  \code{1/rate}).
}
\usage{
mexp(order, rate = 1)
levexp(limit, rate = 1, order = 1)
mgfexp(x, rate = 1, log = FALSE)
}
\arguments{
  \item{order}{order of the moment.}
  \item{limit}{limit of the loss variable.}
  \item{rate}{vector of rates.}
  \item{x}{numeric vector.}
  \item{log}{logical; if \code{TRUE}, the cumulant generating function
    is returned.}
}
\details{
  The \eqn{k}th raw moment of the random variable \eqn{X} is
  \eqn{E[X^k]}{E[X^k]}, the \eqn{k}th limited moment at some limit
  \eqn{d} is \eqn{E[\min(X, d)^k]}{E[min(X, d)^k]} and the moment
  generating function is \eqn{E[e^{xX}]}.
}
\value{
  \code{mexp} gives the \eqn{k}th raw moment,
  \code{levexp} gives the \eqn{k}th moment of the limited loss
  variable, and
  \code{mgfexp} gives the moment generating function in \code{x}.
  
  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\seealso{
  \code{\link{Exponential}}
}
\references{
  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2004),
  \emph{Loss Models, From Data to Decisions, Second Edition}, Wiley.

  Johnson, N. L. and Kotz, S. (1970), \emph{Continuous univariate
    distributions, Volume 1}, Wiley.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca},
  Christophe Dutang and Mathieu Pigeon.
}
\examples{
mexp(2, 3) - mexp(1, 3)^2
levexp(10, 3, order = 2)
mgfexp(1,2)
}
\keyword{distribution}
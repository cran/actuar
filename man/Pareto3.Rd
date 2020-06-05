\name{Pareto3}
\alias{Pareto3}
\alias{dpareto3}
\alias{ppareto3}
\alias{qpareto3}
\alias{rpareto3}
\alias{mpareto3}
\alias{levpareto3}
\title{The Pareto III Distribution}
\description{
  Density function, distribution function, quantile function, random generation,
  raw moments and limited moments for the Pareto III distribution with
  parameters \code{min}, \code{shape} and \code{scale}.
}
\usage{
dpareto3(x, min, shape, rate = 1, scale = 1/rate,
         log = FALSE)
ppareto3(q, min, shape, rate = 1, scale = 1/rate,
         lower.tail = TRUE, log.p = FALSE)
qpareto3(p, min, shape, rate = 1, scale = 1/rate,
         lower.tail = TRUE, log.p = FALSE)
rpareto3(n, min, shape, rate = 1, scale = 1/rate)
mpareto3(order, min, shape, rate = 1, scale = 1/rate)
levpareto3(limit, min, shape, rate = 1, scale = 1/rate,
           order = 1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{min}{lower bound of the support of the distribution.}
  \item{shape, scale}{parameters. Must be strictly positive.}
  \item{rate}{an alternative way to specify the scale.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities/densities
    \eqn{p} are returned as \eqn{\log(p)}{log(p)}.}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}.}
  \item{order}{order of the moment.}
  \item{limit}{limit of the loss variable.}
}
\details{
  The Pareto III (or \dQuote{type III}) distribution with parameters
  \code{min} \eqn{= \mu}{= m},
  \code{shape} \eqn{= \gamma}{= b} and
  \code{scale} \eqn{= \theta}{= s} has density:
  \deqn{f(x) = \frac{\gamma ((x - \mu)/\theta)^{\gamma - 1}}{%
      \theta [1 + ((x - \mu)/\theta)^\gamma]^2}}{%
    f(x) = (b ((x - m)/s)^(b - 1))/(s [1 + ((x - m)/s)^b]^2)}
  for \eqn{x > \mu}{x > m}, \eqn{-\infty < \mu < \infty}{-Inf < m < Inf},
  \eqn{\gamma > 0}{b > 0} and \eqn{\theta > 0}{s > 0}.

  The Pareto III is the distribution of the random variable
  \deqn{\mu + \theta \left(\frac{X}{1 - X}\right)^{1/\gamma},}{%
    m + s (X/(1 - X))^(1/b),}
  where \eqn{X} has a uniform distribution on \eqn{(0, 1)}. It derives
  from the \link[=dfpareto]{Feller-Pareto} 
  distribution with \eqn{\alpha = \tau = 1}{shape1 = shape3 = 1}.
  Setting \eqn{\mu = 0}{min = 0} yields the \link[=dllogis]{loglogistic}
  distribution.

  The \eqn{k}th raw moment of the random variable \eqn{X} is
  \eqn{E[X^k]}{E[X^k]} for nonnegative integer values of \eqn{k <
    \gamma}{k < shape}.

  The \eqn{k}th limited moment at some limit \eqn{d} is \eqn{E[\min(X,
  d)^k]}{E[min(X, d)^k]} for nonnegative integer values of \eqn{k}
  and \eqn{1 - j/\gamma}{1 - j/shape}, \eqn{j = 1, \dots, k}
  not a negative integer.
}
\value{
  \code{dpareto3} gives the density,
  \code{ppareto3} gives the distribution function,
  \code{qpareto3} gives the quantile function,
  \code{rpareto3} generates random deviates,
  \code{mpareto3} gives the \eqn{k}th raw moment, and
  \code{levpareto3} gives the \eqn{k}th moment of the limited loss
  variable.

  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\note{
  \code{levpareto3} computes the limited expected value using
  \code{\link{betaint}}.

  For Pareto distributions, we use the classification of Arnold (2015)
  with the parametrization of Klugman et al. (2012).

  The \code{"distributions"} package vignette provides the
  interrelations between the continuous size distributions in
  \pkg{actuar} and the complete formulas underlying the above functions.
}
\references{
  Arnold, B.C. (2015), \emph{Pareto Distributions}, Second Edition, CRC
  Press.

  Kleiber, C. and Kotz, S. (2003), \emph{Statistical Size Distributions
  in Economics and Actuarial Sciences}, Wiley.

  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2012),
  \emph{Loss Models, From Data to Decisions, Fourth Edition}, Wiley.
}
\seealso{
  \code{\link{dllogis}} for the loglogistic distribution.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca}
}
\examples{
exp(dpareto3(1, min = 10, 3, 4, log = TRUE))
p <- (1:10)/10
ppareto3(qpareto3(p, min = 10, 2, 3), min = 10, 2, 3)

## mean
mpareto3(1, min = 10, 2, 3)

## case with 1 - order/shape > 0
levpareto3(20, min = 10, 2, 3, order = 1)

## case with 1 - order/shape < 0
levpareto3(20, min = 10, 2/3, 3, order = 1)
}
\keyword{distribution}

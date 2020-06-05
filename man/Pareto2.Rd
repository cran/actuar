\name{Pareto2}
\alias{Pareto2}
\alias{dpareto2}
\alias{ppareto2}
\alias{qpareto2}
\alias{rpareto2}
\alias{mpareto2}
\alias{levpareto2}
\title{The Pareto II Distribution}
\description{
  Density function, distribution function, quantile function, random generation,
  raw moments and limited moments for the Pareto II distribution with
  parameters \code{min}, \code{shape} and \code{scale}.
}
\usage{
dpareto2(x, min, shape, rate = 1, scale = 1/rate,
         log = FALSE)
ppareto2(q, min, shape, rate = 1, scale = 1/rate,
         lower.tail = TRUE, log.p = FALSE)
qpareto2(p, min, shape, rate = 1, scale = 1/rate,
         lower.tail = TRUE, log.p = FALSE)
rpareto2(n, min, shape, rate = 1, scale = 1/rate)
mpareto2(order, min, shape, rate = 1, scale = 1/rate)
levpareto2(limit, min, shape, rate = 1, scale = 1/rate,
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
  The Pareto II (or \dQuote{type II}) distribution with parameters
  \code{min} \eqn{= \mu}{= m},
  \code{shape} \eqn{= \alpha}{= a} and
  \code{scale} \eqn{= \theta}{= s} has density:
  \deqn{f(x) = \frac{\alpha}{%
      \theta [1 + (x - \mu)/\theta]^{\alpha + 1}}}{%
    f(x) = a/(s [1 + (x - m)/s]^(a + 1))}
  for \eqn{x > \mu}{x > m}, \eqn{-\infty < \mu < \infty}{-Inf < m < Inf},
  \eqn{\alpha > 0}{a > 0} and \eqn{\theta > 0}{s > 0}.

  The Pareto II is the distribution of the random variable
  \deqn{\mu + \theta \left(\frac{X}{1 - X}\right),}{%
    m + s X/(1 - X),}
  where \eqn{X} has a beta distribution with parameters \eqn{1} and
  \eqn{\alpha}{a}. It derives from the \link[=dfpareto]{Feller-Pareto}
  distribution with \eqn{\tau = \gamma = 1}{shape2 = shape3 = 1}.
  Setting \eqn{\mu = 0}{min = 0} yields the familiar
  \link[=dpareto]{Pareto} distribution.

  The \link[=dpareto1]{Pareto I} (or Single parameter Pareto)
  distribution is a special case of the Pareto II with \code{min ==
  scale}.

  The \eqn{k}th raw moment of the random variable \eqn{X} is
  \eqn{E[X^k]}{E[X^k]} for nonnegative integer values of \eqn{k <
  \alpha}{k < shape}.

  The \eqn{k}th limited moment at some limit \eqn{d} is \eqn{E[\min(X,
  d)^k]}{E[min(X, d)^k]} for nonnegative integer values of \eqn{k}
  and \eqn{\alpha - j}{shape1 - j}, \eqn{j = 1, \dots, k}
  not a negative integer.
}
\value{
  \code{dpareto2} gives the density,
  \code{ppareto2} gives the distribution function,
  \code{qpareto2} gives the quantile function,
  \code{rpareto2} generates random deviates,
  \code{mpareto2} gives the \eqn{k}th raw moment, and
  \code{levpareto2} gives the \eqn{k}th moment of the limited loss
  variable.

  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\note{
  \code{levpareto2} computes the limited expected value using
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
  \code{\link{dpareto}} for the Pareto distribution without a location
  parameter.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca}
}
\examples{
exp(dpareto2(1, min = 10, 3, 4, log = TRUE))
p <- (1:10)/10
ppareto2(qpareto2(p, min = 10, 2, 3), min = 10, 2, 3)

## variance
mpareto2(2, min = 10, 4, 1) - mpareto2(1, min = 10, 4, 1)^2

## case with shape - order > 0
levpareto2(10, min = 10, 3, scale = 1, order = 2)

## case with shape - order < 0
levpareto2(10, min = 10, 1.5, scale = 1, order = 2)
}
\keyword{distribution}

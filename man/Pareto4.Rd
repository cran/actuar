\name{Pareto4}
\alias{Pareto4}
\alias{dpareto4}
\alias{ppareto4}
\alias{qpareto4}
\alias{rpareto4}
\alias{mpareto4}
\alias{levpareto4}
\title{The Pareto IV Distribution}
\description{
  Density function, distribution function, quantile function, random generation,
  raw moments and limited moments for the Pareto IV distribution with
  parameters \code{min}, \code{shape1}, \code{shape2} and \code{scale}.
}
\usage{
dpareto4(x, min, shape1, shape2, rate = 1, scale = 1/rate,
         log = FALSE)
ppareto4(q, min, shape1, shape2, rate = 1, scale = 1/rate,
         lower.tail = TRUE, log.p = FALSE)
qpareto4(p, min, shape1, shape2, rate = 1, scale = 1/rate,
         lower.tail = TRUE, log.p = FALSE)
rpareto4(n, min, shape1, shape2, rate = 1, scale = 1/rate)
mpareto4(order, min, shape1, shape2, rate = 1, scale = 1/rate)
levpareto4(limit, min, shape1, shape2, rate = 1, scale = 1/rate,
           order = 1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{min}{lower bound of the support of the distribution.}
  \item{shape1, shape2, scale}{parameters. Must be strictly positive.}
  \item{rate}{an alternative way to specify the scale.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities/densities
    \eqn{p} are returned as \eqn{\log(p)}{log(p)}.}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}.}
  \item{order}{order of the moment.}
  \item{limit}{limit of the loss variable.}
}
\details{
  The Pareto IV (or \dQuote{type IV}) distribution with parameters
  \code{min} \eqn{= \mu}{= m},
  \code{shape1} \eqn{= \alpha}{= a},
  \code{shape2} \eqn{= \gamma}{= b} and
  \code{scale} \eqn{= \theta}{= s} has density:
  \deqn{f(x) = \frac{\alpha \gamma ((x - \mu)/\theta)^{\gamma - 1}}{%
      \theta [1 + ((x - \mu)/\theta)^\gamma]^{\alpha + 1}}}{%
    f(x) = (a b ((x - m)/s)^(b - 1))/(s [1 + ((x - m)/s)^b]^(a + 1))}
  for \eqn{x > \mu}{x > m}, \eqn{-\infty < \mu < \infty}{-Inf < m < Inf},
  \eqn{\alpha > 0}{a > 0}, \eqn{\gamma > 0}{b > 0}
  and \eqn{\theta > 0}{s > 0}.

  The Pareto IV is the distribution of the random variable
  \deqn{\mu + \theta \left(\frac{X}{1 - X}\right)^{1/\gamma},}{%
    m + s (X/(1 - X))^(1/b),}
  where \eqn{X} has a beta distribution with parameters \eqn{1} and
  \eqn{\alpha}{a}. It derives from the \link[=dfpareto]{Feller-Pareto}
  distribution with \eqn{\tau = 1}{shape3 = 1}. Setting \eqn{\mu = 0}{min
  = 0} yields the \link[=dburr]{Burr} distribution.

  The Pareto IV distribution also has the following direct special
  cases:
  \itemize{
    \item A \link[=dpareto3]{Pareto III} distribution when \code{shape1
      == 1};
    \item A \link[=dpareto2]{Pareto II} distribution when \code{shape1
      == 1}.
  }

  The \eqn{k}th raw moment of the random variable \eqn{X} is
  \eqn{E[X^k]}{E[X^k]} for nonnegative integer values of \eqn{k <
  \alpha\gamma}{k < shape1 * shape2}.

  The \eqn{k}th limited moment at some limit \eqn{d} is \eqn{E[\min(X,
  d)^k]}{E[min(X, d)^k]} for nonnegative integer values of \eqn{k}
  and \eqn{\alpha - j/\gamma}{shape1 - j/shape2}, \eqn{j = 1, \dots, k}
  not a negative integer.
}
\value{
  \code{dpareto4} gives the density,
  \code{ppareto4} gives the distribution function,
  \code{qpareto4} gives the quantile function,
  \code{rpareto4} generates random deviates,
  \code{mpareto4} gives the \eqn{k}th raw moment, and
  \code{levpareto4} gives the \eqn{k}th moment of the limited loss
  variable.

  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\note{
  \code{levpareto4} computes the limited expected value using
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
  \code{\link{dburr}} for the Burr distribution.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca}
}
\examples{
exp(dpareto4(1, min = 10, 2, 3, log = TRUE))
p <- (1:10)/10
ppareto4(qpareto4(p, min = 10, 2, 3, 2), min = 10, 2, 3, 2)

## variance
mpareto4(2, min = 10, 2, 3, 1) - mpareto4(1, min = 10, 2, 3, 1) ^ 2

## case with shape1 - order/shape2 > 0
levpareto4(10, min = 10, 2, 3, 1, order = 2)

## case with shape1 - order/shape2 < 0
levpareto4(10, min = 10, 1.5, 0.5, 1, order = 2)
}
\keyword{distribution}

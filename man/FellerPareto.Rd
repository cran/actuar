\name{FellerPareto}
\alias{FellerPareto}
\alias{dfpareto}
\alias{pfpareto}
\alias{qfpareto}
\alias{rfpareto}
\alias{mfpareto}
\alias{levfpareto}
\title{The Feller Pareto Distribution}
\description{
  Density function, distribution function, quantile function, random generation,
  raw moments and limited moments for the Feller Pareto distribution
  with parameters \code{min}, \code{shape1}, \code{shape2}, \code{shape3} and
  \code{scale}.
}
\usage{
dfpareto(x, min, shape1, shape2, shape3, rate = 1, scale = 1/rate,
        log = FALSE)
pfpareto(q, min, shape1, shape2, shape3, rate = 1, scale = 1/rate,
        lower.tail = TRUE, log.p = FALSE)
qfpareto(p, min, shape1, shape2, shape3, rate = 1, scale = 1/rate,
        lower.tail = TRUE, log.p = FALSE)
rfpareto(n, min, shape1, shape2, shape3, rate = 1, scale = 1/rate)
mfpareto(order, min, shape1, shape2, shape3, rate = 1, scale = 1/rate)
levfpareto(limit, min, shape1, shape2, shape3, rate = 1, scale = 1/rate,
          order = 1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{min}{lower bound of the support of the distribution.}
  \item{shape1, shape2, shape3, scale}{parameters. Must be strictly
    positive.}
  \item{rate}{an alternative way to specify the scale.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities/densities
    \eqn{p} are returned as \eqn{\log(p)}{log(p)}.}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}.}
  \item{order}{order of the moment.}
  \item{limit}{limit of the loss variable.}
}
\details{
  The Feller-Pareto distribution with parameters \code{min} \eqn{= \mu}{= m},
  \code{shape1} \eqn{= \alpha}{= a}, \code{shape2} \eqn{= \gamma}{= b},
  \code{shape3} \eqn{= \tau}{= c} and \code{scale} \eqn{= \theta}{= s}, has
  density:
  \deqn{f(x) = \frac{\Gamma(\alpha + \tau)}{\Gamma(\alpha)\Gamma(\tau)}
    \frac{\gamma ((x - \mu)/\theta)^{\gamma \tau - 1}}{%
      \theta [1 + ((x - \mu)/\theta)^\gamma]^{\alpha + \tau}}}{%
    f(x) = Gamma(a + c)/(Gamma(a) * Gamma(c)) (b ((x - m)/s)^(bc - 1))/%
    (s [1 + ((x - m)/s)^b]^(a + c))}
  for \eqn{x > \mu}{x > m}, \eqn{-\infty < \mu < \infty}{-Inf < m < Inf},
  \eqn{\alpha > 0}{a > 0}, \eqn{\gamma > 0}{b > 0},
  \eqn{\tau > 0}{c > 0} and \eqn{\theta > 0}{s > 0}.
  (Here \eqn{\Gamma(\alpha)}{Gamma(a)} is the function implemented
  by \R's \code{\link{gamma}()} and defined in its help.)

  The Feller-Pareto is the distribution of the random variable
  \deqn{\mu + \theta \left(\frac{1 - X}{X}\right)^{1/\gamma},}{%
   m + s ((1 - X)/X)^(1/b),}
  where \eqn{X} has a beta distribution with parameters \eqn{\alpha}{a}
  and \eqn{\tau}{c}.

  The Feller-Pareto defines a large family of distributions encompassing
  the transformed beta family and many variants of the Pareto
  distribution. Setting \eqn{\mu = 0}{min = 0} yields the
  \link[=dtrbeta]{transformed beta} distribution.

  The Feller-Pareto distribution also has the following direct special
  cases:
  \itemize{
    \item A \link[=dpareto4]{Pareto IV} distribution when \code{shape3
      == 1};
    \item A \link[=dpareto3]{Pareto III} distribution when \code{shape1
      shape3 == 1};
    \item A \link[=dpareto2]{Pareto II} distribution when \code{shape1
      shape2 == 1};
    \item A \link[=dpareto1]{Pareto I} distribution when \code{shape1
      shape2 == 1} and \code{min = scale}.
  }

  The \eqn{k}th raw moment of the random variable \eqn{X} is
  \eqn{E[X^k]}{E[X^k]} for nonnegative integer values of \eqn{k <
  \alpha\gamma}{k < shape1 * shape2}.

  The \eqn{k}th limited moment at some limit \eqn{d} is \eqn{E[\min(X,
  d)^k]}{E[min(X, d)^k]} for nonnegative integer values of \eqn{k}
  and \eqn{\alpha - j/\gamma}{shape1 - j/shape2}, \eqn{j = 1, \dots, k}
  not a negative integer.

  Note that the range of admissible values for \eqn{k} in raw and
  limited moments is larger when \eqn{\mu = 0}{min == 0}.
}
\value{
  \code{dfpareto} gives the density,
  \code{pfpareto} gives the distribution function,
  \code{qfpareto} gives the quantile function,
  \code{rfpareto} generates random deviates,
  \code{mfpareto} gives the \eqn{k}th raw moment, and
  \code{levfpareto} gives the \eqn{k}th moment of the limited loss
  variable.

  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\note{
  \code{levfpareto} computes the limited expected value using
  \code{\link{betaint}}.

  For the Feller-Pareto and other Pareto distributions, we use the
  classification of Arnold (2015) with the parametrization of Klugman et
  al. (2012).  

  The \code{"distributions"} package vignette provides the
  interrelations between the continuous size distributions in
  \pkg{actuar} and the complete formulas underlying the above functions.
}
\references{
  Dutang, C., Goulet, V., Langevin, N. (2022).
  Feller-Pareto and Related Distributions: Numerical Implementation and Actuarial
  Applications.
  \emph{Journal of Statistical Software}, \bold{103}(6), 1--22.
  \doi{10.18637/jss.v103.i06}.

  Abramowitz, M. and Stegun, I. A. (1972), \emph{Handbook of 
  Mathematical Functions}, Dover.
  
  Arnold, B. C. (2015), \emph{Pareto Distributions}, Second Edition, CRC
  Press.

  Kleiber, C. and Kotz, S. (2003), \emph{Statistical Size Distributions
  in Economics and Actuarial Sciences}, Wiley.

  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2012),
  \emph{Loss Models, From Data to Decisions, Fourth Edition}, Wiley.
}
\seealso{
  \code{\link{dtrbeta}} for the transformed beta distribution.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Nicholas Langevin
}
\examples{
exp(dfpareto(2, 1, 2, 3, 4, 5, log = TRUE))
p <- (1:10)/10
pfpareto(qfpareto(p, 1, 2, 3, 4, 5), 1, 2, 3, 4, 5)

## variance
mfpareto(2, 1, 2, 3, 4, 5) - mfpareto(1, 1, 2, 3, 4, 5)^2

## case with shape1 - order/shape2 > 0
levfpareto(10, 1, 2, 3, 4, scale = 1, order = 2)

## case with shape1 - order/shape2 < 0
levfpareto(20, 10, 0.1, 14, 2, scale = 1.5, order = 2)
}
\keyword{distribution}

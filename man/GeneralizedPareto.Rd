\name{GeneralizedPareto}
\alias{GeneralizedPareto}
\alias{dgenpareto}
\alias{pgenpareto}
\alias{qgenpareto}
\alias{rgenpareto}
\alias{mgenpareto}
\alias{levgenpareto}
\title{The Generalized Pareto Distribution}
\description{
  Density function, distribution function, quantile function, random generation,
  raw moments and limited moments for the Generalized Pareto
  distribution with parameters \code{shape1}, \code{shape2} and
  \code{scale}.
}
\usage{
dgenpareto(x, shape1, shape2, rate = 1, scale = 1/rate,
           log = FALSE)
pgenpareto(q, shape1, shape2, rate = 1, scale = 1/rate,
           lower.tail = TRUE, log.p = FALSE)
qgenpareto(p, shape1, shape2, rate = 1, scale = 1/rate,
           lower.tail = TRUE, log.p = FALSE)
rgenpareto(n, shape1, shape2, rate = 1, scale = 1/rate)
mgenpareto(order, shape1, shape2, rate = 1, scale = 1/rate)
levgenpareto(limit, shape1, shape2, rate = 1, scale = 1/rate,
             order = 1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
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
  The Generalized Pareto distribution with parameters \code{shape1}
  \eqn{= \alpha}{= a}, \code{shape2} \eqn{= \tau}{= b} and \code{scale}
  \eqn{= \theta}{= s} has density:
  \deqn{f(x) = \frac{\Gamma(\alpha + \tau)}{\Gamma(\alpha)\Gamma(\tau)}
    \frac{\theta^\alpha x^{\tau - 1}}{%
      (x + \theta)^{\alpha + \tau}}}{%
    f(x) = Gamma(a + b)/(Gamma(a) * Gamma(b)) *
           (s^a x^(b - 1))/(x + s)^(a + b)}
  for \eqn{x > 0}, \eqn{\alpha > 0}{a > 0}, \eqn{\tau > 0}{b > 0} and
  \eqn{\theta > 0}{s > 0}.
  (Here \eqn{\Gamma(\alpha)} is the function implemented
  by \R's \code{\link{gamma}()} and defined in its help.)

  The Generalized Pareto is the distribution of the random variable
  \deqn{\theta \left(\frac{X}{1 - X}\right),}{\theta (X/(1 - X)),}
  where \eqn{X} has a beta distribution with parameters \eqn{\alpha}
  and \eqn{\tau}.

  The Generalized Pareto distribution has the following special cases:
  \itemize{
    \item A \link[=dpareto]{Pareto} distribution when \code{shape2 ==
      1};
    \item An \link[=dinvpareto]{Inverse Pareto} distribution when
    \code{shape1 == 1}.
  }

  The \eqn{k}th raw moment of the random variable \eqn{X} is
  \eqn{E[X^k]}, \eqn{-\tau < k < \alpha}{-shape2 < k < shape1}.

  The \eqn{k}th limited moment at some limit
  \eqn{d} is \eqn{E[\min(X, d)^k]}{E[min(X, d)^k]},
  \eqn{k > -\tau}{k > -shape2} and \eqn{\alpha - k}{shape1 - k} not a
  negative integer.
}
\value{
  \code{dgenpareto} gives the density,
  \code{pgenpareto} gives the distribution function,
  \code{qgenpareto} gives the quantile function,
  \code{rgenpareto} generates random deviates,
  \code{mgenpareto} gives the \eqn{k}th raw moment, and
  \code{levgenpareto} gives the \eqn{k}th moment of the limited loss
  variable.

  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\note{
  \code{levgenpareto} computes the limited expected value using
  \code{\link{betaint}}.

  Distribution also known as the Beta of the Second Kind. See also
  Kleiber and Kotz (2003) for alternative names and parametrizations.

  The Generalized Pareto distribution defined here is different from the
  one in Embrechts et al. (1997) and in
  \href{https://en.wikipedia.org/wiki/Generalized_Pareto_distribution}{Wikipedia};
  see also Kleiber and Kotz (2003, section 3.12). One may most likely
  compute quantities for the latter using functions for the
  \link[=dpareto]{Pareto} distribution with the appropriate change of
  parametrization.

  The \code{"distributions"} package vignette provides the
  interrelations between the continuous size distributions in
  \pkg{actuar} and the complete formulas underlying the above functions.
}
\references{
  Embrechts, P., Klüppelberg, C. and Mikisch, T. (1997), \emph{Modelling
    Extremal Events for Insurance and Finance}, Springer.
  
  Kleiber, C. and Kotz, S. (2003), \emph{Statistical Size Distributions
  in Economics and Actuarial Sciences}, Wiley.

  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2012),
  \emph{Loss Models, From Data to Decisions, Fourth Edition}, Wiley.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Mathieu Pigeon
}
\examples{
exp(dgenpareto(3, 3, 4, 4, log = TRUE))
p <- (1:10)/10
pgenpareto(qgenpareto(p, 3, 3, 1), 3, 3, 1)
qgenpareto(.3, 3, 4, 4, lower.tail = FALSE)

## variance
mgenpareto(2, 3, 2, 1) - mgenpareto(1, 3, 2, 1)^2

## case with shape1 - order > 0
levgenpareto(10, 3, 3, scale = 1, order = 2)

## case with shape1 - order < 0
levgenpareto(10, 1.5, 3, scale = 1, order = 2)
}
\keyword{distribution}

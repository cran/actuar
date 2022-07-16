\name{actuar-package}
\alias{actuar-package}
\alias{actuar}
\docType{package}
\title{\packageTitle{actuar}}
\description{\packageDescription{actuar}}
\details{
  \pkg{actuar} provides additional actuarial science
  functionality and support for heavy tailed distributions to the \R
  statistical system.

  The current feature set of the package can be split into five main
  categories.
  \enumerate{
    \item{Additional probability distributions: 23 continuous heavy
      tailed distributions from the Feller-Pareto and Transformed
      Gamma families, the loggamma, the Gumbel, the inverse Gaussian
      and the generalized beta; phase-type distributions; the
      Poisson-inverse Gaussian discrete distribution; zero-truncated and
      zero-modified extensions of the standard discrete distributions;
      computation of raw moments, limited moments and the moment
      generating function (when it exists) of continuous distributions.
      See the \dQuote{distributions} package vignette for details.}
    \item{Loss distributions modeling: extensive support of grouped
      data; functions to compute empirical raw and limited moments;
      support for minimum distance estimation using three different
      measures; treatment of coverage modifications (deductibles,
      limits, inflation, coinsurance). See the \dQuote{modeling} and
      \dQuote{coverage} package vignettes for details.}
    \item{Risk and ruin theory: discretization of the claim amount
      distribution; calculation of the aggregate claim amount
      distribution; calculation of the adjustment coefficient;
      calculation of the probability of ruin, including using phase-type
      distributions. See the \dQuote{risk} package vignette for details.}
    \item{Simulation of discrete mixtures, compound models (including
      the compound Poisson), and compound hierarchical models. See the
      \dQuote{simulation} package vignette for details.}
    \item{Credibility theory: function \code{cm} fits hierarchical
      (including Bühlmann, Bühlmann-Straub), regression and linear Bayes
      credibility models. See the \dQuote{credibility} package vignette
      for details.}
  }
}
\author{
  Christophe Dutang, Vincent Goulet, Mathieu Pigeon and many other
  contributors; use \code{packageDescription("actuar")} for the complete
  list.

  Maintainer: Vincent Goulet.
}
\references{
  Dutang, C., Goulet, V. and Pigeon, M. (2008).
  actuar: An R Package for Actuarial Science.
  \emph{Journal of Statistical Software}, \bold{25}(7), 1--37.
  \doi{10.18637/jss.v025.i07}.
  
  Dutang, C., Goulet, V., Langevin, N. (2022).
  Feller-Pareto and Related Distributions: Numerical Implementation and Actuarial
  Applications.
  \emph{Journal of Statistical Software}, \bold{103}(6), 1--22.
  \doi{10.18637/jss.v103.i06}.

  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2012),
  \emph{Loss Models, From Data to Decisions, Fourth Edition}, Wiley.
}
\seealso{
  For probability distributions support functions, use as starting
  points:
  \code{\link{FellerPareto}}, 
  \code{\link{TransformedGamma}}, 
  \code{\link{Loggamma}}, 
  \code{\link{Gumbel}}, 
  \code{\link{InverseGaussian}}, 
  \code{\link{PhaseType}}, 
  \code{\link{PoissonInverseGaussian}} and, e.g.,
  \code{\link{ZeroTruncatedPoisson}}, 
  \code{\link{ZeroModifiedPoisson}}.

  For loss modeling support functions:
  \code{\link{grouped.data}},
  \code{\link{ogive}},
  \code{\link{emm}},
  \code{\link{elev}},
  \code{\link{mde}},
  \code{\link{coverage}}.

  For risk and ruin theory functions:
  \code{\link{discretize}},
  \code{\link{aggregateDist}},
  \code{\link{adjCoef}},
  \code{\link{ruin}}.

  For credibility theory functions and datasets:
  \code{\link{cm}},
  \code{\link{hachemeister}}.
}
\examples{
## The package comes with extensive demonstration scripts;
## use the following command to obtain the list.
\dontrun{demo(package = "actuar")}
}
\keyword{package}
\keyword{distribution}
\keyword{models}
\keyword{univar}

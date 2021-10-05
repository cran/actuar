\name{var}
\alias{var}
\alias{var.default}
\alias{var.grouped.data}
\alias{sd}
\alias{sd.default}
\alias{sd.grouped.data}
\title{Variance and Standard Deviation}
\description{
  Generic functions for the variance and standard deviation, and methods
  for individual and grouped data.

  The default methods for individual data are the functions from the
  \pkg{stats} package.
}
\usage{
var(x, \dots)

\method{var}{default}(x, y = NULL, na.rm = FALSE, use, \dots)

\method{var}{grouped.data}(x, \dots)

sd(x, \dots)

\method{sd}{default}(x, na.rm = FALSE, \dots)

\method{sd}{grouped.data}(x, \dots)
}
\arguments{
  \item{x}{a vector or matrix of individual data, or an object of class
    \code{"grouped data"}.}
  \item{y}{see \code{\link[stats:var]{stats::var}}.}
  \item{na.rm}{see \code{\link[stats:var]{stats::var}}.}
  \item{use}{see \code{\link[stats:var]{stats::var}}.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  This page documents variance and standard deviation computations for
  grouped data. For individual data, see \code{\link[stats]{var}} and
  \code{\link[stats]{sd}} from the \pkg{stats} package.

  For grouped data with group boundaries \eqn{c_0, c_1, \dots,
  c_r}{c[0], c[1], \dots, c[r]} and group frequencies \eqn{n_1, \dots,
  n_r}{n[1], \dots, n[r]}, \code{var} computes the sample variance
  \deqn{\frac{1}{n - 1} \sum_{j = 1}^r n_j (a_j - m_1)^2,}{%
    (1/(n - 1)) * sum(j; n[j] * (a[j] - m)^2,}
  where
  \eqn{a_j = (c_{j - 1} + c_j)/2}{a[j] = (c[j - 1] + c[j])/2}
  is the midpoint of the \eqn{j}th interval,
  \eqn{m_1}{m} is the sample mean (or sample first moment) of the data,
  and
  \eqn{n = \sum_{j = 1}^r n_j}{n = sum(j; n[j])}.
  The sample sample standard deviation is the square root of the sample
  variance.

  The sample variance for grouped data differs from the variance
  computed from the empirical raw moments with \code{\link{emm}} in two
  aspects. First, it takes into account the degrees of freedom. Second,
  it applies Sheppard's correction factor to compensate for the
  overestimation of the true variation in the data. For groups of equal
  width \eqn{k}, Sheppard's correction factor is equal to \eqn{-k^2/12}.
}
\value{
  A named vector of variances or standard deviations.
}
\seealso{
  \code{\link{grouped.data}} to create grouped data objects;
  \code{\link{mean.grouped.data}} for the mean and \code{\link{emm}} for
  higher moments.
}
\references{
  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (1998),
  \emph{Loss Models, From Data to Decisions}, Wiley.

  Heumann, C., Schomaker, M., Shalabh (2016), \emph{Introduction to
  Statistics and Data Analysis}, Springer.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca}. Variance and
  standard deviation methods for grouped data contributed by Walter
  Garcia-Fontes \email{walter.garcia@upf.edu}.
}
\examples{
data(gdental)
var(gdental)
sd(gdental)

## Illustration of Sheppard's correction factor
cj <- c(0, 2, 4, 6, 8)
nj <- c(1, 5,  3,  2)
gd <- grouped.data(Group = cj, Frequency = nj)
(sum(nj) - 1)/sum(nj) * var(gd)
(emm(gd, 2) - emm(gd)^2) - 4/12
}
\keyword{univar}

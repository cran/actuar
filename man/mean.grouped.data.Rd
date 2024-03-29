\name{mean.grouped.data}
\alias{mean.grouped.data}
\title{Arithmetic Mean}
\description{
  Mean of grouped data objects.
}
\usage{
\method{mean}{grouped.data}(x, \dots)
}
\arguments{
  \item{x}{an object of class \code{"grouped.data"}.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  The mean of grouped data with group boundaries \eqn{c_0, c_1, \dots,
    c_r}{c[0], c[1], \dots, c[r]} and group frequencies \eqn{n_1, \dots,
    n_r}{n[1], \dots, n[r]} is
  \deqn{\frac{1}{n} \sum_{j = 1}^r a_j n_j,}{%
    (1/n) * sum(j; a[j] * n[j]),}
  where
  \eqn{a_j = (c_{j - 1} + c_j)/2}{a[j] = (c[j - 1] + c[j])/2}
  is the midpoint of the \eqn{j}th interval, and
  \eqn{n = \sum_{j = 1}^r n_j}{n = sum(j; n[j])}.
}
\value{
  A named vector of means.
}
\seealso{
  \code{\link{grouped.data}} to create grouped data objects;
  \code{\link{emm}} to compute higher moments.
}
\references{
  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (1998),
  \emph{Loss Models, From Data to Decisions}, Wiley.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca}
}
\examples{
data(gdental)
mean(gdental)
}
\keyword{univar}

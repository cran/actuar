\name{NEWS}
\title{\pkg{actuar} News}
\encoding{UTF-8}

\section{CHANGES IN \pkg{actuar} VERSION 3.1-4}{
  \subsection{BUG FIXES}{
    \itemize{
      \item{Due to its use of \code{log1mexp} since the previous
	release, the package depends on R >= 4.1.0.}
    }
  }
}

\section{CHANGES IN \pkg{actuar} VERSION 3.1-3}{
  \subsection{BUG FIXES}{
    \itemize{
      \item{Carry over the new implementation of the Cornish-Fisher
	Expansion of base R used by \code{qlogarithmic} and
	\code{qpoisinvgauss}.}
      \item{Fix computation of \code{[pq]zmpois}, \code{[pq]zmbinom} and
	\code{[pq]zmnbinom} following fixes to the underlying base R
	functions introduced in r80271 of R sources. With thanks to
	B.D. Ripley and Martin Maechler.}
    }
  }
}

\section{CHANGES IN \pkg{actuar} VERSION 3.1-2}{
  \subsection{BUG FIXES}{
    \itemize{
      \item{\code{qinvgauss} now returns a finite value when
	\eqn{\code{1.5/shape} > 1000}. Thanks to Bettina Grün
	\email{bettina.gruen@wu.ac.at} for the fix.}
      \item{A protection against rounding errors now ensures that
	\code{qzmlogarithmic(1 - pzmlogarithmic(x), lower.tail = FALSE)
	== x} is always \code{TRUE}.}
      \item{In \code{?dburr}, the scale parameter appeared in the
	denominator of the density instead of \eqn{x}. Thanks to Etienne
	Guy for the heads up.}
      \item{The package tests now correctly use \code{stopifnot} with
	argument \code{exprs} explicitly named.}
      \item{The formula for the moment of order \eqn{k} for grouped data
	in \code{?emm} fixed in version 2.3-3 for the LaTeX version is
	now also fixed for the text version. Thanks (again) to Walter
	Garcia-Fontes.}
    }
  }
}

\section{CHANGES IN \pkg{actuar} VERSION 3.1-1}{
  \subsection{BUG FIXES}{
    \itemize{
      \item{\code{rcompound} and \code{rmixture} now correctly find
	objects defined higher in the call stack.}
    }
  }
}

\section{CHANGES IN \pkg{actuar} VERSION 3.1-0}{
  \subsection{BUG FIXES}{
    \itemize{
      \item{\code{rmixture} now randomly shuffles the variates by
	default and gains an argument \code{shuffle} (\code{TRUE} by
	default). Using \code{shuffle = FALSE} restores the previous
	behaviour where the output vector contains all the random
	variates from the first model, then all the random variates from
	the second model, and so on. When the order of the random
	variates is irrelevant, this cuts execution time roughly in
	half. Thanks to Adam Kałdus \email{akaldus@wp.pl} for the
	stimulating comments on this matter.}
    }
  }
  \subsection{USER VISIBLE CHANGES}{
    \itemize{
      \item{The number of variates returned by \code{rmixture} is now
	the length of argument \code{n} if larger than 1, like other
	\code{r<dist>} functions.}
      \item{\code{rmixture} now checks the validity of its arguments.}
    }
  }
}

\section{CHANGES IN \pkg{actuar} VERSION 3.0-0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item{Support functions \code{[dpqrm,lev]fpareto} for the
	Feller-Pareto distribution and related Pareto distributions with
	a location parameter. The Feller-Pareto defines a large family
	of distributions encompassing the transformed beta family and
	many variants of the Pareto distribution. Using the nomenclature
	of Arnold (2015), the following distributions are now supported
	by \pkg{actuar}: Feller-Pareto, Pareto IV, Pareto III, and
	Pareto II. The Pareto I was already supported under the name
	Single Parameter Pareto. Contributed by Christophe Dutang,
	Vincent Goulet and Nicholas Langevin.}
      \item{The package now exposes through an API its 200+ C routines
	for probability functions and the beta integral. This is
	documented in a new section of the \dQuote{distributions}
	package vignette. See file \file{include/actuarAPI.h} in the
	package installation directory for the complete list of exported
	routines.}
      \item{Improvements to the accuracy in the right tail of the
	\code{p<dist>} and \code{lev<dist>} functions for most
	probability distributions of the transformed beta family.
	Achieved by replacing \code{pbeta(u, a, b, lower.tail)} for
	\eqn{u > 0.5} with \code{pbeta(1 - u, b, a, !lower.tail)} and an
	accurate computation of \code{u}. Contributed by Nicholas
	Langevin.}
      \item{The C workhorse \code{betaint_raw} behind \code{betaint}
	gains an additional argument to receive an accurate value of
	\eqn{1 - x}. Used extensively to improve accuracy of the
	\code{lev<dist>} functions for the transformed beta family.
	Contributed by Nicholas Langevin.}
      \item{The \dQuote{distributions} package vignette now regroups
	distributions of the transformed beta families and the
	single parameter Pareto under the umbrella of the Feller-Pareto
	family of distributions. The vignette now also includes diagrams
	showing the interrelations between the members of this family,
	as well as between the members of the transformed gamma and
	inverse transformed gamma families.}
      \item{Exhaustive regression tests for probability functions.}
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item{Improvements to the simulation algorithm for zero-modified
	discrete distributions in the \eqn{p_0^M < p_0}{p0m < p0} case.
	Contributed by Nicholas Langevin.}
      \item{\code{dpoisinvgauss} no longer returns \code{NaN} for large
	values of \code{x}. Solved by computing probabilities
	recursively instead of by calling \code{bessel_k} (the latter
	would overflow for large \code{nu} and propagate \code{NaN}).
	Computations are actually about twice as fast.}
      \item{\code{ppoisinvgauss} now honors argument \code{lower_tail}.}
      \item{\code{qpoisinvgauss} no longer fails with \code{mu = Inf}
	and \code{log.p = TRUE}.}
      \item{\code{betaint(x, Inf, b)} now returns \code{Inf} instead of
	\code{NaN}.}
      \item{\code{betaint(.Machine$double.xmin, a, b)}, with \eqn{b <
	  0}, now returns 0 instead of \code{NaN}.}
      \item{\code{d<dist>} and \code{p<dist>} functions for all
	continuous size distributions now handle limiting cases for
	infinite scale parameter, or for zero non-scale parameters,
	consistently with functions of base R.

	Affected functions are:
	\code{[dp]trbeta},
	\code{[dp]burr},
	\code{[dp]llogis},
	\code{[dp]paralogis},
	\code{[dp]genpareto},
	\code{[dp]pareto},
	\code{[dp]invburr},
	\code{[dp]invpareto},
	\code{[dp]invparalogis} in the Transformed Beta family;
	\code{[dp]trgamma},
	\code{[dp]invtrgamma},
	\code{[dp]invgamma},
	\code{[dp]invweibull},
	\code{[dp]invexp} in the Transformed Gamma family;
	\code{[dp]lgamma},
	\code{[dp]gumbel},
	\code{[dp]invgauss},
	\code{[dp]genbeta}.}
      \item{\code{levinvexp} no longer returns \code{NaN} for finite
	order.}
    }
  }
  \subsection{BREAKING CHANGE}{
    \itemize{
      \item{Support for the Pareto II distributions comes from functions
	\code{[dpqrm,lev]pareto2}. These functions were \emph{aliases} to
	\code{[dpqrm,lev]pareto} in previous version of \pkg{actuar}.
	The new functions are \emph{not} backward compatible. Therefore,
	calls to the \code{*pareto2} functions of previous versions of
	\pkg{actuar} will return wrong results and should be replaced by
	calls to \code{*pareto} functions.}
    }
  }
  \subsection{DEFUNCT}{
    \itemize{
      \item{Functions \code{[m,lev,mgf]invGauss} that were deprecated in
	version 2.0-0.}
    }
  }
}

\note{
  \itemize{Older news can be found in files \file{NEWS.2.Rd} (2.x series),
    \file{NEWS.1.Rd} (1.x series) and \file{NEWS.0.Rd} (0.x series).}
}

### actuar: Actuarial Functions and Heavy Tailed Distributions
###
### Tests of functions for continuous and discrete probability
### distributions.
###
### Despite the name of the file, the tests are for [dpqrm,lev]<dist>
### functions (for continuous distributions):
###
###   d<dist>: density or probability mass
###   p<dist>: cumulative distribution
###   q<dist>: quantile
###   r<dist>: random number generation
###   m<dist>: moment
###   lev<dist>: limited moment
###
### Distributions are classified and sorted as in appendix A and
### appendix B of the 'distributions' package vignette.
###
### Inspired by (and some parts taken from) `tests/d-p-q-r-tests.R` in
### R sources.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca> with
###         indirect help from the R Core Team

## Load the package
library(actuar)
library(expint)                         # for gammainc

## Define a "local" version of the otherwise non-exported function
## 'betaint'.
betaint <- actuar:::betaint

## No warnings, unless explicitly asserted via tools::assertWarning.
options(warn = 2)
assertWarning <- tools::assertWarning

## Special values and utilities. Taken from `tests/d-p-q-r-tests.R`.
Meps <- .Machine$double.eps
xMax <- .Machine$double.xmax
xMin <- .Machine$double.xmin
All.eq <- function(x, y)
{
    all.equal.numeric(x, y, tolerance = 64 * .Machine$double.eps,
                      scale = max(0, mean(abs(x), na.rm = TRUE)))
}

if(!interactive())
    set.seed(123)

###
### CONTINUOUS DISTRIBUTIONS
###

##
## FELLER-PARETO AND PARETO II, III, IV DISTRIBUTIONS
##

## When reasonable, we also test consistency with the special cases
## min = 0:
##
##   Feller-Pareto -> Transformated beta
##   Pareto IV     -> Burr
##   Pareto III    -> Loglogistic
##   Pareto II     -> Pareto

## Density: first check that functions return 0 when scale = Inf, and
## when x = scale = Inf.
stopifnot(exprs = {
    dfpareto(c(42, Inf), min = 1, shape1 = 2, shape2 = 3, shape3 = 4, scale = Inf) == c(0, 0)
    dpareto4(c(42, Inf), min = 1, shape1 = 2, shape2 = 3, scale = Inf) == c(0, 0)
    dpareto3(c(42, Inf), min = 1, shape  = 3, scale = Inf) == c(0, 0)
    dpareto2(c(42, Inf), min = 1, shape  = 2, scale = Inf) == c(0, 0)
})

## Next test density functions for an array of standard values.
nshpar <- 3                     # (maximum) number of shape parameters
min <- round(rnorm(30, 2), 2)
shpar <- replicate(30, rlnorm(nshpar, 2), simplify = FALSE)
scpar <- rlnorm(30, 2)          # scale parameters
for (i in seq_along(min))
{
    m <- min[i]
    a <- shpar[[c(i, 1)]]; g <- shpar[[c(i, 2)]]; t <- shpar[[c(i, 3)]]
    Be <- beta(a, t)
    for (s in scpar)
    {
        x <- rfpareto(100, min = m, shape1 = a, shape2 = g, shape3 = t, scale = s)
        y <- (x - m)/s
        u <- 1/(1 + y^(-g))
        stopifnot(exprs = {
            all.equal(d1 <- dfpareto(x, min = m,
                                     shape1 = a, shape2 = g, shape3 = t,
                                     scale = s),
                      d2 <- dfpareto(y, min = 0,
                                     shape1 = a, shape2 = g, shape3 = t,
                                     scale = 1)/s,
                      tolerance = 1e-10)
            all.equal(d2,
                      dtrbeta(y,
                              shape1 = a, shape2 = g, shape3 = t,
                              scale = 1)/s,
                      tolerance = 1e-10)
            all.equal(d1,
                      g * y^(g*t - 1)/(s * Be * (1 + y^g)^(a + t)),
                      tolerance = 1e-10)
            all.equal(d1,
                      g * u^t * (1 - u)^a/((x - m) * Be),
                      tolerance = 1e-10)
        })
        x <- rpareto4(100, min = m, shape1 = a, shape2 = g, scale = s)
        y <- (x - m)/s
        u <- 1/(1 + y^g)
        stopifnot(exprs = {
            all.equal(d1 <- dpareto4(x, min = m,
                                     shape1 = a, shape2 = g,
                                     scale = s),
                      d2 <- dpareto4(y, min = 0,
                                     shape1 = a, shape2 = g,
                                     scale = 1)/s,
                      tolerance = 1e-10)
            all.equal(d2,
                      dburr(y,
                            shape1 = a, shape2 = g,
                            scale = 1)/s,
                      tolerance = 1e-10)
            all.equal(d1,
                      a * g * y^(g - 1)/(s * (1 + y^g)^(a + 1)),
                      tolerance = 1e-10)
            all.equal(d1,
                      a * g * u^a * (1 - u)/(x - m),
                      tolerance = 1e-10)
        })
        x <- rpareto3(100, min = m, shape = g, scale = s)
        y <- (x - m)/s
        u <- 1/(1 + y^(-g))
        stopifnot(exprs = {
            all.equal(d1 <- dpareto3(x, min = m,
                                     shape = g,
                                     scale = s),
                      d2 <- dpareto3(y, min = 0,
                                     shape = g,
                                     scale = 1)/s,
                      tolerance = 1e-10)
            all.equal(d2,
                      dllogis(y,
                              shape = g,
                              scale = 1)/s,
                      tolerance = 1e-10)
            all.equal(d1,
                      g * y^(g - 1)/(s * (1 + y^g)^2),
                      tolerance = 1e-10)
            all.equal(d1,
                      g * u * (1 - u)/(x - m),
                      tolerance = 1e-10)
        })
        x <- rpareto2(100, min = m, shape = a, scale = s)
        y <- (x - m)/s
        u <- 1/(1 + y)
        stopifnot(exprs = {
            all.equal(d1 <- dpareto2(x, min = m,
                                     shape = a,
                                     scale = s),
                      d2 <- dpareto2(y, min = 0,
                                     shape = a,
                                     scale = 1)/s,
                      tolerance = 1e-10)
            all.equal(d2,
                      dpareto(y,
                              shape = a,
                              scale = 1)/s,
                      tolerance = 1e-10)
            all.equal(d1,
                      a/(s * (1 + y)^(a + 1)),
                      tolerance = 1e-10)
            all.equal(d1,
                      a * u^a * (1 - u)/(x - m),
                      tolerance = 1e-10)
        })
    }
}

## Tests on the cumulative distribution function.
##
## Note: when shape1 = shape3 = 1, the underlying beta distribution is
## a uniform. Therefore, pfpareto(x, min, 1, shape2, 1, scale) should
## return the value of u = v/(1 + v), v = ((x - min)/scale)^shape2.
##
## x = 2/Meps = 2^53 (with min = 0, shape2 = scale = 1) is the value
## where the cdf would jump to 1 if we weren't using the trick to
## compute the cdf with pbeta(1 - u, ..., lower = FALSE).
scLrg <- 1e300 * c(0.5, 1, 2)
m <- rnorm(1)
stopifnot(exprs = {
    pfpareto(Inf,         min = 10,   1, 2, 3, scale = xMax) == 1
    pfpareto(2^53,        min = 0,    1, 1, 1, scale = 1) != 1
    pfpareto(2^53 + xMax, min = xMax, 1, 1, 1, scale = 1) != 1
    all.equal(pfpareto(xMin + m, min = m, 1, 1, 1, scale = 1), xMin)
    all.equal(y <- pfpareto(1e300 + m, min = m,
                            shape1 = 3, shape2 = rep(c(1, 2), each = length(scLrg)),
                            shape3 = 1,
                            scale = scLrg, log = TRUE),
              ptrbeta(1e300,
                      shape1 = 3, shape2 = rep(c(1, 2), each = length(scLrg)),
                      shape3 = 1,
                      scale = scLrg, log = TRUE))
    all.equal(y,
              c(pbeta(c(2/3, 1/2), 1, 3, lower.tail = TRUE,  log = TRUE),
                pbeta(2/3, 3, 1, lower.tail = FALSE, log = TRUE),
                pbeta(c(4/5, 1/2), 1, 3, lower.tail = TRUE,  log = TRUE),
                pbeta(4/5, 3, 1, lower.tail = FALSE, log = TRUE)))
})
stopifnot(exprs = {
    ppareto4(Inf,         min = 10,   1, 3, scale = xMax) == 1
    ppareto4(2^53,        min = 0,    1, 1, scale = 1) != 1
    ppareto4(2^53 + xMax, min = xMax, 1, 1, scale = 1) != 1
    all.equal(ppareto4(xMin + m, min = m, 1, 1, scale = 1), xMin)
    all.equal(y <- ppareto4(1e300 + m, min = m,
                            shape1 = 3, shape2 = rep(c(1, 2), each = length(scLrg)),
                            scale = scLrg, log = TRUE),
              pburr(1e300,
                    shape1 = 3, shape2 = rep(c(1, 2), each = length(scLrg)),
                    scale = scLrg, log = TRUE))
    all.equal(y,
              c(log(1 - c(1/3, 1/2, 2/3)^3),
                log(1 - c(1/5, 1/2, 4/5)^3)))
})
stopifnot(exprs = {
    ppareto3(Inf,         min = 10,   3, scale = xMax) == 1
    ppareto3(2^53,        min = 0,    1, scale = 1) != 1
    ppareto3(2^53 + xMax, min = xMax, 1, scale = 1) != 1
    all.equal(ppareto3(xMin + m, min = m, 1, scale = 1), xMin)
    all.equal(y <- ppareto3(1e300 + m, min = m,
                            shape = rep(c(1, 2), each = length(scLrg)),
                            scale = scLrg, log = TRUE),
              pllogis (1e300,
                       shape = rep(c(1, 2), each = length(scLrg)),
                       scale = scLrg, log = TRUE))
    all.equal(y,
              c(log(c(2/3, 1/2, 1/3)),
                log(c(4/5, 1/2, 1/5))))
})
stopifnot(exprs = {
    ppareto2(Inf,         min = 10,   3, scale = xMax) == 1
    ppareto2(2^53,        min = 0,    1, scale = 1) != 1
    ppareto2(2^53 + xMax, min = xMax, 1, scale = 1) != 1
    all.equal(ppareto2(xMin + m, min = m, 1, scale = 1), xMin)
    all.equal(y <- ppareto2(1e300 + m, min = m,
                            shape = 3,
                            scale = scLrg, log = TRUE),
              ppareto (1e300,
                       shape = 3,
                       scale = scLrg, log = TRUE))
    all.equal(y,
              c(log(1 - c(1/3, 1/2, 2/3)^3)))
})

## Also check that distribution functions return 0 when scale = Inf.
stopifnot(exprs = {
    pfpareto(x, min = m, shape1 = a, shape2 = g, shape3 = t, scale = Inf) == 0
    ppareto4(x, min = m, shape1 = a, shape2 = g, scale = Inf) == 0
    ppareto3(x, min = m, shape  = g, scale = Inf) == 0
    ppareto2(x, min = m, shape  = a, scale = Inf) == 0
})

## Tests for first three (positive) moments
##
## Simulation of new parameters ensuring that the first three moments
## exist.
set.seed(123)                   # reset the seed
nshpar <- 3                     # (maximum) number of shape parameters
min <- round(rnorm(30, 2), 2)
shpar <- replicate(30, c(3, 3, 0) + rlnorm(nshpar, 2), simplify = FALSE)
scpar <- rlnorm(30, 2)          # scale parameters
for (i in seq_along(min))
{
    m <- min[i]
    a <- shpar[[c(i, 1)]]; g <- shpar[[c(i, 2)]]; t <- shpar[[c(i, 3)]]
    Be <- beta(a, t)
    Ga <- gamma(a)
    for (s in scpar)
    {
        stopifnot(exprs = {
            All.eq(mfpareto(1, min = m,
                            shape1 = a, shape2 = g, shape3 = t,
                            scale = s),
                   m * (Be + (s/m) * beta(t + 1/g, a - 1/g))/Be)
            All.eq(mfpareto(2, min = m,
                            shape1 = a, shape2 = g, shape3 = t,
                            scale = s),
                   m^2 * (Be + 2 * (s/m) * beta(t + 1/g, a - 1/g)
                       + (s/m)^2 * beta(t + 2/g, a - 2/g))/Be)
            All.eq(mfpareto(3, min = m,
                            shape1 = a, shape2 = g, shape3 = t,
                            scale = s),
                   m^3 * (Be + 3 * (s/m) * beta(t + 1/g, a - 1/g)
                       + 3 * (s/m)^2 * beta(t + 2/g, a - 2/g)
                       + (s/m)^3 * beta(t + 3/g, a - 3/g))/Be)
        })
        stopifnot(exprs = {
            All.eq(mpareto4(1, min = m,
                            shape1 = a, shape2 = g,
                            scale = s),
                   m * (Ga + (s/m) * gamma(1 + 1/g) * gamma(a - 1/g))/Ga)
            All.eq(mpareto4(2, min = m,
                            shape1 = a, shape2 = g,
                            scale = s),
                   m^2 * (Ga + 2 * (s/m) * gamma(1 + 1/g) * gamma(a - 1/g)
                       + (s/m)^2 * gamma(1 + 2/g) * gamma(a - 2/g))/Ga)
            All.eq(mpareto4(3, min = m,
                            shape1 = a, shape2 = g,
                            scale = s),
                   m^3 * (Ga + 3 * (s/m) * gamma(1 + 1/g) * gamma(a - 1/g)
                       + 3 * (s/m)^2 * gamma(1 + 2/g) * gamma(a - 2/g)
                       + (s/m)^3 * gamma(1 + 3/g) * gamma(a - 3/g))/Ga)
        })
        stopifnot(exprs = {
            All.eq(mpareto3(1, min = m,
                            shape = g,
                            scale = s),
                   m * (1 + (s/m) * gamma(1 + 1/g) * gamma(1 - 1/g)))
            All.eq(mpareto3(2, min = m,
                            shape = g,
                            scale = s),
                   m^2 * (1 + 2 * (s/m) * gamma(1 + 1/g) * gamma(1 - 1/g)
                       + (s/m)^2 * gamma(1 + 2/g) * gamma(1 - 2/g)))
            All.eq(mpareto3(3, min = m,
                            shape = g,
                            scale = s),
                   m^3 * (1 + 3 * (s/m) * gamma(1 + 1/g) * gamma(1 - 1/g)
                       + 3 * (s/m)^2 * gamma(1 + 2/g) * gamma(1 - 2/g)
                       + (s/m)^3 * gamma(1 + 3/g) * gamma(1 - 3/g)))
        })
        stopifnot(exprs = {
            All.eq(mpareto2(1, min = m,
                            shape = a,
                            scale = s),
                   m * (Ga + (s/m) * gamma(1 + 1) * gamma(a - 1))/Ga)
            All.eq(mpareto2(2, min = m,
                            shape = a,
                            scale = s),
                   m^2 * (Ga + 2 * (s/m) * gamma(1 + 1) * gamma(a - 1)
                       + (s/m)^2 * gamma(1 + 2) * gamma(a - 2))/Ga)
            All.eq(mpareto2(3, min = m,
                            shape = a,
                            scale = s),
                   m^3 * (Ga + 3 * (s/m) * gamma(1 + 1) * gamma(a - 1)
                       + 3 * (s/m)^2 * gamma(1 + 2) * gamma(a - 2)
                       + (s/m)^3 * gamma(1 + 3) * gamma(a - 3))/Ga)
        })
    }
}

## Tests for first three limited moments
##
## Limits are taken from quantiles of each distribution.
q <- c(0.25, 0.50, 0.75, 0.9, 0.95)
for (i in seq_along(min))
{
    m <- min[i]
    a <- shpar[[c(i, 1)]]; g <- shpar[[c(i, 2)]]; t <- shpar[[c(i, 3)]]
    Ga <- gamma(a)
    Gt <- gamma(t)
    for (s in scpar)
    {
        limit <- qfpareto(q, min = m,
                          shape1 = a, shape2 = g, shape3 = t,
                          scale = s)
        y <- (limit - m)/s
        u <- 1/(1 + y^(-g))
        stopifnot(exprs = {
            All.eq(levfpareto(limit, order = 1, min = m,
                              shape1 = a, shape2 = g, shape3 = t,
                              scale = s),
                   m * (betaint(u, t, a) + (s/m) * betaint(u, t + 1/g, a - 1/g))/(Ga * Gt) +
                   limit * pbeta(u, t, a, lower = FALSE))
            All.eq(levfpareto(limit, order = 2, min = m,
                              shape1 = a, shape2 = g, shape3 = t,
                              scale = s),
                   m^2 * (betaint(u, t, a) + 2 * (s/m) * betaint(u, t + 1/g, a - 1/g)
                       + (s/m)^2 * betaint(u, t + 2/g, a - 2/g))/(Ga * Gt) +
                   limit^2 * pbeta(u, t, a, lower = FALSE))
            All.eq(levfpareto(limit, order = 3, min = m,
                              shape1 = a, shape2 = g, shape3 = t,
                              scale = s),
                   m^3 * (betaint(u, t, a) + 3 * (s/m) * betaint(u, t + 1/g, a - 1/g)
                       + 3 * (s/m)^2 * betaint(u, t + 2/g, a - 2/g)
                       + (s/m)^3 * betaint(u, t + 3/g, a - 3/g))/(Ga * Gt) +
                   limit^3 * pbeta(u, t, a, lower = FALSE))
        })
        limit <- qpareto4(q, min = m,
                          shape1 = a, shape2 = g,
                          scale = s)
        y <- (limit - m)/s
        u <- 1/(1 + y^g)
        u1m <- 1/(1 + y^(-g))
        stopifnot(exprs = {
            All.eq(levpareto4(limit, order = 1, min = m,
                              shape1 = a, shape2 = g,
                              scale = s),
                   m * (betaint(u1m, 1, a) + (s/m) * betaint(u1m, 1 + 1/g, a - 1/g))/Ga +
                   limit * u^a)
            All.eq(levpareto4(limit, order = 2, min = m,
                              shape1 = a, shape2 = g,
                              scale = s),
                   m^2 * (betaint(u1m, 1, a) + 2 * (s/m) * betaint(u1m, 1 + 1/g, a - 1/g)
                       + (s/m)^2 * betaint(u1m, 1 + 2/g, a - 2/g))/Ga +
                   limit^2 * u^a)
            All.eq(levpareto4(limit, order = 3, min = m,
                              shape1 = a, shape2 = g,
                              scale = s),
                   m^3 * (betaint(u1m, 1, a) + 3 * (s/m) * betaint(u1m, 1 + 1/g, a - 1/g)
                       + 3 * (s/m)^2 * betaint(u1m, 1 + 2/g, a - 2/g)
                       + (s/m)^3 * betaint(u1m, 1 + 3/g, a - 3/g))/Ga +
                   limit^3 * u^a)
        })
        limit <- qpareto3(q, min = m,
                          shape = g,
                          scale = s)
        y <- (limit - m)/s
        u <- 1/(1 + y^(-g))
        u1m <- 1/(1 + y^g)
        stopifnot(exprs = {
            All.eq(levpareto3(limit, order = 1, min = m,
                              shape = g,
                              scale = s),
                   m * (u + (s/m) * betaint(u, 1 + 1/g, 1 - 1/g)) +
                   limit * u1m)
            All.eq(levpareto3(limit, order = 2, min = m,
                              shape = g,
                              scale = s),
                   m^2 * (u + 2 * (s/m) * betaint(u, 1 + 1/g, 1 - 1/g)
                       + (s/m)^2 * betaint(u, 1 + 2/g, 1 - 2/g)) +
                   limit^2 * u1m)
            All.eq(levpareto3(limit, order = 3, min = m,
                              shape = g,
                              scale = s),
                   m^3 * (u + 3 * (s/m) * betaint(u, 1 + 1/g, 1 - 1/g)
                       + 3 * (s/m)^2 * betaint(u, 1 + 2/g, 1 - 2/g)
                       + (s/m)^3 * betaint(u, 1 + 3/g, 1 - 3/g)) +
                   limit^3 * u1m)
        })
        limit <- qpareto2(q, min = m,
                          shape = a,
                          scale = s)
        y <- (limit - m)/s
        u <- 1/(1 + y)
        u1m <- 1/(1 + y^(-1))
        stopifnot(exprs = {
            All.eq(levpareto2(limit, order = 1, min = m,
                              shape = a,
                              scale = s),
                   m * (betaint(u1m, 1, a) + (s/m) * betaint(u1m, 1 + 1, a - 1))/Ga +
                   limit * u^a)
            All.eq(levpareto2(limit, order = 2, min = m,
                              shape = a,
                              scale = s),
                   m^2 * (betaint(u1m, 1, a) + 2 * (s/m) * betaint(u1m, 1 + 1, a - 1)
                       + (s/m)^2 * betaint(u1m, 1 + 2, a - 2))/Ga +
                   limit^2 * u^a)
            All.eq(levpareto2(limit, order = 3, min = m,
                              shape = a,
                              scale = s),
                   m^3 * (betaint(u1m, 1, a) + 3 * (s/m) * betaint(u1m, 1 + 1, a - 1)
                       + 3 * (s/m)^2 * betaint(u1m, 1 + 2, a - 2)
                       + (s/m)^3 * betaint(u1m, 1 + 3, a - 3))/Ga +
                   limit^3 * u^a)
        })
    }
}

##
## TRANSFORMED BETA FAMILY
##

## Density: first check that functions return 0 when scale = Inf, and
## when x = scale = Inf.
stopifnot(exprs = {
    dtrbeta      (c(42, Inf), shape1 = 2, shape2 = 3, shape3 = 4, scale = Inf) == c(0, 0)
    dburr        (c(42, Inf), shape1 = 2, shape2 = 3, scale = Inf) == c(0, 0)
    dllogis      (c(42, Inf), shape  = 3, scale = Inf) == c(0, 0)
    dparalogis   (c(42, Inf), shape  = 2, scale = Inf) == c(0, 0)
    dgenpareto   (c(42, Inf), shape1 = 2, shape2 = 4, scale = Inf) == c(0, 0)
    dpareto      (c(42, Inf), shape  = 2, scale = Inf) == c(0, 0)
    dinvburr     (c(42, Inf), shape1 = 4, shape2 = 3, scale = Inf) == c(0, 0)
    dinvpareto   (c(42, Inf), shape  = 4, scale = Inf) == c(0, 0)
    dinvparalogis(c(42, Inf), shape  = 4, scale = Inf) == c(0, 0)
})

## Next test density functions for an array of standard values.
set.seed(123)                   # reset the seed
nshpar <- 3                     # (maximum) number of shape parameters
shpar <- replicate(30, rlnorm(nshpar, 2), simplify = FALSE)
scpar <- rlnorm(30, 2)          # scale parameters
for (i in seq_along(shpar))
{
    a <- shpar[[c(i, 1)]]; g <- shpar[[c(i, 2)]]; t <- shpar[[c(i, 3)]]
    Be <- beta(a, t)
    for (s in scpar)
    {
        x <- rtrbeta(100, shape1 = a, shape2 = g, shape3 = t, scale = s)
        y <- x/s
        u <- 1/(1 + y^(-g))
        stopifnot(exprs = {
            all.equal(d1 <- dtrbeta(x, shape1 = a, shape2 = g, shape3 = t,
                                    scale = s),
                      d2 <- dtrbeta(y, shape1 = a, shape2 = g, shape3 = t,
                                    scale = 1)/s,
                      tolerance = 1e-10)
            all.equal(d1,
                      g * y^(g*t - 1)/(s * Be * (1 + y^g)^(a + t)),
                      tolerance = 1e-10)
            all.equal(d1,
                      g * u^t * (1 - u)^a/(x * Be),
                      tolerance = 1e-10)
        })
        x <- rburr(100, shape1 = a, shape2 = g, scale = s)
        y <- x/s
        u <- 1/(1 + y^g)
        stopifnot(exprs = {
            all.equal(d1 <- dburr(x, shape1 = a, shape2 = g,
                                  scale = s),
                      d2 <- dburr(y, shape1 = a, shape2 = g,
                                  scale = 1)/s,
                      tolerance = 1e-10)
            all.equal(d1,
                      a * g * y^(g - 1)/(s * (1 + y^g)^(a + 1)),
                      tolerance = 1e-10)
            all.equal(d1,
                      a * g * u^a * (1 - u)/x,
                      tolerance = 1e-10)
        })
        x <- rllogis(100, shape = g, scale = s)
        y <- x/s
        u <- 1/(1 + y^(-g))
        stopifnot(exprs = {
            all.equal(d1 <- dllogis(x, shape = g,
                                    scale = s),
                      d2 <- dllogis(y, shape = g,
                                    scale = 1)/s,
                      tolerance = 1e-10)
            all.equal(d1,
                      g * y^(g - 1)/(s * (1 + y^g)^2),
                      tolerance = 1e-10)
            all.equal(d1,
                      g * u * (1 - u)/x,
                      tolerance = 1e-10)
        })
        x <- rparalogis(100, shape = a, scale = s)
        y <- x/s
        u <- 1/(1 + y^a)
        stopifnot(exprs = {
            all.equal(d1 <- dparalogis(x, shape = a,
                                       scale = s),
                      d2 <- dparalogis(y, shape = a,
                                       scale = 1)/s,
                      tolerance = 1e-10)
            all.equal(d1,
                      a^2 * y^(a - 1)/(s * (1 + y^a)^(a + 1)),
                      tolerance = 1e-10)
            all.equal(d1,
                      a^2 * u^a * (1 - u)/x,
                      tolerance = 1e-10)
        })
        x <- rgenpareto(100, shape1 = a, shape2 = t, scale = s)
        y <- x/s
        u <- 1/(1 + y^(-1))
        stopifnot(exprs = {
            all.equal(d1 <- dgenpareto(x, shape1 = a, shape2 = t,
                                       scale = s),
                      d2 <- dgenpareto(y, shape1 = a, shape2 = t,
                                       scale = 1)/s,
                      tolerance = 1e-10)
            all.equal(d1,
                      y^(t - 1)/(s * Be * (1 + y)^(a + t)),
                      tolerance = 1e-10)
            all.equal(d1,
                      u^t * (1 - u)^a/(x * Be),
                      tolerance = 1e-10)
        })
        x <- rpareto(100, shape = a, scale = s)
        y <- x/s
        u <- 1/(1 + y)
        stopifnot(exprs = {
            all.equal(d1 <- dpareto(x, shape = a,
                                    scale = s),
                      d2 <- dpareto(y, shape = a,
                                    scale = 1)/s,
                      tolerance = 1e-10)
            all.equal(d1,
                      a/(s * (1 + y)^(a + 1)),
                      tolerance = 1e-10)
            all.equal(d1,
                      a * u^a * (1 - u)/x,
                      tolerance = 1e-10)
        })
        x <- rpareto1(100, min = s, shape = a)
        stopifnot(exprs = {
            all.equal(d1 <- dpareto1(x, min = s, shape = a),
                      a * s^a/(x^(a + 1)),
                      tolerance = 1e-10)
        })
        x <- rinvburr(100, shape1 = t, shape2 = g, scale = s)
        y <- x/s
        u <- 1/(1 + y^(-g))
        stopifnot(exprs = {
            all.equal(d1 <- dinvburr(x, shape1 = t, shape2 = g,
                                     scale = s),
                      d2 <- dinvburr(y, shape1 = t, shape2 = g,
                                     scale = 1)/s,
                      tolerance = 1e-10)
            all.equal(d1,
                      t * g * y^(g*t - 1)/(s * (1 + y^g)^(t + 1)),
                      tolerance = 1e-10)
            all.equal(d1,
                      t * g * u^t * (1 - u)/x,
                      tolerance = 1e-10)
        })
        x <- rinvpareto(100, shape = t, scale = s)
        y <- x/s
        u <- 1/(1 + y^(-1))
        stopifnot(exprs = {
            all.equal(d1 <- dinvpareto(x, shape = t,
                                 scale = s),
                      d2 <- dinvpareto(y, shape = t,
                                       scale = 1)/s,
                      tolerance = 1e-10)
            all.equal(d1,
                      t * y^(t - 1)/(s * (1 + y)^(t + 1)),
                      tolerance = 1e-10)
            all.equal(d1,
                      t * u^t * (1 - u)/x,
                      tolerance = 1e-10)
        })
        x <- rinvparalogis(100, shape = t, scale = s)
        y <- x/s
        u <- 1/(1 + y^(-t))
        stopifnot(exprs = {
            all.equal(d1 <- dinvparalogis(x, shape = t,
                                          scale = s),
                      d2 <- dinvparalogis(y, shape = t,
                                          scale = 1)/s,
                      tolerance = 1e-10)
            all.equal(d1,
                      t^2 * y^(t^2 - 1)/(s * (1 + y^t)^(t + 1)),
                      tolerance = 1e-10)
            all.equal(d1,
                      t^2 * u^t * (1 - u)/x,
                      tolerance = 1e-10)
        })
    }
}

## Tests on the cumulative distribution function.
##
## Note: when shape1 = shape3 = 1, the underlying beta distribution is
## a uniform. Therefore, ptrbeta(x, 1, shape2, 1, scale) should return
## the value of u = v/(1 + v), v = (x/scale)^shape2.
##
## x = 2/Meps = 2^53 (with, shape2 = scale = 1) is the value where the
## cdf would jump to 1 if we weren't using the trick to compute the
## cdf with pbeta(1 - u, ..., lower = FALSE).
scLrg <- 1e300 * c(0.5, 1, 2)
stopifnot(exprs = {
    ptrbeta(Inf,  1, 2, 3, scale = xMax) == 1
    ptrbeta(2^53, 1, 1, 1, scale = 1) != 1
    all.equal(ptrbeta(xMin, 1, 1, 1, scale = 1), xMin)
    all.equal(ptrbeta(1e300,
                      shape1 = 3, shape2 = rep(c(1, 2), each = length(scLrg)),
                      shape3 = 1,
                      scale = scLrg, log = TRUE),
              c(pbeta(c(2/3, 1/2), 1, 3, lower.tail = TRUE,  log = TRUE),
                pbeta(2/3, 3, 1, lower.tail = FALSE, log = TRUE),
                pbeta(c(4/5, 1/2), 1, 3, lower.tail = TRUE,  log = TRUE),
                pbeta(4/5, 3, 1, lower.tail = FALSE, log = TRUE)))
})
stopifnot(exprs = {
    pburr(Inf,  1, 3, scale = xMax) == 1
    pburr(2^53, 1, 1, scale = 1) != 1
    all.equal(pburr(xMin, 1, 1, scale = 1), xMin)
    all.equal(pburr(1e300,
                    shape1 = 3, shape2 = rep(c(1, 2), each = length(scLrg)),
                    scale = scLrg, log = TRUE),
              c(log(1 - c(1/3, 1/2, 2/3)^3),
                log(1 - c(1/5, 1/2, 4/5)^3)))
})
stopifnot(exprs = {
    pllogis(Inf,  3, scale = xMax) == 1
    pllogis(2^53, 1, scale = 1) != 1
    all.equal(pllogis(xMin, 1, scale = 1), xMin)
    all.equal(pllogis(1e300,
                      shape = rep(c(1, 2), each = length(scLrg)),
                      scale = scLrg, log = TRUE),
              c(log(c(2/3, 1/2, 1/3)),
                log(c(4/5, 1/2, 1/5))))
})
stopifnot(exprs = {
    pparalogis(Inf,  3, scale = xMax) == 1
    pparalogis(2^53, 1, scale = 1) != 1
    all.equal(pparalogis(xMin, 1, scale = 1), xMin)
    all.equal(pparalogis(1e300,
                         shape = rep(c(2, 3), each = length(scLrg)),
                         scale = scLrg, log = TRUE),
              c(log(1 - c(1/5, 1/2, 4/5)^2),
                log(1 - c(1/9, 1/2, 8/9)^3)))
})
stopifnot(exprs = {
    pgenpareto(Inf,  1, 3, scale = xMax) == 1
    pgenpareto(2^53, 1, 1, scale = 1) != 1
    all.equal(pgenpareto(xMin, 1, 1, scale = 1), xMin)
    all.equal(pgenpareto(1e300,
                         shape1 = 3, shape2 = 1,
                         scale = scLrg, log = TRUE),
              c(pbeta(c(2/3, 1/2), 1, 3, lower.tail = TRUE,  log = TRUE),
                pbeta(2/3, 3, 1, lower.tail = FALSE, log = TRUE)))
})
stopifnot(exprs = {
    ppareto(Inf,  3, scale = xMax) == 1
    ppareto(2^53, 1, scale = 1) != 1
    all.equal(ppareto(xMin, 1, scale = 1), xMin)
    all.equal(ppareto(1e300,
                      shape = 3,
                      scale = scLrg, log = TRUE),
              c(log(1 - c(1/3, 1/2, 2/3)^3)))
})
stopifnot(exprs = {
    ppareto1(Inf,  3, min = xMax) == 1
    ppareto1(2^53, 1, min = 1) != 1
    all.equal(ppareto1(xMin, 1, min = 1), xMin)
    all.equal(ppareto1(1e300,
                      shape = 3,
                      min = 1e300 * c(0.001, 0.1, 0.5), log = TRUE),
              c(log(1 - c(0.001, 0.1, 0.5)^3)))
})
stopifnot(exprs = {
    pinvburr(Inf,  1, 3, scale = xMax) == 1
    pinvburr(2^53, 1, 1, scale = 1) != 1
    all.equal(pinvburr(xMin, 1, 1, scale = 1), xMin)
    all.equal(pinvburr(1e300,
                       shape1 = 3, shape2 = rep(c(1, 2), each = length(scLrg)),
                       scale = scLrg, log = TRUE),
              c(log(c(2/3, 1/2, 1/3)^3),
                log(c(4/5, 1/2, 1/5)^3)))
})
stopifnot(exprs = {
    pinvpareto(Inf,  3, scale = xMax) == 1
    pinvpareto(2^53, 1, scale = 1) != 1
    all.equal(pinvpareto(xMin, 1, scale = 1), xMin)
    all.equal(pinvpareto(1e300,
                         shape = 3,
                         scale = scLrg, log = TRUE),
              c(log(c(2/3, 1/2, 1/3)^3)))
})
stopifnot(exprs = {
    pinvparalogis(Inf,  3, scale = xMax) == 1
    pinvparalogis(2^53, 1, scale = 1) != 1
    all.equal(pinvparalogis(xMin, 1, scale = 1), xMin)
    all.equal(pinvparalogis(1e300,
                            shape = rep(c(2, 3), each = length(scLrg)),
                            scale = scLrg, log = TRUE),
              c(log(c(4/5, 1/2, 1/5)^2),
                log(c(8/9, 1/2, 1/9)^3)))
})

## Also check that distribution functions return 0 when scale = Inf.
stopifnot(exprs = {
    ptrbeta      (x, shape1 = a, shape2 = g, shape3 = t, scale = Inf) == 0
    pburr        (x, shape1 = a, shape2 = g, scale = Inf) == 0
    pllogis      (x, shape  = g, scale = Inf) == 0
    pparalogis   (x, shape  = a, scale = Inf) == 0
    pgenpareto   (x, shape1 = a, shape2 = t, scale = Inf) == 0
    ppareto      (x, shape  = a, scale = Inf) == 0
    pinvburr     (x, shape1 = t, shape2 = g, scale = Inf) == 0
    pinvpareto   (x, shape  = t, scale = Inf) == 0
    pinvparalogis(x, shape  = t, scale = Inf) == 0
})

## Tests for first three positive moments and first two negative
## moments.
##
## Simulation of new parameters ensuring that said moments exist.
set.seed(123)                   # reset the seed
nshpar <- 3                     # (maximum) number of shape parameters
shpar <- replicate(30, c(3, 3, 3) + rlnorm(nshpar, 2), simplify = FALSE)
scpar <- rlnorm(30, 2)          # scale parameters
k <- c(-2, -1, 1, 2, 3)         # orders
for (i in seq_along(shpar))
{
    a <- shpar[[c(i, 1)]]; g <- shpar[[c(i, 2)]]; t <- shpar[[c(i, 3)]]
    Be <- beta(a, t)
    Ga <- gamma(a)
    for (s in scpar)
    {
        stopifnot(exprs = {
            All.eq(mtrbeta(k, shape1 = a, shape2 = g, shape3 = t, scale = s),
                   s^k * beta(t + k/g, a - k/g)/Be)
            All.eq(mburr(k, shape1 = a, shape2 = g, scale = s),
                   s^k * gamma(1 + k/g) * gamma(a - k/g)/Ga)
            All.eq(mllogis(k, shape = g, scale = s),
                   s^k * gamma(1 + k/g) * gamma(1 - k/g))
            All.eq(mparalogis(k, shape = a, scale = s),
                   s^k * gamma(1 + k/a) * gamma(a - k/a)/Ga)
            All.eq(mgenpareto(k, shape1 = a, shape2 = t, scale = s),
                   s^k * beta(t + k, a - k)/Be)
            All.eq(mpareto(k[k > -1], shape = a, scale = s),
                   s^k[k > -1] * gamma(1 + k[k > -1]) * gamma(a - k[k > -1])/Ga)
            All.eq(mpareto1(k, shape = a, min = s),
                   s^k * a/(a - k))
            All.eq(minvburr(k, shape1 = a, shape2 = g, scale = s),
                   s^k * gamma(a + k/g) * gamma(1 - k/g)/Ga)
            All.eq(minvpareto(k[k < 1], shape = a, scale = s),
                s^k[k < 1] * gamma(a + k[k < 1]) * gamma(1 - k[k < 1])/Ga)
            All.eq(minvparalogis(k, shape = a, scale = s),
                   s^k * gamma(a + k/a) * gamma(1 - k/a)/Ga)
        })
    }
}

## Tests for first three positive limited moments and first two
## negative limited moments.
##
## Limits are taken from quantiles of each distribution.
order <- c(-2, -1, 1, 2, 3)             # orders
q <- c(0.25, 0.50, 0.75, 0.9, 0.95)     # quantiles
for (i in seq_along(shpar))
{
    a <- shpar[[c(i, 1)]]; g <- shpar[[c(i, 2)]]; t <- shpar[[c(i, 3)]]
    Ga <- gamma(a)
    Gt <- gamma(t)
    for (s in scpar)
    {
        limit <- qtrbeta(q, shape1 = a, shape2 = g, shape3 = t, scale = s)
        y <- limit/s
        u <- 1/(1 + y^(-g))
        for (k in order)
            stopifnot(exprs = {
                All.eq(levtrbeta(limit, order = k, shape1 = a, shape2 = g, shape3 = t, scale = s),
                       s^k * betaint(u, t + k/g, a - k/g)/(Ga * Gt) +
                       limit^k * pbeta(u, t, a, lower = FALSE))
            })
        limit <- qburr(q, shape1 = a, shape2 = g, scale = s)
        y <- limit/s
        u <- 1/(1 + y^g)
        for (k in order)
            stopifnot(exprs = {
                All.eq(levburr(limit, order = k, shape1 = a, shape2 = g, scale = s),
                       s^k * betaint(1 - u, 1 + k/g, a - k/g)/Ga +
                       limit^k * u^a)
            })
        limit <- qllogis(q, shape = g, scale = s)
        y <- limit/s
        u <- 1/(1 + y^(-g))
        for (k in order)
            stopifnot(exprs = {
                All.eq(levllogis(limit, order = k, shape = g, scale = s),
                       s^k * betaint(u, 1 + k/g, 1 - k/g) +
                       limit^k * (1 - u))
            })
        limit <- qparalogis(q, shape = a, scale = s)
        y <- limit/s
        u <- 1/(1 + y^a)
        for (k in order)
            stopifnot(exprs = {
                All.eq(levparalogis(limit, order = k, shape = a, scale = s),
                       s^k * betaint(1 - u, 1 + k/a, a - k/a)/Ga +
                       limit^k * u^a)
            })
        limit <- qgenpareto(q, shape1 = a, shape2 = t, scale = s)
        y <- limit/s
        u <- 1/(1 + y^(-1))
        for (k in order)
            stopifnot(exprs = {
                All.eq(levgenpareto(limit, order = k, shape1 = a, shape2 = t, scale = s),
                       s^k * betaint(u, t + k, a - k)/(Ga * Gt) +
                       limit^k * pbeta(u, t, a, lower = FALSE))
            })
        limit <- qpareto(q, shape = a, scale = s)
        y <- limit/s
        u <- 1/(1 + y)
        for (k in order[order > -1])
            stopifnot(exprs = {
                All.eq(levpareto(limit, order = k, shape = a, scale = s),
                       s^k * betaint(1 - u, 1 + k, a - k)/Ga +
                       limit^k * u^a)
        })
        limit <- qpareto1(q, shape = a, min = s)
        for (k in order)
            stopifnot(exprs = {
                All.eq(levpareto1(limit, order = k, shape = a, min = s),
                       s^k * a/(a - k) - k * s^a/((a - k) * limit^(a - k)))
            })
        limit <- qinvburr(q, shape1 = a, shape2 = g, scale = s)
        y <- limit/s
        u <- 1/(1 + y^(-g))
        for (k in order)
            stopifnot(exprs = {
                All.eq(levinvburr(limit, order = k, shape1 = a, shape2 = g, scale = s),
                       s^k * betaint(u, a + k/g, 1 - k/g)/Ga +
                       limit^k * (1 - u^a))
            })
        limit <- qinvpareto(q, shape = a, scale = s)
        y <- limit/s
        u <- 1/(1 + y^(-1))
        for (k in order[order < 1])
            stopifnot(exprs = {
                All.eq(levinvpareto(limit, order = k, shape = a, scale = s),
                       s^k * a *
                       sapply(u,
                              function(upper)
                                  integrate(function(x) x^(a+k-1) * (1-x)^(-k),
                                            lower = 0, upper = upper)$value) +
                       limit^k * (1 - u^a))
            })
        limit <- qinvparalogis(q, shape = a, scale = s)
        y <- limit/s
        u <- 1/(1 + y^(-a))
        for (k in order)
            stopifnot(exprs = {
                All.eq(levinvparalogis(limit, order = k, shape = a, scale = s),
                       s^k * betaint(u, a + k/a, 1 - k/a)/Ga +
                       limit^k * (1 - u^a))
        })
    }
}

##
## TRANSFORMED GAMMA AND INVERSE TRANSFORMED GAMMA FAMILIES
##

## Density: first check that functions return 0 when scale = Inf, and
## when x = scale = Inf (transformed gamma), or when scale = 0 and
## when x = scale = 0 (inverse distributions).
stopifnot(exprs = {
    dtrgamma   (c(42, Inf), shape1 = 2, shape2 = 3, scale = Inf) == c(0, 0)
    dinvtrgamma(c(42, 0), shape1 = 2, shape2 = 3, scale = 0) == c(0, 0)
    dinvgamma  (c(42, 0), shape  = 2, scale = 0) == c(0, 0)
    dinvweibull(c(42, 0), shape = 3, scale = 0) == c(0, 0)
    dinvexp    (c(42, 0), scale = 0) == c(0, 0)
})

## Tests on the density
set.seed(123)                   # reset the seed
nshpar <- 2                     # (maximum) number of shape parameters
shpar <- replicate(30, rgamma(nshpar, 5), simplify = FALSE)
scpar <- rlnorm(30, 2)          # scale parameters
for (i in seq_along(shpar))
{
    a <- shpar[[c(i, 1)]]; t <- shpar[[c(i, 2)]]
    Ga <- gamma(a)
    for (s in scpar)
    {
        x <- rtrgamma(100, shape1 = a, shape2 = t, scale = s)
        y <- x/s
        u <- y^t
        stopifnot(exprs = {
            all.equal(d1 <- dtrgamma(x, shape1 = a, shape2 = t,
                                    scale = s),
                      d2 <- dtrgamma(y, shape1 = a, shape2 = t,
                                     scale = 1)/s,
                      tolerance = 1e-10)
            all.equal(d2,
                      t/(Ga * s^(a * t)) * x^(a * t - 1) * exp(-u),
                      tolerance = 1e-10)
            all.equal(d1,
                      t/(Ga * x) * u^a * exp(-u),
                      tolerance = 1e-10)
        })
        x <- rinvtrgamma(100, shape1 = a, shape2 = t, scale = s)
        y <- x/s
        u <- y^(-t)
        stopifnot(exprs = {
            all.equal(d1 <- dinvtrgamma(x, shape1 = a, shape2 = t,
                                        scale = s),
                      d2 <- dinvtrgamma(y, shape1 = a, shape2 = t,
                                        scale = 1)/s,
                      tolerance = 1e-10)
            all.equal(d2,
                      t * s^(a * t)/(Ga * x^(a * t + 1)) * exp(-u),
                      tolerance = 1e-10)
            all.equal(d1,
                      t/(Ga * x) * u^a * exp(-u),
                      tolerance = 1e-10)
        })
        x <- rinvgamma(100, shape = a, scale = s)
        y <- x/s
        u <- y^(-1)
        stopifnot(exprs = {
            all.equal(d1 <- dinvgamma(x, shape = a, scale = s),
                      d2 <- dinvgamma(y, shape = a, scale = 1)/s,
                      tolerance = 1e-10)
            all.equal(d2,
                      s^a/(Ga * x^(a + 1)) * exp(-u),
                      tolerance = 1e-10)
            all.equal(d1,
                      1/(Ga * x) * u^a * exp(-u),
                      tolerance = 1e-10)
        })
        x <- rinvweibull(100, shape = t, scale = s)
        y <- x/s
        u <- y^(-t)
        stopifnot(exprs = {
            all.equal(d1 <- dinvweibull(x, shape = t, scale = s),
                      d2 <- dinvweibull(y, shape = t, scale = 1)/s,
                      tolerance = 1e-10)
            all.equal(d2,
                      t * s^t/x^(t + 1) * exp(-u),
                      tolerance = 1e-10)
            all.equal(d1,
                      t/x * u * exp(-u),
                      tolerance = 1e-10)
        })
        x <- rinvexp(100, scale = s)
        y <- x/s
        u <- y^(-1)
        stopifnot(exprs = {
            all.equal(d1 <- dinvexp(x, scale = s),
                      d2 <- dinvexp(y, scale = 1)/s,
                      tolerance = 1e-10)
            all.equal(d2,
                      s/x^2 * exp(-u),
                      tolerance = 1e-10)
            all.equal(d1,
                      1/x * u * exp(-u),
                      tolerance = 1e-10)
        })
    }
}

## Tests on the cumulative distribution function.
scLrg <- c(2, 100, 1e300 * c(0.1, 1, 10, 100), 1e307, xMax, Inf)
stopifnot(exprs = {
    ptrgamma(Inf,  2, 3, scale = xMax) == 1
    ptrgamma(xMax, 2, 3, scale = xMax) == pgamma(1, 2, 1)
    ptrgamma(xMin, 2, 1, scale = 1)    == pgamma(xMin, 2, 1)
    all.equal(ptrgamma(1e300, shape1 = 2, shape2 = 1, scale = scLrg, log = TRUE),
              pgamma(c(5e299, 1e+298, 10, 1, 0.1, 0.01, 1e-7, 1e+300/xMax, 0),
                     2, 1, log = TRUE))
})
scLrg <- c(2, 100, 1e300 * c(0.1, 1, 10, 100), 1e307, xMax, 0)
stopifnot(exprs = {
    pinvtrgamma(Inf,  2, 3, scale = xMax) == 1
    pinvtrgamma(xMax, 2, 3, scale = xMax) == pgamma(1, 2, 1, lower = FALSE)
    pinvtrgamma(xMin, 2, 1, scale = 1)    == pgamma(1/xMin, 2, 1, lower = FALSE)
    all.equal(pinvtrgamma(1e300, shape1 = 2, shape2 = 1, scale = scLrg, log = TRUE),
              pgamma(c(2e-300, 1e-298, 0.1, 1, 10, 100, 1e+7, xMax/1e+300, 0),
                     2, 1, lower = FALSE, log = TRUE))
})
stopifnot(exprs = {
    pinvgamma(Inf,  2, scale = xMax) == 1
    pinvgamma(xMax, 2, scale = xMax) == pgamma(1, 2, 1, lower = FALSE)
    pinvgamma(xMin, 2, scale = 1)    == pgamma(1/xMin, 2, 1, lower = FALSE)
    all.equal(pinvgamma(1e300, shape = 2, scale = scLrg, log = TRUE),
              pgamma(c(2e-300, 1e-298, 0.1, 1, 10, 100, 1e+7, xMax/1e+300, 0),
                     2, 1, lower = FALSE, log = TRUE))
})
stopifnot(exprs = {
    pinvweibull(Inf,  3, scale = xMax) == 1
    pinvweibull(xMax, 3, scale = xMax) == exp(-1)
    pinvweibull(xMin, 1, scale = 1)    == exp(-1/xMin)
    all.equal(pinvweibull(1e300, shape = 1, scale = scLrg, log = TRUE),
              -c(2e-300, 1e-298, 0.1, 1, 10, 100, 1e+7, xMax/1e+300, 0))
})
stopifnot(exprs = {
    pinvexp(Inf,  3, scale = xMax) == 1
    pinvexp(xMax, 3, scale = xMax) == exp(-1)
    pinvexp(xMin, 1, scale = 1)    == exp(-1/xMin)
    all.equal(pinvexp(1e300, scale = scLrg, log = TRUE),
              -c(2e-300, 1e-298, 0.1, 1, 10, 100, 1e+7, xMax/1e+300, 0))
})

## Tests for first three positive moments and first two negative
## moments. (Including for the Gamma, Weibull and Exponential
## distributions of base R.)
##
## Simulation of new parameters ensuring that said moments exist.
set.seed(123)                   # reset the seed
nshpar <- 2                     # (maximum) number of shape parameters
shpar <- replicate(30, c(3, 3) + rlnorm(nshpar, 2), simplify = FALSE)
scpar <- rlnorm(30, 2)          # scale parameters
k <- c(-2, -1, 1, 2, 3)         # orders
for (i in seq_along(shpar))
{
    a <- shpar[[c(i, 1)]]; t <- shpar[[c(i, 2)]]
    Ga <- gamma(a)
    for (s in scpar)
    {
        stopifnot(exprs = {
            All.eq(mtrgamma(k, shape1 = a, shape2 = t, scale = s),
                   s^k * gamma(a + k/t)/Ga)
            All.eq(mgamma(k, shape = a, scale = s),
                   s^k * gamma(a + k)/Ga)
            All.eq(mweibull(k, shape = t, scale = s),
                   s^k * gamma(1 + k/t))
            All.eq(mexp(k[k > -1], rate = 1/s),
                   s^k[k > -1] * gamma(1 + k[k > -1]))
            All.eq(minvtrgamma(k, shape1 = a, shape2 = t, scale = s),
                   s^k * gamma(a - k/t)/Ga)
            All.eq(minvgamma(k, shape = a, scale = s),
                   s^k * gamma(a - k)/Ga)
            All.eq(minvweibull(k, shape = t, scale = s),
                   s^k * gamma(1 - k/t))
            All.eq(minvexp(k[k < 1], scale = s),
                   s^k[k < 1] * gamma(1 - k[k < 1]))
        })
    }
}

## Tests for first three positive limited moments and first two
## negative limited moments. (Including for the Gamma, Weibull and
## Exponential distributions of base R.)
##
## Limits are taken from quantiles of each distribution.
order <- c(-2, -1, 1, 2, 3)             # orders
q <- c(0.25, 0.50, 0.75, 0.9, 0.95)     # quantiles
for (i in seq_along(shpar))
{
    a <- shpar[[c(i, 1)]]; t <- shpar[[c(i, 2)]]
    Ga <- gamma(a)
    for (s in scpar)
    {
        limit <- qtrgamma(q, shape1 = a, shape2 = t, scale = s)
        y <- limit/s
        u <- y^t
        for (k in order)
            stopifnot(exprs = {
                All.eq(levtrgamma(limit, order = k, shape1 = a, shape2 = t, scale = s),
                       s^k * gamma(a + k/t)/Ga * pgamma(u, a + k/t, scale = 1) +
                       limit^k * pgamma(u, a, scale = 1, lower = FALSE))
            })
        limit <- qgamma(q, shape = a, scale = s)
        y <- limit/s
        for (k in order)
            stopifnot(exprs = {
                All.eq(levgamma(limit, order = k, shape = a, scale = s),
                       s^k * gamma(a + k)/Ga * pgamma(y, a + k, scale = 1) +
                       limit^k * pgamma(y, a, scale = 1, lower = FALSE))
            })
        limit <- qweibull(q, shape = t, scale = s)
        y <- limit/s
        u <- y^t
        for (k in order)
            stopifnot(exprs = {
                All.eq(levweibull(limit, order = k, shape = t, scale = s),
                       s^k * gamma(1 + k/t) * pgamma(u, 1 + k/t, scale = 1) +
                       limit^k * pgamma(u, 1, scale = 1, lower = FALSE))
            })
        limit <- qexp(q, rate = 1/s)
        y <- limit/s
        for (k in order[order > -1])
            stopifnot(exprs = {
                All.eq(levexp(limit, order = k, rate = 1/s),
                       s^k * gamma(1 + k) * pgamma(y, 1 + k, scale = 1) +
                       limit^k * pgamma(y, 1, scale = 1, lower = FALSE))
            })
        limit <- qinvtrgamma(q, shape1 = a, shape2 = t, scale = s)
        y <- limit/s
        u <- y^(-t)
        for (k in order)
            stopifnot(exprs = {
                All.eq(levinvtrgamma(limit, order = k, shape1 = a, shape2 = t, scale = s),
                       s^k * (gammainc(a - k/t, u)/Ga) +
                       limit^k * pgamma(u, a, scale = 1))
            })
        limit <- qinvgamma(q, shape = a, scale = s)
        y <- limit/s
        u <- y^(-1)
        for (k in order)
            stopifnot(exprs = {
                All.eq(levinvgamma(limit, order = k, shape = a, scale = s),
                       s^k * (gammainc(a - k, u)/Ga) +
                       limit^k * pgamma(u, a, scale = 1))
            })
        limit <- qinvweibull(q, shape = t, scale = s)
        y <- limit/s
        u <- y^(-t)
        for (k in order)
            stopifnot(exprs = {
                All.eq(levinvweibull(limit, order = k, shape = t, scale = s),
                       s^k * gammainc(1 - k/t, u) +
                       limit^k * (-expm1(-u)))
            })
        limit <- qinvexp(q, scale = s)
        y <- limit/s
        u <- y^(-1)
        for (k in order)
            stopifnot(exprs = {
                All.eq(levinvexp(limit, order = k, scale = s),
                       s^k * gammainc(1 - k, u) +
                       limit^k * (-expm1(-u)))
            })
    }
}

##
## OTHER DISTRIBUTIONS
##

## Distributions in this category are quite different, so let's treat
## them separately.

## LOGGAMMA

## Tests on the density.
stopifnot(exprs = {
    dlgamma(c(42, Inf), shapelog = 2, ratelog = 0) == c(0, 0)
})
assertWarning(stopifnot(exprs = {
    is.nan(dlgamma(c(0, 42, Inf), shapelog = 2, ratelog = Inf))
}))
x <- rlgamma(100, shapelog = 2, ratelog = 1)
for(a in round(rlnorm(30), 2))
{
    Ga <- gamma(a)
    for(r in round(rlnorm(30), 2))
	stopifnot(exprs = {
            All.eq(dlgamma(x, shapelog = a, ratelog = r),
                   r^a * (log(x))^(a - 1)/(Ga * x^(r + 1)))
        })
}

## Tests on the cumulative distribution function.
assertWarning(stopifnot(exprs = {
    is.nan(plgamma(Inf,   1, ratelog = Inf))
    is.nan(plgamma(Inf, Inf, ratelog = Inf))
}))
scLrg <- log(c(2, 100, 1e300 * c(0.1, 1, 10, 100), 1e307, xMax, Inf))
stopifnot(exprs = {
    plgamma(Inf,  2, ratelog = xMax) == 1
    plgamma(xMax, 2, ratelog = 0) == 0
    all.equal(plgamma(1e300, 2, ratelog = 1/scLrg, log = TRUE),
              pgamma(log(1e300), 2, scale = scLrg, log = TRUE))
})

## Tests for first three positive moments and first two negative
## moments.
k <- c(-2, -1, 1, 2, 3)         # orders
for(a in round(rlnorm(30), 2))
{
    Ga <- gamma(a)
    for(r in 3 + round(rlnorm(30), 2))
	stopifnot(exprs = {
            All.eq(mlgamma(k, shapelog = a, ratelog = r),
                   (1 - k/r)^(-a))
        })
}

## Tests for first three positive limited moments and first two
## negative limited moments.
order <- c(-2, -1, 1, 2, 3)             # orders
q <- c(0.25, 0.50, 0.75, 0.9, 0.95)     # quantiles
for(a in round(rlnorm(30), 2))
{
    Ga <- gamma(a)
    for(r in 3 + round(rlnorm(30), 2))
    {
        limit <- qlgamma(q, shapelog = a, ratelog = r)
        for (k in order)
        {
            u <- log(limit)
        stopifnot(exprs = {
            All.eq(levlgamma(limit, order = k, shapelog = a, ratelog = r),
                   (1 - k/r)^(-a) * pgamma((r - k) * u, a, scale = 1) +
                   limit^k * pgamma(r * u, a, scale = 1,lower = FALSE))
            })
        }
    }
}

## GUMBEL

## Tests on the density.
stopifnot(exprs = {
    dgumbel(c(1, 3, Inf), alpha = 2,   scale = Inf) == c(0, 0, 0)
    dgumbel(c(1, 2, 3),   alpha = 2,   scale = 0) == c(0, Inf, 0)
    dgumbel(c(-Inf, Inf), alpha = 1,   scale = 1) == c(0, 0)
    dgumbel(1,            alpha = Inf, scale = 1) == 0
})
assertWarning(stopifnot(exprs = {
    is.nan(dgumbel(Inf,  alpha = Inf,  scale =  1))
    is.nan(dgumbel(-Inf, alpha = -Inf, scale =  1))
    is.nan(dgumbel(Inf,  alpha = 1,    scale = -1))
    is.nan(dgumbel(1,    alpha = 1,    scale = -1))
    is.nan(dgumbel(1,    alpha = Inf,  scale = -1))
}))
x <- rgumbel(100, alpha = 2, scale = 5)
for(a in round(rlnorm(30), 2))
{
    Ga <- gamma(a)
    for(s in round(rlnorm(30), 2))
    {
        u <- (x - a)/s
	stopifnot(exprs = {
            All.eq(dgumbel(x, alpha = a, scale = s),
                   exp(-(u + exp(-u)))/s)
        })
    }
}

## Tests on the cumulative distribution function.
assertWarning(stopifnot(exprs = {
    is.nan(pgumbel(Inf,  alpha = Inf,  scale =  1))
    is.nan(pgumbel(-Inf, alpha = -Inf, scale =  1))
    is.nan(pgumbel(Inf,  alpha = 1,    scale = -1))
    is.nan(pgumbel(1,    alpha = 1,    scale = -1))
    is.nan(pgumbel(1,    alpha = Inf,  scale = -1))
}))
scLrg <- c(2, 100, 1e300 * c(0.1, 1, 10, 100), 1e307, xMax, Inf)
stopifnot(exprs = {
    pgumbel(c(-Inf, Inf),  2, scale = xMax) == c(0, 1)
    pgumbel(c(xMin, xMax), 2, scale = 0) == c(0, 1)
    all.equal(pgumbel(1e300, 0, scale = scLrg, log = TRUE),
              -exp(-c(5e299, 1e+298, 10, 1, 0.1, 0.01, 1e-7, 1e+300/xMax, 0)))
})

## Test the first two moments, the only ones implemented.
assertWarning(stopifnot(exprs = {
    is.nan(mgumbel(c(-2, -1, 3, 4), alpha = 2, scale = 5))
}))
stopifnot(exprs = {
    All.eq(mgumbel(1, alpha = 2, scale = 5),
           2 + 5 * 0.577215664901532860606512090082)
    All.eq(mgumbel(2, alpha = 2, scale = 5),
           pi^2 * 25/6 + (2 + 5 * 0.577215664901532860606512090082)^2)
})

## INVERSE GAUSSIAN

## Tests on the density.
stopifnot(exprs = {
    dinvgauss(c(1, 3, Inf),  mean = 2,   dispersion = Inf) == c(0, 0, 0)
    dinvgauss(c(0, 42, Inf), mean = 2,   dispersion = 0) == c(Inf, 0, 0)
    dinvgauss(c(0, Inf),     mean = 1,   dispersion = 1) == c(0, 0)
    dinvgauss(1,             mean = Inf, dispersion = 2) == dinvgamma(1, 0.5, scale = 0.25)
})
assertWarning(stopifnot(exprs = {
    is.nan(dinvgauss(-Inf, mean = -1,  dispersion =  1))
    is.nan(dinvgauss(Inf,  mean = 1,   dispersion = -1))
    is.nan(dinvgauss(1,    mean = 1,   dispersion = -1))
    is.nan(dinvgauss(1,    mean = Inf, dispersion = -1))
}))
x <- rinvgauss(100, mean = 2, dispersion = 5)
for(mu in round(rlnorm(30), 2))
{
    for(phi in round(rlnorm(30), 2))
	stopifnot(exprs = {
            All.eq(dinvgauss(x, mean = mu, dispersion = phi),
                   1/sqrt(2*pi*phi*x^3) * exp(-((x/mu - 1)^2)/(2*phi*x)))
        })
}

## Tests on the cumulative distribution function.
assertWarning(stopifnot(exprs = {
    is.nan(pinvgauss(-Inf, mean = -Inf, dispersion =  1))
    is.nan(pinvgauss(Inf,  mean = 1,    dispersion = -1))
    is.nan(pinvgauss(1,    mean = Inf,  dispersion = -1))
}))
x <- c(1:50, 10^c(3:10, 20, 50, 150, 250))
sqx <- sqrt(x)
stopifnot(exprs = {
    pinvgauss(c(0, Inf),  mean = 2, dispersion = xMax) == c(0, 1)
    pinvgauss(c(0, xMax), mean = xMax, dispersion = 0) == c(0, 1)
    all.equal(pinvgauss(x, 1, dispersion = 1, log = TRUE),
              log(pnorm(sqx - 1/sqx) + exp(2) * pnorm(-sqx - 1/sqx)))
})

## Tests for small value of 'shape'. Added for the patch in 4294e9c.
q <- runif(100)
stopifnot(exprs = {
    all.equal(q,
              pinvgauss(qinvgauss(q, 0.1, 1e-2), 0.1, 1e-2))
    all.equal(q,
              pinvgauss(qinvgauss(q, 0.1, 1e-6), 0.1, 1e-6))
})

## Tests for first three positive, integer moments.
k <- 1:3
for(mu in round(rlnorm(30), 2))
{
    for(phi in round(rlnorm(30), 2))
	stopifnot(exprs = {
            All.eq(minvgauss(k, mean = mu, dispersion = phi),
                   c(mu,
                     mu^2 * (1 + phi * mu),
                     mu^3 * (1 + 3 * phi * mu + 3 * (phi * mu)^2)))
        })
}

## Tests for limited expected value.
q <- c(0.25, 0.50, 0.75, 0.9, 0.95)     # quantiles
for(mu in round(rlnorm(30), 2))
{
    for(phi in round(rlnorm(30), 2))
    {
        limit <- qinvgauss(q, mean = mu, dispersion = phi)
        stopifnot(exprs = {
            All.eq(levinvgauss(limit, mean = mu, dispersion = phi),
                   mu * (pnorm((limit/mu - 1)/sqrt(phi * limit)) -
                         exp(2/phi/mu) * pnorm(-(limit/mu + 1)/sqrt(phi * limit))) +
                   limit * pinvgauss(limit, mean = mu, dispersion = phi, lower = FALSE))
        })
    }
}

## GENERALIZED BETA
stopifnot(exprs = {
    dgenbeta(c(0, 2.5, 5), shape1 = 0,   shape2 = 0,   shape3 = 3,   scale = 5) == c(Inf, 0, Inf)
    dgenbeta(c(0, 2.5, 5), shape1 = 0,   shape2 = 0,   shape3 = 0,   scale = 5) == c(Inf, 0, Inf)
    dgenbeta(c(0, 2.5, 5), shape1 = 0,   shape2 = 2,   shape3 = 0,   scale = 5) == c(Inf, 0, 0)
    dgenbeta(c(0, 2.5, 5), shape1 = 0,   shape2 = Inf, shape3 = 3,   scale = 5) == c(Inf, 0, 0)
    dgenbeta(c(0, 2.5, 5), shape1 = 1,   shape2 = Inf, shape3 = 3,   scale = 5) == c(Inf, 0, 0)
    dgenbeta(c(0, 2.5, 5), shape1 = Inf, shape2 = Inf, shape3 = 3,   scale = 5) == c(0, Inf, 0)
    dgenbeta(c(0, 2.5, 5), shape1 = Inf, shape2 = Inf, shape3 = Inf, scale = 5) == c(0, 0, Inf)
})
nshpar <- 3                     # number of shape parameters
shpar <- replicate(30, rlnorm(nshpar, 2), simplify = FALSE)
scpar <- rlnorm(30, 2)          # scale parameters
for (i in seq_along(shpar))
{
    a <- shpar[[c(i, 1)]]; b <- shpar[[c(i, 2)]]; t <- shpar[[c(i, 3)]]
    Be <- beta(a, b)
    for (s in scpar)
    {
        u <- rbeta(100, a, b)
        y <- u^(1/t)
        x <- s * y
        stopifnot(exprs = {
            all.equal(d1 <- dgenbeta(x, shape1 = a, shape2 = b, shape3 = t,
                                     scale = s),
                      d2 <- dgenbeta(y, shape1 = a, shape2 = b, shape3 = t,
                                     scale = 1)/s,
                      tolerance = 1e-10)
            all.equal(d1,
                      t * y^(a*t - 1) * (1 - y^t)^(b - 1)/(s * Be),
                      tolerance = 1e-10)
            all.equal(d1,
                      t * u^a * (1 - u)^(b - 1)/(x * Be),
                      tolerance = 1e-10)
        })
    }
}

## Tests on the cumulative distribution function.
scLrg <- 1e300 * c(0.5, 1, 2, 4)
stopifnot(exprs = {
    all.equal(pgenbeta(1e300,
                       shape1 = 3, shape2 = 1,
                       shape3 = rep(c(1, 2), each = length(scLrg)),
                       scale = scLrg, log = TRUE),
              c(0, pbeta(c(1, 1/2, 1/4), 3, 1, log = TRUE),
                0, pbeta(c(1, 1/4, 1/16), 3, 1, log = TRUE)))
})

## Tests for first three positive moments and first two negative
## moments.
##
## Simulation of new parameters ensuring that said moments exist.
set.seed(123)                   # reset the seed
nshpar <- 3                     # number of shape parameters
shpar <- replicate(30, sqrt(c(3, 0, 3)) + rlnorm(nshpar, 2), simplify = FALSE)
scpar <- rlnorm(30, 2)          # scale parameters
k <- c(-2, -1, 1, 2, 3)         # orders
for (i in seq_along(shpar))
{
    a <- shpar[[c(i, 1)]]; b <- shpar[[c(i, 2)]]; t <- shpar[[c(i, 3)]]
    Be <- beta(a, b)
    for (s in scpar)
        stopifnot(exprs = {
            All.eq(mgenbeta(k, shape1 = a, shape2 = b, shape3 = t, scale = s),
                   s^k * beta(a + k/t, b)/Be)
        })
}

## Tests for first three positive limited moments and first two
## negative limited moments.
##
## Simulation of new parameters ensuring that said moments exist.
order <- c(-2, -1, 1, 2, 3)             # orders
q <- c(0.25, 0.50, 0.75, 0.9, 0.95)     # quantiles
for (i in seq_along(shpar))
{
    a <- shpar[[c(i, 1)]]; g <- shpar[[c(i, 2)]]; t <- shpar[[c(i, 3)]]
    Be <- beta(a, b)
    for (s in scpar)
    {
        limit <- qgenbeta(q, shape1 = a, shape2 = b, shape3 = t, scale = s)
        u <- (limit/s)^t
        for (k in order)
            stopifnot(exprs = {
                All.eq(levgenbeta(limit, order = k, shape1 = a, shape2 = b, shape3 = t, scale = s),
                       s^k * beta(a + k/t, b)/Be * pbeta(u, a + k/t, b) +
                       limit^k * pbeta(u, a, b, lower = FALSE))
            })
    }
}

##
## RANDOM NUMBERS (all continuous distributions)
##
set.seed(123)
n <- 20
m <- rnorm(1)

## Generate variates
Rfpareto <- rfpareto(n, min = m, shape1 = 0.8, shape2 = 1.5, shape3 = 2, scale = 2)
Rpareto4 <- rpareto4(n, min = m, shape1 = 0.8, shape2 = 1.5,             scale = 2)
Rpareto3 <- rpareto3(n, min = m,               shape  = 1.5,             scale = 2)
Rpareto2 <- rpareto2(n, min = m, shape  = 0.8,                           scale = 2)
Rtrbeta        <- rtrbeta      (n, shape1 = 0.8, shape2 = 1.5, shape3 = 2, scale = 2)
Rburr          <- rburr        (n, shape1 = 0.8, shape2 = 1.5,             scale = 2)
Rllogis        <- rllogis      (n, shape  = 1.5, scale = 2)
Rparalogis     <- rparalogis   (n, shape  = 0.8, scale = 2)
Rgenpareto     <- rgenpareto   (n, shape1 = 0.8, shape2 = 2, scale = 2)
Rpareto        <- rpareto      (n, shape  = 0.8, scale = 2)
Rpareto1       <- rpareto1     (n, shape  = 0.8, min = 2)
Rinvburr       <- rinvburr     (n, shape1 = 1.5, shape2 = 2, scale = 2)
Rinvpareto     <- rinvpareto   (n, shape  = 2,   scale = 2)
Rinvparalogis  <- rinvparalogis(n, shape  = 2,   scale = 2)
Rtrgamma       <- rtrgamma     (n, shape1 = 2, shape2 = 3, scale = 5)
Rinvtrgamma    <- rinvtrgamma  (n, shape1 = 2, shape2 = 3, scale = 5)
Rinvgamma      <- rinvgamma    (n, shape  = 2,             scale = 5)
Rinvweibull    <- rinvweibull  (n,             shape  = 3, scale = 5)
Rinvexp        <- rinvexp      (n,                         scale = 5)
Rlgamma <- rlgamma(n, shapelog = 1.5, ratelog = 5)
Rgumbel <- rgumbel(n, alpha = 2, scale = 5)
Rinvgauss <- rinvgauss(n, mean = 2, dispersion = 5)
Rgenbeta <- rgenbeta(n, shape1 = 0.8, shape2 = 1.5, shape3 = 2, scale = 2)

## Compute quantiles
Pfpareto <- pfpareto(Rfpareto, min = m, shape1 = 0.8, shape2 = 1.5, shape3 = 2, scale = 2)
Ppareto4 <- ppareto4(Rpareto4, min = m, shape1 = 0.8, shape2 = 1.5,             scale = 2)
Ppareto3 <- ppareto3(Rpareto3, min = m,               shape  = 1.5,             scale = 2)
Ppareto2 <- ppareto2(Rpareto2, min = m, shape  = 0.8,                           scale = 2)
Ptrbeta        <- ptrbeta      (Rtrbeta,       shape1 = 0.8, shape2 = 1.5, shape3 = 2, scale = 2)
Pburr          <- pburr        (Rburr,         shape1 = 0.8, shape2 = 1.5,             scale = 2)
Pllogis        <- pllogis      (Rllogis,       shape  = 1.5, scale = 2)
Pparalogis     <- pparalogis   (Rparalogis,    shape  = 0.8, scale = 2)
Pgenpareto     <- pgenpareto   (Rgenpareto,    shape1 = 0.8, shape2 = 2, scale = 2)
Ppareto        <- ppareto      (Rpareto,       shape  = 0.8, scale = 2)
Ppareto1       <- ppareto1     (Rpareto1,      shape  = 0.8, min = 2)
Pinvburr       <- pinvburr     (Rinvburr,      shape1 = 1.5, shape2 = 2, scale = 2)
Pinvpareto     <- pinvpareto   (Rinvpareto,    shape  = 2,   scale = 2)
Pinvparalogis  <- pinvparalogis(Rinvparalogis, shape  = 2,   scale = 2)
Ptrgamma       <- ptrgamma     (Rtrgamma,      shape1 = 2, shape2 = 3, scale = 5)
Pinvtrgamma    <- pinvtrgamma  (Rinvtrgamma,   shape1 = 2, shape2 = 3, scale = 5)
Pinvgamma      <- pinvgamma    (Rinvgamma,     shape  = 2,             scale = 5)
Pinvweibull    <- pinvweibull  (Rinvweibull,               shape  = 3, scale = 5)
Pinvexp        <- pinvexp      (Rinvexp,                               scale = 5)
Plgamma <- plgamma(Rlgamma, shapelog = 1.5, ratelog = 5)
Pgumbel <- pgumbel(Rgumbel, alpha = 2, scale = 5)
Pinvgauss <- pinvgauss(Rinvgauss, mean = 2, dispersion = 5)
Pgenbeta <- pgenbeta(Rgenbeta, shape1 = 0.8, shape2 = 1.5, shape3 = 2, scale = 2)

## Just compute pdf
Dfpareto <- dfpareto(Rfpareto, min = m, shape1 = 0.8, shape2 = 1.5, shape3 = 2, scale = 2)
Dpareto4 <- dpareto4(Rpareto4, min = m, shape1 = 0.8, shape2 = 1.5,             scale = 2)
Dpareto3 <- dpareto3(Rpareto3, min = m,               shape  = 1.5,             scale = 2)
Dpareto2 <- dpareto2(Rpareto2, min = m, shape  = 0.8,                           scale = 2)
Dtrbeta        <- dtrbeta      (Rtrbeta,       shape1 = 0.8, shape2 = 1.5, shape3 = 2, scale = 2)
Dburr          <- dburr        (Rburr,         shape1 = 0.8, shape2 = 1.5,             scale = 2)
Dllogis        <- dllogis      (Rllogis,       shape  = 1.5, scale = 2)
Dparalogis     <- dparalogis   (Rparalogis,    shape  = 0.8, scale = 2)
Dgenpareto     <- dgenpareto   (Rgenpareto,    shape1 = 0.8, shape2 = 2, scale = 2)
Dpareto        <- dpareto      (Rpareto,       shape  = 0.8, scale = 2)
Dpareto1       <- dpareto1     (Rpareto1,      shape  = 0.8, min = 2)
Dinvburr       <- dinvburr     (Rinvburr,      shape1 = 1.5, shape2 = 2, scale = 2)
Dinvpareto     <- dinvpareto   (Rinvpareto,    shape  = 2,   scale = 2)
Dinvparalogis  <- dinvparalogis(Rinvparalogis, shape  = 2,   scale = 2)
Dtrgamma       <- dtrgamma     (Rtrgamma,      shape1 = 2, shape2 = 3, scale = 5)
Dinvtrgamma    <- dinvtrgamma  (Rinvtrgamma,   shape1 = 2, shape2 = 3, scale = 5)
Dinvgamma      <- dinvgamma    (Rinvtrgamma,   shape  = 2,             scale = 5)
Dinvweibull    <- dinvweibull  (Rinvweibull,               shape  = 3, scale = 5)
Dinvexp        <- dinvexp      (Rinvexp,                               scale = 5)
Dlgamma <- dlgamma(Rlgamma, shapelog = 1.5, ratelog = 5)
Dgumbel <- dgumbel(Rgumbel, alpha = 2, scale = 5)
Dinvgauss <- dinvgauss(Rinvgauss, mean = 2, dispersion = 5)
Dgenbeta <- dgenbeta(Rgenbeta, shape1 = 0.8, shape2 = 1.5, shape3 = 2, scale = 2)

## Check q<dist>(p<dist>(.)) identity
stopifnot(exprs = {
    All.eq(Rfpareto, qfpareto(Pfpareto, min = m, shape1 = 0.8, shape2 = 1.5, shape3 = 2, scale = 2))
    All.eq(Rpareto4, qpareto4(Ppareto4, min = m, shape1 = 0.8, shape2 = 1.5,             scale = 2))
    All.eq(Rpareto3, qpareto3(Ppareto3, min = m,               shape  = 1.5,             scale = 2))
    All.eq(Rpareto2, qpareto2(Ppareto2, min = m, shape  = 0.8,                           scale = 2))
    All.eq(Rtrbeta,       qtrbeta      (Ptrbeta,       shape1 = 0.8, shape2 = 1.5, shape3 = 2, scale = 2))
    All.eq(Rburr,         qburr        (Pburr,         shape1 = 0.8, shape2 = 1.5,             scale = 2))
    All.eq(Rllogis,       qllogis      (Pllogis,       shape  = 1.5, scale = 2))
    All.eq(Rparalogis,    qparalogis   (Pparalogis,    shape  = 0.8, scale = 2))
    All.eq(Rgenpareto,    qgenpareto   (Pgenpareto,    shape1 = 0.8, shape2 = 2, scale = 2))
    All.eq(Rpareto,       qpareto      (Ppareto,       shape  = 0.8, scale = 2))
    All.eq(Rpareto1,      qpareto1     (Ppareto1,      shape  = 0.8, min = 2))
    All.eq(Rinvburr,      qinvburr     (Pinvburr,      shape1 = 1.5, shape2 = 2, scale = 2))
    All.eq(Rinvpareto,    qinvpareto   (Pinvpareto,    shape  = 2,   scale = 2))
    All.eq(Rinvparalogis, qinvparalogis(Pinvparalogis, shape  = 2,   scale = 2))
    All.eq(Rtrgamma,      qtrgamma     (Ptrgamma,      shape1 = 2, shape2 = 3, scale = 5))
    All.eq(Rinvtrgamma,   qinvtrgamma  (Pinvtrgamma,   shape1 = 2, shape2 = 3, scale = 5))
    All.eq(Rinvgamma,     qinvgamma    (Pinvgamma,     shape  = 2,             scale = 5))
    All.eq(Rinvweibull,   qinvweibull  (Pinvweibull,               shape  = 3, scale = 5))
    All.eq(Rinvexp,       qinvexp      (Pinvexp,                               scale = 5))
    All.eq(Rlgamma, qlgamma(Plgamma, shapelog = 1.5, ratelog = 5))
    All.eq(Rgumbel, qgumbel(Pgumbel, alpha = 2, scale = 5))
    All.eq(Rinvgauss, qinvgauss(Pinvgauss, mean = 2, dispersion = 5))
    All.eq(Rgenbeta, qgenbeta(Pgenbeta, shape1 = 0.8, shape2 = 1.5, shape3 = 2, scale = 2))
})

## Check q<dist>(p<dist>(.)) identity for special cases
stopifnot(exprs = {
    All.eq(Rfpareto - m, qtrbeta(Pfpareto, shape1 = 0.8, shape2 = 1.5, shape3 = 2, scale = 2))
    All.eq(Rpareto4 - m, qburr  (Ppareto4, shape1 = 0.8, shape2 = 1.5,             scale = 2))
    All.eq(Rpareto3 - m, qllogis(Ppareto3,               shape  = 1.5,             scale = 2))
    All.eq(Rpareto2 - m, qpareto(Ppareto2, shape  = 0.8,                           scale = 2))
})

## Check q<dist>(p<dist>(.)) identity with upper tail
stopifnot(exprs = {
    All.eq(Rfpareto, qfpareto(1 - Pfpareto, min = m, shape1 = 0.8, shape2 = 1.5, shape3 = 2, scale = 2, lower = FALSE))
    All.eq(Rpareto4, qpareto4(1 - Ppareto4, min = m, shape1 = 0.8, shape2 = 1.5,             scale = 2, lower = FALSE))
    All.eq(Rpareto3, qpareto3(1 - Ppareto3, min = m,               shape  = 1.5,             scale = 2, lower = FALSE))
    All.eq(Rpareto2, qpareto2(1 - Ppareto2, min = m, shape  = 0.8,                           scale = 2, lower = FALSE))
    All.eq(Rtrbeta,       qtrbeta      (1 - Ptrbeta,       shape1 = 0.8, shape2 = 1.5, shape3 = 2, scale = 2, lower = FALSE))
    All.eq(Rburr,         qburr        (1 - Pburr,         shape1 = 0.8, shape2 = 1.5,             scale = 2, lower = FALSE))
    All.eq(Rllogis,       qllogis      (1 - Pllogis,       shape  = 1.5, scale = 2, lower = FALSE))
    All.eq(Rparalogis,    qparalogis   (1 - Pparalogis,    shape  = 0.8, scale = 2, lower = FALSE))
    All.eq(Rgenpareto,    qgenpareto   (1 - Pgenpareto,    shape1 = 0.8, shape2 = 2, scale = 2, lower = FALSE))
    All.eq(Rpareto,       qpareto      (1 - Ppareto,       shape  = 0.8, scale = 2, lower = FALSE))
    All.eq(Rpareto1,      qpareto1     (1 - Ppareto1,      shape  = 0.8, min = 2, lower = FALSE))
    All.eq(Rinvburr,      qinvburr     (1 - Pinvburr,      shape1 = 1.5, shape2 = 2, scale = 2, lower = FALSE))
    All.eq(Rinvpareto,    qinvpareto   (1 - Pinvpareto,    shape  = 2,   scale = 2, lower = FALSE))
    All.eq(Rinvparalogis, qinvparalogis(1 - Pinvparalogis, shape  = 2,   scale = 2, lower = FALSE))
    All.eq(Rtrgamma,      qtrgamma     (1 - Ptrgamma,      shape1 = 2, shape2 = 3, scale = 5, lower = FALSE))
    All.eq(Rinvtrgamma,   qinvtrgamma  (1 - Pinvtrgamma,   shape1 = 2, shape2 = 3, scale = 5, lower = FALSE))
    All.eq(Rinvgamma,     qinvgamma    (1 - Pinvgamma,     shape  = 2,             scale = 5, lower = FALSE))
    All.eq(Rinvweibull,   qinvweibull  (1 - Pinvweibull,               shape  = 3, scale = 5, lower = FALSE))
    All.eq(Rinvexp,       qinvexp      (1 - Pinvexp,                               scale = 5, lower = FALSE))
    All.eq(Rlgamma, qlgamma(1 - Plgamma, shapelog = 1.5, ratelog = 5, lower = FALSE))
    All.eq(Rgumbel, qgumbel(1 - Pgumbel, alpha = 2, scale = 5, lower = FALSE))
    All.eq(Rinvgauss, qinvgauss(1 - Pinvgauss, mean = 2, dispersion = 5, lower = FALSE))
    All.eq(Rgenbeta, qgenbeta(1 - Pgenbeta, shape1 = 0.8, shape2 = 1.5, shape3 = 2, scale = 2, lower = FALSE))
})

## Check q<dist>(p<dist>(., log), log) identity
stopifnot(exprs = {
    All.eq(Rfpareto, qfpareto(log(Pfpareto), min = m, shape1 = 0.8, shape2 = 1.5, shape3 = 2, scale = 2, log = TRUE))
    All.eq(Rpareto4, qpareto4(log(Ppareto4), min = m, shape1 = 0.8, shape2 = 1.5,             scale = 2, log = TRUE))
    All.eq(Rpareto3, qpareto3(log(Ppareto3), min = m,               shape  = 1.5,             scale = 2, log = TRUE))
    All.eq(Rpareto2, qpareto2(log(Ppareto2), min = m, shape  = 0.8,                           scale = 2, log = TRUE))
    All.eq(Rtrbeta,       qtrbeta      (log(Ptrbeta),       shape1 = 0.8, shape2 = 1.5, shape3 = 2, scale = 2, log = TRUE))
    All.eq(Rburr,         qburr        (log(Pburr),         shape1 = 0.8, shape2 = 1.5,             scale = 2, log = TRUE))
    All.eq(Rllogis,       qllogis      (log(Pllogis),       shape  = 1.5, scale = 2, log = TRUE))
    All.eq(Rparalogis,    qparalogis   (log(Pparalogis),    shape  = 0.8, scale = 2, log = TRUE))
    All.eq(Rgenpareto,    qgenpareto   (log(Pgenpareto),    shape1 = 0.8, shape2 = 2, scale = 2, log = TRUE))
    All.eq(Rpareto,       qpareto      (log(Ppareto),       shape  = 0.8, scale = 2, log = TRUE))
    All.eq(Rpareto1,      qpareto1     (log(Ppareto1),      shape  = 0.8, min = 2, log = TRUE))
    All.eq(Rinvburr,      qinvburr     (log(Pinvburr),      shape1 = 1.5, shape2 = 2, scale = 2, log = TRUE))
    All.eq(Rinvpareto,    qinvpareto   (log(Pinvpareto),    shape  = 2,   scale = 2, log = TRUE))
    All.eq(Rinvparalogis, qinvparalogis(log(Pinvparalogis), shape  = 2,   scale = 2, log = TRUE))
    All.eq(Rtrgamma,      qtrgamma     (log(Ptrgamma),      shape1 = 2, shape2 = 3, scale = 5, log = TRUE))
    All.eq(Rinvtrgamma,   qinvtrgamma  (log(Pinvtrgamma),   shape1 = 2, shape2 = 3, scale = 5, log = TRUE))
    All.eq(Rinvgamma,     qinvgamma    (log(Pinvgamma),     shape  = 2,             scale = 5, log = TRUE))
    All.eq(Rinvweibull,   qinvweibull  (log(Pinvweibull),               shape  = 3, scale = 5, log = TRUE))
    All.eq(Rinvexp,       qinvexp      (log(Pinvexp),                               scale = 5, log = TRUE))
    All.eq(Rlgamma, qlgamma(log(Plgamma), shapelog = 1.5, ratelog = 5, log = TRUE))
    All.eq(Rgumbel, qgumbel(log(Pgumbel), alpha = 2, scale = 5, log = TRUE))
    All.eq(Rinvgauss, qinvgauss(log(Pinvgauss), mean = 2, dispersion = 5, log = TRUE))
    All.eq(Rgenbeta, qgenbeta(log(Pgenbeta), shape1 = 0.8, shape2 = 1.5, shape3 = 2, scale = 2, log = TRUE))
})

## Check q<dist>(p<dist>(., log), log) identity with upper tail
stopifnot(exprs = {
    All.eq(Rfpareto, qfpareto(log1p(-Pfpareto), min = m, shape1 = 0.8, shape2 = 1.5, shape3 = 2, scale = 2, lower = FALSE, log = TRUE))
    All.eq(Rpareto4, qpareto4(log1p(-Ppareto4), min = m, shape1 = 0.8, shape2 = 1.5,             scale = 2, lower = FALSE, log = TRUE))
    All.eq(Rpareto3, qpareto3(log1p(-Ppareto3), min = m,               shape  = 1.5,             scale = 2, lower = FALSE, log = TRUE))
    All.eq(Rpareto2, qpareto2(log1p(-Ppareto2), min = m, shape  = 0.8,                           scale = 2, lower = FALSE, log = TRUE))
    All.eq(Rtrbeta,       qtrbeta      (log1p(-Ptrbeta),       shape1 = 0.8, shape2 = 1.5, shape3 = 2, scale = 2, lower = FALSE, log = TRUE))
    All.eq(Rburr,         qburr        (log1p(-Pburr),         shape1 = 0.8, shape2 = 1.5,             scale = 2, lower = FALSE, log = TRUE))
    All.eq(Rllogis,       qllogis      (log1p(-Pllogis),       shape  = 1.5, scale = 2, lower = FALSE, log = TRUE))
    All.eq(Rparalogis,    qparalogis   (log1p(-Pparalogis),    shape  = 0.8, scale = 2, lower = FALSE, log = TRUE))
    All.eq(Rgenpareto,    qgenpareto   (log1p(-Pgenpareto),    shape1 = 0.8, shape2 = 2, scale = 2, lower = FALSE, log = TRUE))
    All.eq(Rpareto,       qpareto      (log1p(-Ppareto),       shape  = 0.8, scale = 2, lower = FALSE, log = TRUE))
    All.eq(Rpareto1,      qpareto1     (log1p(-Ppareto1),      shape  = 0.8, min = 2, lower = FALSE, log = TRUE))
    All.eq(Rinvburr,      qinvburr     (log1p(-Pinvburr),      shape1 = 1.5, shape2 = 2, scale = 2, lower = FALSE, log = TRUE))
    All.eq(Rinvpareto,    qinvpareto   (log1p(-Pinvpareto),    shape  = 2,   scale = 2, lower = FALSE, log = TRUE))
    All.eq(Rinvparalogis, qinvparalogis(log1p(-Pinvparalogis), shape  = 2,   scale = 2, lower = FALSE, log = TRUE))
    All.eq(Rtrgamma,      qtrgamma     (log1p(-Ptrgamma),      shape1 = 2, shape2 = 3, scale = 5, lower = FALSE, log = TRUE))
    All.eq(Rinvtrgamma,   qinvtrgamma  (log1p(-Pinvtrgamma),   shape1 = 2, shape2 = 3, scale = 5, lower = FALSE, log = TRUE))
    All.eq(Rinvgamma,     qinvgamma    (log1p(-Pinvgamma),     shape  = 2,             scale = 5, lower = FALSE, log = TRUE))
    All.eq(Rinvweibull,   qinvweibull  (log1p(-Pinvweibull),               shape  = 3, scale = 5, lower = FALSE, log = TRUE))
    All.eq(Rinvexp,       qinvexp      (log1p(-Pinvexp),                               scale = 5, lower = FALSE, log = TRUE))
    All.eq(Rlgamma, qlgamma(log1p(-Plgamma), shapelog = 1.5, ratelog = 5, lower = FALSE, log = TRUE))
    All.eq(Rgumbel, qgumbel(log1p(-Pgumbel), alpha = 2, scale = 5, lower = FALSE, log = TRUE))
    All.eq(Rinvgauss, qinvgauss(log1p(-Pinvgauss), mean = 2, dispersion = 5, lower = FALSE, log = TRUE))
    All.eq(Rgenbeta, qgenbeta(log1p(-Pgenbeta), shape1 = 0.8, shape2 = 1.5, shape3 = 2, scale = 2, lower = FALSE, log = TRUE))
})


###
### DISCRETE DISTRIBUTIONS
###

## Reset seed
set.seed(123)

## Define a small function to compute probabilities for the (a, b, 1)
## family of discrete distributions using the recursive relation
##
##   p[k] = (a + b/k)p[k - 1], k = 2, 3, ...
##
## for a, b and p[1] given.
dab1 <- function(x, a, b, p1)
{
    x <- floor(x)
    if (x < 1)
        stop("recursive computations possible for x >= 2 only")
    for (k in seq(2, length.out = x - 1))
    {
        p2 <- (a + b/k) * p1
        p1 <- p2
    }
    p1
}

## ZERO-TRUNCATED (a, b, 1) CLASS

## Tests on the probability mass function:
##
## 1. probability is 0 at x = 0;
## 2. pmf satisfies the recursive relation
lambda <- rlnorm(30, 2)                 # Poisson parameters
r <- lambda                             # size for negative binomial
prob <- runif(30)                       # probs
size <- round(lambda)                   # size for binomial
stopifnot(exprs = {
    dztpois(0, lambda) == 0
    dztnbinom(0, r, prob) == 0
    dztgeom(0, prob) == 0
    dztbinom(0, size, prob) == 0
    dlogarithmic(0, prob) == 0
})

x <- sapply(size, sample, size = 1)
stopifnot(exprs = {
    All.eq(dztpois(x, lambda),
           mapply(dab1, x,
                  a = 0,
                  b = lambda,
                  p1 = lambda/(exp(lambda) - 1)))
    All.eq(dztnbinom(x, r, prob),
           mapply(dab1, x,
                  a = 1 - prob,
                  b = (r - 1) * (1 - prob),
                  p1 = r * prob^r * (1 - prob)/(1 - prob^r)))
    All.eq(dztgeom(x, prob),
           mapply(dab1, x,
                  a = 1 - prob,
                  b = 0,
                  p1 = prob))
    All.eq(dztbinom(x, size, prob),
           mapply(dab1, x,
                  a = -prob/(1 - prob),
                  b = (size + 1) * prob/(1 - prob),
                  p1 = size * prob * (1 - prob)^(size - 1)/(1 - (1 - prob)^size)))
    All.eq(dlogarithmic(x, prob),
           mapply(dab1, x,
                  a = prob,
                  b = -prob,
                  p1 = -prob/log1p(-prob)))
})

## Tests on cumulative distribution function.
for (l in lambda)
    stopifnot(exprs = {
        all.equal(cumsum(dztpois(0:20, l)),
                  pztpois(0:20, l),
                  tolerance = 1e-8)
    })
for (i in seq_along(r))
    stopifnot(exprs = {
        all.equal(cumsum(dztnbinom(0:20, r[i], prob[i])),
                  pztnbinom(0:20, r[i], prob[i]),
                  tolerance = 1e-8)
    })
for (i in seq_along(r))
    stopifnot(exprs = {
        all.equal(cumsum(dztgeom(0:20, prob[i])),
                  pztgeom(0:20, prob[i]),
                  tolerance = 1e-8)
    })
for (i in seq_along(size))
    stopifnot(exprs = {
        all.equal(cumsum(dztbinom(0:20, size[i], prob[i])),
                  pztbinom(0:20, size[i], prob[i]),
                  tolerance = 1e-8)
    })
for (p in prob)
    stopifnot(exprs = {
        all.equal(cumsum(dlogarithmic(0:20, p)),
                  plogarithmic(0:20, p),
                  tolerance = 1e-8)
    })

## ZERO-MODIFIED (a, b, 1) CLASS

## Tests on the probability mass function:
##
## 1. probability is p0 at x = 0 (trivial, but...);
## 2. pmf satisfies the recursive relation
lambda <- rlnorm(30, 2)                 # Poisson parameters
r <- lambda                             # size for negative binomial
prob <- runif(30)                       # probs
size <- round(lambda)                   # size for binomial
p0 <- runif(30)                         # probs at 0
stopifnot(exprs = {
    dzmpois(0, lambda, p0) == p0
    dzmnbinom(0, r, prob, p0) == p0
    dzmgeom(0, prob, p0) == p0
    dzmbinom(0, size, prob, p0) == p0
    dzmlogarithmic(0, prob, p0) == p0
})

x <- sapply(size, sample, size = 1)
stopifnot(exprs = {
    All.eq(dzmpois(x, lambda, p0),
           mapply(dab1, x,
                  a = 0,
                  b = lambda,
                  p1 = (1 - p0) *lambda/(exp(lambda) - 1)))
    All.eq(dzmnbinom(x, r, prob, p0),
           mapply(dab1, x,
                  a = 1 - prob,
                  b = (r - 1) * (1 - prob),
                  p1 = (1 - p0) * r * prob^r * (1 - prob)/(1 - prob^r)))
    All.eq(dzmgeom(x, prob, p0),
           mapply(dab1, x,
                  a = 1 - prob,
                  b = 0,
                  p1 = (1 - p0) * prob))
    All.eq(dzmbinom(x, size, prob, p0),
           mapply(dab1, x,
                  a = -prob/(1 - prob),
                  b = (size + 1) * prob/(1 - prob),
                  p1 = (1 - p0) * size * prob * (1 - prob)^(size - 1)/(1 - (1 - prob)^size)))
    All.eq(dzmlogarithmic(x, prob, p0),
           mapply(dab1, x,
                  a = prob,
                  b = -prob,
                  p1 = -(1 - p0) * prob/log1p(-prob)))
})

## Tests on cumulative distribution function.
for (i in seq_along(lambda))
    stopifnot(exprs = {
        all.equal(cumsum(dzmpois(0:20, lambda[i], p0 = p0[i])),
                  pzmpois(0:20, lambda[i], p0 = p0[i]),
                  tolerance = 1e-8)
    })
for (i in seq_along(r))
    stopifnot(exprs = {
        all.equal(cumsum(dzmnbinom(0:20, r[i], prob[i], p0[i])),
                  pzmnbinom(0:20, r[i], prob[i], p0[i]),
                  tolerance = 1e-8)
    })
for (i in seq_along(r))
    stopifnot(exprs = {
        all.equal(cumsum(dzmgeom(0:20, prob[i], p0[i])),
                  pzmgeom(0:20, prob[i], p0[i]),
                  tolerance = 1e-8)
    })
for (i in seq_along(size))
    stopifnot(exprs = {
        all.equal(cumsum(dzmbinom(0:20, size[i], prob[i], p0[i])),
                  pzmbinom(0:20, size[i], prob[i], p0[i]),
                  tolerance = 1e-8)
    })
for (i in seq_along(prob))
    stopifnot(exprs = {
        all.equal(cumsum(dzmlogarithmic(0:20, prob[i], p0[i])),
                  pzmlogarithmic(0:20, prob[i], p0[i]),
                  tolerance = 1e-8)
    })

## POISSON-INVERSE GAUSSIAN

## Reset seed
set.seed(123)

## Define a small function to compute probabilities for the PIG
## directly using the Bessel function.
dpigBK <- function(x, mu, phi)
{
    M_LN2 <- 0.693147180559945309417232121458
    M_SQRT_2dPI <- 0.225791352644727432363097614947

    phimu  <- phi * mu
    lphi <- log(phi)
    y <- x - 0.5

    logA = -lphi/2 - M_SQRT_2dPI
    logB = (M_LN2 + lphi + log1p(1/(2 * phimu * mu)))/2;

    exp(logA + 1/phimu - lfactorial(x) - y * logB) *
        besselK(exp(logB - lphi), y)
}

## Tests on the probability mass function.
mu <- rlnorm(30, 2)
phi <- rlnorm(30, 2)
x <- 0:100
for (i in seq_along(phi))
{
    stopifnot(exprs = {
        all.equal(dpoisinvgauss(x, mean = mu[i], dispersion = phi[i]),
                  dpigBK(x, mu[i], phi[i]))
        all.equal(dpoisinvgauss(x, mean = Inf, dispersion = phi[i]),
                  dpigBK(x, Inf, phi[i]))
    })
}

## Tests on cumulative distribution function.
for (i in seq_along(phi))
    stopifnot(exprs = {
        all.equal(cumsum(dpoisinvgauss(0:20, mu[i], phi[i])),
                  ppoisinvgauss(0:20, mu[i], phi[i]),
                  tolerance = 1e-8)
        all.equal(cumsum(dpoisinvgauss(0:20, Inf, phi[i])),
                  ppoisinvgauss(0:20, Inf, phi[i]),
                  tolerance = 1e-8)
    })

##
## RANDOM NUMBERS (all discrete distributions)
##
set.seed(123)
n <- 20

## Generate variates.
##
## For zero-modified distributions, we simulate two sets of values:
## one with p0m < p0 (suffix 'p0lt') and one with p0m > p0 (suffix
## 'p0gt').
Rztpois      <- rztpois     (n, lambda = 12)
Rztnbinom    <- rztnbinom   (n, size = 7, prob = 0.01)
Rztgeom      <- rztgeom     (n, prob = pi/16)
Rztbinom     <- rztbinom    (n, size = 55, prob = pi/16)
Rlogarithmic <- rlogarithmic(n, prob = 0.99)
Rzmpoisp0lt        <- rzmpois       (n, lambda = 6, p0 = 0.001)
Rzmpoisp0gt        <- rzmpois       (n, lambda = 6, p0 = 0.010)
Rzmnbinomp0lt      <- rzmnbinom     (n, size = 7, prob = 0.8, p0 = 0.01)
Rzmnbinomp0gt      <- rzmnbinom     (n, size = 7, prob = 0.8, p0 = 0.40)
Rzmgeomp0lt        <- rzmgeom       (n, prob = pi/16, p0 = 0.01)
Rzmgeomp0gt        <- rzmgeom       (n, prob = pi/16, p0 = 0.40)
Rzmbinomp0lt       <- rzmbinom      (n, size = 12, prob = pi/16, p0 = 0.01)
Rzmbinomp0gt       <- rzmbinom      (n, size = 12, prob = pi/16, p0 = 0.12)
Rzmlogarithmicp0lt <- rzmlogarithmic(n, prob = 0.99, p0 = 0.05)
Rzmlogarithmicp0gt <- rzmlogarithmic(n, prob = 0.99, p0 = 0.55)
Rpoisinvgauss    <- rpoisinvgauss(n, mean = 12,  dispersion = 0.1)
RpoisinvgaussInf <- rpoisinvgauss(n, mean = Inf, dispersion = 1.1)

## Compute quantiles
Pztpois      <- pztpois     (Rztpois,      lambda = 12)
Pztnbinom    <- pztnbinom   (Rztnbinom,    size = 7, prob = 0.01)
Pztgeom      <- pztgeom     (Rztgeom,      prob = pi/16)
Pztbinom     <- pztbinom    (Rztbinom,     size = 55, prob = pi/16)
Plogarithmic <- plogarithmic(Rlogarithmic, prob = 0.99)
Pzmpoisp0lt        <- pzmpois       (Rzmpoisp0lt,        lambda = 6, p0 = 0.001)
Pzmpoisp0gt        <- pzmpois       (Rzmpoisp0gt,        lambda = 6, p0 = 0.010)
Pzmnbinomp0lt      <- pzmnbinom     (Rzmnbinomp0lt,      size = 7, prob = 0.8, p0 = 0.01)
Pzmnbinomp0gt      <- pzmnbinom     (Rzmnbinomp0gt,      size = 7, prob = 0.8, p0 = 0.40)
Pzmgeomp0lt        <- pzmgeom       (Rzmgeomp0lt,        prob = pi/16, p0 = 0.01)
Pzmgeomp0gt        <- pzmgeom       (Rzmgeomp0gt,        prob = pi/16, p0 = 0.40)
Pzmbinomp0lt       <- pzmbinom      (Rzmbinomp0lt,       size = 12, prob = pi/16, p0 = 0.01)
Pzmbinomp0gt       <- pzmbinom      (Rzmbinomp0gt,       size = 12, prob = pi/16, p0 = 0.12)
Pzmlogarithmicp0lt <- pzmlogarithmic(Rzmlogarithmicp0lt, prob = 0.99, p0 = 0.05)
Pzmlogarithmicp0gt <- pzmlogarithmic(Rzmlogarithmicp0gt, prob = 0.99, p0 = 0.55)
Ppoisinvgauss    <- ppoisinvgauss(Rpoisinvgauss,    mean = 12,  dispersion = 0.1)
PpoisinvgaussInf <- ppoisinvgauss(RpoisinvgaussInf, mean = Inf, dispersion = 1.1)

## Just compute pmf
Dztpois      <- dztpois     (Rztpois,      lambda = 12)
Dztnbinom    <- dztnbinom   (Rztnbinom,    size = 7, prob = 0.01)
Dztgeom      <- dztgeom     (Rztgeom,      prob = pi/16)
Dztbinom     <- dztbinom    (Rztbinom,     size = 55, prob = pi/16)
Dlogarithmic <- dlogarithmic(Rlogarithmic, prob = pi/16)
Dzmpoisp0lt        <- dzmpois       (Rzmpoisp0lt,        lambda = 6, p0 = 0.001)
Dzmpoisp0gt        <- dzmpois       (Rzmpoisp0gt,        lambda = 6, p0 = 0.010)
Dzmnbinomp0lt      <- dzmnbinom     (Rzmnbinomp0lt,      size = 7, prob = 0.8, p0 = 0.01)
Dzmnbinomp0gt      <- dzmnbinom     (Rzmnbinomp0gt,      size = 7, prob = 0.8, p0 = 0.40)
Dzmgeomp0lt        <- dzmgeom       (Rzmgeomp0lt,        prob = pi/16, p0 = 0.01)
Dzmgeomp0gt        <- dzmgeom       (Rzmgeomp0gt,        prob = pi/16, p0 = 0.40)
Dzmbinomp0lt       <- dzmbinom      (Rzmbinomp0lt,       size = 12, prob = pi/16, p0 = 0.01)
Dzmbinomp0gt       <- dzmbinom      (Rzmbinomp0gt,       size = 12, prob = pi/16, p0 = 0.12)
Dzmlogarithmicp0lt <- dzmlogarithmic(Rzmlogarithmicp0lt, prob = 0.99, p0 = 0.05)
Dzmlogarithmicp0gt <- dzmlogarithmic(Rzmlogarithmicp0gt, prob = 0.99, p0 = 0.55)
Dpoisinvgauss    <- dpoisinvgauss(Rpoisinvgauss,    mean = 12,  dispersion = 0.1)
DpoisinvgaussInf <- dpoisinvgauss(RpoisinvgaussInf, mean = Inf, dispersion = 1.1)

## Check q<dist>(p<dist>(.)) identity
stopifnot(exprs = {
    Rztpois      == qztpois     (Pztpois,      lambda = 12)
    Rztnbinom    == qztnbinom   (Pztnbinom,    size = 7, prob = 0.01)
    Rztgeom      == qztgeom     (Pztgeom,      prob = pi/16)
    Rztbinom     == qztbinom    (Pztbinom,     size = 55, prob = pi/16)
    Rlogarithmic == qlogarithmic(Plogarithmic, prob = 0.99)
    Rzmpoisp0lt        == qzmpois       (Pzmpoisp0lt,        lambda = 6, p0 = 0.001)
    Rzmpoisp0gt        == qzmpois       (Pzmpoisp0gt,        lambda = 6, p0 = 0.010)
    Rzmnbinomp0lt      == qzmnbinom     (Pzmnbinomp0lt,      size = 7, prob = 0.8, p0 = 0.01)
    Rzmnbinomp0gt      == qzmnbinom     (Pzmnbinomp0gt,      size = 7, prob = 0.8, p0 = 0.40)
    Rzmgeomp0lt        == qzmgeom       (Pzmgeomp0lt,        prob = pi/16, p0 = 0.01)
    Rzmgeomp0gt        == qzmgeom       (Pzmgeomp0gt,        prob = pi/16, p0 = 0.40)
    Rzmbinomp0lt       == qzmbinom      (Pzmbinomp0lt,       size = 12, prob = pi/16, p0 = 0.01)
    Rzmbinomp0gt       == qzmbinom      (Pzmbinomp0gt,       size = 12, prob = pi/16, p0 = 0.12)
    Rzmlogarithmicp0lt == qzmlogarithmic(Pzmlogarithmicp0lt, prob = 0.99, p0 = 0.05)
    Rzmlogarithmicp0gt == qzmlogarithmic(Pzmlogarithmicp0gt, prob = 0.99, p0 = 0.55)
    Rpoisinvgauss    == qpoisinvgauss(Ppoisinvgauss,    mean = 12,  dispersion = 0.1)
    RpoisinvgaussInf == qpoisinvgauss(PpoisinvgaussInf, mean = Inf, dispersion = 1.1)
})

## Check q<dist>(p<dist>(.)) identity with upper tail
stopifnot(exprs = {
    Rztpois      == qztpois     (1 - Pztpois,      lambda = 12, lower = FALSE)
    Rztnbinom    == qztnbinom   (1 - Pztnbinom,    size = 7, prob = 0.01, lower = FALSE)
    Rztgeom      == qztgeom     (1 - Pztgeom,      prob = pi/16, lower = FALSE)
    Rztbinom     == qztbinom    (1 - Pztbinom,     size = 55, prob = pi/16, lower = FALSE)
    Rlogarithmic == qlogarithmic(1 - Plogarithmic, prob = 0.99, lower = FALSE)
    Rzmpoisp0lt        == qzmpois       (1 - Pzmpoisp0lt,        lambda = 6, p0 = 0.001, lower = FALSE)
    Rzmpoisp0gt        == qzmpois       (1 - Pzmpoisp0gt,        lambda = 6, p0 = 0.010, lower = FALSE)
    Rzmnbinomp0lt      == qzmnbinom     (1 - Pzmnbinomp0lt,      size = 7, prob = 0.8, p0 = 0.01, lower = FALSE)
    Rzmnbinomp0gt      == qzmnbinom     (1 - Pzmnbinomp0gt,      size = 7, prob = 0.8, p0 = 0.40, lower = FALSE)
    Rzmgeomp0lt        == qzmgeom       (1 - Pzmgeomp0lt,        prob = pi/16, p0 = 0.01, lower = FALSE)
    Rzmgeomp0gt        == qzmgeom       (1 - Pzmgeomp0gt,        prob = pi/16, p0 = 0.40, lower = FALSE)
    Rzmbinomp0lt       == qzmbinom      (1 - Pzmbinomp0lt,       size = 12, prob = pi/16, p0 = 0.01, lower = FALSE)
    Rzmbinomp0gt       == qzmbinom      (1 - Pzmbinomp0gt,       size = 12, prob = pi/16, p0 = 0.12, lower = FALSE)
    Rzmlogarithmicp0lt == qzmlogarithmic(1 - Pzmlogarithmicp0lt, prob = 0.99, p0 = 0.05, lower = FALSE)
    Rzmlogarithmicp0gt == qzmlogarithmic(1 - Pzmlogarithmicp0gt, prob = 0.99, p0 = 0.55, lower = FALSE)
    Rpoisinvgauss    == qpoisinvgauss(1 - Ppoisinvgauss,    mean = 12,  dispersion = 0.1, lower = FALSE)
    RpoisinvgaussInf == qpoisinvgauss(1 - PpoisinvgaussInf, mean = Inf, dispersion = 1.1, lower = FALSE)
})

## Check q<dist>(p<dist>(., log), log) identity
stopifnot(exprs = {
    Rztpois      == qztpois     (log(Pztpois),      lambda = 12, log = TRUE)
    Rztnbinom    == qztnbinom   (log(Pztnbinom),    size = 7, prob = 0.01, log = TRUE)
    Rztgeom      == qztgeom     (log(Pztgeom),      prob = pi/16, log = TRUE)
    Rztbinom     == qztbinom    (log(Pztbinom),     size = 55, prob = pi/16, log = TRUE)
    Rlogarithmic == qlogarithmic(log(Plogarithmic), prob = 0.99, log = TRUE)
    Rzmpoisp0lt        == qzmpois       (log(Pzmpoisp0lt),        lambda = 6, p0 = 0.001, log = TRUE)
    Rzmpoisp0gt        == qzmpois       (log(Pzmpoisp0gt),        lambda = 6, p0 = 0.010, log = TRUE)
    Rzmnbinomp0lt      == qzmnbinom     (log(Pzmnbinomp0lt),      size = 7, prob = 0.8, p0 = 0.01, log = TRUE)
    Rzmnbinomp0gt      == qzmnbinom     (log(Pzmnbinomp0gt),      size = 7, prob = 0.8, p0 = 0.40, log = TRUE)
    Rzmgeomp0lt        == qzmgeom       (log(Pzmgeomp0lt),        prob = pi/16, p0 = 0.01, log = TRUE)
    Rzmgeomp0gt        == qzmgeom       (log(Pzmgeomp0gt),        prob = pi/16, p0 = 0.40, log = TRUE)
    Rzmbinomp0lt       == qzmbinom      (log(Pzmbinomp0lt),       size = 12, prob = pi/16, p0 = 0.01, log = TRUE)
    Rzmbinomp0gt       == qzmbinom      (log(Pzmbinomp0gt),       size = 12, prob = pi/16, p0 = 0.12, log = TRUE)
    Rzmlogarithmicp0lt == qzmlogarithmic(log(Pzmlogarithmicp0lt), prob = 0.99, p0 = 0.05, log = TRUE)
    Rzmlogarithmicp0gt == qzmlogarithmic(log(Pzmlogarithmicp0gt), prob = 0.99, p0 = 0.55, log = TRUE)
    Rpoisinvgauss    == qpoisinvgauss(log(Ppoisinvgauss),    mean = 12,  dispersion = 0.1, log = TRUE)
    RpoisinvgaussInf == qpoisinvgauss(log(PpoisinvgaussInf), mean = Inf, dispersion = 1.1, log = TRUE)
})

## Check q<dist>(p<dist>(., log), log) identity with upper tail
stopifnot(exprs = {
    Rztpois      == qztpois     (log1p(-Pztpois),      lambda = 12, lower = FALSE, log = TRUE)
    Rztnbinom    == qztnbinom   (log1p(-Pztnbinom),    size = 7, prob = 0.01, lower = FALSE, log = TRUE)
    Rztgeom      == qztgeom     (log1p(-Pztgeom),      prob = pi/16, lower = FALSE, log = TRUE)
    Rztbinom     == qztbinom    (log1p(-Pztbinom),     size = 55, prob = pi/16, lower = FALSE, log = TRUE)
    Rlogarithmic == qlogarithmic(log1p(-Plogarithmic), prob = 0.99, lower = FALSE, log = TRUE)
    Rzmpoisp0lt        == qzmpois       (log1p(-Pzmpoisp0lt),        lambda = 6, p0 = 0.001, lower = FALSE, log = TRUE)
    Rzmpoisp0gt        == qzmpois       (log1p(-Pzmpoisp0gt),        lambda = 6, p0 = 0.010, lower = FALSE, log = TRUE)
    Rzmnbinomp0lt      == qzmnbinom     (log1p(-Pzmnbinomp0lt),      size = 7, prob = 0.8, p0 = 0.01, lower = FALSE, log = TRUE)
    Rzmnbinomp0gt      == qzmnbinom     (log1p(-Pzmnbinomp0gt),      size = 7, prob = 0.8, p0 = 0.40, lower = FALSE, log = TRUE)
    Rzmgeomp0lt        == qzmgeom       (log1p(-Pzmgeomp0lt),        prob = pi/16, p0 = 0.01, lower = FALSE, log = TRUE)
    Rzmgeomp0gt        == qzmgeom       (log1p(-Pzmgeomp0gt),        prob = pi/16, p0 = 0.40, lower = FALSE, log = TRUE)
    Rzmbinomp0lt       == qzmbinom      (log1p(-Pzmbinomp0lt),       size = 12, prob = pi/16, p0 = 0.01, lower = FALSE, log = TRUE)
    Rzmbinomp0gt       == qzmbinom      (log1p(-Pzmbinomp0gt),       size = 12, prob = pi/16, p0 = 0.12, lower = FALSE, log = TRUE)
    Rzmlogarithmicp0lt == qzmlogarithmic(log1p(-Pzmlogarithmicp0lt), prob = 0.99, p0 = 0.05, lower = FALSE, log = TRUE)
    Rzmlogarithmicp0gt == qzmlogarithmic(log1p(-Pzmlogarithmicp0gt), prob = 0.99, p0 = 0.55, lower = FALSE, log = TRUE)
    Rpoisinvgauss    == qpoisinvgauss(log1p(-Ppoisinvgauss),    mean = 12,  dispersion = 0.1, lower = FALSE, log = TRUE)
    RpoisinvgaussInf == qpoisinvgauss(log1p(-PpoisinvgaussInf), mean = Inf, dispersion = 1.1, lower = FALSE, log = TRUE)
})

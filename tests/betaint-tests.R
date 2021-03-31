### actuar: Actuarial Functions and Heavy Tailed Distributions
###
### Tests for the "beta integral"
###
###    B(a, b; x) = Gamma(a + b) int_0^x t^(a-1) (1 - t)^(b-1) dt
###
### Inspired by (and some parts taken from) `tests/d-p-q-r-tests.R` in
### R sources.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca> with
###         indirect help from the R Core Team

## Load the package
library(actuar)

## Define a "local" version of the otherwise non-exported function
## 'betaint'.
betaint <- actuar:::betaint

## Special values and utilities. Taken from `tests/d-p-q-r-tests.R`.
xMax <- 1 - .Machine$double.eps
xMin <- .Machine$double.xmin
All.eq <- function(x, y)
{
    all.equal.numeric(x, y, tolerance = 64 * .Machine$double.eps,
                      scale = max(0, mean(abs(x), na.rm = TRUE)))
}

if(!interactive())
    set.seed(123)

## Limiting cases
stopifnot(exprs = {
    !is.finite(betaint(0.3, Inf, 2))
    !is.finite(betaint(0.3, Inf, -2.2))
    is.nan    (betaint(0.3,   0, 2))
    !is.finite(betaint(0.3,   2, Inf))
    is.nan    (betaint(0.3,   2, -2.2)) # a <= 1 + floor(-b)
    is.nan    (betaint(0.3,   2, 0))
})

## Tests for cases with b > 0
x <- c(xMin, runif(10), xMax)
b <- 2
for (a in rlnorm(5, 2))
    stopifnot(exprs = {
        All.eq(betaint(x, a, b),
               gamma(a) * gamma(b) * pbeta(x, a, b))
    })

## Tests for cases with b < 0
b <- -2.2
r <- floor(-b)        # r = 2
for (a in 1 + r + rlnorm(5, 2))
{
    s <- (x^(a-1) * (1-x)^b)/b +
        ((a-1) * x^(a-2) * (1-x)^(b+1))/(b * (b+1)) +
        ((a-1) * (a-2) * x^(a-3) * (1-x)^(b+2))/(b * (b+1) * (b+2))
    stopifnot(exprs = {
        all.equal(betaint(x, a, b),
                  -gamma(a+b) * s +
                  (a-1)*(a-2)*(a-3) * gamma(a-r-1)/(b*(b+1)*(b+2)) *
                  gamma(b+r+1)*pbeta(x, a-r-1, b+r+1))
    })
}

### actuar: Actuarial Functions and Heavy Tailed Distributions
###
### Tests for the 'var' and 'sd' methods for individual and grouped
### data.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

## Load the package
library(actuar)

###
### Individual data
###

## Check that results are identical to stats::var and stats::sd, as it
## should be.
stopifnot(exprs = {
    identical(var(dental), stats::var(dental))
    identical(sd(dental), stats::sd(dental))
})

## Check correct handling of missing data (issue #5 fixed with
## 62b92eff).
x <- c(dental, NA)
stopifnot(exprs = {
    identical(var(x, na.rm = TRUE), stats::var(x, na.rm = TRUE))
    identical(sd(x, na.rm = TRUE), stats::sd(x, na.rm = TRUE))
})

###
### Grouped data
###

## Extract group boundaries and frequencies from a grouped data
## object.
cj <- gdental[, 1]
nj <- gdental[, 2]

## Compute variance and standard deviation by hand.
midpoints <- cj[-length(cj)] + diff(cj)/2
means <- drop(crossprod(nj, midpoints)/sum(nj))
v <- drop(crossprod(nj, (midpoints - means)^2)/(sum(nj) - 1))
s <- sqrt(v)
stopifnot(exprs = {
    all.equal(v, unname(var(gdental)))
    all.equal(s, unname(sd(gdental)))
})

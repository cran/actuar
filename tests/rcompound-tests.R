### actuar: Actuarial Functions and Heavy Tailed Distributions
###
### Tests for the simulation of compound models with 'rcompound' and
### 'rcomppois'.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

## Load the package
library(actuar)

## Copy of tools::assertError.
assertError <- tools::assertError

###
### Tests for rcompound
###

## Test the function itself with various types of arguments.
n <- 20
fmodel <- expression(rnbinom(2, 0.8))
smodel <- expression(rgamma(2, 1))
set.seed(123)
x <- numeric(n)
N <- rnbinom(n, 2, 0.8)
y <- rgamma(sum(N), 2, 1)
x[which(N != 0)] <- tapply(y, rep(seq_len(n), N), sum)
stopifnot(exprs = {
    identical(x, {
        set.seed(123)
        rcompound(n, rnbinom(2, 0.8), rgamma(2, 1))
    })
    identical(x, {
        set.seed(123)
        rcompound(n, rnbinom(2, 0.8), expression(rgamma(2, 1)))
    })
    identical(x, {
        set.seed(123)
        rcompound(n, expression(rnbinom(2, 0.8)), rgamma(2, 1))
    })
    identical(x, {
        set.seed(123)
        rcompound(n, fmodel, smodel)
    })
})

## Test the calling environment, that is that arguments are correctly
## identified when 'rcompound' is called inside another function.
n <- 20
lambda <- 2
smodel <- expression(rgamma(2, 1))
set.seed(123)
x <- rcompound(n, rpois(2), rgamma(2, 1))
f <- function(n, p, model.sev)
{
    ## safe way to pass down the arguments
    model.freq <- substitute(rpois(p), list(p = p))
    model.sev <- substitute(model.sev)
    if (is.name(model.sev))
        model.sev <- eval.parent(model.sev)
    rcompound(n, model.freq, model.sev)
}
g1 <- function(n, p, s, r)
    rcompound(n, rpois(p), rgamma(s, r))
g2 <- function(n, p, s, r)
    rcompound(n, expression(rpois(p)), expression(rgamma(s, r)))
h <- function(n, p, model.sev)
{
    ## safe way to pass down the arguments
    model.sev <- substitute(model.sev)
    if (is.name(model.sev))
        model.sev <- eval.parent(model.sev)
    f(n, p, model.sev)
}
stopifnot(exprs = {
    identical(x, {
        set.seed(123)
        f(n, 2, rgamma(2, 1))
    })
    identical(x, {
        set.seed(123)
        f(n, lambda, expression(rgamma(2, 1)))
    })
    identical(x, {
        set.seed(123)
        f(n, lambda, smodel)
    })
    identical(x, {
        set.seed(123)
        g1(n, lambda, 2, 1)
    })
    identical(x, {
        set.seed(123)
        g2(n, lambda, 2, 1)
    })
    identical(x, {
        set.seed(123)
        h(n, 2, rgamma(2, 1))
    })
    identical(x, {
        set.seed(123)
        h(n, lambda, smodel)
    })
})

## Test invalid arguments.
assertError(rcompound(-1, rpois(2), rgamma(2, 1)))


###
### Tests for rcomppois
###

## Test the function itself with various types of arguments.
n <- 20
lambda <- 2
smodel <- expression(rgamma(2, 1))
set.seed(123)
x <- numeric(n)
N <- rpois(n, 2)
y <- rgamma(sum(N), 2, 1)
x[which(N != 0)] <- tapply(y, rep(seq_len(n), N), sum)
stopifnot(exprs = {
    identical(x, {
        set.seed(123)
        rcomppois(n, 2, rgamma(2, 1))
    })
    identical(x, {
        set.seed(123)
        rcomppois(n, lambda, expression(rgamma(2, 1)))
    })
    identical(x, {
        set.seed(123)
        rcomppois(n, lambda, smodel)
    })
})

## Test the calling environment, that is that arguments are correctly
## identified when 'rcomppois' is called inside another function.
n <- 20
lambda <- 2
smodel <- expression(rgamma(2, 1))
set.seed(123)
x <- rcomppois(n, lambda, smodel)
f <- function(n, p, model)
{
    ## safe way to pass down all sorts of 'model' objects
    model <- substitute(model)
    if (is.name(model))
        model <- eval.parent(model)
    rcomppois(n, p, model)
}
g1 <- function(n, p, s, r)
    rcomppois(n, p, rgamma(s, r))
g2 <- function(n, p, s, r)
    rcomppois(n, p, expression(rgamma(s, r)))
h <- function(n, p, model)
{
    ## safe way to pass down all sorts of 'model' objects
    model <- substitute(model)
    if (is.name(model))
        model <- eval.parent(model)
    f(n, p, model)
}
stopifnot(exprs = {
    identical(x, {
        set.seed(123)
        f(n, 2, rgamma(2, 1))
    })
    identical(x, {
        set.seed(123)
        f(n, lambda, expression(rgamma(2, 1)))
    })
    identical(x, {
        set.seed(123)
        f(n, lambda, smodel)
    })
    identical(x, {
        set.seed(123)
        g1(n, 2, 2, 1)
    })
    identical(x, {
        set.seed(123)
        g2(n, 2, 2, 1)
    })
    identical(x, {
        set.seed(123)
        h(n, 2, rgamma(2, 1))
    })
    identical(x, {
        set.seed(123)
        h(n, lambda, smodel)
    })
})

## Test invalid arguments.
assertError(rcomppois(-1, lambda, smodel))
assertError(rcomppois(n, -1, smodel))
assertError(rcomppois(n, c(3, -1), smodel))

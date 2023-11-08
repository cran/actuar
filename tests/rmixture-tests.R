### actuar: Actuarial Functions and Heavy Tailed Distributions
###
### Tests for the simulation of discrete mixtures with 'rmixture'.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

## Load the package
library(actuar)

## Copy of tools::assertError.
assertError <- tools::assertError

## Set common values for the tests
n <- 20
bmodels <- expression(rexp(1/20),
                      rlnorm(3.6, 0.6),
                      rpareto(shape = 4, scale = 240))

## Function to inject the number of variates in an expression and
## evaluate it.
f <- function(n, expr)
{
    expr$n <- n
    eval(expr)
}

## Test a "normal" case (with data that is not reshuffled).
set.seed(123)
probs <- c(2, 3, 5)/10
nj <- rmultinom(1, n, prob = probs)
x <- c(f(nj[1], bmodels[[1]]), f(nj[2], bmodels[[2]]), f(nj[3], bmodels[[3]]))
set.seed(123)
stopifnot(exprs = {
    identical(x, rmixture(n, probs, bmodels, shuffle = FALSE))
})

## Test recycling of the probability vector.
set.seed(123)
probs <- 1
nj <- rmultinom(1, n, prob = rep_len(probs, 3))
x <- c(f(nj[1], bmodels[[1]]), f(nj[2], bmodels[[2]]), f(nj[3], bmodels[[3]]))
set.seed(123)
stopifnot(exprs = {
    identical(x, rmixture(n, probs, bmodels, shuffle = FALSE))
})

## Test recycling of the models vector.
set.seed(123)
probs <- c(2, 3, 5)
nj <- rmultinom(1, n, prob = probs)
x <- f(n, bmodels[[1]])
set.seed(123)
stopifnot(exprs = {
    identical(x, rmixture(n, probs, bmodels[1], shuffle = FALSE))
})

## Test special cases.
stopifnot(exprs = {
    identical(numeric(0), rmixture(0, probs, bmodels))
    identical(2L, length(rmixture(c(n, n), probs, bmodels)))
})

## Test the calling environment, that is that arguments are correctly
## identified when 'rmixture' is called inside another function.
set.seed(123)
probs <- c(2, 3, 5)/10
x <- rmixture(n, probs, bmodels)
f <- function(n, p, model)
    rmixture(n, p, model)
g <- function(n, p, m, q)
    rmixture(n, p, expression(rexp(m[1]), rlnorm(m[2], q[2]), rpareto(m[3], q[3])))
h <- function(n, p, model)
    f(n, c(p[1], p[2], p[3]),
             c(model[1], model[2], model[3]))
k <- function(n, p, m, q)
{
    ## Pathological case where the models expression does not evaluate
    ## in the frame of 'rmixture' as 'm' and 'q' will not be bound.
    ## The fix is to substitute variables by their values.
    models <- substitute(expression(rexp(m[1]), rlnorm(m[2], q[2]), rpareto(m[3], q[3])),
                         list(m = m, q = q))

    f(n, p, eval(models))
}
stopifnot(exprs = {
    identical(x, {
        set.seed(123)
        f(n, probs, bmodels)
    })
    identical(x, {
        set.seed(123)
        f(n, c(probs[1], probs[2], probs[3]),
          c(bmodels[1], bmodels[2], bmodels[3]))
    })
    identical(x, {
        set.seed(123)
        g(n, p = probs,
          m = c(eval(bmodels[[c(1, 2)]]), eval(bmodels[[c(2, 2)]]), eval(bmodels[[c(3, 2)]])),
          q = c(NA,                       eval(bmodels[[c(2, 3)]]), eval(bmodels[[c(3, 3)]])))
    })
    identical(x, {
        set.seed(123)
        h(n, probs,
          expression(rexp(eval(bmodels[[c(1, 2)]])),
                     rlnorm(eval(bmodels[[c(2, 2)]]), eval(bmodels[[c(2, 3)]])),
                     rpareto(shape = eval(bmodels[[c(3, 2)]]), scale = eval(bmodels[[c(3, 3)]]))))
    })
    identical(x, {
        set.seed(123)
        k(n, p = probs,
          m = c(eval(bmodels[[c(1, 2)]]), eval(bmodels[[c(2, 2)]]), eval(bmodels[[c(3, 2)]])),
          q = c(NA,                       eval(bmodels[[c(2, 3)]]), eval(bmodels[[c(3, 3)]])))
    })
})

## Finally, test invalid arguments.
assertError(rmixture(-1, probs, bmodels))
assertError(rmixture(c(3, -1), probs, bmodels))
assertError(rmixture(n, numeric(0), bmodels))
assertError(rmixture(n, 0, bmodels))
assertError(rmixture(n, c(0, 0), bmodels))
assertError(rmixture(n, probs, c(rexp(2), rexp(7))))

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
models <- expression(rexp(1/20),
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
x <- c(f(nj[1], models[[1]]), f(nj[2], models[[2]]), f(nj[3], models[[3]]))
set.seed(123)
stopifnot({
    identical(x, rmixture(n, probs, models, shuffle = FALSE))
})

## Test recycling of the probabilty vector.
set.seed(123)
probs <- 1
nj <- rmultinom(1, n, prob = rep_len(probs, 3))
x <- c(f(nj[1], models[[1]]), f(nj[2], models[[2]]), f(nj[3], models[[3]]))
set.seed(123)
stopifnot({
    identical(x, rmixture(n, probs, models, shuffle = FALSE))
})

## Test recycling of the models vector.
set.seed(123)
probs <- c(2, 3, 5)
nj <- rmultinom(1, n, prob = probs)
x <- f(n, models[[1]])
set.seed(123)
stopifnot({
    identical(x, rmixture(n, probs, models[1], shuffle = FALSE))
})

## Test special cases.
stopifnot({
    identical(numeric(0), rmixture(0, probs, models))
    identical(2L, length(rmixture(c(n, n), probs, models)))
})

## Finally, test invalid arguments.
assertError(rmixture(-1, probs, models))
assertError(rmixture(c(3, -1), probs, models))
assertError(rmixture(n, numeric(0), models))
assertError(rmixture(n, 0, models))
assertError(rmixture(n, c(0, 0), models))
assertError(rmixture(n, probs, c(rexp(2), rexp(7))))

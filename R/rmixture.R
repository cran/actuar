### actuar: Actuarial Functions and Heavy Tailed Distributions
###
### Simulation of discrete mixtures
###
###    f(x) = p_1 f_1(x) + ... + p_n f_n(x).
###
### Uses the syntax of rcomphierarc() for model specfification.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

rmixture <- function(n, probs, models, shuffle = TRUE)
{
    ## Validity checks (similar to other r<dist> functions and to
    ## rmultinom).
    if (any(is.na(n)) || any(n < 0))
        stop(sprintf("invalid first argument %s", sQuote("n")))
    if (all(probs <= 0))
        stop("no positive probabilities")
    if ((!is.expression(models)) || (length(models) == 0L))
        stop(sprintf("invalid third argument %s", sQuote("models")))

    ## Number of models in the mixture.
    m <- max(length(probs), length(models))

    ## Number of variates to generate: 'length(n)' if length of 'n' is
    ## > 1, like other 'r<dist>' functions.
    if (length(n) > 1L)
        n <- length(n)

    ## Number of variates from each model. By definition of the
    ## multinomial distribution, sum(nj) == n.
    ##
    ## Note that 'rmultinom' will normalize probabilities to sum 1.
    nj <- rmultinom(1, size = n, prob = rep_len(probs, m))

    ## Auxiliary function to generate 'n' variates from the model
    ## given in 'expr'. The expressions end up being evaluated three
    ## frames below the current one.
    f <- function(n, expr)
    {
        expr$n <- n
        eval.parent(expr, n = 3)
    }

    ## Simulate from each model the appropriate number of times and
    ## return result as an atomic vector. Variates are ordered by
    ## model: all random variates from model 1, then all random
    ## variates from model 2, and so on.
    x <- unlist(mapply(f, n = nj, expr = rep_len(models, m),
                       SIMPLIFY = FALSE))

    ## Return variates reshuffled or in the order above as per
    ## argument 'shuffle'.
    if (shuffle)
        x[sample.int(n)]
    else
        x
}

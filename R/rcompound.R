### actuar: Actuarial Functions and Heavy Tailed Distributions
###
### Simulation of standard, non hierarchical, compound models. Uses a
### simplified version of the syntax of 'rcomphierarc' for model
### specfification.
###
### Where 'rcomphierarc' was developed for flexibility, the functions
### therein aim at execution speed. Various algorithms were tested.
### No argument validity checks.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

rcompound <- function(n, model.freq, model.sev, SIMPLIFY = TRUE)
{
    ## Validity checks.
    if (any(is.na(n)) || any(n < 0))
        stop(sprintf("invalid first argument %s", sQuote("n")))

    ## Convert model expressions into language objects.
    cl.freq <- substitute(model.freq)
    cl.sev <- substitute(model.sev)

    ## If a model expression was actually an object containing the
    ## model, we need to evaluate the object to retrieve the model.
    ## If the resulting object is an expression object, its first
    ## element is the language object we are after.
    if (is.name(cl.freq))
    {
        cl.freq <- eval.parent(cl.freq)
        if (is.expression(cl.freq))
            cl.freq <- cl.freq[[1L]]
    }
    if (is.name(cl.sev))
    {
        cl.sev <- eval.parent(cl.sev)
        if (is.expression(cl.sev))
            cl.sev <- cl.sev[[1L]]
    }

    ## If a model expression is wrapped into 'expression' (as in
    ## 'rcomphierarc'), get rid of the call.
    if (cl.freq[[1L]] == "expression")
        cl.freq <- cl.freq[[-1L]]
    if (cl.sev[[1L]] == "expression")
        cl.sev <- cl.sev[[-1L]]

    ## Initialize the output vector. We will use the fact that 'res'
    ## is filled with zeros later.
    res <- numeric(n)

    ## Add the number of variates to the 'model.freq' call.
    cl.freq$n <- n

    ## Generate frequencies.
    N <- eval.parent(cl.freq)

    ## Add the number of variates to the 'model.sev' call.
    cl.sev$n <- sum(N)

    ## Generate all severities.
    x <- eval.parent(cl.sev)

    ## Create a vector that will be used as a factor to regroup
    ## severities for the computation of aggregate values. Idea:
    ## assign one integer to each frequency and repeat that integer a
    ## number of times equal to the frequency. For example, if the
    ## frequencies are (2, 0, 1, 3), then the vector will be (1, 1, 3,
    ## 4, 4, 4).
    f <- rep.int(seq_len(n), N)

    ## Compute aggregate values and put them in the appropriate
    ## positions in the output vector. The positions corresponding to
    ## zero frequencies are already initialized with zeros.
    res[which(N != 0)] <- tapply(x, f, sum)

    if (SIMPLIFY)
        res
    else
        list(aggregate = res,
             frequency = N,
             severity = x)
}

rcomppois <- function(n, lambda, model.sev, SIMPLIFY = TRUE)
{
    ## Validity checks.
    if (any(is.na(n)) || any(n < 0))
        stop(sprintf("invalid first argument %s", sQuote("n")))
    if (any(lambda < 0))
        stop(sprintf("invalid values in %s", sQuote("lambda")))

    ## Convert model expression into language object.
    cl.sev <- substitute(model.sev)

    ## If the model expression was actually an object containing the
    ## model, we need to evaluate the object to retrieve the model.
    ## If the resulting object is an expression object, its first
    ## element is the language object we are after.
    if (is.name(cl.sev))
    {
        cl.sev <- eval.parent(cl.sev)
        if (is.expression(cl.sev))
            cl.sev <- cl.sev[[1L]]
    }

    ## Get rid of the eventual 'expression' call in the language
    ## object.
    if (cl.sev[[1L]] == "expression")
        cl.sev <- cl.sev[[-1L]]

    ## Initialize the output vector.
    res <- numeric(n)

    ## Generate frequencies from Poisson distribution.
    N <- rpois(n, lambda)

    ## Add the number of variates to the 'model.sev' call.
    cl.sev$n <- sum(N)

    ## Generate all severities.
    x <- eval.parent(cl.sev)

    ## Create a vector that will be used as a factor to regroup
    ## severities for the computation of aggregate values. (See
    ## comments in 'rcompound' for details.)
    f <- rep.int(seq_len(n), N)

    ## Compute aggregate values and put them in the appropriate
    ## positions in the output vector.
    res[which(N != 0)] <- tapply(x, f, sum)

    if (SIMPLIFY)
        res
    else
        list(aggregate = res,
             frequency = N,
             severity = x)
}

### actuar: Actuarial Functions and Heavy Tailed Distributions
###
### Empirical moments for individual and grouped data.
###
### See Klugman, Panjer & Willmot, Loss Models, Wiley, 1998.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

emm <- function(x, order = 1, ...) UseMethod("emm")

emm.default <- function(x, order = 1, ...)
{
    if (any(order < 0))
      stop(sprintf("%s must be positive", sQuote("order")))

    colMeans(outer(x, order, "^"), ...)
}

emm.grouped.data <- function(x, order = 1, ...)
{
    ## Function does not work for negative moments
    if (any(order < 0))
      stop(sprintf("%s must be positive", sQuote("order")))

    ## Extract group boundaries
    cj <- eval(expression(cj), envir = environment(x))

    ## Compute the factor
    ##
    ## f_j = (c_j^{k + 1} - c_{j-1}^{k+1})/((k+1) * (c_j - c_{j-1}))
    ##
    ## for all values of 'j' and 'k' == 'order'.
    y <- diff(outer(cj, order + 1, "^")) / outer(diff(cj), order + 1)

    ## Drop the group boundaries column
    x <- as.matrix(x[-1L])

    ## Compute sum(n_j * f_j)/sum(nj) for all values of 'order'.
    drop(crossprod(x, y)) / colSums(x, ...)
}

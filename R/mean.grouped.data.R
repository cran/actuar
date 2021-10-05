### actuar: Actuarial Functions and Heavy Tailed Distributions
###
### Mean of grouped data objects
###
### See Klugman, Panjer & Willmot, Loss Models, Wiley, 1998.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

## New method of base::mean generic for grouped data
mean.grouped.data <- function(x, ...)
{
    ## Get group boundaries
    cj <- eval(expression(cj), envir = environment(x))

    ## Compute group midpoints
    midpoints <- cj[-length(cj)] + diff(cj)/2

    ## Extract frequencies columns by dropping the boundaries column;
    ## convert to matrix for use in crossprod()
    x <- as.matrix(x[-1L])

    ## Compute mean per column
    drop(crossprod(x, midpoints))/colSums(x)
}

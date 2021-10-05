### actuar: Actuarial Functions and Heavy Tailed Distributions
###
### Variance and standard deviation
###
### See Klugman, Panjer & Willmot, Loss Models, Wiley, 1998.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>
### Walter Garcia-Fontes

## New generics for functions of the stats package
var <- function(x, ...) UseMethod("var")
sd <- function(x, ...) UseMethod("sd")

## Default methods are stats::var and stats:sd
var.default <- function(x, y = NULL, na.rm = FALSE, use, ...)
    stats::var(x, y = NULL, na.rm = FALSE, use)
sd.default <- function(x, na.rm = FALSE, ...)
    stats::sd(x, na.rm = FALSE)

## Methods for grouped data
var.grouped.data <- function(x, ...)
{
    ## Get group boundaries
    cj <- eval(expression(cj), envir = environment(x))

    ## Compute group midpoints
    midpoints <- cj[-length(cj)] + diff(cj)/2

    ## Compute midpoints minus mean and square it
    midsquare <- (midpoints - mean(x))^2

    ## Extract frequencies columns by dropping the boundaries column;
    ## convert to matrix for use in crossprod()
    x <- as.matrix(x[-1L])

    ## Compute mean per column
    drop(crossprod(x, midsquare))/(colSums(x) - 1)
}

sd.grouped.data <- function(x, ...)
{
    ## Square root of variance
    drop(sqrt(var.grouped.data(x)))
}

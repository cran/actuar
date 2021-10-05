### actuar: Actuarial Functions and Heavy Tailed Distributions
###
### Quantiles (inverse of the ogive) for grouped data
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>
### Walter Garcia-Fontes

quantile.grouped.data <- function(x, probs = seq(0, 1, 0.25),
                                  names = TRUE, ...)
{
    ## We keep the first frequencies column only; group boundaries are
    ## in the environment of 'x'
    y <- x[, 2L]
    x <- eval(expression(cj), envir = environment(x))

    ## Inverse of the ogive
    fun <- approxfun(c(0, cumsum(y))/sum(y), x,
                     yleft = min(x), yright = max(x),
                     method = "linear", ties = "ordered")

    ## Quantiles
    res <- fun(probs)

    if (names)
    {
        dig <- max(2, getOption("digits"))
        names(res) <- formatC(paste(100 * probs, "%", sep = ""),
                              format = "fg", width = 1, digits = dig)
    }
    res
}

summary.grouped.data <- function(object, ...)
{
    ## Keep only the first frequencies column
    object <- object[1L:2L]

    res <- quantile(object)
    res <- c(res[1L:3L], mean(object), res[4L:5L])
    names(res) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
    class(res) <- c("summaryDefault", "table")
    res
}

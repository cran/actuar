### actuar: Actuarial Functions and Heavy Tailed Distributions
###
### Histogram for grouped data
###
### See Klugman, Panjer & Willmot, Loss Models, Wiley.
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>, Mathieu Pigeon

hist.grouped.data <-
    function(x, freq = NULL, probability = !freq,
             density = NULL, angle = 45, col = NULL, border = NULL,
             main = paste("Histogram of", xname), xlim = range(x),
             ylim = NULL, xlab = xname, ylab, axes = TRUE,
             plot = TRUE, labels = FALSE, ...)
{
    ## We keep the first frequencies column only; group boundaries are
    ## in the environment of 'x'
    y <- x[, 2L]
    x <- eval(expression(cj), envir = environment(x))

    ## If any frequency is non finite, omit the group
    keep <- which(is.finite(y))
    y <- y[keep]
    x <- x[c(1L, keep + 1L)]

    ## Some useful values
    n <- sum(y)                         # total number of observations
    h <- diff(x)                        # group widths
    dens <- y/(n * h)                   # group "densities"

    ## Cannot plot histogram with infinite group
    if (any(is.infinite(x)))
        stop("infinite group boundaries")

    ## The rest is taken from hist.default()
    xname <- paste(deparse(substitute(x), 500), collapse = "\n")
    equidist <- diff(range(h)) < 1e-07 * mean(h)
    if (is.null(freq))
    {
        freq <- if (!missing(probability))
                    !as.logical(probability)
                else equidist
    }
    else if (!missing(probability) && any(probability == freq))
        stop(sprintf("%s is an alias for %s, however they differ.",
                     sQuote("probability"), sQuote("!freq")))
    mids <- 0.5 * (x[-1L] + x[-length(x)])
    r <- structure(list(breaks = x, counts = y, intensities = dens,
                        density = dens, mids = mids, xname = xname,
                        equidist = equidist),
                   class = "histogram")
    if (plot)
    {
        plot(r, freq = freq, col = col, border = border, angle = angle,
             density = density, main = main, xlim = xlim, ylim = ylim,
             xlab = xlab, ylab = ylab, axes = axes, labels = labels, ...)
        invisible(r)
    }
    else
        r
}

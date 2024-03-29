### actuar: Actuarial Functions and Heavy Tailed Distributions
###
### Extraction and replacement methods for grouped data
### objects
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>,
### Mathieu Pigeon, Louis-Philippe Pouliot

"[.grouped.data" <- function(x, i, j)
{
    ## Only columns to extract are specified.
    if (nargs() < 3L)
    {
        if (missing(i))
            return(x)
        if (is.matrix(i))
            return(as.matrix(x)[i])
        res <- as.data.frame(NextMethod())
        if (length(i) > 1 && 1 %in% seq(ncol(x))[i])
        {
            environment(res) <- environment(x)
            class(res) <- c("grouped.data", class(res))
        }
        return(res)
    }

    ## Convert row and column indexes to strictly positive integers.
    ii <- if (missing(i)) seq.int(nrow(x)) else seq.int(nrow(x))[i]
    ij <- if (missing(j)) integer(0) else seq.int(ncol(x))[j]

    ## Extraction of at least the group boundaries (the complicated case).
    if (!length(ij) || 1L %in% ij)
    {
        ## Extraction of group boundaries in increasing order only
        ## (untractable otherwise).
        if (is.unsorted(ii))
        {
            warning("rows extracted in increasing order")
            ii <- sort(ii)
        }

        ## Fetch the appropriate group boundaries.
        cj <- eval(expression(cj), envir = environment(x))
        cj <- cj[sort(unique(c(ii, ii + 1L)))]

        ## Extraction of the first column only: return the vector of group
        ## boundaries.
        if (identical(ij, 1L))
            return(cj)

        ## Return a modified 'grouped.data' object.
        res <- NextMethod()
        environment(res) <- new.env()
        assign("cj", cj, environment(res))
        return(res)
    }

    ## All other cases handled like a regular data frame.
    NextMethod()
}

"[<-.grouped.data" <- function(x, i, j, value)
{
    nA <- nargs()
    if (nA == 4L)
    {
        ii <- if (missing(i)) NULL else i
        ij <- if (missing(j)) NULL else j
    }
    else if (nA == 3L)
    {
        ## No arguments inside [ ]: only replacing by NULL is supported.
        if (missing(i) && missing(j))
        {
            if (is.null(value))
                return(x[logical(0)])
            stop("impossible to replace boundaries and frequencies simultaneously")
        }
        ## Indexing by a logical matrix is supported, but only two
        ## types of replacement are allowed: replacing in the
        ## first column only, or replacing in any column but the
        ## first.
        if (is.logical(i) && is.matrix(i) && all(dim(i) == dim(x)))
        {
            ij <- apply(i, 2, any)      # columns with replacements
            if (match(TRUE, ij) == 1)   # boundaries to replace
            {
                if (length(ij) > 1)     # boundaries and frequencies
                    stop("impossible to replace boundaries and frequencies simultaneously")
                ii <- i[, ij]           # boundaries only
            }
            return(NextMethod())        # frequencies only
        }
        ## Indexing by a non logical matrix is not supported.
        if (is.matrix(i))
            stop("only logical matrix subscripts are allowed in replacement")
        ## Indexing by a vector: the argument specifies columns to
        ## replace.
        ij <- i
        ii <- NULL
    }
    else
        stop("need 0, 1, or 2 subscripts")

    ## Convert row and column indexes to integers.
    ii <- if (is.null(ii)) seq.int(nrow(x)) else seq.int(nrow(x))[ii]
    ij <- if (is.null(ij)) integer(0) else seq.int(ncol(x))[ij]

    ## Replacement at least in the group boundaries column.
    if (!length(ij) || 1L %in% ij)
    {
        ## supported: replacement of group boundaries only
        if (identical(ij, 1L))
        {
            cj <- eval(expression(cj), envir = environment(x))
            cj[sort(unique(c(ii, ii + 1L)))] <- value
            res <- grouped.data(cj, x[, -1L])
            names(res) <- names(x)
            return(res)
        }
        ## not supported (untractable): replacement in the column of
        ## boundaries and any other column
        stop("impossible to replace boundaries and frequencies simultaneously")
    }

    ## All other cases handled like a regular data frame.
    NextMethod()
}

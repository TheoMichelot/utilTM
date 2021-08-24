
#' Block-wise autocorrelation function
#'
#' This original version of this function was written by Stacy DeRuiter.
#' It computes the ACF of a vector passed as input, ignoring time lags that
#' overlap blocks (i.e., time series).
#'
#' @param x Numeric vector
#' @param blocks Vector of indices of blocks, of same length as x
#' @param lag_max Maximum lag at which to calculate the ACF. Defaults to the
#' length of the shortest block.
#' @param make_plot Logical: should the ACF be plotted? Defaults to TRUE.
#'
#' @return Vector of ACF at each time lag
#'
#' @export
acf_by_block <- function(x, blocks, lag_max = NULL, make_plot = TRUE, ...) {
    # Check arguments
    blocks <- factor(blocks, levels = unique(blocks))
    if(length(blocks) != length(x)) {
        stop("x and blocks should be the same length.")
    }
    if(is.null(lag_max)) {
        # By default, the max lag is the length of the shortest block
        lag_max <- min(tapply(blocks, blocks, length))
    }
    n <- length(x)

    # Get indices of last element of each block excluding the last
    i1 <- which(blocks[-1] != blocks[-n])
    x_with_NAs <- x
    block_acf <- matrix(1, nrow = lag_max + 1, ncol = 1)
    for(lag in 1:lag_max) {
        for(i_block in 1:length(i1)){
            x_with_NAs <- append(x_with_NAs, NA, after = i1[i_block])
        }
        # Adjust for the growing x_with_NAs
        i1 <- i1 + head(c(0:(nlevels(blocks) - 1)), -1)
        this_acf <- acf(x_with_NAs, lag.max = lag_max,
                        plot = FALSE, na.action = na.pass)
        block_acf[lag + 1] <- this_acf$acf[lag + 1, 1, 1]
    }

    if(make_plot) {
        # This is just to get an acf object...
        A <- acf(x, lag.max = lag_max, plot = FALSE, na.action = na.pass)
        # Into which we'll insert our coefficients
        A$acf[,1,1] <- block_acf
        # And then plot OUR results
        plot(A, ...)
    }

    return(block_acf)
}

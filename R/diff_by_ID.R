
#' First-order difference for grouped data
#'
#' @param x Vector of values
#' @param ID Vector of group IDs, of same length as x
#' @param ... Optional parameters passed to difftime, if x is a vector
#' of dates
#'
#' @return Vector of same length as x of first-order differences of x
#' (and NAs for last value of each group ID)
#'
#' @export
diff_by_ID <- function(x, ID, ...) {
    if(length(x) != length(ID)) {
        stop("'x' and 'ID' must be vectors of the same length")
    }
    n <- length(ID)

    # Indices of first and last value of each group
    i0 <- which(ID[-1] != ID[-n])
    i_first <- c(1, i0 + 1)
    i_last <- c(i0, n)

    # First-order differences
    dx <- rep(NA, n)
    if(!inherits(x, "POSIXct")) {
        dx[-i_last] <- x[-i_first] - x[-i_last]
    } else {
        dx[-i_last] <- difftime(time1 = x[-i_first], time2 = x[-i_last], ...)
    }


    return(dx)
}

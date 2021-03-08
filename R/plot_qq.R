
#' Quantile-quantile plot
#'
#' @param x Sample to compare to the standard normal distribution
#'
#' @importFrom ggplot2 geom_abline stat_qq
#' @export
plot_qq <- function(x) {
    # qqnorm gives the coordinates of the points, used for axis limits
    qq_pts <- qqnorm(x, plot = FALSE)
    lim <- range(qq_pts$x, qq_pts$y, na.rm = TRUE)

    df <- data.frame(x = x)
    p <- ggplot(df, aes(sample = x)) +
        stat_qq() +
        coord_equal(xlim = lim, ylim = lim) +
        geom_abline(slope = 1, intercept = 0) +
        theme_light()

    return(p)
}

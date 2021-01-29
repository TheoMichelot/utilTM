
#' Raster to dataframe
#'
#' @param rast Raster object
#'
#' @return Data frame with columns x and y (coordinates), and z (raster value)
#'
#' @importFrom raster coordinates values
#'
#' @export
rast2df <- function(rast) {
    data.frame(coordinates(rast), z = values(rast))
}

#' Plot raster
#'
#' @param rast Raster object
#' @param contour Logical. If TRUE, contour lines are added. (Default: FALSE)
#' @param guide Argument passed to scale_fill. Defaults to "colourbar" for
#' continuous colour scale, and can be set to "none" to remove the scale.
#'
#' @return A ggplot object
#'
#' @importFrom ggplot2 ggplot coord_equal theme_light geom_raster geom_contour
#' scale_x_continuous scale_y_continuous aes
#' @importFrom scico scale_fill_scico
#'
#' @export
plot_rast <- function(rast, contour = FALSE, guide = "colourbar") {
    # Create data frame for ggplot
    rast_df <- rast2df(rast)

    # Create plot
    p <- ggplot(rast_df, aes(x, y)) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        coord_equal() +
        theme_light() +
        geom_raster(aes(fill = z)) +
        scale_fill_scico(guide = guide, palette = "bamako", direction = -1)

    if(contour) {
        p <- p + geom_contour(mapping = aes(z = z), colour = "black", size = 0.3)
    }

    return(p)
}

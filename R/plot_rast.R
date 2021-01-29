
#' Raster to dataframe
#'
#' @param rast Raster object
#'
#' @return Data frame with columns x and y (coordinates), and z (raster value)
#'
#' @importFrom raster coordinates values
rast2df <- function(rast) {
    data.frame(coordinates(rast), z = values(rast))
}

#' Plot raster
#'
#' @param rast Raster object
#'
#' @return A ggplot object
#'
#' @importFrom ggplot2 ggplot coord_equal theme_light geom_raster geom_contour
#' scale_x_continuous scale_y_continuous
#' @importFrom scico scale_fill_scico
plot_rast <- function(rast) {
    # Create data frame for ggplot
    rast_df <- rast2df(rast)

    # Create plot
    p <- ggplot(rast_df, aes(x, y)) +
        coord_equal() +
        theme_light() +
        geom_raster(aes(fill = z)) +
        geom_contour(mapping = aes(z = z), colour = "black", size = 0.3) +
        scale_fill_scico(name = "", palette = "bamako", direction = -1) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0))

    return(p)
}

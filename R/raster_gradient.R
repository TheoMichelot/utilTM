
#' Gradient of raster using bilinear interpolation
#'
#' @details The gradient has a closed form solution if bilinear interpolation
#' is used, as described in Appendix C of Michelot (2019, PhD thesis).
#' https://etheses.whiterose.ac.uk/23688/
#'
#' @param rast Raster object
#' @param pts Matrix with two columns (one for each coordinate)
#'
#' @return Matrix with same dimensions as pts, where the first column
#' gives the gradient along the first dimension, and the second column
#' along the second dimension.
#'
#' @export
grad_raster <- function(rast, pts) {
    # Loop over points
    grads <- apply(pts, 1, function(pt) {
        x <- pt[1]
        y <- pt[2]

        # Get coordinates of four adjacent points
        x1 <- max(coordinates(rast)[which(coordinates(rast)[,1] < x), 1])
        x2 <- min(coordinates(rast)[which(coordinates(rast)[,1] > x), 1])
        y1 <- max(coordinates(rast)[which(coordinates(rast)[,2] < y), 2])
        y2 <- min(coordinates(rast)[which(coordinates(rast)[,2] > y), 2])

        # Evaluate raster at four adjacent points
        f11 <- rast[cellFromXY(rast, c(x1, y1))]
        f12 <- rast[cellFromXY(rast, c(x1, y2))]
        f21 <- rast[cellFromXY(rast, c(x2, y1))]
        f22 <- rast[cellFromXY(rast, c(x2, y2))]

        # Gradient
        dfdx <- ((y2 - y) * (f21 - f11) + (y - y1) * (f22 - f12)) /
            ((y2 - y1) * (x2 - x1))
        dfdy <- ((x2 - x) * (f12 - f11) + (x - x1) * (f22 - f21)) /
            ((y2 - y1) * (x2 - x1))
        return(c(dfdx, dfdy))
    })
    data.frame(gradx = grads[1,], grady = grads[2,])
}

#' Transform SpatRast to convenient format for grad_terra
#'
#' @param rast SpatRast object
#'
#' @return List with elements x (grid of x values), y (grid of y values),
#' and z (matrix of raster values), as expected by grad_terra()
#'
#' @export
rast_to_xyz <- function(rast) {
    lim <- as.vector(ext(rast))
    res <- res(rast)
    xgrid <- seq(lim[1] + res[1]/2, lim[2] - res[1]/2, by = res[1])
    ygrid <- seq(lim[3] + res[2]/2, lim[4] - res[2]/2, by = res[2])
    # Matrix needs to be reversed vertically so that increasing
    # row index corresponds to increasing y coordinate
    zmat <- t(apply(as.matrix(rast, wide = TRUE), 2, rev))
    return(list(x = xgrid, y = ygrid, z = zmat))
}

#' Gradient of terra raster using bilinear interpolation
#'
#' @details The gradient has a closed form solution if bilinear interpolation
#' is used, as described in Appendix C of Michelot (2019, PhD thesis).
#' https://etheses.whiterose.ac.uk/23688/
#'
#' @param rast SpatRast object
#' @param pts Matrix with two columns (one for each coordinate)
#' @param xgrid Grid of x values (if rast not provided)
#' @param ygrid Grid of y values (if rast not provided)
#' @param zmat Matrix of raster values (if rast not provided)
#'
#' @return Matrix with same dimensions as pts, where the first column
#' gives the gradient along the first dimension, and the second column
#' along the second dimension.
#'
#' @export
grad_terra <- function(rast = NULL, pts, xgrid = NULL, ygrid = NULL, zmat = NULL) {
    if(!is.null(rast)) {
        xyz <- rast_to_xyz(rast = rast)
        xgrid <- xyz$x
        ygrid <- xyz$y
        zmat <- xyz$z
    }
    if(!is.matrix(pts)) {
        # "pts" could be either a vector with two elements (i.e.,
        # a single point) or a data frame
        pts <- matrix(as.matrix(pts), ncol = 2)
    }
    x <- pts[,1]
    y <- pts[,2]
    x_ind <- findInterval(x = x, vec = xgrid)
    y_ind <- findInterval(x = y, vec = ygrid)
    x1 <- xgrid[x_ind]
    x2 <- xgrid[x_ind + 1]
    y1 <- ygrid[y_ind]
    y2 <- ygrid[y_ind + 1]
    f11 <- zmat[cbind(x_ind, y_ind)]
    f12 <- zmat[cbind(x_ind, y_ind + 1)]
    f21 <- zmat[cbind(x_ind + 1, y_ind)]
    f22 <- zmat[cbind(x_ind + 1, y_ind + 1)]
    dfdx <- ((y2 - y) * (f21 - f11) + (y - y1) * (f22 - f12)) /
        ((y2 - y1) * (x2 - x1))
    dfdy <- ((x2 - x) * (f12 - f11) + (x - x1) * (f22 - f21)) /
        ((y2 - y1) * (x2 - x1))
    return(cbind(dfdx, dfdy))
}

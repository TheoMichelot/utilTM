
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



#' Gradient of terra raster using bilinear interpolation
#'
#' @details The gradient has a closed form solution if bilinear interpolation
#' is used, as described in Appendix C of Michelot (2019, PhD thesis).
#' https://etheses.whiterose.ac.uk/23688/
#'
#' @param rast SpatRast object
#' @param pts Matrix with two columns (one for each coordinate)
#'
#' @return Matrix with same dimensions as pts, where the first column
#' gives the gradient along the first dimension, and the second column
#' along the second dimension.
#'
#' @export
grad_terra <- function(rast, pts) {
    if(is.null(dim(pts))) {
        pts <- matrix(pts, nrow = 1)
    }
    # Loop over points
    grads <- apply(pts, 1, function(pt) {
        x <- pt[1]
        y <- pt[2]

        # Get coordinates of four adjacent points
        x1 <- max(crds(rast)[which(crds(rast)[,1] < x), 1])
        x2 <- min(crds(rast)[which(crds(rast)[,1] > x), 1])
        y1 <- max(crds(rast)[which(crds(rast)[,2] < y), 2])
        y2 <- min(crds(rast)[which(crds(rast)[,2] > y), 2])

        # Evaluate raster at four adjacent points
        f11 <- as.numeric(rast[cellFromXY(rast, matrix(c(x1, y1), nrow = 1))])
        f12 <- as.numeric(rast[cellFromXY(rast, matrix(c(x1, y2), nrow = 1))])
        f21 <- as.numeric(rast[cellFromXY(rast, matrix(c(x2, y1), nrow = 1))])
        f22 <- as.numeric(rast[cellFromXY(rast, matrix(c(x2, y2), nrow = 1))])

        # Gradient
        dfdx <- ((y2 - y) * (f21 - f11) + (y - y1) * (f22 - f12)) /
            ((y2 - y1) * (x2 - x1))
        dfdy <- ((x2 - x) * (f12 - f11) + (x - x1) * (f22 - f21)) /
            ((y2 - y1) * (x2 - x1))
        return(c(dfdx, dfdy))
    })
    cbind(gradx = grads[1,], grady = grads[2,])
}

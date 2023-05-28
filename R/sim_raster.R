
#' Simulate random covariate field
#' 
#' @param rho Coefficient of spatial autocorrelation
#' @param lim Vector of limits of the map (xmin, xmax, ymin, ymax)
#' @param res Grid resolution (defaults to 1)
#' 
#' @return Raster object
#' 
#' @importFrom terra focal crop rast values values<-
#' @importFrom stats runif
#' @export
sim_raster <- function(rho, lim, res = 1) {
    # rho must be odd
    if(rho %% 2 == 0) {
        rho <- rho + 1        
    }
    if(rho <= 1) {
        # Not sure why, but terra::focal doesn't like rho = 1
        stop("'rho' should be at least 3")
    }
    
    # Generate uniform random weights on larger map 
    # (it will be cropped by moving average)
    xgrid <- seq(lim[1] - rho/2, lim[2] + rho/2, by = res)
    ygrid <- seq(lim[3] - rho/2, lim[4] + rho/2, by = res)
    xygrid <- expand.grid(xgrid, ygrid)
    xyz <- data.frame(x = xygrid[,1], y = xygrid[,2], z = runif(nrow(xygrid)))
    rand <- rast(xyz, type = "xyz")
    
    # Define the matrix of weights for smoothing 
    # (ones inside disc, zeros outside)
    foox <- (-(rho - 1)/2) : ((rho - 1)/2)
    fooy <- (-(rho - 1)/2) : ((rho - 1)/2)
    fooxy <- expand.grid(foox, fooy)
    dist <- sqrt(fooxy[,1]^2 + fooxy[,2]^2)
    food <- rep(0, nrow(fooxy))
    food[which(dist <= rho/2)] <- 1
    m <- matrix(food, rho, rho)
    
    # Smooth random field using 2-d moving average over disc
    smooth <- focal(x = rand, w = m)
    smooth <- crop(smooth, lim)
    
    # Scale to [0,1]
    values(smooth) <- (values(smooth) - min(values(smooth)))/
        (max(values(smooth)) - min(values(smooth)))
    
    return(smooth)
}


#' Transform longitude-latitude to UTM
#'
#' @param lon Vector of longitudes
#' @param lat Vector of latitudes
#' @param zone UTM zone for projection
#'
#' @return Data frame with columns x (Easting) and y (Northing)
#'
#' @export
LL_to_UTM <- function(lon, lat, zone) {
    # @importFrom sp coordinates proj4string CRS
    # @importFrom rgdal spTransform

    xy <- data.frame(x = lon, y = lat)
    coordinates(xy) <- c("x", "y")
    proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")
    res <- spTransform(xy, CRS(paste("+proj=utm +zone=", zone,
                                     " ellps=WGS84", sep='')))
    return(as.data.frame(res))
}

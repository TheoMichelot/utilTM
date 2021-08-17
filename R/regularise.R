
#' Interpolate time series to regular time intervals
#' 
#' This is a wrapper for the function stats::spline, which is applied
#' to each variable to interpolate.
#' 
#' @param data Data frame with columns for ID (time series ID),
#' time, and variables to interpolate
#' @param new_interval Time interval of interpolation, in the syntax
#' expected by 'seq'. For example, '5 sec', '1 min', '2 hour', etc.
#' @param var_names Vector of names of variables to interpolate
#' 
#'  @return Data frame with columns for ID, time, and interpolated
#'  variables 
regularise <- function(data, new_interval, var_names) {
    # If there is no ID column, assume that there is only one time series
    if(is.null(data$ID)) {
        data$ID <- "ts1"
    }
    IDs <- unique(data$ID)
    new_data <- NULL
    
    # Loop over time series
    for(id in IDs) {
        sub_data <- subset(data, ID == id)
        
        # Regular time sequence for interpolation
        new_times <- seq(sub_data$time[1], 
                         sub_data$time[nrow(sub_data)], 
                         by = new_interval)
        
        # Interpolate to regular time intervals
        new_sub_data <- as.data.frame(sapply(var_names, function(var) {
            spline(x = sub_data$time, y = sub_data[,var], 
                   xout = new_times)$y
        }))
        new_sub_data$ID <- id
        new_sub_data$time <- new_times
        
        new_data <- rbind(new_data, new_sub_data)
    }
    
    return(new_data)
}


#' Split time series data at gaps
#'
#' @param data Data frame with columns 'ID' and 'time', and any other
#' data variables
#' @param max_interval Longest permissible time interval within a
#' time series, in minutes The time series will be split where
#' there are gaps longer than \code{max_interval}.
#' @param min_ts_length Minimum permissible length for a time series,
#' in minutes. Time series shorter than this will be excluded from
#' the output data frame. Defaults to 0, i.e., all time series are
#' kept.
#' @param units Time unit (default: "min")
#' @param tz Time zone
#'
#' @return Data frame with same columns as input, plus a column 'ID_split'
#' for the ID of the time series segments after splitting.
#'
#' @export
split_ts <- function(data, max_interval, min_ts_length = 0,
                     units = "min", tz = NULL) {
    dt <- diff_by_ID(x = data$time, ID = data$ID, units = units, tz = tz)

    # Indices where ID changes
    ind_change <- which(data$ID[-1] != data$ID[-nrow(data)])
    ind_start <- c(1, ind_change + 1)
    ind_end <- c(ind_change, nrow(data))

    # Indices where ID should change: either new track or long interval
    ind_end_new <- sort(c(ind_end, which(dt >= max_interval)))
    ind_start_new <- c(1, ind_end_new[-length(ind_end_new)] + 1)

    # Define ID_split to change each time there is a long gap
    data$ID_split <- NA
    for(i in 1:length(ind_end_new)) {
        ind_sub_track <- ind_start_new[i]:ind_end_new[i]
        data$ID_split[ind_sub_track] <- paste0(data$ID[ind_sub_track], "-", i)
    }

    # Only keep sub-tracks longer than some duration
    ts_lengths <- sapply(unique(data$ID_split), function(id) {
        ind <- which(data$ID_split == id)
        difftime(data$time[ind[length(ind)]], data$time[ind[1]],
                 units = units, tz = tz)
    })
    ID_keep <- names(ts_lengths)[which(ts_lengths >= min_ts_length)]
    data <- subset(data, ID_split %in% ID_keep)

    return(data)
}

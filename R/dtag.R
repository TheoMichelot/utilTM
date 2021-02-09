
#' Prepare DTag data
#'
#' @param file File name
#' @param start Time for start of data in 'file' (as POSIX)
#' @param freq Required frequency of output
#' @param keep_interval Vector of two times: start and end of time interval
#' over which to keep the data (e.g., to exclude periods where the tag was not
#' on the animal). If not specified, the whole data set is kept.
#'
#' @return Data frame with columns:
#' \itemize{
#' \item{}{}
#' \item{}{}
#' \item{}{}
#' }
#'
#' @importFrom R.matlab readMat
#' @importFrom tagtools find_dives
#' @importFrom lubridate ymd_hms %within%
#' @export
prep_dtag <- function(file, start_time, out_freq, keep_interval = NULL) {
    # Read file
    dtag <- readMat(file)
    # Get tag ID from file name
    tag_name <- tools::file_path_sans_ext(basename(file))

    # Get data variables
    n <- length(dtag$p)
    pitch <- dtag$pitch
    roll <- dtag$roll
    dhead <- atan2(sin(dtag$head[-n] - dtag$head[-1]),
                   cos(dtag$head[-n] - dtag$head[-1]))
    head <- cumsum(c(0, dhead))
    depth <- dtag$p

    # Sequence of times
    dt <- 1/as.numeric(dtag$fs)
    start_time <- ymd_hms(start_time)
    time <- seq(start_time, by = dt, length = n)
    time_num <- as.numeric(time - time[1])

    # Identify dives using tagtools
    dives <- find_dives(p = matrix(depth, ncol = 1), mindepth = 10, sampling_rate = 1/dt)
    diveID <- rep(NA, length = length(pitch))
    dive_type <- rep(NA, length = length(pitch))
    for(d in 1:nrow(dives)) {
        if(dives$max[d] > 100) {
            ind <- which(time_num >= dives$start[d] & time_num <= dives$end[d])
            diveID[ind] <- paste0(tag_name, "-", d)
            dive_type[ind] <- ifelse(dives$max[d] > 700, "deep", "shallow")
        }
    }

    # Create data frame for this DTag
    data_all <- data.frame(tagID = tag_name,
                           ID = diveID,
                           Ax = dtag$Aw[,1],
                           Ay = dtag$Aw[,2],
                           Az = dtag$Aw[,3],
                           pitch = pitch,
                           roll = roll,
                           head = head,
                           depth = depth,
                           time = time,
                           time_num = time_num,
                           type = dive_type)

    # Thin to required frequency
    thin <- as.numeric(dtag$fs/out_freq)
    data <- data_all[seq(1, nrow(data_all), by = thin),]

    # Remove surface time and non-complete dives
    data <- subset(data, !is.na(ID))

    # Only keep times in keep_times interval
    if(!is.null(keep_interval)) {
        keep_times <- interval(ymd_hms(keep_interval[1]), ymd_hms(keep_interval[2]))
        data <- subset(data, time %within% keep_times)
    }

    # Add dive proportion column
    data$diveprop <- NA
    for(id in unique(data$ID)) {
        ind <- which(data$ID == id)
        data$diveprop[ind] <- seq(0, 1, length = length(ind))
    }

    return(data)
}

#' Prepare data from several DTags
#'
#' @param files Vector of file names
#' @param start_times Vector of simes for start of data in 'file' (as POSIX)
#' @param out_freq Required frequency of output
#' @param keep_intervals List where each element is a vector of two times: start
#' and end of time interval over which to keep the data (e.g., to exclude periods
#' where the tag was not on the animal). If not specified, the whole data set is
#' kept.
#'
#' @return Same data frame as prep_dtag
#'
#' @importFrom pbmcapply pbmclapply
#' @export
prep_dtags <- function(files, start_times, out_freq, keep_intervals = NULL,
                       n_cores = 1) {
    # Call prep_dtag on each data file (in parallel if possible)
    data_list <- pbmclapply(seq_along(files), function(i) {
        prep_dtag(file = files[i],
                  start_time = start_times[i],
                  out_freq = out_freq,
                  keep_interval = keep_intervals[[i]])
    }, mc.cores = n_cores)

    return(do.call("rbind", data_list))
}

#' Add exposure columns to DTag data frame
#'
#' @param data Data frame, at least with columns 'tagID', 'ID' (dive ID),
#' and 'time'
#' @param expo Data frame with one row for each exposure event, and one column
#' each for 'tagID', 'start' and 'end'
#'
#' @return Input data frame with additional columns for
#' \itemize{
#' \item{expo}{either 'before', 'during' or 'after' (an exposure event)}
#' \item{expo_dive}{'yes' for exposure dives, 'no' for other dives}
#' \item{expo_time}{time since start of exposure (0 before)}
#' }
#'
#' @importFrom lubridate interval %within%
#' @export
add_expo <- function(data, expo) {
    # Create output data frame
    out <- data

    # Make sure the data formats are okay
    expo$tagID <- as.character(expo$tagID)
    expo$start <- ymd_hms(expo$start)
    expo$end <- ymd_hms(expo$end)
    data$tagID <- as.character(data$tagID)
    data$ID <- as.character(data$ID)
    data$time <- ymd_hms(data$time)

    # Add 'expo' column ('before', 'during' or 'after')
    out$expo <- "before"
    for(i in seq_along(unique(data$tagID))) {
        # This tag ID (as string rather than factor)
        thisID <- as.character(unique(data$tagID)[i])
        ind_thisID <- which(data$tagID == thisID)

        # Exposure information for this tag ID
        expo_thisID <- subset(expo, tagID == thisID)

        # Loop over exposure events for this tag ID
        for(k in which(expo$tagID == thisID)) {
            # Start and end of exposure
            start <- expo$start[k]
            end <- expo$end[k]
            # Indices of exposure
            ind_expo <- ind_thisID[which(out$time[ind_thisID] %within% interval(start, end))]
            out$expo[ind_expo] <- "during"
            # Indices after exposure
            ind_after <- ind_thisID[which(out$time[ind_thisID] > end)]
            out$expo[ind_after] <- "after"
        }
    }

    # Add 'expo_dive' column ('yes' if exposed dive, 'no' otherwise)
    out$expo_dive <- "no"
    for(i in seq_along(unique(data$ID))) {
        # This dive ID
        ID <- unique(data$ID)[i]
        ind_thisID <- which(data$ID == ID)

        # If there is any exposed time during the dive, mark as exposed
        if(any(out$expo[ind_thisID] == "during")) {
            out$expo_dive[ind_thisID] <- "yes"
        }
    }

    # Add 'expo_time' column (time since start of exposure)
    out$expo_time <- 0
    for(i in seq_along(unique(data$tagID))) {
        # This tag ID
        ID <- unique(data$tagID)[i]
        ind_thisID <- which(data$tagID == ID)

        ind_expo <- ind_thisID[which(out$expo[ind_thisID] != "before")]
        if(length(ind_expo) > 0) {
            t0 <- out$time[ind_expo[1]]
            out$expo_time[ind_expo] <- out$time[ind_expo] - t0
        }
    }

    return(out)
}

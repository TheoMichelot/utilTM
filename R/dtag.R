
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
#' \item{tagID}{Identifier for tag (i.e., time series)}
#' \item{Ax, Ay, Az}{Acceleration in x, y, and z, obtained directly from the tag data}
#' \item{pitch, roll, head}{Eulerian angles}
#' \item{depth}{Depth in metres}
#' \item{time}{Time as POSIX}
#' \item{time_num}{Time as numeric}
#' \item{type}{Dive type (either "shallow" or "deep")}
#' }
#'
#' @importFrom R.matlab readMat
#' @importFrom tagtools find_dives njerk
#' @importFrom lubridate ymd_hms %within%
#' @export
prep_dtag <- function(file, start_time = NULL, out_freq, downsample = "thin_regular",
                      keep_interval = NULL) {
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
    njerk <- njerk(A = dtag$Aw, sampling_rate = as.numeric(dtag$fs))

    # Sequence of times
    dt <- 1/as.numeric(dtag$fs)
    if(!is.null(start_time)) {
        start_time <- ymd_hms(start_time)
        time <- seq(start_time, by = dt, length = n)
        time_num <- as.numeric(time - time[1])
    } else {
        time <- rep(NA, n)
        time_num <- seq(0, by = dt, length = n)
    }

    # Identify dives using tagtools
    dives <- find_dives(p = matrix(depth, ncol = 1), mindepth = 100, surface = 2,
                        sampling_rate = 1/dt, findall = 1)
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
                           dhead = c(dhead, NA),
                           depth = depth,
                           njerk = njerk,
                           time = time,
                           time_num = time_num,
                           type = dive_type)

    # Remove surface time and non-complete dives
    data_all <- subset(data_all, !is.na(ID))

    # Thin to required frequency
    thin <- round(as.numeric(dtag$fs/out_freq))
    if(downsample == "thin_regular") {
        # Keep observations on regular time grid
        wh_keep <- seq(1, nrow(data_all), by = thin)
        data <- data_all[wh_keep,]
    } else if(downsample == "thin_random") {
        # Keep observations on irregular time grid (at random)
        wh_keep <- sort(sample(1:nrow(data_all), size = nrow(data_all)/thin))
        data <- data_all[wh_keep,]
    } else if(downsample == "running_mean") {
        # Running mean of observations
        data <- NULL
        ids <- na.omit(unique(data_all$ID))
        for(i in seq_along(ids)) {
            cat("\nDive:", ids[i], "\n")
            sub_data_all <- subset(data_all, ID == ids[i])
            keep <- seq(1, by = thin, length = floor(nrow(sub_data_all)/thin))
            sub_data <- sub_data_all[keep,]

            for(var in c("pitch", "roll", "head")) {
                sub_data[[var]] <- gtools::running(sub_data_all[[var]], width = thin,
                                                   by = thin, fun = mean, trim = 0,
                                                   na.rm = TRUE)
            }

            data <- rbind(data, sub_data)
        }
    } else {
        stop(paste0("'downsample' must be one of: 'thin_regular',",
                    "'thin_random', or 'running_mean'"))
    }

    # Only keep times in keep_times interval
    if(!is.null(start_time) & !is.null(keep_interval)) {
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
prep_dtags <- function(files, start_times = NULL, out_freq, downsample = "thin_regular",
                       keep_intervals = NULL, n_cores = 1) {

    data_list <- NULL
    # Call prep_dtag on each data file (in parallel if possible)
    for(i in seq_along(files)){
        cat("\n\nFile:", files[i], "\n")
        # data_list <- pbmclapply(seq_along(files), function(i) {
        data_list[[i]] <- prep_dtag(file = files[i],
                                    start_time = start_times[i],
                                    out_freq = out_freq,
                                    downsample = downsample,
                                    keep_interval = keep_intervals[[i]])
        # }, mc.cores = n_cores)
    }


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
    # and 'after_expo' ('no' or name of exposure event)
    out$expo <- "before"
    out$after_expo <- "no"
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
            out$after_expo[c(ind_expo, ind_after)] <- paste0("expo_", thisID, "_", k)
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

#' Add received level columns to DTag data frame
#'
#' @param data Data frame, at least with columns 'tagID' and 'time'
#' @param RL Data frame of received level, with with columns 'tagID',
#' 'time', and 'received_level'
#'
#' @return Input data frame with additional column for received level
#'
#' @export
add_RL <- function(data, RL) {
    # Create output data frame
    out <- data
    out$received_level <- 0

    # Make sure the data formats are okay
    RL$tagID <- as.character(RL$tagID)
    RL$time <- ymd_hms(RL$time)
    data$tagID <- as.character(data$tagID)
    data$ID <- as.character(data$ID)
    data$time <- ymd_hms(data$time)

    # TODO: fix case where several exposures per tagID
    # Loop over IDs
    for(id in unique(RL$tagID)) {
        subRL <- subset(RL, tagID == id)

        # Loop over lines of data and add RL value with matching time
        k <- 1
        i0 <- which(data$tagID == id & data$time >= subRL$time[1])[1]
        i1 <- which(data$tagID == id & data$time >= subRL$time[nrow(subRL)])[1]
        for(i in i0:i1) {
            out$received_level[i] <- subRL$received_level[k]
            while(k < nrow(subRL) & subRL$time[k] < data$time[i + 1]) {
                k <- k + 1
            }
        }
    }

    return(out)
}

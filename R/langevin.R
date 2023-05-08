
#' Gradient of log(pi(x)) in Langevin model
#'
#' @param pt 2d point
#' @param beta Vector of selection coefficients
#' @param cov_list List of rasters
#'
#' @return 2d gradient
#' @export
grad_log_pi <- function(pt, beta, cov_list) {
    pts <- matrix(pt, nrow = 1)
    grad_list <- lapply(cov_list, function(cov)
        grad_raster(rast = cov, pts = pts))
    sum_grad <- rowSums(sapply(seq_along(grad_list), function(i)
        unlist(beta[i] * grad_list[[i]])))
    return(sum_grad)
}

#' Simulate from Langevin process
#'
#' @param times Vector of (numeric) times
#' @param beta Vector of selection coefficients
#' @param cov_list List of rasters
#' @param gamma Speed parameter
#' @param x0 Initial point
#'
#' @return Matrix of simulated locations
#' @export
sim_langevin <- function(times, beta, cov_list, gamma = 1, x0 = c(0, 0)) {
    n <- length(times)
    dt <- diff(times)
    x <- matrix(0, nrow = n, ncol = 2)
    x[1,] <- x0
    for(i in 2:n) {
        if(i%%(n/20) == 0) cat("\rIteration", i, "/", n, "\n")
        grad <- grad_log_pi(x[i-1,], beta = beta, cov_list = cov_list)
        x[i,] <- x[i-1,] + 0.5 * gamma^2 * dt[i-1] * grad +
            rnorm(2, 0, gamma * sqrt(dt[i-1]))
    }
    return(x)
}

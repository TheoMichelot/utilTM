
#' Simulate Ornstein-Uhlenbeck process
#'
#' @param times Times of observation
#' @param mu Centre of attraction
#' @param beta Reversion parameter
#' @param sigma Diffusion parameter
#' @param z0 Initial value for simulation (defaults to mu)
#'
#' @return Data frame with columns z (simulated process) and time
#'
#' @export
sim_OU <- function(times, mu = 0, beta = 1, sigma = 1, z0 = NULL) {
    # Number of observations
    n <- length(times)
    # Time intervals
    dt <- diff(times)

    # Check input
    if(length(mu) == 1) {
        mu <- rep(mu, n)
    } else if(length(mu) != n) {
        stop("'mu' should be of length 1 or", n)
    }
    if(length(beta) == 1) {
        beta <- rep(beta, n)
    } else if(length(beta) != n) {
        stop("'beta' should be of length 1 or", n)
    }
    if(length(sigma) == 1) {
        sigma <- rep(sigma, n)
    } else if(length(sigma) != n) {
        stop("'sigma' should be of length 1 or", n)
    }
    if(is.null(z0)) {
        z0 <- mu[1]
    }

    # Initialise process
    z <- rep(z0, n)
    # Iterate through observations
    for(i in 1:(n-1)) {
        mean <- exp(-beta[i] * dt[i]) * z[i] +
            (1 - exp(-beta[i] * dt[i])) * mu[i]
        sd <- sigma[i]/sqrt(2*beta[i]) * sqrt(1 - exp(-2 * beta[i] * dt[i]))
        z[i+1] <- rnorm(1, mean = mean, sd = sd)
    }

    return(data.frame(time = times, z = z))
}

#' Simulate Brownian motion
#'
#' @param times Times of observation
#' @param mu Drift
#' @param sigma Diffusion
#' @param z0 Initial value for simulation (defaults to 0)
#'
#' @return Data frame with columns z (simulated process) and time
#'
#' @export
sim_BM <- function(times, mu = 0, sigma = 1, z0 = 0) {
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
    if(length(sigma) == 1) {
        sigma <- rep(sigma, n)
    } else if(length(sigma) != n) {
        stop("'sigma' should be of length 1 or", n)
    }

    # Generate increments and derive process
    dz <- rnorm(n - 1, mean = mu[-n] * dt, sd = sigma[-n] * sqrt(dt))
    z <- cumsum(c(z0, dz))

    return(data.frame(time = times, z = z))
}

#' Simulate Ornstein-Uhlenbeck process
#'
#' @param times Times of observation
#' @param mu Centre of attraction
#' @param beta Reversion parameter
#' @param sigma Diffusion parameter
#' @param tau Time scale of autocorrelation parameter
#' @param kappa Long-term variance parameter
#' @param z0 Initial value for simulation (defaults to mu)
#'
#' @return Data frame with columns z (simulated process) and time
#'
#' @export
sim_OU <- function(times, mu = 0, beta = 1, sigma = 1,
                   tau = NULL, kappa = NULL, z0 = mu[1]) {
    # Number of observations
    n <- length(times)
    # Time intervals
    dt <- diff(times)

    # Check input
    if(!is.null(tau) & !is.null(kappa)) {
        beta <- 1/tau
        sigma <- sqrt(2 * kappa / tau)
    } else if(!is.null(tau) | !is.null(kappa)) {
        stop("You need to specify both 'tau' and 'kappa'")
    }
    if(length(mu) == 1) {
        mu <- rep(mu, n)
    } else if(length(mu) != n) {
        stop("'mu' should be of length 1 or", n)
    }
    if(length(beta) == 1) {
        beta <- rep(beta, n)
    } else if(length(beta) != n) {
        stop("'beta' (or 'tau') should be of length 1 or", n)
    }
    if(length(sigma) == 1) {
        sigma <- rep(sigma, n)
    } else if(length(sigma) != n) {
        stop("'sigma' (or 'kappa') should be of length 1 or", n)
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

#' Make covariance matrix for CTCRW simulation
#'
#' @param beta Parameter beta of the OU process
#' @param sigma Parameter sigma of the OU process
#' @param dt Time interval
make_cov <- function(beta, sigma, dt) {
    Q <- matrix(0, 2, 2)
    Q[1,1] <- sigma^2/(2*beta) * (1-exp(-2*beta*dt))
    Q[2,2] <- (sigma/beta)^2 * (dt + (1-exp(-2*beta*dt))/(2*beta) -
                                    2*(1-exp(-beta*dt))/beta)
    Q[1,2] <- sigma^2/(2*beta^2) * (1 - 2*exp(-beta*dt) + exp(-2*beta*dt))
    Q[2,1] <- Q[1,2]
    return(Q)
}

#' Simulate from CTCRW process
#'
#' @param times Vector of times of observations
#' @param mu Mean velocity parameter
#' @param beta Reversion parameter of velocity
#' @param sigma Variance parameter of velocity
#' @param tau Time scale of autocorrelation parameter
#' @param nu Mean speed parameter
#' @param z0 Initial value for simulation (default: 0)
#'
#' @return Simulated track
#'
#' @export
sim_CTCRW <- function(times, mu = 0, beta = 1, sigma = 1,
                      tau = NULL, nu = NULL, z0 = 0) {
    # Number of observations
    n <- length(times)
    # Time intervals
    dt <- diff(times)

    # Check input
    if(!is.null(tau) & !is.null(nu)) {
        beta <- 1/tau
        sigma <- 2 * nu / sqrt(tau * pi)
    } else if(!is.null(tau) | !is.null(nu)) {
        stop("You need to specify both 'tau' and 'nu'")
    }
    if(length(mu) == 1) {
        mu <- rep(mu, n)
    } else if(length(mu) != n) {
        stop("'mu' should be of length 1 or", n)
    }
    if(length(beta) == 1) {
        beta <- rep(beta, n)
    } else if(length(beta) != n) {
        stop("'beta' (or 'tau') should be of length 1 or", n)
    }
    if(length(sigma) == 1) {
        sigma <- rep(sigma, n)
    } else if(length(sigma) != n) {
        stop("'sigma' (or 'nu') should be of length 1 or", n)
    }

    data <- matrix(0, nrow = n, ncol = 2)
    mean <- rep(NA, 2)
    colnames(data) <- c("v", "z")
    data[1,2] <- z0
    for(i in 2:n) {
        # Mean of next state vector (V, Z)
        p <- exp(-beta[i-1] * dt[i-1])
        mean[1] <- p * data[i-1, "v"] + (1 - p) * mu[i-1]
        mean[2] <- data[i-1, "z"] + mu[i-1] * dt[i-1] +
            (data[i-1, "v"] - mu[i-1]) / beta[i-1] * (1 - p)

        # Covariance of next state vector
        V <- make_cov(beta = beta[i-1], sigma = sigma[i-1], dt = dt[i-1])

        data[i,] <- mgcv::rmvn(1, mu = mean, V = V)
    }

    return(data.frame(z = data[,"z"], v = data[,"v"], time = times))
}

#' Simulate Cox-Ingersoll-Ross model
#'
#' @param times Times of observation
#' @param mu Mean (> 0)
#' @param beta Reversion (> 0)
#' @param sigma Diffusion (> 0)
#' @param z0 Initial value for simulation (> 0)
#'
#' @return Data frame with columns z (simulated process) and time
#'
#' @export
sim_CIR <- function(times, mu = 1, beta = 1, sigma = 1, z0 = 1) {
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

    # Generate increments and derive process
    z <- rep(z0, n)
    for(i in 2:n) {
        c <- 2 * beta[i-1] / ((1 - exp(-beta[i-1] * dt[i-1])) * sigma[i-1]^2)
        df <- 4 * beta[i-1] * mu[i-1] / sigma[i-1]^2
        ncp <- 2 * c * z[i-1] * exp(- beta * dt[i-1])
        Y <- rchisq(n = 1, df = df, ncp = ncp)
        z[i] <- Y/(2 * c)
    }

    return(data.frame(time = times, z = z))
}

#' Simulate from geometric Brownian motion
#'
#' @param times Times of observation
#' @param mu Percentage drift (> 0)
#' @param sigma Percentage volatility (> 0)
#' @param z0 Initial value for simulation (> 0)
#'
#' @return Data frame with columns z (simulated process) and time
#'
#' @export
sim_GBM <- function(times, mu = 1, sigma = 1, z0 = 1) {
    n <- length(times)
    z <- rep(NA, n)
    z[1] <- z0
    dt <- diff(times)

    # Check input
    if(length(mu) == 1) {
        mu <- rep(mu, n)
    } else if(length(mu) != n) {
        stop("'mu' should be of length 1 or", n)
    }
    if(length(sigma) == 1) {
        sigma <- rep(sigma, n)
    } else if(length(sigma) != n) {
        stop("'sigma' should be of length 1 or", n)
    }

    # Loop over time steps
    for(i in 1:(n-1)) {
        w <- rnorm(1, 0, sigma)
        z[i+1] <- z[i] * exp((mu[i] - sigma[i]^2/2) * dt[i] + w)
    }

    return(data.frame(time = times, z = z))
}

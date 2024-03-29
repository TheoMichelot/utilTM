% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sde.R
\name{sim_OU}
\alias{sim_OU}
\title{Simulate Ornstein-Uhlenbeck process}
\usage{
sim_OU(
  times,
  mu = 0,
  beta = 1,
  sigma = 1,
  tau = NULL,
  kappa = NULL,
  z0 = mu[1]
)
}
\arguments{
\item{times}{Times of observation}

\item{mu}{Centre of attraction}

\item{beta}{Reversion parameter}

\item{sigma}{Diffusion parameter}

\item{tau}{Time scale of autocorrelation parameter}

\item{kappa}{Long-term variance parameter}

\item{z0}{Initial value for simulation (defaults to mu)}
}
\value{
Data frame with columns z (simulated process) and time
}
\description{
Simulate Ornstein-Uhlenbeck process
}

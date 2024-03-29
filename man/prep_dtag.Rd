% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dtag.R
\name{prep_dtag}
\alias{prep_dtag}
\title{Prepare DTag data}
\usage{
prep_dtag(
  file,
  start_time = NULL,
  out_freq,
  downsample = "thin_regular",
  keep_interval = NULL
)
}
\arguments{
\item{file}{File name}

\item{keep_interval}{Vector of two times: start and end of time interval
over which to keep the data (e.g., to exclude periods where the tag was not
on the animal). If not specified, the whole data set is kept.}

\item{start}{Time for start of data in 'file' (as POSIX)}

\item{freq}{Required frequency of output}
}
\value{
Data frame with columns:
\itemize{
\item{tagID}{Identifier for tag (i.e., time series)}
\item{Ax, Ay, Az}{Acceleration in x, y, and z, obtained directly from the tag data}
\item{pitch, roll, head}{Eulerian angles}
\item{depth}{Depth in metres}
\item{time}{Time as POSIX}
\item{time_num}{Time as numeric}
\item{type}{Dive type (either "shallow" or "deep")}
}
}
\description{
Prepare DTag data
}

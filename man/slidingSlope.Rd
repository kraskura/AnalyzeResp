% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/slidingSlope.R
\name{slidingSlope}
\alias{slidingSlope}
\title{steepest slope analysis only}
\usage{
slidingSlope(data, Ch, local_path, N_Ch)
}
\arguments{
\item{data}{Dataframe with each measurement present (O2 observed each second )}

\item{Ch}{Channel that needs to be analyzed (use format: c(1, 2), c(3), etc. up to 4 channels)}

\item{local_path}{Logical. If TRUE (no default) all returned files will be saved in the local working directory. Can also provide a path if this function is run independently}

\item{N_Ch}{Number of channels for the oxygen meter}
}
\value{
The output from \code{\link{print}}
}
\description{
Estimates slopes at every time iteration (often second), and pulls out the steepest 30 s, 60 s, 90 s, and 180 s slopes, as well as the regression slope for the entire recorded duration.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bounds_functions_public.R
\name{plot.ste_cdf}
\alias{plot.ste_cdf}
\title{Plot CDFs of Spectral of Treatment Effects}
\usage{
\method{plot}{ste_cdf}(ste_cdf_obj)
}
\arguments{
\item{ste_cdf_obj}{Object of class \code{ste_cdf} computed by the function \code{ste_cdf()}}
}
\value{
A \code{ggplot2} plot of the STE.
}
\description{
This function plots CDFs of spectral treatment effects computed by \code{ste_cdf()}.
}
\examples{
# Generate Data
temp <- matrix(rnorm(100), 10,10)
temp <- (temp + t(temp))/2
Treat.Mat <- temp*1
diag(Treat.Mat) <- 0

temp <- matrix(rnorm(400), 20,20)
temp <- (temp + t(temp))/2
Control.Mat <- temp*1
diag(Control.Mat) <- 0
# Compute STE
effects <- ste(Treat.Mat, Control.Mat)
# Compute CDF of STE and plot them
effects_cdf <- ste_cdf(effects)
plot(effects_cdf$stt)
plot(effects_cdf$stc)
}

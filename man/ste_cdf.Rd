% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bounds_functions_public.R
\name{ste_cdf}
\alias{ste_cdf}
\title{Construct CDFs of Spectral Treatment Effects}
\usage{
ste_cdf(
  ste_obj,
  smooth = FALSE,
  bandwidth = NULL,
  start.ext = 1,
  end.ext = 1,
  n.grid = c(NULL, NULL)
)
}
\arguments{
\item{ste_obj}{Object of clas 'ste' computed by the function \code{ste()}.}

\item{smooth}{Logical indicating whether the smoothing should be performed. Smoothing is recommended if potential outcomes are estimated.}

\item{bandwidth}{2x1 Array indicating bandwidth to be used for smoothing STT and STC respectively.}

\item{start.ext}{Number of standard deviations by which to extend the start point of the CDF}

\item{end.ext}{Number of standard deviations by which to extend the end point of the CDF}

\item{n.grid}{2x1 Array indicating the number of grids used to compute the CDF for treatment and control respectively.}
}
\value{
A list \code{ste} containing \tabular{ll}{
\code{stt_cdf} \tab An \code{ste_cdf} object containing the CDF of the STT \cr
\tab \cr
\code{stc_cdf} \tab An \code{ste_cdf} object containing the CDF of the STC \cr
}
}
\description{
This function computes CDFs of spectral treatment effects computed by \code{ste()}. The output of this function can be plotted directly using \code{plot()}.
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

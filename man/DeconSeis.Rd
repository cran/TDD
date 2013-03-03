\name{DeconSeis}
\alias{DeconSeis}
\title{Deconvolve discrete instrument response from many traces}
\description{
  Deconvolves instrument responses from a seismogram structure from the
  RSEIS package.
}
\usage{
DeconSeis(GH, inst, L, fl = 0.1, fh = NaN, bitweight = NULL, dec =
rep(1, length(GH$JSTR)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GH}{Seismogram structure}
  \item{inst}{Vector of indices of instrument responses within L to
    deconvolve from each trace in GH}
  \item{L}{List in which each element is a discrete instrument response
    (from MakeDPZ, for example)}
  \item{fl}{Low corner of filter (NaN for no high-pass filtering) (Hz)}
  \item{fh}{High corner of filter (NaN for no low-pass filtering) (Hz)}
  \item{bitweight}{Vector of optional counts-to-volts factors for data in counts
    (volts/counts)--NULL if data are already in volts}
  \item{dec}{Oversampling/decimation factor (optional); vector equal to
    the number of traces considered}
}
\details{
  Discrete instrument responses are specific to a given sampling rate.
  If the response you give has a different sample rate (given by DPZ$dt)
  from the trace x, you will get incorrect results.  DPZ$dt * dec should
  be equal to the sample intervals of the traces.
}
\value{
  GH, with the instrument response removed from every trace.
}
\author{Jake Anderson}

\examples{
data(COLOC)
swig(COLOC)
L = GetDPZ(c(4, 14), c(0.01, 0.01)) # get responses for 3T and 40T-1s
inst = c(1,1,1,2,2,2) # deconvolve 3T response from channel 1-3,
                      # 40T-1 response from channel 4-6
D = DeconSeis(COLOC, inst, L)
swig(D)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }


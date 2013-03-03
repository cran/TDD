\name{DeconTrace}
\alias{DeconTrace}
\title{Deconvolve Instrument Response (Single Trace)}
\description{
  Deconvolves a discrete instrument response from a seismic trace.  
}
\usage{
DeconTrace(x, DPZ, fl = 0.05, fh = NaN, bitweight = NULL, dec = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Trace from which instrument response is deconvolved}
  \item{DPZ}{Discrete instrument response list (from MakeDPZ, for example)}
  \item{fl}{Low corner of filter (NaN for no high-pass filtering) (Hz)}
  \item{fh}{High corner of filter (NaN for no low-pass filtering) (Hz)}
  \item{bitweight}{Optional counts-to-volts factor for data in counts
    (volts/counts)--NULL if data are already in volts}
  \item{dec}{Oversampling/decimation factor (optional)}
}
\details{
  Discrete instrument responses are specific to a given sampling rate.
  If the response you give has a different sample rate (given by DPZ$dt)
  from the trace x, you will get incorrect results.  DPZ$dt * dec should
  be equal to the trace's sample interval.
}
\value{
  Deconvolved velocity trace (vector).
}

\author{Jake Anderson}

\seealso{
  ConvolveTrace, DeconSeis
}
\examples{
# Response of Guralp CMG-3T
DPZ = GetDPZ(4, 0.01)[[1]]

data(COLOC)
x = COLOC$JSTR[[1]]

DeconTrace(x, DPZ)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }


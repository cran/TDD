\name{ConvolveTrace}
\alias{ConvolveTrace}
\title{Convolve Trace with Instrument Response}
\description{Convolves a single velocity trace (m/s) with a discrete instrument
  response to get the voltage signal it returns.}
\usage{
ConvolveTrace(x, DPZ, dec = 1)
}

\arguments{
  \item{x}{Velocity trace (m/s)}
  \item{DPZ}{Discrete instrument response (from MakeDPZ, for example)}
  \item{dec}{Oversampling/decimation factor (optional)}
}
\details{
  Discrete instrument responses are specific to a given sampling rate.
  If the response you give has a different sample rate (given by DPZ$dt)
  from the trace x, you will get incorrect results.  DPZ$dt times dec
  should be equal to the sample interval of the trace.
}
\value{
  Convolved trace in volts (vector).
}
\author{Jake Anderson}

\seealso{
  DeconTrace, DeconSeis
}
\examples{
# Response of Guralp CMG-40T
DPZ = GetDPZ(12, 1)[[1]]

x = rnorm(1000)
ConvolveTrace(x, DPZ)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }


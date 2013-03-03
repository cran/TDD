\name{PZ2Resp}
\alias{PZ2Resp}
\title{Calculate Instrument Response from Poles/Zeros}
\description{Returns complex instrument response for a vector of
  frequencies and set of analog poles and zeros.  Optionally plots the
  magnitude of the instrument response.}
\usage{
PZ2Resp(PZ, f, PLOT = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PZ}{Analog poles and zeros}
  \item{f}{Vector of frequencies for which response is calculated}
  \item{PLOT}{Logical: whether to plot magnitude of response}
}
\details{
  The response is calculated by the following equation:
  R(s) = PZ$Sense * PZ$Knorm * prod(s - PZ$zeros)/prod(s - PZ$poles)

  where s = 2 * pi * 1i * f.
  
  PZ requires the following elements:
  poles: Vector of poles
  np: number of poles
  zeros: Vector of zeros
  nz: Number of zeros
  Knorm: Normalization constant
  Sense: Instrument sensitivity (V/(m/s))
}
\value{Vector of instrument response values corresponding to the
  frequencies in f.
}

\author{Jake Anderson}

\examples{
# Response of Guralp CMG-40T

PZ = GetPZ(12)[[1]]

f = (1:10000 - 1)/1000

PZ2Resp(PZ, f)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }


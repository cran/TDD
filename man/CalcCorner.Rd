\name{CalcCorners}
\alias{CalcCorners}
\title{Calculate Corner Frequencies}
\description{Inputs a continuous instrument response list and returns a
  vector of its cutoff frequencies (defined here as the -3dB point).
  Optionally plots results.}
\usage{
CalcCorners(PZ, f = 1:1000000/1000, PLOT = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PZ}{Continuous instrument response list (from GetPZ, for
    example)}
  \item{f}{Vector of frequency values to test (Hz)}
  \item{PLOT}{Logical: plot response spectrum and mark corner
    frequencies?}
}
\value{
  Vector of corner frequencies in response spectrum.
}

\author{Jake Anderson}

\examples{
# Response of CMG-40T
PZ = GetPZ(12)[[1]]
CalcCorners(PZ, PLOT = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }


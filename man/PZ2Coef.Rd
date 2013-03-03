\name{PZ2Coef}
\alias{PZ2Coef}
\title{Calculate Recursive Filter Coefficients}
\description{Returns coefficients of recursive filter approximating
  instrument response, given poles/zeros and sample interval.  
}
\usage{
PZ2Coef(PZ, dt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PZ}{Poles/zeros list}
  \item{dt}{Sample interval}
}
\details{
  PZ requires the following elements:
  poles: Vector of poles
  np: number of poles
  zeros: Vector of zeros
  nz: Number of zeros
  Knorm: Normalization constant
  Sense: Instrument sensitivity (V/(m/s))

  Output recursive filter is of the form a[1] * y_i + a[2] * y_(i-1) +
  a[3] * y_(i-2) + ... = b[1] * x_i + b[2] * x_(i-1)
  + b[3] * x_(i-2) + ..., where x is ground motion velocity and y is the
  recorded voltage.

}
\value{List including the following elements:
  \item{b}{Coefficients of filter input terms}
  \item{a}{Coefficients of filter output terms}
}
\author{Jake Anderson}
\note{
  Due to effects of discretization, the spectrum of the recursive filter
  DOES NOT match that of the poles/zeros.  So, poles and zeroes must be
  adjusted in order to make them match, either by inversion or by the
  bilinear transform.
}


\examples{
PZ_40T = list(poles = c(-0.149 + 0.149i, -0.149 - 0.149i, -503, -1010,
-1130), zeros = c(0, 0), Knorm = 574095649, Sense = 800)

dt = 0.01
PZ2Coef(PZ_40T, dt)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }


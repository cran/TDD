\name{ButterPZ}
\alias{ButterPZ}
\title{Calculate Butterworth Filter Poles and Zeros}
\description{Calculates continuous poles and zeros (of the Laplace
  transform) of Butterworth filters.}
\usage{
ButterPZ(fl = NaN, nl = NaN, fh = NaN, nh = NaN, g = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fl}{Low corner frequency (Hz)}
  \item{nl}{Order of high-pass filter}
  \item{fh}{High corner frequency (Hz)}
  \item{nh}{Order of low-pass filter}
  \item{g}{gain (unitless)}
}
\details{
  For a bandpass filter, all inputs should be non-NaN.  For high pass,
  fh and nh should be NaN; for low pass, fl and nl should be NaN.  Input
  corner frequencies are in cycles/second, not radians/second.
}
\value{List including the following elements:
  \item{poles}{Vector of poles (rad/s)}
  \item{zeros}{Vector of zeros (rad/s)}
  \item{np}{Number of poles}
  \item{nz}{Number of zeros}
  \item{Knorm}{Normalization constant}
  \item{Sense}{Sensitivity (V * s/m)}
}
\author{Jake Anderson}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
# Calculate poles and zeros of Butterworth filter with a second-order
# high-pass above 1 Hz, and a fourth-order low-pass below 10 Hz
ButterPZ(1, 2, 10, 4) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }


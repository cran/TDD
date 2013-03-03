\name{Oversample}
\alias{Oversample}
\title{Oversample by Nearest-Neighbor Interpolation}
\description{In order to maintain digital filter fidelity at very high
  frequencies, it is sometimes necessary to oversample a signal.  This
  function oversamples a signal by nearest-neighbor interpolation.}
\usage{
Oversample(x, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Signal to be oversampled}
  \item{n}{Factor by which it should be oversampled}
}
\details{
  The output probably needs to be decimated after deconvolution.
}
\value{Vector of oversampled data.
}
\author{Jake Anderson}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
# Oversample a random trace by a factor of 10
x = rnorm(100)
Oversample(x, 10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }


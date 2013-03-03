\name{PlotResp}
\alias{PlotResp}
\title{Plot Instrument Responses}
\description{Plots responses of analog and digital poles/zeros lists.}
\usage{
PlotResp(PZ, DPZ, fmin = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PZ}{Poles and zeros of continuous instrument response}
  \item{DPZ}{Digital poles and zeros of discrete instrument response}
  \item{fmin}{Minimum frequency to plot}
}
\details{
  PZ and DPZ must contain the elements poles, zeros, np, nz, Knorm, and
  Sense.  Additionally, DPZ must contain the element dt.
}
\value{
  Graphical side effects only.
}
\author{Jake Anderson}

\seealso{
PZ2Resp
}
\examples{
# Response of Guralp CMG-40T

PZ = GetPZ(12)[[1]]
DPZ = GetDPZ(12, 1)[[1]]
PlotResp(PZ, DPZ, fmin = 1/50)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }


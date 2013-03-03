\name{GetPZ}
\alias{GetPZ}
\title{Retrieve Included Continuous Instrument Response}
\description{Continuous responses for common seismometers are included
  in this package.  GetPZ retrieves them.}
\usage{
GetPZ(w)
}
\arguments{
  \item{w}{Vector of indices of seismometers used (see Details)}
}
\details{Seismometers are numbered as follows:

  Broadband Seismometers:
  1.   Streckeisen STS-1 (360 s)
  2.   Trillium 240 (generation 1)
  3.   Trillium 240 (generation 2)
  4.   Guralp CMG-3T
  5.   Streckeisen STS-2 (generation 1)
  6.   Streckeisen STS-2 (generation 2)
  7.   Streckeisen STS-2 (generation 3)
  8.   Trillium 120
  9.   Compact Trillium

  Intermediate Seismometers:
  10.  Trillium 40
  11.  Guralp CMG-3ESP
  12.  Guralp CMG-40T (30 s)
  13.  Streckeisen STS-1 (20 s)

  Short-period Seismometers:
  14:  Guralp CMG-40T (1 s)
}
\value{
  List of instrument responses corresponding to instruments given in w.
}
\references{
  Sources for all instrument responses are given in the comments of
  GetPZ; the reference list is too long to include here.
}
\author{Jake Anderson}
\note{
  The STS-2 and Trillium 240 come in multiple generations, each with a
  slightly different response.  For the Trillium 240,
  serial numbers less than 400 belong to generation 1, while serial
  numbers greater than or equal to 400 are in generation 2.  To
  determine which generation an STS-2 is, see
  http://www.iris.edu/NRL/sensors/streckeisen/streckeisen_sts2_sensors.htm .  
  
  Certain short-period instruments were intentionally omitted
  because their responses depend on installation-specific
  parameters and are therefore not completely standardized.
  Given that discrete responses of short-period instruments
  can be calculated quickly, I consider the convenience of having
  pre-calculated responses for these instruments to not be worth
  the risk of the user selecting the wrong response and getting
  inaccurate results.
}
\seealso{
  ReadInstr
  MakeRespSP
  GetDPZ
  MakeDPZ
}
\examples{
# responses for 3T sampling at 1 Hz and 40T (30 s) 
PZLIST = GetPZ(c(4, 12))
}
\keyword{ misc }


\name{PZLIST}
\alias{PZLIST}
\docType{data}
\title{List of Continuous Instrument Responses}
\description{
List of continuous instrument responses of 14 common seismometers.
}
\usage{ data(PZLIST) }
\format{
  List of 14 lists, each consisting of:
  \describe{
    \item{Sense}{Instrument passband sensitivity (V * s/m)}
    \item{Knorm}{Normalization constant}
    \item{poles}{Poles of Laplace transform of instrument response
      (rad/s)}
    \item{np}{Number of poles}
    \item{zeros}{Zeros of Laplace transform of instrument response
      (rad/s)}
    \item{nz}{Number of zeros}
  }
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

\references{
  Each of these responses was drawn from either a manufacturer document
  or from IRIS; sources are noted in the comments of GetPZ.
}
\examples{
# 40T response:
data(PZLIST)
PZ = PZLIST[[12]]
}
\keyword{datasets}

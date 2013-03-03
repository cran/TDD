\name{DPZLIST}
\alias{DPZLIST}
\docType{data}
\title{List of Pre-Calculated Discrete Instrument Responses}
\description{
List of discrete instrument responses of 14 common seismometers for 6
common sample intervals.
}
\usage{ data(DPZLIST) }
\format{
  List of 14 lists (each corresponding to a different seismometer), each
  including 6 lists (each corresponding to a different sample interval),
  each containing the following elements:
  \describe{
    \item{Sense}{Instrument passband sensitivity (V * s/m)}
    \item{Knorm}{Normalization constant}
    \item{poles}{Poles of Laplace transform of instrument response
      (rad/s)}
    \item{np}{Number of poles}
    \item{zeros}{Zeros of Laplace transform of instrument response
      (rad/s)}
    \item{nz}{Number of zeros}
    \item{dt}{Sample interval for this response (s)}
    \item{fmax}{Maximum frequency for which this digital response
      matches the true analog response of the sensor within 1\%}
    \item{Zpg}{Zpg-class element (from package 'signal') giving the
      digital response of the filter in terms of its zeros and poles (in
      Z-transform space) and gain.}
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

  Digital responses are provided for the following six sample intervals
  (in seconds): 1, 0.1, 0.05, 0.025, 0.02, 0.01.
}
\note{
  The STS-2 and Trillium 240 come in multiple generations, each with a
  slightly different response.  For the Trillium 240,
  serial numbers less than 400 belong to generation 1, while serial
  numbers greater than or equal to 400 are in generation 2.  To
  determine which generation an STS-2 is, see
  http://www.iris.edu/NRL/sensors/streckeisen/streckeisen_sts2_sensors.htm .  
  
  Many short-period instruments were intentionally omitted
  because their responses depend on installation-specific
  parameters and are therefore not completely standardized.
  Given that discrete responses of short-period instruments
  can be calculated quickly, I consider the convenience of having
  pre-calculated responses for these instruments to not be worth
  the risk of the user selecting the wrong response and getting
  inaccurate results.
}

\examples{
# 40T, 0.01-s response:
data(DPZLIST)
DPZ = DPZLIST[[12]][[6]]
}
\keyword{datasets}

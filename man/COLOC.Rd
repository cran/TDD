\name{COLOC}
\alias{COLOC}
\docType{data}
\title{Seismic Data From Co-Located Seismometers}
\description{
  Example of seismic data structure. Channels 1-3 are from a broadband
  Guralp CMG-3T; channels 4-6 are from a short-period Guralp CMG-40T-1.
  Data were logged using a RefTek RT130 six-channel logger and are in
  volts.
  Deconvolution of instrument responses will demonstrate that the two
  sensors experience the same ground motion (see DeconSeis).
}
\usage{ data(COLOC) }
\format{
  List, consisting of:
  \describe{
    \item{JSTR}{list of digital seismic data traces}
    \item{STNS}{vector of stations}
    \item{dir}{directory}
    \item{ifile}{original file names}
    \item{COMPS}{Component names, V N E, e.g.}
    \item{OCOMPS}{Old Component names}
    \item{dt}{vector of delta-t, sampling time intervals}
    \item{KNOTES}{Notes for plotting on panels}
    \item{info}{List, detailed information about traces, including}
    \item{dat}{not used}
    \item{nn}{Number of traces}
    \item{ex}{time axis for plotting}
    \item{pcol}{colors for plotting}
    \item{ok}{which traces are okay}
    \item{wintim}{window span time, seconds}
    \item{ftime}{alphanumeric time stamp}
    \item{pickfile}{pickfile, see below}
    \item{velfile}{velocity model list}
    \item{stafile}{station information list including lat, lon, z}
    \item{aname}{source name for loading}
    \item{UWFILEID}{event ID number}
  }
  
  The info list consists of:
  \describe{
    \item{fn}{file name}
    \item{name}{identification name}
    \item{yr}{start year}
    \item{jd}{start julianday}
    \item{mo}{month}
    \item{dom}{day of month}
    \item{hr}{hour}
    \item{mi}{minute}
    \item{sec}{second}
    \item{msec}{millisecond}
    \item{dt}{delta-t}
    \item{t1}{time 1}
    \item{t2}{time 2}
    \item{off}{offset}
    \item{n1}{number of samples}
    \item{n2}{not used}
    \item{n3}{not used}
    \item{n}{number of samples}
  }
}
\seealso{
  DeconSeis
}
\examples{
library(RSEIS)
data(COLOC)
swig(COLOC)
}
\keyword{datasets}

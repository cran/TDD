\name{ReadInstr}
\alias{ReadInstr}
\title{Read IRIS Instrument Response File}
\description{
  Scans an instrument response file from IRIS and returns a list with poles, zeros, normalization constant, and sensitivity.
}
\usage{
  ReadInstr(fn)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fn}{List of filenames of instrument response files to read}
}

\value{List including the following elements:
  \item{Sense}{Sensitivity of instrument (V * s/m)}
  \item{Knorm}{Normalization constant}
  \item{poles}{Poles of Laplace transform of instrument impulse response
    (rad/s)}
  \item{np}{Number of poles}
  \item{zeros}{Zeros of Laplace transform of instrument impulse response
    (rad/s)}
  \item{nz}{Number of zeros}
}

\author{Jake Anderson}


\seealso{
  GetPZ
  MakeRespSP
}
\examples{
# not run:
# ReadInstr('SAC_PZs_IU_OTAV_BHZ_00_2009.091.00.00.00.0000_2010.136.22.12.60.99999')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }


\name{MakeRespSP}
\alias{MakeRespSP}
\title{Calculate Poles and Zeros of Short-Period Sensor}
\description{
  Many mechanical short-period seismometers are characterized by three
  parameters: the natural angular frequency (omega naught), damping coefficient
  (h), and sensitivity.  The response of these sensors to a velocity
  impulse has two poles and two zeros; the response is zero at the
  origin, increases roughly proportionately with f^2 up to a low corner,
  and is flat above the low corner.

  The differential equation describing this, where y is the output
  voltage and v is the velocity of the ground, is

  y'' + 2*h*omega_0*y' + omega_0^2*y = Sense * v''
  Some short-period sensors are customizable, so it is very important to
  make sure you use the correct parameters for your installation here.
}
\usage{
MakeRespSP(h, omega0, Sense, f_Sense = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{h}{Damping coefficient (unitless, 1 for critical damping)}
  \item{omega0}{Natural angular frequency (rad/s)}
  \item{Sense}{Sensitivity in passband (V * s/m)}
  \item{f_Sense}{If given, the frequency (Hz) at which Sense is valid.  If
  NULL (which should ordinarily be the case), Sense is assumed to be
  valid at frequencies much higher than the low corner.}
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
  ReadInstr
}
\examples{
# L4C3D
omega0 = 2*pi # 1 Hz natural frequency * 2pi
h = 0.707
Sense = 171
MakeRespSP(h, omega0, Sense)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }


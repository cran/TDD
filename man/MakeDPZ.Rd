\name{MakeDPZ}
\alias{MakeDPZ}
\title{Calculate Find Digital Match to Analog Poles/Zeros}
\description{Calculates digital poles and zeros to match a continuous
  instrument response given in poles, zeros, and sensitivity.
  Discretization effects mean that the analog poles and zeros do not work
  for finite sample rates, with discrepancies increasing as the time
  interval increases.

  This function uses two methods to match the responses.  The first 
  uses a finite difference approximation and Markov Chain Monte Carlo
  (MCMC) routine to optimize the response.  The second uses the bilinear
  transform to approximate the analog response.  Whichever of these
  responses best matches the analog response is used; usually, the
  first method provides a better fit.}
\usage{
MakeDPZ(PZ, dt, fmin = 1/360, niter = 50000, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PZ}{List including poles and zeros of instrument response}
  \item{dt}{Sample interval (s)}
  \item{fmin}{Lowest frequency to match (Hz)}
  \item{niter}{Number of iterations in Markov Chain Monte Carlo}
  \item{...}{Additional arguments for MatchCoefDPZ}
}
\details{
  Large N allow it to match very low frequencies, but take longer to
  calculate.  Large niter means longer calculation time, but probably a
  closer match.  The burn-in period should be set to zero unless you
  want the posterior distribution of the poles and zeros.  Large sigfac
  means that the MCMC makes smaller jumps, meaning it explores the
  sample space more slowly, but is less likely to make large jumps away
  from the interesting region.  Note that the standard deviations of
  proposal distributions of the model parameters are proportional to the
  magnitude of the "guess" model--meaning that model parameters
  identically equal to zero (such as zeros at the origin) are fixed.

  Discretization effects often make it difficult to match higher
  frequencies.  Close match of somewhat high frequencies is done at the
  expense of poor match of very high frequencies.  If very high
  frequencies are not interesting, fh should be left at its default
  value.  Otherwise, it should be set to the highest interesting frequency.

}
\value{List including the following elements:
  \item{poles}{Vector of "analog poles" (rad/s)}
  \item{zeros}{Vector of "analog zeros" (rad/s)}
  \item{np}{Number of poles}
  \item{nz}{Number of zeros}
  \item{Knorm}{Normalization constant}
  \item{Sense}{Sensitivity (V * s/m)}
  \item{dt}{Sample interval (s)}
  \item{fmax}{Maximum frequency for which this digital response
    matches the true analog response of the sensor within 1\%}
  \item{Zpg}{Zpg-class element (from package 'signal') giving the
    digital response of the filter in terms of its zeros and poles (in
    Z-transform space) and gain.}
}
\author{Jake Anderson}
\note{
  This is a wrapper function for MatchCoefDPZ, and should be used for
  most applications.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  MatchCoefDPZ
}
\examples{
# Response of Guralp CMG-40T

PZ = list(poles = c(-0.149 + 0.149i, -0.149 - 0.149i, -503, -1010,
-1130), zeros = c(0, 0), Knorm = 574095649, Sense = 800)
# MakeDPZ(PZ, dt = 0.01, fmin = 1/60) # takes minutes to run

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }


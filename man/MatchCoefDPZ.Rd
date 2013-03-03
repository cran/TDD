\name{MatchCoefDPZ}
\alias{MatchCoefDPZ}
\title{Find Digital Match to Analog Poles/Zeros}
\description{Calculates digital poles and zeros to match a continuous
  instrument response given in poles, zeros, and sensitivity.
  Discretization effects mean that the given poles and zeros do not work
  for finite sample rates, with discrepancies increasing as the time
  interval increases.  This function uses a Markov Chain Monte Carlo (MCMC)
  routine to match the responses.
}
\usage{
MatchCoefDPZ(PZ, dt, N, niter = 50000, burn = 0, sigfac = 1, fh =
0.25/dt, k = 0.001, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PZ}{List including poles and zeros of instrument response}
  \item{dt}{Sample interval (s)}
  \item{N}{Number of samples to use when matching response (higher to
    match lower frequencies)}
  \item{niter}{Number of iterations in Markov Chain Monte Carlo}
  \item{burn}{Burn-in period of MCMC}
  \item{sigfac}{Factor by which standard deviations are reduced in MCMC}
  \item{fh}{Highest frequency to try to match (default 0.25 * sampling
    rate)}
  \item{k}{Weight to give to misfit for frequencies over fh--should be
    low to prevent high frequencies from being matched at the expense of
    low frequencies}
  \item{verbose}{Logical: if TRUE, progress updates are printed to the
    screen}
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
  \item{b}{Moving Average polynomial coefficients}
  \item{a}{Autoregressive polynomial coefficients}
  \item{analogresp}{Continuous "analog" response}
  \item{digitalresp}{Response of digital filter}
  \item{inv}{Detailed MCMC results}
  \item{error}{Geometric root-mean-square error between digital and
    analog response}
  \item{DPZ}{Digital Poles and Zeros}
}
\author{Jake Anderson}
\note{
  MakeDPZ is a higher-level routine and should be used for most applications.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  MakeDPZ
}
\examples{
# Response of Guralp CMG-40T

PZ = list(poles = c(-0.149 + 0.149i, -0.149 - 0.149i, -503, -1010,
-1130), zeros = c(0, 0), Knorm = 574095649, Sense = 800)
# MatchCoefDPZ(PZ, dt = 0.01, N = 10000) # takes minutes to run
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }


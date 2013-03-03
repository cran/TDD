\name{ConvDiffCoef}
\alias{ConvDiffCoef}
\title{Convert Differential to Difference Eq. Coefficients}
\description{Returns difference equation coefficients corresponding to
  input differential equation coefficients and sample interval.
}
\usage{
ConvDiffCoef(db, dt)
}

\arguments{
  \item{db}{Vector of differential equation coefficients}
  \item{dt}{Time interval (s)}
}
\details{
  Input differential equation is of the form db[1] * f(t) + db[2] *
  f'(t) + db[3] * f''(t)...

  Output difference equation is of the form b[1] * x_i + b[2] * x_(i-1)
  + b[3] * x_(i-2)...
}
\value{
  Coefficients of difference equation.
}

\author{Jake Anderson}

\examples{
db = c(0, 0, 1) # represents f''(t)
dt = 0.1

ConvDiffCoef(db, dt)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }


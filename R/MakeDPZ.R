MakeDPZ = function(PZ, dt, fmin = 1/360, niter = 50000, ...){
  N = 2^ceiling(log(1/(fmin * dt), 2)) # to make sure the fft is fast
  X = MatchCoefDPZ(PZ, dt, N, niter = niter, ...)
  return(X$DPZ)
}

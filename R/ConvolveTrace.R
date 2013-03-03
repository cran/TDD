ConvolveTrace = function(x, DPZ){
  dt = DPZ$dt
  # obtain coefficients
  ba = PZ2Coef(DPZ, dt)

  # convolve to voltage
  volt = filter(ba$b, ba$a, x)
  return(volt)
}

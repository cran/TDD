DeconTrace = function(x, DPZ, fl = 0.05, fh = NaN, bitweight = NULL){
  # bitweight is optional A-to-D constant (when data are in counts) [V/count]
  dt = DPZ$dt
  
  x = x - mean(x)

  # set up an "inverse response" iDPZ
  iDPZ = DPZ
  iDPZ$poles = DPZ$zeros
  iDPZ$np = DPZ$nz
  iDPZ$zeros = DPZ$poles
  iDPZ$nz = DPZ$np
  iDPZ$Knorm = 1/DPZ$Knorm
  iDPZ$Sense = 1/DPZ$Sense

  # convert counts to volts, if necessary
  if(!is.null(bitweight)){
    iDPZ$Sense = iDPZ$Sense * bitweight
  }
  
  # obtain coefficients
  ba = PZ2Coef(iDPZ, dt)

  # deconvolve to velocity
  v = filter(ba$b, ba$a, x)
  
  # filter it
  if(!is.na(fl) && (fl > 0 && fl < 1/(2*dt))){
    ba = butter(2, fl*2*dt, 'high')
    v = filter(ba$b, ba$a, v)
  }
  if(!is.na(fh) && (fh > 0 && fh < 1/(2*dt))){
    ba = butter(2, fh*2*dt, 'low')
    v = filter(ba$b, ba$a, v)
  }
  
  return(v)
}

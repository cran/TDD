DeconSeis = function(GH, inst, L, fl = 0.1, fh = NaN, bitweight = NULL){
  # GH is seismogram
  # inst is vector of indices in L corresponding to traces in GH (0 for no deconvolution)
  # L is list of DPZs of instruments used
  for(i in which(inst != 0)){
    GH$JSTR[[i]] = DeconTrace(GH$JSTR[[i]], L[[inst[i]]], fl = fl, fh = fh, bitweight = bitweight[i])
    GH$units[i] = 'm/s'
  }
  GH$process = c(GH$process, paste('DeconRec', fl, fh))
  invisible(GH)
}

#########################################################################
#      $Log: Cwt_Morlet.S,v $
#########################################################################
#
#               (c) Copyright  1997                             
#                          by                                   
#      Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   
#                  Princeton University
#                  All right reserved                           
#########################################################################

cwt <- function(input, noctave, nvoice = 1, w0 = 2*pi, twoD = TRUE,
                plot = TRUE)
#########################################################################
#      cwt:
#      ---
#       continuous wavelet transform function
# 	compute the continuous wavelet transform with (complex-valued)
#			Morlet wavelet
#
#       Input:
#       ------
# 	 input: input signal (possibly complex-valued)
#	 noctave: number of powers of 2 for the scale variable
#	 nvoice: number of scales between 2 consecutive powers of 2
#        w0: central frequency of Morlet wavelet
#	 twoD: if set to TRUE, organizes the output as a 2D array 
#			(signal_size X nb_scales)
#		      if not: 3D array (signal_size X noctave X nvoice)
#	 plot: if set to TRUE, displays the modulus of cwt on the graphic
#		device.
#
#       Output:
#       -------
#        tmp: continuous (complex) wavelet transform
#
#########################################################################
{
  oldinput <- input
  isize <- length(oldinput)
  
  tmp <- adjust.length(oldinput)
  input <- tmp$signal
  newsize <- length(input)
  
  pp <- noctave * nvoice
  Routput <- matrix(0,newsize,pp)
  Ioutput <- matrix(0,newsize,pp)
  output <- matrix(0,newsize,pp)
  dim(Routput) <- c(pp * newsize,1)
  dim(Ioutput) <- c(pp * newsize,1)
  dim(input) <- c(newsize,1)
  
  z <- .C("Scwt_morlet",
          as.double(Re(input)),
          as.double(Im(input)),
          Rtmp = as.double(Routput),
          Itmp = as.double(Ioutput),
          as.integer(noctave),
          as.integer(nvoice),
          as.integer(newsize),
          as.double(w0),
          PACKAGE="Rcwt")
  
  Routput <- z$Rtmp
  Ioutput <- z$Itmp
  dim(Routput) <- c(newsize,pp)
  dim(Ioutput) <- c(newsize,pp)
  if(twoD) {
    output <- Routput[1:isize,] + 1i*Ioutput[1:isize,]
    if(plot) {
      image(Mod(output), xlab="Time", ylab="log(scale)",
            main="Wavelet Transform Modulus")
    }
    output
  } 
  else {
    Rtmp <- array(0,c(isize,noctave,nvoice))
    Itmp <- array(0,c(isize,noctave,nvoice))
    for(i in 1:noctave)
      for(j in 1:nvoice) {
        Rtmp[,i,j] <- Routput[1:isize,(i-1)*nvoice+j]
        Itmp[,i,j] <- Ioutput[1:isize,(i-1)*nvoice+j]
      }
    Rtmp + 1i*Itmp
  }
}


morlet <- function(sigsize, location, scale, w0 = 2*pi)
#########################################################################
#       morlet: Morlet's wavelet  
#       -------
#        Generates a Morlet wavelet for given location and scale
#
#       input:
#       ------
#        sigsize: signal size (dimension of the array)
#        location: location of the wavelet
#        scale: value of the scale at which the transform is computed
#        w0: central frequency of Morlet wavelet
#
#       output:
#       -------
#        z$wavelet.r + z$wavelet.i * i: wavelet (complex 1D array
#          of size sigsize)
#
#########################################################################
{
  wavelet.r <- numeric(sigsize)
  wavelet.i <- numeric(sigsize)
  
  z <- .C("morlet_time",
          as.double(w0),
          as.double(scale),
          as.integer(location),
          wavelet.r = as.double(wavelet.r),
          wavelet.i = as.double(wavelet.i),
          as.integer(sigsize),
          PACKAGE="Rcwt")
  
  z$wavelet.r + 1i*z$wavelet.i
}

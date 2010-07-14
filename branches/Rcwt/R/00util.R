#########################################################################
#       $Log: Util.S,v $
#
#               (c) Copyright  1997
#                          by                                   
#      Author: Rene Carmona,Bruno Torresani,Wen-Liang Hwang,Andra Wang   
#                  Princeton University
#                  All right reserved                           
#########################################################################

adjust.length <- function(inputdata)
#########################################################################
# adjust.length adds zeros to the end of the data if necessary so that 
# its length is a power of 2.  It returns the data with zeros added if
# nessary and the length of the adjusted data.
#
# Arguments: 
#	<inputdata> is either a text file or an S object containing
# 		    data.
#########################################################################
{
  if (is.character(inputdata))
    s <- scan(inputdata)
  else
    s <- inputdata
  np <- length(s)
  
  pow <- 1
  while ( 2*pow < np )
    pow <- 2*pow
  new.np <- 2*pow
  
  if ( np == new.np )
    list( signal=s, length=np )
  else {
    new.s <- 1:new.np
    new.s[1:new.np] <- 0
    new.s[1:np] <- s
    list( signal=new.s, length=new.np )
  }
}



cleanph <- function(tfrep, thresh = .01, plot = F)
## moved from TF_Maxima.R, a nice example of a clean, 1-purpose function
#########################################################################
#  cleanph:
#  --------
# 	sets to zero the phase of time-frequency transform when
#        modulus is below a certain value.
#
#      input:
#      ------
# 	tfrep: continuous time-frequency transform (2D array)
#       thresh: (relative) threshold.
#	plot: if set to TRUE, displays the maxima of cwt on the graphic
#          device.
#
#      output:
#      -------
#       output: thresholded phase (2D array)
#
#########################################################################
{
  thrmod1 <- thrmod2 <- Mod(tfrep)
  limit <- range(thrmod1)[2] * thresh
  thrmod1 <- (Mod(tfrep) > limit)
  thrmod2 <- (Mod(tfrep) <= limit)

  output <- thrmod1 * Arg(tfrep) - pi * thrmod2

  if(plot) image(output)
  output
}




























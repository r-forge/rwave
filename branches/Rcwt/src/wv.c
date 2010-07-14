#include <stdlib.h>

/****************************************************************
*               (c) Copyright  1998                             *
*                          by                                   *
*      Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   *
*                  Princeton University                         *
*                  All right reserved                           *
****************************************************************/

#include "Swave.h"
#include "denoise.h"

/***********************************************************
*  Function: WV_mult
*  ---------
*          Manipulation for Wigner-Ville
*
*   frequency: value of the frequency parameter (double)
*   Ri,Ii: input signal (real and imag. parts)
*   Ro,Io: output signal (real and imag. parts)
*   isize: signal size
*
************************************************************/
void WV_mult(int n, double *Ri,double *Ii,
  double *Ro, double *Io,int isize)
{
  int p,n2;
  int iplus,iminus;
  
  n2 = 2*n;

  for(p = 0; p < isize; p++) {
    iplus = (n2 + p + isize/2)%isize;
    iminus = (n2 - p + 3*isize/2)%isize;
    Ro[p] = Ri[iplus] * Ri[iminus] + Ii[iplus] * Ii[iminus];
    Io[p] = -Ri[iplus] * Ii[iminus] + Ii[iplus] * Ri[iminus];
  }
  return;
}





/***********************************************************
*  Function: WV:
*  ---------
*      Continuous Wigner-Ville transform.
*
*   input: input signal
*   Oreal, Oimage: Wigner Ville transform
*                  (real and imaginary parts)
*   pinputsize: signal size
*   pnbfreq: Number of values for the frequency
*   pfreqstep: frequency step
*
***********************************************************/
void WV(double *input,double *Oreal,double *Oimage,int *pnbfreq,
  double *pfreqstep,int *pinputsize)
{
  int nbfreq, i, p, k, ii, inputsize, locsize;
  double freqstep, frequency;
  double *Ri1, *Ii1, *Ii, *Ri, *tmpreal, *tmpimage;


  /* Initialization of S variables */
  nbfreq = *pnbfreq;
  freqstep = *pfreqstep;
  inputsize = *pinputsize;
  locsize = 2*inputsize;

  /* Memory allocation */
  if(!(Ri = (double *)S_alloc(locsize, sizeof(double))))
    error("Memory allocation failed for Ri in WV.c \n");
  if(!(Ii = (double *)S_alloc(locsize, sizeof(double))))
    error("Memory allocation failed for Ii in WV.c \n");

  if(!(Ri1 = (double *)S_alloc(locsize, sizeof(double))))
    error("Memory allocation failed for Ri1 in WV.c \n");
  if(!(Ii1 = (double *)S_alloc(locsize, sizeof(double))))
    error("Memory allocation failed for Ii1 in WV.c \n");

  if(!(tmpreal = (double *)S_alloc(locsize, sizeof(double))))
    error("Memory allocation failed for tmpreal in WV.c \n");
  if(!(tmpimage = (double *)S_alloc(locsize, sizeof(double))))
    error("Memory allocation failed for tmpimage in WV.c \n");

  /* Load signal for FFT */
  for(i = 0; i < inputsize; i++){
    *Ri = (double)(*input);
    Ri++; input++;
  }
  Ri -= inputsize;
  input -= inputsize;
  
  /* Compute short FFT of the signal */
  double_fft(Ri1,Ii1,Ri,Ii,inputsize,-1);   
  
  /* Interpolate and analytize */
  for(i=1+3*inputsize/2;i<locsize;i++){
    /* Ri1[i] = 2. * Ri1[i-inputsize];*/
    Ri1[i] = 0.0;
    /* Ii1[i] = 2. * Ii1[i-inputsize]; */
    Ii1[i] = 0.0;
  }
  for(i=1+inputsize/2;i<locsize;i++){
    Ri1[i] = 0.0;
    Ii1[i] = 0.0;
  }
  Ri1[3*inputsize/2] = 0.0;
  Ii1[3*inputsize/2] = 0.0;
  
  /* compute long inverse FFT */
  double_fft(Ri,Ii,Ri1,Ii1,locsize,1);

  /* Compute Wigner transform */

  for(p = 0; p < inputsize; p++) {
    WV_mult(p,Ri,Ii,tmpreal,tmpimage,locsize);
    double_fft(tmpreal,tmpimage,tmpreal,tmpimage,locsize,-1); 

    for(k=0;k<inputsize;k++){
      Oreal[p+k*inputsize] = tmpreal[2*k];
      Oimage[p+k*inputsize] = tmpimage[2*k];
    }

  }


  return;
}


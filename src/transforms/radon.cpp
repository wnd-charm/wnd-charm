//---------------------------------------------------------------------------


#include <math.h>
#include "radon.h"

//---------------------------------------------------------------------------
/*
pPtr -array of double- output column. a pre-allocated vector of the size 2*ceil(norm(size(I)-floor((size(I)-1)/2)-1))+3
iPtr -double *- input pixels
thetaPtr -array of double- array of the size numAngles (degrees)
numAngles -int- the number of theta angles to compute
*/
void radon(double *pPtr, double *iPtr, double *thetaPtr, int M, int N,
      int xOrigin, int yOrigin, int numAngles, int rFirst, int rSize)
{
  int k, m, n, p;           /* loop counters */
  double angle;             /* radian angle value */
  double cosine, sine;      /* cosine and sine of current angle */
  double *pr;               /* points inside output array */
  double *pixelPtr;         /* points inside input array */
  double pixel;             /* current pixel value */
  double rIdx;              /* r value offset from initial array element */
  int rLow;                 /* (int) rIdx */
  double pixelLow;          /* amount of pixel's mass to be assigned to */
                            /* the bin below Idx */
  double *yTable, *xTable;  /* x- and y-coordinate tables */
  double *ySinTable, *xCosTable;
                            /* tables for x*cos(angle) and y*sin(angle) */


  /* Allocate space for the lookup tables */
  yTable = new double[2*M];
  xTable = new double[2*N];
  xCosTable = new double[2*N];
  ySinTable = new double[2*M];

  /* x- and y-coordinates are offset from pixel locations by 0.25 */
  /* spaced by intervals of 0.5. */

  /* We want bottom-to-top to be the positive y direction */
  yTable[2*M-1] = -yOrigin - 0.25;
  for (k = 2*M-2; k >=0; k--)       
    yTable[k] = yTable[k+1] + 0.5;

  xTable[0] = -xOrigin - 0.25;
  for (k = 1; k < 2*N; k++)
    xTable[k] = xTable[k-1] + 0.5;

  for (k = 0; k < numAngles; k++) {
    angle = thetaPtr[k];
    angle = (angle*M_PI)/180;
    pr = pPtr + k*rSize;  /* pointer to the top of the output column */
    cosine = cos(angle);
    sine = sin(angle);   

    /* Radon impulse response locus:  R = X*cos(angle) + Y*sin(angle) */
    /* Fill the X*cos table and the Y*sin table.  Incorporate the */
    /* origin offset into the X*cos table to save some adds later. */
    for (p = 0; p < 2*N; p++)
      xCosTable[p] = xTable[p] * cosine - rFirst;
    for (p = 0; p < 2*M; p++)
      ySinTable[p] = yTable[p] * sine;

    /* Remember that n and m will each change twice as fast as the */
    /* pixel pointer should change. */
    for (n = 0; n < 2*N; n++) {
      pixelPtr = iPtr + (n/2)*M;
      for (m = 0; m < 2*M; m++)
      {
	pixel = *pixelPtr;
	if (pixel)
        {
	  pixel *= 0.25;                         /* 1 flop/pixel */
	  rIdx = (xCosTable[n] + ySinTable[m]);  /* 1 flop/pixel */
	  rLow = (int) rIdx;                     /* 1 flop/pixel */
	  pixelLow = pixel*(1 - rIdx + rLow);    /* 3 flops/pixel */
	  pr[rLow++] += pixelLow;                /* 1 flop/pixel */
	  pr[rLow] += pixel - pixelLow;          /* 2 flops/pixel */
        }
	if (m%2)
	  pixelPtr++;   
      }
    }
  }

  delete [] yTable;
  delete [] xTable;
  delete [] xCosTable;
  delete [] ySinTable;
}

/* vd_RadonTextures
   just change the order of the vector
   vec -pointer to double- a pre-allocated vector with 12 enteries.
//out = [in(1) in(2) in(3) in(10) in(11) in(12) in(4) in(5) in(6) in(7) in(8) in(9)];
*/
void vd_RadonTextures(double *vec)
{  double temp[12];
   int a;
   temp[0]=vec[0];
   temp[1]=vec[1];
   temp[2]=vec[2];
   temp[3]=vec[9];
   temp[4]=vec[10];
   temp[5]=vec[11];
   temp[6]=vec[3];
   temp[7]=vec[4];
   temp[8]=vec[5];
   temp[9]=vec[6];
   temp[10]=vec[7];
   temp[11]=vec[8];
   for (a=0;a<12;a++)
     vec[a]=temp[a];
}

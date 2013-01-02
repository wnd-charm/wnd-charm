/****************************************************************************/
/*                                                                          */
/*                                                                          */
/*                              mb_Znl.c                                    */
/*                                                                          */
/*                                                                          */
/*                           Michael Boland                                 */
/*                            09 Dec 1998                                   */
/*                                                                          */     
/*  Revisions:                                                              */
/*  9-1-04 Tom Macura <tmacura@nih.gov> modified to make the code ANSI C    */
/*         and work with included complex arithmetic library from           */
/*         Numerical Recepies in C instead of using the system's C++ STL    */
/*         Libraries.                                                       */
/*                                                                          */
/*  1-29-06 Lior Shamir <shamirl (-at-) mail.nih.gov> modified "factorial"  */
/*  to a loop, replaced input structure with ImageMatrix class.             */
/*  2011-04-25 Ilya Goldberg. Optimization due to this function accounting  */
/*    for 35% of the total wndchrm run-time.  Now ~4x faster.               */
/*  2012-12-13 Ilya Goldberg. Added 10x faster mb_zernike2D_2               */
/*    the feature values this algorithm produces are not the same as before */
/*    however, the weights assigned to these features in classification     */
/*    are as good or better than mb_zernike2D                               */
/****************************************************************************/


//---------------------------------------------------------------------------
#include <complex>
#include <cmath>
#include <cfloat> // Has definition of DBL_EPSILON
#include <assert.h>
#include <stdio.h>
#include "gsl/specfunc.h"

#include "cmatrix.h"

#include "zernike.h"
#define PI 3.14159265358979323846264338328
#define MAX_L 32

// This sets the maximum D parameter (15)
// The D parameter has to match MAX_D. See mb_Znl() below.
#define MAX_D 15
// This is based on the maximum D parameter being 15, which sets the number of returned values.
#define MAX_Z 72
// This is also based on the maximum D parameter - contains pre-computed factorials
#define MAX_LUT 240


//---------------------------------------------------------------------------

/* mb_Znl
  Zernike moment generating function.  The moment of degree n and
  angular dependence l for the pixels defined by coordinate vectors
  X and Y and intensity vector P.  X, Y, and P must have the same
  length
*/
void mb_Znl(double *X, double *Y, double *P, int size, double D, double m10_m00, double m01_m00, double R, double psum, double *zvalues, long *output_size) {
	static double LUT[MAX_LUT];
	static int n_s[MAX_Z], l_s[MAX_Z];
	static char init_lut=0;

	double x, y, p ;   /* individual values of X, Y, P */
	int i,m, theZ, theLUT, numZ=0;
	int n=0,l=0;
	using namespace std;

	complex<double> sum [MAX_Z];
	complex<double> Vnl;


// The LUT indexes don't work unless D == MAX_D
// To make it more flexible, store the LUT by [m][n][l].  Needs [(D+1)/2][D+1][D+1] of storage.
// Other hard-coded D values should just need changing MAX_D, MAX_Z and MAX_LUT above.
	assert (D == MAX_D);

	if (!init_lut) {
		theZ=0;
		theLUT=0;
		for (n = 0; n <= MAX_D; n++) {
			for (l = 0; l <= n; l++) {
				if ( (n-l) % 2 == 0 ) {
					for (m = 0; m <= (n-l)/2; m++) {
						LUT[theLUT] = pow((double)-1.0,(double)m) * ( (long double) gsl_sf_fact(n-m) / ( (long double)gsl_sf_fact(m) * (long double)gsl_sf_fact((n - 2*m + l) / 2) *
							(long double)gsl_sf_fact((n - 2*m - l) / 2) ) );
						theLUT++;
					}
					n_s[theZ] = n;
					l_s[theZ] = l;
					theZ++;
				}
			}
		}
		init_lut = 1;
	}

// Get the number of Z values, and clear the sums.
	for (n = 0; n <= D; n++) {
		for (l = 0; l <= n; l++) {
			if ( (n-l) % 2 == 0 ) {
				sum [numZ] = complex<double>(0.0,0.0);
				numZ++;
			}
		}
	}

// int nfoo=0;
	for(i = 0 ; i < size ; i++) {
		x = (X[i] - m10_m00)/R;
		y = (Y[i] - m01_m00)/R;
		double sqr_x2y2 = sqrt (x*x + y*y);
		if (sqr_x2y2 > 1.0) continue;

		p = P[i] / psum;
// if (nfoo < 50) printf ("(%g,%g) = %g\n",x,y,p);
// nfoo++;
		double atan2yx = atan2(y,x);
		theLUT = 0;
		for (theZ = 0; theZ < numZ; theZ++) {
			n = n_s[theZ];
			l = l_s[theZ];
			Vnl = complex<double>(0.0,0.0);
			for( m = 0; m <= (n-l)/2; m++ ) {
				Vnl += ( polar (1.0, l*atan2yx) * LUT[theLUT] * pow( sqr_x2y2, (double)(n - 2*m)) );
				theLUT++;
			}
			sum [theZ] += (conj(Vnl) * p);
		}
	}

	double preal, pimag;
	for (theZ = 0; theZ < numZ; theZ++) {
		sum [theZ] *= ((n_s[theZ]+1)/PI);
		preal = real ( sum [theZ] );
		pimag = imag ( sum [theZ] );
//printf ("[%d][%d] = (%g,%g)\n",n_s[theZ],l_s[theZ],real ( sum [theZ] ),imag ( sum [theZ] ));
		zvalues[theZ] = fabs(sqrt(preal*preal+pimag*pimag));
	}

	*output_size = numZ;
}


void mb_zernike2D_OLD(ImageMatrix *I, double D, double R, double *zvalues, long *output_size) {
	double *Y,*X,*P,psum;
	double intensity;
	int x,y,size;

	*output_size=0;
	int rows = I->height,cols = I->width;
	if (D<=0) D=15;
	if (R<=0) R=rows/2;

	Y=new double[rows*cols];
	X=new double[rows*cols];
	P=new double[rows*cols];

	readOnlyPixels I_pix_plane = I->ReadablePixels();
   /* Find all non-zero pixel coordinates and values */
	size=0;
	psum=0;
	double moment10 = 0.0, moment00 = 0.0, moment01 = 0.0;
	for (y=0;y<rows;y++)
		for (x=0;x<cols;x++) {
			intensity = I_pix_plane(y,x);
			if (intensity != 0) {
				Y[size] = y+1;
				X[size] = x+1;
				P[size] = intensity;
				psum   += intensity;
				size++;
			}
		// moments
			moment10 += (x+1) * intensity;
			moment00 += intensity;
			moment01 += (y+1) * intensity;
// if (x < 10 && y < 10) printf ("(%d,%d) = %g\n",x,y,intensity);
		}

   /* Normalize the coordinates to the center of mass and normalize
      pixel distances using the maximum radius argument (R) */

	double m10_m00 = moment10/moment00;
	double m01_m00 = moment01/moment00;
	mb_Znl (X,Y,P,size,D,m10_m00,m01_m00,R,psum,zvalues,output_size);

	delete Y;
	delete X;
	delete P;

}

/*
  Algorithms for fast computation of Zernike moments and their numerical stability
  Chandan Singh and Ekta Walia, Image and Vision Computing 29 (2011) 251–259

  Implemented from pseudo-code by Ilya Goldberg 2011-04-27
  This code is 10x faster than the previous code, and 50x faster than previous unoptimized code.
  200 x 200 image 2 GHz CoreDuo (ms):
  3813 previous unoptimized
   877 previous optimized
    75 this implementation
  The values of the returned features are different, but their Fischer scores are slightly
  better on average than the previous version, and they produce better classification in problems
  where zernike features are useful.
*/
void mb_zernike2D (ImageMatrix *I, double order, double rad, double *zvalues, long *output_size) {
	int L, N, D;

// N is the smaller of I->width and I->height
	N = I->width < I->height ? I->width : I->height;
	if (order > 0) L = (int)order;
	else L = 15;
	assert (L < MAX_L);

	if (! rad > 0.0) rad = N;
	D = (int)(rad * 2);

	static double H1[MAX_L][MAX_L];
	static double H2[MAX_L][MAX_L];
	static double H3[MAX_L][MAX_L];
	static char init=1;

	double COST[MAX_L], SINT[MAX_L], R[MAX_L];
	double Rn, Rnm, Rnm2, Rnnm2, Rnmp2, Rnmp4;

	double a,b,x, y, area, r, r2, f, const_t;
	int n,m,i,j;
	
	double AR[MAX_L][MAX_L], AI[MAX_L][MAX_L];
	
	double sum = 0;
	int cols = I->width;
	int rows = I->height;
	readOnlyPixels I_pix_plane = I->ReadablePixels();

// compute x/0, y/0 and 0/0 moments to center the unit circle on the centroid
	double moment10 = 0.0, moment00 = 0.0, moment01 = 0.0;
	double intensity;
	for (i = 0; i < cols; i++)
		for (j = 0; j < rows; j++) {
			intensity = I_pix_plane(j,i);
			sum += intensity;
			moment10 += (i+1) * intensity;
			moment00 += intensity;
			moment01 += (j+1) * intensity;
		}
	double m10_m00 = moment10/moment00;
	double m01_m00 = moment01/moment00;
			

// Pre-initialization of statics
	if (init) {
		for (n = 0; n < MAX_L; n++) {
			for (m = 0; m <= n; m++) {
				if (n != m) {
					H3[n][m] = -(double)(4.0 * (m+2.0) * (m + 1.0) ) / (double)( (n+m+2.0) * (n - m) ) ;
					H2[n][m] = ( (double)(H3[n][m] * (n+m+4.0)*(n-m-2.0)) / (double)(4.0 * (m+3.0)) ) + (m+2.0);
					H1[n][m] = ( (double)((m+4.0)*(m+3.0))/2.0) - ( (m+4.0)*H2[n][m] ) + ( (double)(H3[n][m]*(n+m+6.0)*(n-m-4.0)) / 8.0 );
				}
			}
		}
		init = 0;
	}

// Zero-out the Zernike moment accumulators
	for (n = 0; n <= L; n++) {
		for (m = 0; m <= n; m++) {
			AR[n][m] = AI[n][m] = 0.0;
		}
	}

	area = PI * rad * rad;
	for (i = 0; i < cols; i++) {
	// In the paper, the center of the unit circle was the center of the image
	//	x = (double)(2*i+1-N)/(double)D;
		x = (i+1 - m10_m00) / rad;
		for (j = 0; j < rows; j++) {
		// In the paper, the center of the unit circle was the center of the image
		//	y = (double)(2*j+1-N)/(double)D;
			y = (j+1 - m01_m00) / rad;
			r2 = x*x + y*y;
			r = sqrt (r2);
			if ( r < DBL_EPSILON || r > 1.0) continue;
			/*compute all powers of r and save in a table */
			R[0] = 1;
			for (n=1; n <= L; n++) R[n] = r*R[n-1];
			/* compute COST SINT and save in tables */
			a = COST[0] = x/r;
			b = SINT[0] = y/r;
			for (m = 1; m <= L; m++) {
				COST[m] = a * COST[m-1] - b * SINT[m-1];
				SINT[m] = a * SINT[m-1] + b * COST[m-1];
			}

		// compute contribution to Zernike moments for all 
		// orders and repetitions by the pixel at (i,j)
		// In the paper, the intensity was the raw image intensity
			f = I_pix_plane(j,i) / sum;

			Rnmp2 = Rnm2 = 0;
			for (n = 0; n <= L; n++) {
			// In the paper, this was divided by the area in pixels
			// seemed that pi was supposed to be the area of a unit circle.
				const_t = (n+1) * f/PI;
				Rn = R[n];
				if (n >= 2) Rnm2 = R[n-2];
				for (m = n; m >= 0; m -= 2) {
					if (m == n) {
						Rnm = Rn;
						Rnmp4 = Rn;
					} else if (m == n-2) {
						Rnnm2 = n*Rn - (n-1)*Rnm2;
						Rnm = Rnnm2;
						Rnmp2 = Rnnm2;
					} else {
						Rnm = H1[n][m] * Rnmp4 + ( H2[n][m] + (H3[n][m]/r2) ) * Rnmp2;
						Rnmp4 = Rnmp2;
						Rnmp2 = Rnm;
					}
					AR[n][m] += const_t * Rnm * COST[m];
					AI[n][m] -= const_t * Rnm * SINT[m];
				}
			}
		}
	}

	int numZ=0;
	for (n = 0; n <= L; n++) {
		for (m = 0; m <= n; m++) {
			if ( (n-m) % 2 == 0 ) {
				AR[n][m] *= AR[n][m];
				AI[n][m] *= AI[n][m];
				zvalues[numZ] = fabs (sqrt ( AR[n][m] + AI[n][m] ));
				numZ++;
			}
		}
	}
	*output_size = numZ;

}




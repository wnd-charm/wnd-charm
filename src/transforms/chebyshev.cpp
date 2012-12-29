/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                               */
/*    Copyright (C) 2003 Open Microscopy Environment                             */
/*         Massachusetts Institue of Technology,                                 */
/*         National Institutes of Health,                                        */
/*         University of Dundee                                                  */
/*                                                                               */
/*                                                                               */
/*                                                                               */
/*    This library is free software; you can redistribute it and/or              */
/*    modify it under the terms of the GNU Lesser General Public                 */
/*    License as published by the Free Software Foundation; either               */
/*    version 2.1 of the License, or (at your option) any later version.         */
/*                                                                               */
/*    This library is distributed in the hope that it will be useful,            */
/*    but WITHOUT ANY WARRANTY; without even the implied warranty of             */
/*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU          */
/*    Lesser General Public License for more details.                            */
/*                                                                               */
/*    You should have received a copy of the GNU Lesser General Public           */
/*    License along with this library; if not, write to the Free Software        */
/*    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA  */
/*                                                                               */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                               */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Written by:  Lior Shamir <shamirl [at] mail [dot] nih [dot] gov>              */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <stdlib.h>
#include <math.h>

#include "cmatrix.h"
#include "chebyshev.h"

//---------------------------------------------------------------------------

void TNx(double *x, double *out, int N, int height) {
	int ix,iy;
	double *temp,*temp1;
//	if( max(abs(x(:))) > 1 )
//		error(':: Cheb. Polynomials Tn :: abs(arg) > 1');
//	end;
	temp = new double[N*height];
	temp1 = new double[N*height];

// 	T = cos((ones(size(x,2),1)*(0:(N-1))).*acos(x'*ones(1,N)));
//     	T(:,1) = ones(size(x'));

	// x'*ones(1,N)
	for (ix = 0; ix < N; ix++)
		for (iy = 0; iy < height; iy++)
			temp[iy*N+ix] = x[iy];
	// acos
	for (ix = 0; ix < N; ix++)
		for (iy = 0; iy < height; iy++)
			if (fabs(temp[iy*N+ix]) > 1) temp[iy*N+ix] = 0;   /* protect from acos domain error */
			else temp[iy*N+ix] = acos(temp[iy*N+ix]);
	// ones(size(x,2),1)*(0:(N-1))
	for (ix = 0; ix < N; ix++)
		for (iy = 0; iy < height; iy++)
			temp1[iy*N+ix] = ix;
	//.*
	for (ix = 0; ix < N; ix++)
		for (iy = 0; iy < height; iy++)
			out[iy*N+ix] = temp[iy*N+ix]*temp1[iy*N+ix];
	//cos
	for (ix = 0; ix < N; ix++)
		for (iy = 0; iy < height; iy++)
			out[iy*N+ix] = cos(out[iy*N+ix]);

	for (iy = 0; iy < height; iy++)
		out[iy*N+0] = 1;

	delete [] temp1;
	delete [] temp;
}

void getChCoeff1D(double *f,double *out,double *Tj,int N,int width) {
	double *tj;
	int jj,a;

	tj = new double[width];
	for (jj = 0; jj < N; jj++) {
		int jx;
		jx = jj;
		for (a = 0; a < width; a++)
			tj[a] = Tj[a*N+jj];
		if (!jx) {
			for (a = 0; a < width; a++)
				tj[a] = tj[a]/(double)width;
		} else {
			for (a = 0; a < width; a++)
			tj[a] = tj[a]*2/(double)width;
		}
		out[jj] = 0;
		for (a = 0; a < width; a++)
			out[jj] += f[a]*tj[a]/2;
	}
	delete [] tj;
}

void getChCoeff(double *Im, double *out, double *Tj,int N,int width, int height) {
	int iy;

	for (iy = 0; iy < height; iy++) {
		getChCoeff1D(&(Im[iy*width]),&(out[iy*N]),Tj,N,width);
	}
}

/* inputs:
IM - image
N - coefficient
width - width of the image
height - height of the image
*/
void Chebyshev2D(ImageMatrix *Im, double *out, unsigned int N) {
	double *TjIn,*Tj;
	double *in;
	unsigned int a,i,j;

// Make a default value for coeficient order if it was not given as an input
//   if (N< = 0)
//     N = min(Im->width,Im->height);

	TjIn = new double[Im->width];
	for (a = 0; a < Im->width; a++)
		TjIn[a] = 2*(double)(a+1) / (double)Im->width -1;

// Pre-compute Tj on x
	Tj = new double[Im->width*N];
	TNx(TjIn,Tj,N,Im->width);


	in = new double[Im->width*Im->height];
	readOnlyPixels Im_pix_plane = Im->ReadablePixels();

	for (j = 0; j < Im->height; j++)
		for (i = 0; i < Im->width; i++)
			in[j*Im->width+i] = Im_pix_plane(j,i);
	getChCoeff(in,out,Tj,N,Im->width,Im->height);

	/* transpose the matrix "out" into "in" */
	for (j = 0; j < N; j++)
		for (i = 0; i < Im->height/*Im->width*/; i++)
			in[j*Im->height+i] = out[i*N+j];

// If the height is different, re-compute Tj
	if (Im->height != Im->width) {
		delete [] Tj;
		delete [] TjIn;
		TjIn = new double[Im->height];
		for (a = 0; a < Im->height; a++)
			TjIn[a] = 2*(double)(a+1) / (double)Im->height -1;
	// Pre-compute Tj on y
		Tj = new double[Im->height*N];
		TNx(TjIn,Tj,N,Im->height);
	}
	getChCoeff(in,out,Tj,N,Im->height,N);

	delete [] in;
	delete [] TjIn;
	delete [] Tj;
}



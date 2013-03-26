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
/* Modified by:  Ilya G. Goldberg <igg [at] nih [dot] gov> 2012-12-30            */
/*   various rearrangements and optimizations lead to ~9.8x speedup,             */
/*   and much less memory use.                                                   */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

//#include <vcl.h>

#include <math.h>
#include <stdio.h>
#include "ChebyshevFourier.h"

#define min(a, b)  (((a) < (b)) ? (a) : (b))
//---------------------------------------------------------------------------

double *ChebPol(double x,unsigned long N, double *out) {
	unsigned long a;
	if (x < -1) x = -1;      /* protect the acos function */
	if (x > 1) x = 1;
	double acos_x = acos(x);
	for (a = 0; a < N; a++)
		out[a] = cos(a * acos_x);
	out[0] = 0.5;
	return(out);
}



/*
ChebyshevFourier - Chebyshev Fourier transform
"coeff_packed" -array of doubles- a pre-allocated array of 32 doubles
*/
void ChebyshevFourier2D(const ImageMatrix &Im, unsigned long N, double *coeff_packed, unsigned int packingOrder) {
	unsigned long a,m,n,x,y,nLast,NN,Nmax,ind;
	double *sum_r,*sum_i,*f,*r,*img,*Tn,*coeff;
	double min,max;
	long *kk;

	if (N==0) N=11;
	m=Im.height;
	n=Im.width;

	img = new double[m*n];
	f   = new double[m*n];
	r   = new double[m*n];
	kk  = new long[m*n];  /* the required size of kk is equal to nLast */

	readOnlyPixels Im_pix_plane = Im.ReadablePixels();
	double x_ind,x_2, y_ind;
	double two_over_n_minus_1 = (2.0/((double)n-1));
	double two_over_m_minus_1 = (2.0/((double)m-1));
	ind = 0;
	nLast = 0;
	for (x = 0; x < n; x++) {
		x_ind = -1.0 + (double)x * two_over_n_minus_1;
		x_2 = pow (x_ind, 2);
		for (y = 0; y < m; y++) {
			img[ind]=Im_pix_plane(y,x);

			// convert cartesian to polar
			y_ind = -1.0 + (double)y * two_over_m_minus_1;
			r[ind] = sqrt( x_2 + pow (y_ind, 2) );
			f[ind] = -1 * atan2(y_ind, x_ind);
			if (r[ind] < 1)
				kk[nLast++] = ind;
			ind++;
		}
	}

	Nmax=(unsigned long)((min(m,n)-1)/2);
	if (N>Nmax) N=Nmax;
	NN = 2*N + 1;

	Tn=new double[NN];
	sum_r = new double [NN*NN];
	sum_i = new double [NN*NN];
	for (a = 0; a < NN*NN; a++) sum_r [a] = sum_i [a] = 0;

	for (ind = 0; ind < nLast; ind++) {
		double ri,fi,ii;
		unsigned long im;
		ri = r[kk[ind]];
		fi = f[kk[ind]];
		ii = img[kk[ind]];
		if (ri > 1) continue;
		ChebPol(ri*2-1, NN, Tn);
		unsigned long NN_ind = 0;
		for (im = 1; im <= NN; im++) {
			double Ftrm;
			long mf;
			mf = im-1-N;
			if (mf == 0) Ftrm = 0.5;
			else Ftrm = 1.0;
			double Ftrm_cos_mf_fi = Ftrm*cos(mf*fi);
			double Ftrm_sin_mf_fi = Ftrm*sin(-1*mf*fi);
			for (a = 0; a < NN; a++) {
				sum_r[NN_ind] += (Ftrm_cos_mf_fi*Tn[a]) * ii;
				sum_i[NN_ind] += (Ftrm_sin_mf_fi*Tn[a]) * ii;
				NN_ind++;
			}
		}
	}

	min =  INF;
	max = -INF;
	coeff = new double[NN*NN];
	for (a = 0; a < NN*NN; a++) {
		coeff[a]=sqrt( pow (sum_r[a], 2) + pow (sum_i[a], 2) );
		if (coeff[a] < min) min = coeff[a];
		if (coeff[a] > max) max = coeff[a];
	}

	for (a = 0; a < packingOrder; a++)
		coeff_packed[a] = 0;
	if (max != min) {
		for (a = 0; a < NN*NN; a++) {
			if (coeff[a] == max) coeff_packed[31] += 1;
			else coeff_packed[(int)(((coeff[a]-min)/(max-min))*32)] += 1;
		}
	}

	delete [] sum_r;
	delete [] sum_i;
	delete [] coeff;
	delete [] f;
	delete [] r;
	delete [] kk;
	delete [] img;
	delete [] Tn;
}

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


#include <math.h>
#include <cfloat> // DBL_MAX
#include "CombFirst4Moments.h"
#include "Moments.h"

#define N_COMB_SAMPLES 20

// matlab rounds away from 0, just like round() in math.h
// http://pubs.opengroup.org/onlinepubs/000095399/functions/round.html
// double round (double x) {
//   return x < 0 ? floor(x) : ceil(x);
// }
// 

//---------------------------------------------------------------------------

int matr4moments_to_hist(double matr4moments[4][N_COMB_SAMPLES], double *vec, int vec_start) {
	unsigned long a, b, vec_index, bin_index, nbins = 3;
	double bins[3];
	vec_index=vec_start;
	for (a = 0; a < 4; a++) {
		double h_min = INF, h_max = -INF, h_scale;
		// find min and max (for the bins)
		for (b = 0; b < N_COMB_SAMPLES; b++) {
			if (matr4moments[a][b] > h_max) h_max = matr4moments[a][b];
			if (matr4moments[a][b] < h_min) h_min = matr4moments[a][b];
		}

		// get the scale
		if (h_max - h_min > 0) h_scale = (double)nbins / double(h_max - h_min);
		else h_scale = 0;
		// initialize the bins
		memset(bins, 0, nbins * sizeof (double));
		// find the bins
		for (b = 0; b < N_COMB_SAMPLES; b++) {
			bin_index = (unsigned long)(( (matr4moments[a][b] - h_min)*h_scale));
			if (bin_index >= nbins) bin_index = nbins-1;
			bins[bin_index] += 1.0;
		}	   
		/* add the bins to the vector */
		for (bin_index = 0; bin_index < nbins; bin_index++)
			vec[vec_index++] = bins[bin_index];
	}
	return((int)vec_index);
}

//---------------------------------------------------------------------------
int CombFirst4Moments2D(ImageMatrix *Im, double *vec) {
	double **I,**J,**J1,z4[4]={0,0,0,0},z[4];
	double matr4moments[4][N_COMB_SAMPLES];
	long m,n,n2,m2;
	long a,x,y,ii;
	int matr4moments_index;
	int vec_count=0;
	readOnlyPixels pix_plane = Im->ReadablePixels();
	long step;
	Moments tmpMoments;
	for (a = 0; a < 4; a++)    /* initialize */
		for (matr4moments_index = 0; matr4moments_index < N_COMB_SAMPLES; matr4moments_index++)
			matr4moments[a][matr4moments_index] = 0;

	m=Im->height;
	n=Im->width;
	I=new double*[n];
	J=new double*[n];
	J1=new double*[n];
	for (a = 0; a < n; a++) {
		I[a] = new double[m];
		J[a] = new double[m];
		J1[a] = new double[m];
	}

	for (y = 0; y < m; y++) {
		for (x = 0; x < n; x++) {
			I[x][y] = y+1;
			J[x][y] = x+1;
		}
	}

	n2 = (int)(round(n/2));
	m2 = (int)(round(m/2));

	/* major diag -45 degrees */
	matr4moments_index=0;
	step = (int)(round((double)m/10));
	if (step < 1) step = 1;
	for (ii = 1-m; ii <= m; ii = ii+step) {
		for (a = 0; a < 4; a++) matr4moments[a][matr4moments_index]=z4[a];

		tmpMoments.reset();
		for (y = 0; y < m; y++) {
			for (x = 0; x < n; x++) {
				if (fabs(I[x][y] + ii - J[x][y]) < 1)
					tmpMoments.add (pix_plane(y,x));
			}
		}

		tmpMoments.momentVector(z);
		for (a = 0; a < 4; a++) matr4moments[a][matr4moments_index] = z[a];
		matr4moments_index++;
	}
	vec_count=matr4moments_to_hist(matr4moments,vec,vec_count);

	/* major diag +45 degrees */
	/* fliplr J */
	for (y = 0; y < m; y++)
		for (x = 0; x < n; x++)
			J1[x][y] = J[n-1-x][y];

	matr4moments_index=0;
	step = (int)(round((double)m/10));
	if (step < 1) step = 1;
	for (ii = 1-m; ii <= m; ii = ii+step) {
		for (a = 0; a < 4; a++) matr4moments[a][matr4moments_index]=z4[a];

		tmpMoments.reset();
		for (y = 0; y < m; y++) {
			for (x = 0; x < n; x++) {
				if (fabs(I[x][y] + ii - J1[x][y]) < 1)
					tmpMoments.add (pix_plane(y,x));
			}
        }

		tmpMoments.momentVector(z);
		for (a = 0; a < 4; a++) matr4moments[a][matr4moments_index] = z[a];
		matr4moments_index++;
	}
	vec_count=matr4moments_to_hist(matr4moments,vec,vec_count);

	/* vertical comb */
	matr4moments_index=0;
	step = (int)(round((double)n/10));
	if (step < 1) step = 1;
	for (ii = 1-n; ii <= n; ii = ii+step) {
		for (a = 0; a < 4; a++) matr4moments[a][matr4moments_index]=z4[a];

		tmpMoments.reset();
		for (y = 0; y < m; y++) {
			for (x = 0; x < n; x++) {
				if (fabs(J[x][y] + ii - n2) < 1)
					tmpMoments.add (pix_plane(y,x));
			}
		}

		tmpMoments.momentVector(z);
		for (a = 0; a < 4; a++) matr4moments[a][matr4moments_index] = z[a];
		matr4moments_index++;
	}
	vec_count=matr4moments_to_hist(matr4moments,vec,vec_count);

	/* horizontal comb */
	matr4moments_index=0;
	step = (int)(round((double)m/10));
	if (step < 1) step = 1;
	for (ii = 1-m; ii <= m; ii = ii+step) {
		for (a = 0; a < 4; a++) matr4moments[a][matr4moments_index] = z4[a];

		tmpMoments.reset();
		for (y = 0; y < m; y++) {
			for (x = 0; x < n; x++) {
				if (fabs(I[x][y] + ii - m2) < 1)
					tmpMoments.add (pix_plane(y,x));
			}
        }

		tmpMoments.momentVector(z);
		for (a = 0; a < 4; a++) matr4moments[a][matr4moments_index] = z[a];
		matr4moments_index++;
	}
	vec_count=matr4moments_to_hist(matr4moments,vec,vec_count);

	/* free the memory used by the function */
	for (a=0;a<n;a++) {
		delete [] I[a];
		delete [] J[a];
		delete [] J1[a];
	}
	delete [] I;
	delete [] J;
	delete [] J1;

	return(vec_count);
}



void vd_Comb4Moments(double *in) {
	double temp[48];
	int a;
	for (a = 0; a < 48; a++)
		temp[a] = in[a];

	in[0]  = temp[45];
	in[1]  = temp[46];
	in[2]  = temp[47];
	in[3]  = temp[36];
	in[4]  = temp[37];
	in[5]  = temp[38];
	in[6]  = temp[42];
	in[7]  = temp[43];
	in[8]  = temp[44];
	in[9]  = temp[39];
	in[10] = temp[40];
	in[11] = temp[41];
	in[12] = temp[33];
	in[13] = temp[34];
	in[14] = temp[35];
	in[15] = temp[24];
	in[16] = temp[25];
	in[17] = temp[26];
	in[18] = temp[30];
	in[19] = temp[31];
	in[20] = temp[32];
	in[21] = temp[27];
	in[22] = temp[28];
	in[23] = temp[29];
	in[24] = temp[9];
	in[25] = temp[10];
	in[26] = temp[11];
	in[27] = temp[0];
	in[28] = temp[1];
	in[29] = temp[2];
	in[30] = temp[6];
	in[31] = temp[7];
	in[32] = temp[8];
	in[33] = temp[3];
	in[34] = temp[4];
	in[35] = temp[5];
	in[36] = temp[21];
	in[37] = temp[22];
	in[38] = temp[23];
	in[39] = temp[12];
	in[40] = temp[13];
	in[41] = temp[14];
	in[42] = temp[18];
	in[43] = temp[19];
	in[44] = temp[20];
	in[45] = temp[15];
	in[46] = temp[16];
	in[47] = temp[17];
}

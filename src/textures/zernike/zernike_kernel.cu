/*zernike_kernel.cu                                                              */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                               */
/* Copyright (C) 2015                                                            */
/*                                                                               */
/*       eInfochips Ltd.                                                         */
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
/* Written by: Parikshit Dharmale						  */
/*	Email: parikshit.dharmale@einfochips.com                                  */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include <cmatrix.h>
#include <iostream>
#include <common.h>

// Include thrust library header files
#include <thrust/functional.h>
#include <thrust/fill.h>
#include <thrust/count.h>
#include <thrust/device_vector.h>

#define PI 3.14159265358979323846264338328
#define MAX_L 32

struct str_H{
	double *H1;
	double *H2;
	double *H3;
	double *data_ptr;
	double m10_m00;
	double m01_m00;
	double rad;
	double sum;

};

void free_gpu_mb_zernike2D(str_H *ptr)
{
	if(ptr)
		free(ptr);

	cudaDeviceReset();
}

__global__ void
zernike2D_kernel(int width, int height, str_H *H,double *AR, double *AI, int L)
{

	int col_no = blockIdx.x*blockDim.x + threadIdx.x; //width
	int row_no = blockIdx.y*blockDim.y + threadIdx.y; //height
	int m,n;

	//	int total = height * width;
	//	int thread_no = row_no * width + col_no;
	int itr_no = 0;

	double COST[MAX_L], SINT[MAX_L], R[MAX_L];

	double const_t ; /*cost_add, sint_add; */
	double a, b, x, y, r, r2, f;
	double Rn, Rnm, Rnm2, Rnnm2, Rnmp2, Rnmp4;


	if (col_no >= width || row_no >= height){

		return ;
	}

	x = (col_no + 1 - H->m10_m00) / H->rad;
	// In the paper, the center of the unit circle was the center of the image
	//	y = (double)(2*j+1-N)/(double)D;
	y = (row_no + 1 - H->m01_m00) / H->rad;

	r2 = x*x + y*y;
	r = sqrt (r2);

	if ( !(r < DBL_EPSILON  || r > 1.0)){

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
		f = H->data_ptr[row_no * width + col_no] / H->sum;

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
					Rnm = H->H1[n * MAX_L + m] * Rnmp4 +
						( H->H2[n * MAX_L + m] + (H->H3[n * MAX_L + m]/r2) ) * Rnmp2;
					Rnmp4 = Rnmp2;
					Rnmp2 = Rnm;
				}

				AR[itr_no * (height * width) + (row_no * width + col_no)] += const_t * Rnm * COST[m];
				AI[itr_no * (height * width) + (row_no * width + col_no)] += const_t * Rnm * SINT[m];

				itr_no++;

				//d_atomicAdd(&AR[n * MAX_L + m],cost_add);
				//d_atomicAdd(&AI[n * MAX_L + m],sint_add);


			}

		}
	}
}

extern "C" int
gpu_mb_zernike2D (const ImageMatrix &Im, double order, double rad, double *zvalues, long *output_size )
{
	int L, N;

	// N is the smaller of Im.width and Im.height
	N = Im.width < Im.height ? Im.width : Im.height;
	if (order > 0) L = (int)order;
	else L = 15;
	assert (L < MAX_L);

	if (! rad > 0.0) rad = N;

	static double H1[MAX_L][MAX_L];
	static double H2[MAX_L][MAX_L];
	static double H3[MAX_L][MAX_L];
	static char init=1;

	int n,m,i,j;

	double sum = 0;
	int cols = Im.width;
	int rows = Im.height;
	int size = cols * rows;
	int numZ=0;

	double ARR[72] , AII[72];

	double ARR_1[MAX_L][MAX_L], AII_1[MAX_L][MAX_L];


	readOnlyPixels I_pix_plane = Im.ReadablePixels();

	// Error code to check return values for CUDA calls
	cudaError_t err = cudaSuccess;

	double *d_data = NULL; //for pix_plane
	double *d_H1 = NULL;
	double *d_H2 = NULL;
	double *d_H3 = NULL;
	double *d_AR = NULL;
	double *d_AI = NULL;
	double *d_COST = NULL;
	double *d_SINT = NULL;
	double *d_R = NULL;
	double *h_data = NULL;  //host pix data ptr

	dim3 threadsPerBlock;
	dim3 numBlocks;

	struct str_H *H = NULL;
	struct str_H *d_H = NULL;

	thrust::device_ptr<double>dev_ptr;
	thrust::device_ptr<double>dev_ptr_i;


	H = (str_H *) malloc (sizeof(str_H) );

	if(H == NULL)
	{
		fprintf(stderr,"malloc failed","in gpu_mb_zernike2D\n");
		free_gpu_mb_zernike2D(NULL);
	}

	//copy pixel data to h_data
	h_data = (double *) I_pix_plane.data();

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


	//Allocate Device memory

	if (checkCudaErrors(cudaMalloc((void **)&d_H, sizeof(struct str_H))) != cudaSuccess)
	{
		fprintf(stderr, "Failed to allocate device data (error code %s)!\n", cudaGetErrorString(err));
		free_gpu_mb_zernike2D(H);
		return -1;
	}

	if (checkCudaErrors(cudaMalloc((void **)&d_data, size * sizeof(double))) != cudaSuccess)
	{
		fprintf(stderr, "Failed to allocate device data (error code %s)!\n", cudaGetErrorString(err));
		free_gpu_mb_zernike2D(H);
		return -1;
	}

	if (checkCudaErrors(cudaMalloc((void **)&d_H1, (MAX_L*MAX_L) * sizeof(double))) != cudaSuccess)
	{
		fprintf(stderr, "Failed to allocate device data (error code %s)!\n", cudaGetErrorString(err));
		free_gpu_mb_zernike2D(H);
		return -1;
	}

	if (checkCudaErrors(cudaMalloc((void **)&d_H2, (MAX_L*MAX_L) * sizeof(double))) != cudaSuccess)
	{
		fprintf(stderr, "Failed to allocate device data (error code %s)!\n", cudaGetErrorString(err));
		free_gpu_mb_zernike2D(H);
		return -1;
	}

	if (checkCudaErrors(cudaMalloc((void **)&d_H3, (MAX_L*MAX_L) * sizeof(double))) != cudaSuccess)
	{
		fprintf(stderr, "Failed to allocate device data (error code %s)!\n", cudaGetErrorString(err));
		free_gpu_mb_zernike2D(H);
		return -1;
	}

	if (checkCudaErrors(cudaMalloc((void **)&d_AR, (cols * rows * 72) * sizeof(double))) != cudaSuccess)
	{
		fprintf(stderr, "Failed to allocate device data (error code %s)!\n", cudaGetErrorString(err));
		free_gpu_mb_zernike2D(H);
		return -1;
	}

	if (checkCudaErrors(cudaMalloc((void **)&d_AI, (cols * rows * 72) * sizeof(double))) != cudaSuccess)
	{
		fprintf(stderr, "Failed to allocate device data (error code %s)!\n", cudaGetErrorString(err));
		free_gpu_mb_zernike2D(H);
		return -1;
	}

	if (checkCudaErrors(cudaMalloc((void **)&d_COST, (MAX_L) * sizeof(double))) != cudaSuccess)
	{
		fprintf(stderr, "Failed to allocate device data (error code %s)!\n", cudaGetErrorString(err));
		free_gpu_mb_zernike2D(H);
		return -1;
	}

	if (checkCudaErrors(cudaMalloc((void **)&d_SINT, (MAX_L) * sizeof(double))) != cudaSuccess)
	{
		fprintf(stderr, "Failed to allocate device data (error code %s)!\n", cudaGetErrorString(err));
		free_gpu_mb_zernike2D(H);
		return -1;
	}

	if (checkCudaErrors(cudaMalloc((void **)&d_R, (MAX_L) * sizeof(double))) != cudaSuccess)
	{
		fprintf(stderr, "Failed to allocate device data (error code %s)!\n", cudaGetErrorString(err));
		free_gpu_mb_zernike2D(H);
		return -1;
	}

	// Copy data to device

	if (checkCudaErrors(cudaMemcpy(d_data, h_data, size * sizeof(double), cudaMemcpyHostToDevice)) != cudaSuccess)
	{
		fprintf(stderr, "Failed to copy data from host to device (error code %s)!\n", cudaGetErrorString(err));
		free_gpu_mb_zernike2D(H);
		return -1;
	}

	if (checkCudaErrors(cudaMemcpy(d_H1, H1, (MAX_L*MAX_L) * sizeof(double), cudaMemcpyHostToDevice)) != cudaSuccess)
	{
		fprintf(stderr, "Failed to copy data from host to device (error code %s)!\n", cudaGetErrorString(err));
		free_gpu_mb_zernike2D(H);
		return -1;
	}

	if (checkCudaErrors(cudaMemcpy(d_H2, H2, (MAX_L*MAX_L) * sizeof(double), cudaMemcpyHostToDevice)) != cudaSuccess)
	{
		fprintf(stderr, "Failed to copy data from host to device (error code %s)!\n", cudaGetErrorString(err));
		free_gpu_mb_zernike2D(H);
		return -1;
	}

	if (checkCudaErrors(cudaMemcpy(d_H3, H3, (MAX_L*MAX_L) * sizeof(double), cudaMemcpyHostToDevice)) != cudaSuccess)
	{
		fprintf(stderr, "Failed to copy data from host to device (error code %s)!\n", cudaGetErrorString(err));
		free_gpu_mb_zernike2D(H);
		return -1;
	}

	//copy values to structure
	H->H1 = d_H1;
	H->H2 = d_H2;
	H->H3 = d_H3;
	H->data_ptr = d_data;
	H->m01_m00 = m01_m00;
	H->m10_m00 = m10_m00;
	H->rad = rad;
	H->sum = sum;

	if (checkCudaErrors(cudaMemcpy(d_H, H, sizeof(str_H), cudaMemcpyHostToDevice)) != cudaSuccess)
	{
		fprintf(stderr, "Failed to copy data from host to device (error code %s)!\n", cudaGetErrorString(err));
		free_gpu_mb_zernike2D(H);
		return -1;
	}

	if (checkCudaErrors(cudaMemset(d_AR,0,(cols * rows * 72)*sizeof(double))) != cudaSuccess)
	{
		fprintf(stderr, "Failed to copy data from host to device (error code %s)!\n", cudaGetErrorString(err));
		free_gpu_mb_zernike2D(H);
		return -1;
	}

	if (checkCudaErrors(cudaMemset(d_AI,0,(cols * rows * 72)*sizeof(double))) != cudaSuccess)
	{
		fprintf(stderr, "Failed to copy data from host to device (error code %s)!\n", cudaGetErrorString(err));
		free_gpu_mb_zernike2D(H);
		return -1;
	}

	threadsPerBlock.x = 32;
	threadsPerBlock.y = 32;

	numBlocks.x = (int) ceil ((float) cols/threadsPerBlock.x);
	numBlocks.y = (int) ceil ((float) rows/threadsPerBlock.y);


	zernike2D_kernel<<<numBlocks,threadsPerBlock>>>(/*d_data,*/cols,rows,/*d_H1,d_H2,d_H3*/d_H,d_AR,d_AI,/*m10_m00,m01_m00,rad,sum,*/L);
	cudaDeviceSynchronize();

	if (checkCudaErrors(cudaGetLastError()) != cudaSuccess)
	{
		fprintf(stderr, "Failed to launch zernike2D_kernel (error code %s)!\n", cudaGetErrorString(err));
		free_gpu_mb_zernike2D(H);
		return -1;
	}


	dev_ptr = thrust::device_pointer_cast(d_AR);

	dev_ptr_i= thrust::device_pointer_cast(d_AI);


	for(i=0; i<72; i++)
	{
		ARR[i]=thrust::reduce(dev_ptr+(i*(double)size),dev_ptr+(((i+1)*(double)size)),(double) 0,thrust::plus<double>() );
		AII[i]=thrust::reduce(dev_ptr_i+(i*(double)size),dev_ptr_i+(((i+1)*(double)size)),(double) 0,thrust::plus<double>() );

	}


	for (n = 0; n <= L; n++) {
		for (m = n; m >= 0; m -= 2) {

			ARR_1[n][m] = ARR[numZ];
			AII_1[n][m] = AII[numZ];
			numZ++;

		}
	}

	numZ=0;
	for (n = 0; n <= L; n++) {
		for (m = 0; m <= n; m++) {
			if ( (n-m) % 2 == 0 ) {
				ARR_1[n][m] *= ARR_1[n][m];
				AII_1[n][m] *= AII_1[n][m];
				zvalues[numZ] = fabs (sqrt ( ARR_1[n][m] + AII_1[n][m] ));
				numZ++;
			}
		}
	}

	*output_size = numZ;

	//de-allocation of device memory

	if(checkCudaErrors(cudaFree(d_H))!= cudaSuccess)
	{
		free_gpu_mb_zernike2D(H);
		return -1;
	}

	if(checkCudaErrors(cudaFree(d_data))!= cudaSuccess)
	{
		free_gpu_mb_zernike2D(H);
		return -1;
	}

	if(checkCudaErrors(cudaFree(d_H1))!= cudaSuccess)
	{
		free_gpu_mb_zernike2D(H);
		return -1;
	}

	if(checkCudaErrors(cudaFree(d_H2))!= cudaSuccess)
	{
		free_gpu_mb_zernike2D(H);
		return -1;
	}


	if(checkCudaErrors(cudaFree(d_H3))!= cudaSuccess)
	{
		free_gpu_mb_zernike2D(H);
		return -1;
	}


	if(checkCudaErrors(cudaFree(d_AR))!= cudaSuccess)
	{
		free_gpu_mb_zernike2D(H);
		return -1;
	}


	if(checkCudaErrors(cudaFree(d_AI))!= cudaSuccess)
	{
		free_gpu_mb_zernike2D(H);
		return -1;
	}

	if(checkCudaErrors(cudaFree(d_COST))!= cudaSuccess)
	{
		free_gpu_mb_zernike2D(H);
		return -1;
	}


	if(checkCudaErrors(cudaFree(d_SINT))!= cudaSuccess)
	{
		free_gpu_mb_zernike2D(H);
		return -1;
	}

	if(checkCudaErrors(cudaFree(d_R))!= cudaSuccess)
	{
		free_gpu_mb_zernike2D(H);
		return -1;
	}

	if(H)
		free(H);

	return 0;
}

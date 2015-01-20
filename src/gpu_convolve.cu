/*	gpu_convolve.cu	                                                          */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                               */
/* Copyright (C) 2015                                                            */
/*                                                                               */
/*       <eInfochips Ltd.>                                                       */
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
/* Written by:                                                                   */
/*      Pratik Bari                                                              */
/*      pratik.bari@einfochips.com                                               */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include"gpu_convolve.cuh"
#include "common.h"

/************************************************* Function to Free the GPU Memory ****************************************************************/

void free_func_convolve() {

	cudaDeviceReset();
}

/*********************************** Device function to perform atomicAdd operation for double variable ******************************************/

__device__ double atomicAdd_d(double* address, double val)
{
    unsigned long long int* address_as_ull =
                              (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                               __longlong_as_double(assumed)));

    // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);

    return __longlong_as_double(old);
}

/********************************************** The kernel to perform convolve operation *********************************************************/

__global__ void
__launch_bounds__(1024, 1)				// To ensure atleast 1 block runs on the SMX
gpu_convolve_kernel(double *d_pix_img, double *d_filter, double *d_temp, unsigned long height, unsigned long width,
			unsigned long height2_old, unsigned long width2_old) {

	unsigned long i;
	unsigned long thread_cnt = blockDim.x * gridDim.x;
	unsigned long tid = threadIdx.x + blockDim.x * blockIdx.x;

	for(i=tid; i<(height * width); i=i+thread_cnt) {

		unsigned long j, k, rownum, colnum;
		long xx, yy;

		rownum = i/width;
		colnum = i%width;

		unsigned long height2 = height2_old/2;
		unsigned long width2 = width2_old/2;

		for(j=0; j<=(2 * width2); j++) {
			xx = colnum + j - width2;
			if(xx < width && xx >= 0) {
				for(k=0; k<=(2 * height2); k++) {
					yy = rownum + k - height2;
					if(yy < height && yy >= 0) {
						atomicAdd_d(&d_temp[i], (double)((d_filter[j + width2_old * k]) * (d_pix_img[xx + width * yy])));
					}
				}
			}
		}
	}
}

extern "C" void gpu_convolve(double *h_pix_img, double *h_filter, unsigned long height, unsigned long width,
			unsigned long height2, unsigned long width2, double *h_temp) {

	double *d_pix_img = NULL, *d_filter = NULL, *d_temp = NULL;

/************************************************* Device Memory Allocation *******************************************************************/

	if(checkCudaErrors(cudaMalloc((void **)&d_pix_img, sizeof(double) * height * width)) != cudaSuccess){	// Device Memory for Image
		free_func_convolve();
		exit(1);
	}

	if(checkCudaErrors(cudaMalloc((void **)&d_filter, sizeof(double) * height2  * width2)) != cudaSuccess){	// Device Memory for Filter
		free_func_convolve();
		exit(1);
	}

	if(checkCudaErrors(cudaMalloc((void **)&d_temp, sizeof(double) * height * width)) != cudaSuccess) {	// Device Memory to store Results
		free_func_convolve();
		exit(1);
	}

/********************************************************* Host To Device Memcpy Operations ***************************************************/

	if(checkCudaErrors(cudaMemcpy(d_pix_img, h_pix_img, sizeof(double) * height * width, cudaMemcpyHostToDevice)) != cudaSuccess) {
		// Image Memcpy
		free_func_convolve();
		exit(1);
	}

	if(checkCudaErrors(cudaMemcpy(d_filter, h_filter, sizeof(double) * height2 * width2, cudaMemcpyHostToDevice)) != cudaSuccess) {
		// Filter Memcpy
		free_func_convolve();
		exit(1);
	}

	if(checkCudaErrors(cudaMemset(d_temp, 0, sizeof(double) * height * width)) != cudaSuccess) {
		// Memset the Result to 0
		free_func_convolve();
		exit(1);
	}

/*********************************************************** Kernel Call **********************************************************************/

	gpu_convolve_kernel<<<NUMB, NUMT>>>(d_pix_img, d_filter, d_temp, height, width, height2, width2);
	cudaDeviceSynchronize();

/******************************************************** Device to Host Memcpy Operations ****************************************************/

	if(checkCudaErrors(cudaMemcpy(h_temp, d_temp, sizeof(double) * height * width, cudaMemcpyDeviceToHost)) != cudaSuccess) {
		// Result Memcpy
		free_func_convolve();
		exit(1);
	}

	if(checkCudaErrors(cudaFree(d_pix_img)) != cudaSuccess) {
		free_func_convolve();
		exit(1);
	}

	if(checkCudaErrors(cudaFree(d_filter)) != cudaSuccess) {
		free_func_convolve();
		exit(1);
	}

	if(checkCudaErrors(cudaFree(d_temp)) != cudaSuccess) {
		free_func_convolve();
		exit(1);
	}
}

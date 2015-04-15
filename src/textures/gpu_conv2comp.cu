/*     gpu_conv2comp.cu                                                          */
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

#include "gpu_conv2comp.cuh"
#include "common.h"

/********************** Device Function call to multiply the Image pixel with the Kernel Image pixel and perform addition *************************/

__device__ void compute_conv(int row, int col, double2 *d_c, double *d_a, double2 *d_b, int *o_row_vect, int *o_col_vect, int ma,
				int na, int mb, int nb, int mc, int nc) {

	int count_row = o_row_vect[row];
	int count_col = o_col_vect[col];
	int row_idx;
	int col_idx;
	int k_row_idx;
	int k_col_idx;
	int k_col_start_idx;

	int i_row_idx;
	int i_col_idx;
	int i_col_start_idx;

	k_row_idx = row - (ma - 1);
	k_row_idx = k_row_idx < 0 ? 0:k_row_idx;

	k_col_start_idx = col - (na - 1);
	k_col_start_idx = k_col_start_idx < 1? 0: k_col_start_idx;
	k_col_idx = k_col_start_idx;

	i_row_idx = row > (ma - 1) ? (ma - 1) : row;
	i_col_idx = col > (na - 1) ? (na - 1) : col;
	i_col_start_idx = i_col_idx;

	for ( row_idx = 0; row_idx < count_row; row_idx++) {
		for (col_idx = 0; col_idx < count_col; col_idx++) {

			d_c[col + nc * row].x += d_a[i_col_idx + na * i_row_idx] * d_b[k_col_idx + nb * k_row_idx].x;
			d_c[col + nc * row].y += d_a[i_col_idx + na * i_row_idx] * d_b[k_col_idx + nb * k_row_idx].y;

			k_col_idx++;
			i_col_idx--;
		}
		k_row_idx++;
		i_row_idx--;
		k_col_idx = k_col_start_idx;
		i_col_idx = i_col_start_idx;

	}
}

/************************************************** Kernel to perform the convolution *******************************************************/

__global__ void kernel_conv(double2 *d_c, double *d_a, double2 *d_b, int *d_row_vect, int *d_col_vect, int ma, int na, int mb, int nb, int mc, int nc) {

	int i, idx;
	int rownum, colnum, num_threads;

	idx = threadIdx.x + blockIdx.x * blockDim.x;
	num_threads = gridDim.x * blockDim.x;

	for(i=idx; i< (mc *nc); i=i+num_threads){

		rownum = i / nc;
		colnum = i % nc;

		// Device Function call to multiply the Image pixel with the Kernel Image pixel and perform addition
		compute_conv(rownum, colnum, d_c, d_a, d_b, d_row_vect, d_col_vect, ma, na, mb, nb, mc, nc);
	}
}

/***************************************** Function to generate Row vector and Column Vector ***************************************************/

void generate_vect(int *vect, int len, int dim) {

	int i;
        int min = 0;
        int max = len - 1;

        for (i = 0; i < dim; i++) {

                vect[min++] = i+1;
                vect[max--] = i+1;

        }

        while(min <= max) {
                vect[min++] = i;
                vect[max--] = i;
        }

}

/*********************************************** Function to find minimum of the two integers ************************************************/

int min_dim(int a, int b) {

        return (a < b ? a : b);
}

/**************************************************** Function for the Memory Free *********************************************************/

void free_func_conv2comp(int *h_row_vect, int *h_col_vect) {

	if(h_row_vect != NULL){
		free(h_row_vect);
		h_row_vect = NULL;
	}

	if(h_col_vect != NULL){
		free(h_col_vect);
		h_col_vect = NULL;
	}

	cudaDeviceReset();
}

/*
	ma = Height of the Image
	na = Width of the Image
	mb = Height of the Kernel Image
	nb = Width of the Kernel Image
*/

extern "C" void gpu_conv2comp(double *c, double *a, double *b, int na, int ma, int nb, int mb) {

	double *d_a = NULL;
	double2 *d_b = NULL, *d_c = NULL;
	int *h_row_vect = NULL, *h_col_vect = NULL;
	int *d_row_vect = NULL, *d_col_vect = NULL;

	int mc, nc, dim_h, dim_w;

	mc = ma+mb-1;		// Height of the Output
	nc = na+nb-1;		// Width of the Output

/******************************************************* Host Memory Allocations *********************************************************/

	h_row_vect = (int *)malloc(sizeof(int) * mc);		// Host Memory for the row vector
	if(h_row_vect == NULL) {
		printf("Error in malloc : h_row_vect\n");
		free_func_conv2comp(h_row_vect, h_col_vect);
		exit(1);
	}

	h_col_vect = (int *)malloc(sizeof(int) * nc);		// Host Memory for the column vector
	if(h_col_vect == NULL) {
		printf("Error in malloc : h_col_vect\n");
		free_func_conv2comp(h_row_vect, h_col_vect);
		exit(1);
	}

/***************************************************** Device Memory Allocations *********************************************************/

	if(checkCudaErrors(cudaMalloc((void **)&d_a, sizeof(double) * ma * na)) != cudaSuccess) {	// Device Memory for Image
		free_func_conv2comp(h_row_vect, h_col_vect);
		exit(1);
	}

	if(checkCudaErrors(cudaMalloc((void **)&d_b, sizeof(double2) * mb * nb)) != cudaSuccess) {	// Device Memory for Kernel Image
		free_func_conv2comp(h_row_vect, h_col_vect);
		exit(1);
	}

	if(checkCudaErrors(cudaMalloc((void **)&d_c, sizeof(double2) * mc * nc)) != cudaSuccess) {	// Device Memory for Output Array
		free_func_conv2comp(h_row_vect, h_col_vect);
		exit(1);
	}

	if(checkCudaErrors(cudaMalloc((void **)&d_row_vect, sizeof(int) * mc)) != cudaSuccess) {	// Device Memory for Row vector
		free_func_conv2comp(h_row_vect, h_col_vect);
		exit(1);
	}

	if(checkCudaErrors(cudaMalloc((void **)&d_col_vect, sizeof(int) * nc)) != cudaSuccess) {	// Device Memory for Column vector
		free_func_conv2comp(h_row_vect, h_col_vect);
		exit(1);
	}

/************************************************* Host to Device Memcpy operations *********************************************************/

	if(checkCudaErrors(cudaMemcpy(d_b, b, sizeof(double2) * mb * nb, cudaMemcpyHostToDevice)) != cudaSuccess) {	// Kernel Image Memcpy
		free_func_conv2comp(h_row_vect, h_col_vect);
		exit(1);
	}

	if(checkCudaErrors(cudaMemcpy(d_a, a, sizeof(double) * ma * na, cudaMemcpyHostToDevice)) != cudaSuccess) {	// Image Memcpy
		free_func_conv2comp(h_row_vect, h_col_vect);
		exit(1);
	}

	if(checkCudaErrors(cudaMemset(d_c, 0, sizeof(double2) * mc * nc)) != cudaSuccess) {
		// Initializing the initial value of the Output Array to 0
		free_func_conv2comp(h_row_vect, h_col_vect);
		exit(1);
	}

	dim_h = min_dim(ma, mb);	// Function call to calculate minimum of the height of the Image or Kernel Image
	dim_w = min_dim(na, nb);	// Function call to calculate minimum of the width of the Image or Kernel Image

	generate_vect(h_row_vect, mc, dim_h);	// Function call to calculate the Row vector
        generate_vect(h_col_vect, nc, dim_w);	// Function call to calculate the Column vector

	if(checkCudaErrors(cudaMemcpy(d_row_vect, h_row_vect, sizeof(int) * mc, cudaMemcpyHostToDevice)) != cudaSuccess) {
		// Row vector Memcpy
		free_func_conv2comp(h_row_vect, h_col_vect);
		exit(1);
	}

	if(checkCudaErrors(cudaMemcpy(d_col_vect, h_col_vect, sizeof(int) * nc, cudaMemcpyHostToDevice)) != cudaSuccess) {
		// Column vector Memcpy
		free_func_conv2comp(h_row_vect, h_col_vect);
		exit(1);
	}

	kernel_conv<<<NUMB, NUMT>>>(d_c, d_a, d_b, d_row_vect, d_col_vect, ma, na, mb, nb, mc, nc);	// Kernel call
	cudaDeviceSynchronize();

/******************************************************* Device to Host Memcpy operations ***********************************************************/

	if(checkCudaErrors(cudaMemcpy(c, d_c, sizeof(double) * mc * nc * 2, cudaMemcpyDeviceToHost)) != cudaSuccess) {
		// Output Array Memcpy
		free_func_conv2comp(h_row_vect, h_col_vect);
		exit(1);
	}

	if(checkCudaErrors(cudaFree(d_a)) != cudaSuccess) {
		free_func_conv2comp(h_row_vect, h_col_vect);
		exit(1);
	}

	if(checkCudaErrors(cudaFree(d_b)) != cudaSuccess) {
		free_func_conv2comp(h_row_vect, h_col_vect);
		exit(1);
	}

	if(checkCudaErrors(cudaFree(d_c)) != cudaSuccess) {
		free_func_conv2comp(h_row_vect, h_col_vect);
		exit(1);
	}

	if(checkCudaErrors(cudaFree(d_row_vect)) != cudaSuccess) {
		free_func_conv2comp(h_row_vect, h_col_vect);
		exit(1);
	}

	if(checkCudaErrors(cudaFree(d_col_vect)) != cudaSuccess) {
		free_func_conv2comp(h_row_vect, h_col_vect);
		exit(1);
	}

	free(h_row_vect);
	free(h_col_vect);
}

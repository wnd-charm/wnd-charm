/*     fft2.cu		                                                         */
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

#include "fft2.cuh"
#include "common.h"

/**************************************************** Function for the Memory Free **********************************************************/

void free_fft(double *in, cufftDoubleComplex *out) {

	if(in != NULL) {
		free(in);
		in = NULL;
	}

	if(out != NULL) {
		free(out);
		out = NULL;
	}

	cudaDeviceReset();
}

extern "C" void
gpu_fft2 (ImageMatrix *obj, const ImageMatrix &matrix_IN) {

	cufftHandle plan;
	double *in = NULL;
	cufftDoubleComplex *out = NULL;
	cufftDoubleReal *d_in = NULL;
	cufftDoubleComplex *d_out = NULL;

	unsigned int half_height = matrix_IN.height/2+1;
	unsigned int height = matrix_IN.height;
	unsigned int width = matrix_IN.width;

	obj->copyFields(matrix_IN);
	obj->allocate (matrix_IN.width, matrix_IN.height);
	writeablePixels out_plane = obj->WriteablePixels();
	readOnlyPixels in_plane = matrix_IN.ReadablePixels();

/********************************************************* Host Memory Allocations ************************************************************/

	in = (double*)malloc(sizeof(double) * width * height);	 // Host Memory for input
	if(in == NULL){
		printf("Error in malloc : 'in' in fft2");
		free_fft(in, out);
		exit(1);
	}

	out = (cufftDoubleComplex *)malloc(sizeof(cufftDoubleComplex) * width * height); // Host memory for output
	if(out == NULL){
		printf("Error in malloc : 'out' in fft2");
		free_fft(in, out);
		exit(1);
	}

/******************************************************* Device Memory Allocations *************************************************************/

	if(checkCudaErrors(cudaMalloc((void **)&d_in, sizeof(cufftDoubleReal) * width * height)) != cudaSuccess) { // Device Memory for input
		free_fft(in, out);
		exit(1);
	}

	if(checkCudaErrors(cudaMalloc((void **)&d_out, sizeof(cufftDoubleComplex) * width * height)) != cudaSuccess) {	// Device Memory for output
		free_fft(in, out);
		exit(1);
	}

	unsigned int x,y;
	for (x=0;x<width;x++){
		for (y=0;y<height;y++)
			in[height*x+y]=in_plane.coeff(y,x); // Initialization of the data in host input
	}

/***************************************************** Host To Device Memcpy operations ***********************************************************/

	if(checkCudaErrors(cudaMemcpy(d_in, in, sizeof(double) * width * height, cudaMemcpyHostToDevice)) != cudaSuccess) { // Input Memcpy
		free_fft(in, out);
		exit(1);
	}

/*************************************************************** Plan Creation *********************************************************************/

	if (checkCudaErrors(cufftPlan2d(&plan, width, height, CUFFT_D2Z)) != CUFFT_SUCCESS){
		free_fft(in, out);
		exit(1);
	}

/************************************************************ Execution of the Plan ***************************************************************/

	if (checkCudaErrors(cufftExecD2Z(plan, d_in, d_out)) != CUFFT_SUCCESS){
		free_fft(in, out);
		exit(1);
	}

/***************************************************** Device To Host Memcpy operations **********************************************************/

	if(checkCudaErrors(cudaMemcpy(out, d_out, sizeof(cufftDoubleComplex) * width * height, cudaMemcpyDeviceToHost)) != cudaSuccess) {
		free_fft(in, out);
		exit(1);
	}

	unsigned long idx;
	for (x=0;x<width;x++) {
		for (y=0;y<half_height;y++) {
			idx = half_height*x+y;
			double ri = out[idx].x;
			double im = out[idx].y;
			out_plane (y,x) = obj->stats.add (sqrt((ri * ri)+(im * im))); // sqrt(real(X).^2 + imag(X).^2)
		}
	}

	// complete the first column
	for (y=half_height;y<height;y++)
		out_plane (y,0) = obj->stats.add (out_plane (height - y, 0));

	// complete the rest of the columns
	for (y=half_height;y<height;y++)
		for (x=1;x<width;x++)   // 1 because the first column is already completed
			out_plane (y,x) =obj-> stats.add (out_plane (height - y, width - x));

	// clean up
	if(checkCudaErrors(cufftDestroy(plan)) != CUFFT_SUCCESS) {
		free_fft(in, out);
		exit(1);
	}

	if(checkCudaErrors(cudaFree(d_in)) != cudaSuccess) {
		free_fft(in, out);
		exit(1);
	}

	if(checkCudaErrors(cudaFree(d_out)) != cudaSuccess) {
		free_fft(in, out);
		exit(1);
	}

	free(in);
	free(out);
}

/* chebyshev_kernel.cu                                                           */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                               */
/* Copyright (C) 2015                                                            */
/*                                                                               */
/*       eInfochips Limited                                                      */
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
/*      Lalit Chandivade                                                         */
/*      lalit.chandivade@einfochips.com                                          */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <cuda.h>
#include <cublas_v2.h>
#include "common.h"

#include "cmatrix.h"

const dim3 numThreads(32,32);

__global__
void TNx_kernel(double *out, int N, int height, double *temp, double *temp1) {

	int ix,iy;
	int row, col;
	row = threadIdx.y + (blockIdx.y * blockDim.y);
	col = threadIdx.x + (blockIdx.x * blockDim.x);

	if (row < height && col < N) {

		ix = col;
		iy = row;
		// x'*ones(1,N)
		temp[iy*N+ix] = 2*(double)(iy+1) / (double)height -1;
		// acos
		if (fabs(temp[iy*N+ix]) > 1)
			temp[iy*N+ix] = 0;   /* protect from acos domain error */
		else
			temp[iy*N+ix] = acos(temp[iy*N+ix]);
		// ones(size(x,2),1)*(0:(N-1))
		temp1[iy*N+ix] = ix;
		out[iy*N+ix] = temp[iy*N+ix]*temp1[iy*N+ix];
		//cos
		out[iy*N+ix] = cos(out[iy*N+ix]);

		if (col == 0)
			out[iy*N+0] = 1;
	}
}


__global__ void
__launch_bounds__(1024,2)
setTj_kernel(double *d_Tj, double *d_tj, int N, int width)
{
	int row, col;
	int idx;
	double width_d = (double)width;

	row = threadIdx.y + (blockIdx.y * blockDim.y);
	col = threadIdx.x + (blockIdx.x * blockDim.x);

	width_d = col? (width_d): ( width_d * 2) ;

	if (row < width && col < N) {
		idx = row * N + col;
		d_tj[idx] = d_Tj[idx]/width_d;
	}
}

void chebyshev_free(double * ptr1, double * ptr2)
{
	printf("ptr1 = 0x%x ptr2 = 0x%x\n",ptr1,ptr2);
	if(ptr1)
		delete [] ptr1;
	if(ptr2)
		delete [] ptr2;
	/** It Destroy all allocations and reset all state on the current device in the current process*/
	cudaDeviceReset();
}

int TNx_gpu(double *out,int N,int height, cudaStream_t getChStream)
{
	double *d_out;
	double *d_temp;
	double *d_temp1;
	dim3 numBlocks;

        /** device memory allocation **/
	if(checkCudaErrors(cudaMalloc((void **)&d_out, sizeof(double) * height * N))!= cudaSuccess)
        {
           chebyshev_free(NULL,NULL);
           return -1;
        }
	if(checkCudaErrors(cudaMalloc((void **)&d_temp, sizeof(double) * height * N))!= cudaSuccess)
        {
           chebyshev_free(NULL,NULL);
           return -1;
        }
	if(checkCudaErrors(cudaMalloc((void **)&d_temp1, sizeof(double) * height * N))!= cudaSuccess)
        {
           chebyshev_free(NULL,NULL);
           return -1;
        }
        /** kernel configuration **/
        numBlocks.x = (int) ceil ((float) N/numThreads.x);
        numBlocks.y =  (int) ceil ((float) height/numThreads.y);

	TNx_kernel<<<numBlocks, numThreads, 0, getChStream>>>(d_out, N, height, d_temp, d_temp1);

        /** memory transfer Device-> host**/
	if(checkCudaErrors(cudaMemcpy(out, d_out,  sizeof(double) * height * N,  cudaMemcpyDeviceToHost))!= cudaSuccess)
        {
           chebyshev_free(NULL,NULL);
           return -1;
        }

	if(checkCudaErrors(cudaFree(d_out))!= cudaSuccess)
        {
		chebyshev_free(NULL,NULL);
		return -1;
        }

	if(checkCudaErrors(cudaFree(d_temp))!= cudaSuccess)
	{
		chebyshev_free(NULL,NULL);
		return -1;
	}
	if(checkCudaErrors(cudaFree(d_temp1))!= cudaSuccess)
	{
		chebyshev_free(NULL,NULL);
		return -1;
	}

        return 0;
}

/*
	Im = Mat1 (h x w)
	Tj = Mat2 ((h/w) x N)
	out = Mat 3 (h x N)
*/
int getChCoeff_gpu(const ImageMatrix &Im, double *in, double *out, double *Tj,int N,int width, int height, cudaStream_t getChStream) {
	double *d_Im = NULL;
	double *d_out = NULL;
	double *d_Tj = NULL;
	double *d_tj = NULL;
	int im_size, out_size, tj_size;
	const double alpha = 1.0f;
        const double beta  = 0.0f;
        cublasHandle_t handle;

	dim3 numBlocks;

	im_size = height * width * sizeof(double);
	out_size = Im.height * N * sizeof(double);
	tj_size = width * N * sizeof(double);

        /** device memory allocation **/
	if(checkCudaErrors(cudaMalloc((void **)&d_Im, im_size))!= cudaSuccess)
        {
	     chebyshev_free(NULL,NULL);
	     return -1;
        }
	if(checkCudaErrors(cudaMalloc((void **)&d_out, out_size))!= cudaSuccess)
        {
	     chebyshev_free(NULL,NULL);
	     return -1;
        }
	if(checkCudaErrors(cudaMalloc((void **)&d_Tj, tj_size))!= cudaSuccess)
        {
	     chebyshev_free(NULL,NULL);
	     return -1;
        }
	if(checkCudaErrors(cudaMalloc((void **)&d_tj, tj_size))!= cudaSuccess)
        {
	     chebyshev_free(NULL,NULL);
	     return -1;
        }

        /** memory transfer  host -> device**/
	if(checkCudaErrors(cudaMemcpy(d_Im, in, im_size, cudaMemcpyHostToDevice))!= cudaSuccess)
        {
	     chebyshev_free(NULL,NULL);
	     return -1;
        }
	if(checkCudaErrors(cudaMemcpy(d_Tj, Tj, tj_size, cudaMemcpyHostToDevice))!= cudaSuccess)
        {
	     chebyshev_free(NULL,NULL);
	     return -1;
        }

#ifdef GPU_DEBUG
	size_t p_size;
	cudaDeviceGetLimit(&p_size,cudaLimitPrintfFifoSize);
	printf("Current Size of Priint Buffer = %d\n", p_size);
	p_size = 600 * 1024 * 1024;
	cudaDeviceSetLimit(cudaLimitPrintfFifoSize, p_size);
#endif
        // kernel configuration
        numBlocks.x = (int) ceil ((float) N/numThreads.x);
        numBlocks.y =  (int) ceil ((float) width/numThreads.y);

	/* Set TJ once only */
	setTj_kernel<<<numBlocks, numThreads, 0, getChStream>>>(d_Tj, d_tj, N, width);

        /** matrix multiplication using cublass library **/
        if(checkCudaErrors(cublasCreate(&handle))!=CUBLAS_STATUS_SUCCESS)
        {
	     chebyshev_free(NULL,NULL);
	     return -1;
        }

	if(checkCudaErrors(cublasSetStream(handle, getChStream))!=CUBLAS_STATUS_SUCCESS)
        {
	     chebyshev_free(NULL,NULL);
	     return -1;
        }

	if(checkCudaErrors(cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, N,
				    height, width, &alpha, d_tj,
				    N, d_Im, width, &beta,
				    d_out, N))!=CUBLAS_STATUS_SUCCESS)
        {
	     chebyshev_free(NULL,NULL);
	     return -1;
        }


        /** memory transfer  device -> Host **/
	if(checkCudaErrors(cudaMemcpy(out, d_out, out_size, cudaMemcpyDeviceToHost))!= cudaSuccess)
        {
	     chebyshev_free(NULL,NULL);
	     return -1;
        }
        // Destroy the cu-blass handle
        if(checkCudaErrors(cublasDestroy(handle))!=CUBLAS_STATUS_SUCCESS)
        {
	     chebyshev_free(NULL,NULL);
	     return -1;
        }


        //de-allocation of device memory
	if(checkCudaErrors(cudaFree(d_Im))!= cudaSuccess)
        {
           chebyshev_free(NULL,NULL);
           return -1;
        }
	if(checkCudaErrors(cudaFree(d_out))!= cudaSuccess)
        {
           chebyshev_free(NULL,NULL);
           return -1;
        }
	if(checkCudaErrors(cudaFree(d_Tj))!= cudaSuccess)
        {
           chebyshev_free(NULL,NULL);
           return -1;
        }
	if(checkCudaErrors(cudaFree(d_tj))!= cudaSuccess)
        {
           chebyshev_free(NULL,NULL);
           return -1;
        }
        return 0;
}

extern "C"
void Chebyshev2D_gpu(const ImageMatrix &Im, double *out, unsigned int N) {
	double *Tj;
	double *in;
	unsigned int i,j;
	int ret = 0;
	cudaStream_t getChStream;
        /** cuda stream creation **/
	if(checkCudaErrors(cudaStreamCreate(&getChStream))!= cudaSuccess)
        {
	    chebyshev_free(NULL,NULL);
            exit(EXIT_FAILURE);
        }

        // Pre-compute Tj on x
	Tj = new double[Im.width*N];
	ret = TNx_gpu(Tj,N,Im.width, getChStream);
        if(ret == -1)
        {
	    chebyshev_free(NULL,NULL);
            exit(EXIT_FAILURE);
        }

#ifdef GPU_DEBUG
	printf("Dimension of Tj = %d(Im.width, H) x %d (N, Col)", Im.width, N);
	printf("Calling getChCoeff with N=%d, Im.width=%d, Im.height=%d\n", N, Im.width, Im.height);
#endif

	ret = getChCoeff_gpu(Im, (double *)Im.data_ptr(), out, Tj, N, Im.width, Im.height, getChStream);
        if(ret == -1)
        {
	    chebyshev_free(Tj,NULL);
            exit(EXIT_FAILURE);
        }

	in = new double[N*Im.height];
        if(in == NULL)
        {
	    chebyshev_free(Tj,NULL);
            exit(EXIT_FAILURE);
        }
	/* transpose the matrix "out" into "in" */
	for (j = 0; j < N; j++)
		for (i = 0; i < Im.height/*Im.width*/; i++)
			in[j*Im.height+i] = out[i*N+j];

        // If the height is different, re-compute Tj
	if (Im.height != Im.width) {
		delete [] Tj;
	// Pre-compute Tj on y
		Tj = new double[Im.height*N];
		ret = TNx_gpu(Tj,N,Im.height, getChStream);
                if(ret == -1)
                {
			chebyshev_free(Tj,NULL);
			exit(EXIT_FAILURE);
                }
	}
#ifdef GPU_DEBUG
	printf("Dimension of Tj = %d(Im.height, H) x %d (N, Col)", Im.height, N);
	printf("2. Calling getChCoeff with N=%d, width=%d, height=%d\n", N, Im.height, N);
#endif

	ret = getChCoeff_gpu(Im, in, out, Tj, N, Im.height, N, getChStream);
        if(ret == -1)
        {
	    chebyshev_free(Tj,in);
            exit(EXIT_FAILURE);
        }

	if(checkCudaErrors(cudaStreamSynchronize(getChStream))!= cudaSuccess)
        {
	    chebyshev_free(Tj,in);
            exit(EXIT_FAILURE);
        }

	if(checkCudaErrors(cudaStreamDestroy(getChStream))!= cudaSuccess)
        {
	    chebyshev_free(Tj,in);
            exit(EXIT_FAILURE);
        }


	if(in)
	    delete [] in;
        if(Tj)
	    delete [] Tj;

}

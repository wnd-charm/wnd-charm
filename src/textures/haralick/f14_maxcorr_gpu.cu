/* f14_maxcorr_gpu.cu                                                            */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                               */
/* Copyright (C) 2015                                                            */
/*                                                                               */
/*       eInfochips Ltd                                                          */
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
/*      Arifchand Mulani                                                         */
/*      arifchand.mulani@einfochips.com                                          */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/



#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <sys/types.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "CVIPtexture.h"
#include "textures/gpu_conv2comp.cuh"
#include "src/common.h"

void free_f14_maxcorr_gpu_mem(double *cpu_P)
{
	if(cpu_P != NULL){
		free(cpu_P);
		cpu_P=NULL;
	}
}

/* device function for atomicAdd for double values */
__device__ double atomicAdd_dB(double* address, double val)
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

//kennel_matrixQ <<<Ng,Ng>>> (d_P, px, py, d_Q, Ng);
__global__ void kennel_matrixQ(double *d_P, double *d_px, double *d_py, double *d_Q, int Ng)
{
	int rownum,colnum;
	int tid = threadIdx.x + blockDim.x * blockIdx.x;
	double temp;

	rownum = tid/Ng;
	colnum = tid%Ng;
	d_Q[(tid+Ng)]=0;	//Ng is added to point in next row

	for(int k=0; k<Ng; ++k)
	{
		if (d_px[rownum] && d_py[k])  // make sure to protect division by zero
		{
			temp=d_P[k+rownum*Ng]*d_P[k+colnum*Ng]/d_px[rownum]/d_py[k];
			atomicAdd_dB(&d_Q[(tid+Ng)],temp);
		}
	}
}

//f14_maxcorr_gpu(P, px, py, x, iy, Ng, gpu_Q);
extern "C"
int f14_maxcorr_gpu(double **P, double *px, double *py, double *x, double *iy, int Ng, double **gpu_Q)
{
	double *d_px=NULL, *d_py=NULL, *d_Q=NULL, *d_P=NULL, *cpu_P=NULL;
	double *d_x=NULL, *d_iy=NULL;

	int i, j;
	cpu_P=(double *)malloc(sizeof(double)*Ng*Ng);
	if(cpu_P == NULL)
	{
		printf("Error in malloc : at cpu_P in f14_maxcorr_gpu.cu\n");
		free_f14_maxcorr_gpu_mem(cpu_P);
		cudaDeviceReset();
		return -1;
	}

	/* All device allocation */
	if(checkCudaErrors(cudaMalloc((void **)&d_px, sizeof(double)*Ng)) != cudaSuccess){
		free_f14_maxcorr_gpu_mem(cpu_P);
		cudaDeviceReset();
		return -1;
	}

	if(checkCudaErrors(cudaMalloc((void **)&d_py, sizeof(double)*Ng)) != cudaSuccess){
		free_f14_maxcorr_gpu_mem(cpu_P);
		cudaDeviceReset();
		return -1;
	}

	if(checkCudaErrors(cudaMalloc((void **)&d_x, sizeof(double)*(Ng))) != cudaSuccess){
		free_f14_maxcorr_gpu_mem(cpu_P);
		cudaDeviceReset();
		return -1;
	}

	if(checkCudaErrors(cudaMalloc((void **)&d_iy, sizeof(double)*(Ng))) != cudaSuccess){
		free_f14_maxcorr_gpu_mem(cpu_P);
		cudaDeviceReset();
		return -1;
	}

	if(checkCudaErrors(cudaMalloc((void **)&d_Q, sizeof(double)*(Ng+1)*(Ng+1))) != cudaSuccess){
		free_f14_maxcorr_gpu_mem(cpu_P);
		cudaDeviceReset();
		return -1;
	}

	if(checkCudaErrors(cudaMalloc((void **)&d_P, sizeof(double)*Ng*Ng)) != cudaSuccess){
		free_f14_maxcorr_gpu_mem(cpu_P);
		cudaDeviceReset();
		return -1;
	}

	/* All host to device memcpy */
	if(checkCudaErrors(cudaMemcpy(d_px, px, sizeof(double)*Ng, cudaMemcpyHostToDevice)) != cudaSuccess){
		free_f14_maxcorr_gpu_mem(cpu_P);
		cudaDeviceReset();
		return -1;
	}

	if(checkCudaErrors(cudaMemcpy(d_py, py, sizeof(double)*Ng, cudaMemcpyHostToDevice)) != cudaSuccess){
		free_f14_maxcorr_gpu_mem(cpu_P);
		cudaDeviceReset();
		return -1;
	}

	for (i = 0; i < Ng; ++i) {
		for (j = 0; j < Ng; ++j) {
			cpu_P[j+i*Ng]= P[i][j];
		}
	}

	if(checkCudaErrors(cudaMemcpy(d_P, cpu_P, sizeof(double)*Ng*Ng, cudaMemcpyHostToDevice)) != cudaSuccess){
		free_f14_maxcorr_gpu_mem(cpu_P);
		cudaDeviceReset();
		return -1;
	}

	/* Launch kernel this kernel will calculate values for matrix Q*/
	kennel_matrixQ <<<Ng,Ng>>> (d_P, d_px, d_py, d_Q, Ng);
	if(checkCudaErrors(cudaDeviceSynchronize()) != cudaSuccess){
		free_f14_maxcorr_gpu_mem(cpu_P);
		cudaDeviceReset();
		return -1;
	}

	/* Copy back results -- store Q */
	if(checkCudaErrors(cudaMemcpy(*gpu_Q, d_Q, sizeof(double)*(Ng+1)*(Ng+1), cudaMemcpyDeviceToHost)) != cudaSuccess){
		free_f14_maxcorr_gpu_mem(cpu_P);
		cudaDeviceReset();
		return -1;
	}

	/* device ptrs cudaFree and host ptrs free */
	free_f14_maxcorr_gpu_mem(cpu_P);
	if(checkCudaErrors(cudaFree(d_px))!= cudaSuccess){
		cudaDeviceReset();
		return -1;
        }
	if(checkCudaErrors(cudaFree(d_py))!= cudaSuccess){
		cudaDeviceReset();
		return -1;
        }
	if(checkCudaErrors(cudaFree(d_Q))!= cudaSuccess){
		cudaDeviceReset();
		return -1;
        }
	if(checkCudaErrors(cudaFree(d_P))!= cudaSuccess){
		cudaDeviceReset();
		return -1;
        }
	if(checkCudaErrors(cudaFree(d_x))!= cudaSuccess){
		cudaDeviceReset();
		return -1;
        }
	if(checkCudaErrors(cudaFree(d_iy))!= cudaSuccess){
		cudaDeviceReset();
		return -1;
        }
	return 0;
}

/*CombFirst4Moments_gpu.cu                                                       */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                               */
/* Copyright (C) 2015                                                            */
/*                                                                               */
/*       eInfochips Limited                                                      */
/*                                                                               */
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
/* Written by:                                                                   */
/*      <Arifchand Mulani> <arifchand.mulani@einfochips.com>                     */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


/*
CombFirst4Moments_gpu.cu is calculation stastical moments
on GPU.
*/

#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cfloat> // DBL_MAX
#include "CombFirst4Moments.h"
#include "Moments.h"
#include "common.h"
#include <iostream>

#include <thrust/version.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/sequence.h>
#include <thrust/replace.h>
#include <thrust/transform_reduce.h>
#include <thrust/functional.h>

#define N_COMB_SAMPLES 20

void free_CombFirstMoments_Memory(double *h_pix_plane){
	if(h_pix_plane != NULL){
		free(h_pix_plane);
		h_pix_plane=NULL;
	}
}

/* structure holdind all moments */
template <typename T>
struct summary_stats_data
{
	T n;
	T min;
	T max;
	T mean;
	T M2;
	T M3;
	T M4;
	// initialize to the identity element
	void initialize()
	{
		n = mean = M2 = M3 = M4 = 0;
		min = std::numeric_limits<T>::max();
		max = std::numeric_limits<T>::min();
	}
	T variance() { return (n>2 ? sqrt (M2 / (n - 1)) : 0.0) ;}
	T variance_n() { return M2 / n; }
	T skewness() { return (n>3 ? (std::sqrt(n) * M3) / (std::pow(M2, (T) 1.5)) : 0.0); }
	T kurtosis() { return (n>4 ? (n * M4) / (M2 * M2) : 0.0); }
};

template <typename T>
struct summary_stats_unary_op
{
	__host__ __device__
	summary_stats_data<T> operator()(const T& x) const
	{
		summary_stats_data<T> result;
		result.n = 1;
		result.min = x;
		result.max = x;
		result.mean = x;
		result.M2 = 0;
		result.M3 = 0;
		result.M4 = 0;
		return result;
	}
};

/*  Momemts calculation are done here */
template <typename T>
struct summary_stats_binary_op
	: public thrust::binary_function<const summary_stats_data<T>&,
	const summary_stats_data<T>&,
	summary_stats_data<T> >
{
	__host__ __device__
	summary_stats_data<T> operator()(const summary_stats_data<T>& x, const summary_stats_data <T>& y) const
	{
		summary_stats_data<T> result;
		// precompute some common subexpressions
		T n = x.n + y.n;
		T n2 = n * n;
		T n3 = n2 * n;
		T delta = y.mean - x.mean;
		T delta2 = delta * delta;
		T delta3 = delta2 * delta;
		T delta4 = delta3 * delta;
		//Basic number of samples (n), min, and max
		result.n = n;
		result.min = thrust::min(x.min, y.min);
		result.max = thrust::max(x.max, y.max);
		result.mean = x.mean + delta * y.n / n;
		result.M2 = x.M2 + y.M2;
		result.M2 += delta2 * x.n * y.n / n;
		result.M3 = x.M3 + y.M3;
		result.M3 += delta3 * x.n * y.n * (x.n - y.n) / n2;
		result.M3 += (T) 3.0 * delta * (x.n * y.M2 - y.n * x.M2) / n;
		result.M4 = x.M4 + y.M4;
		result.M4 += delta4 * x.n * y.n * (x.n * x.n - x.n * y.n + y.n * y.n) / n3;
		result.M4 += (T) 6.0 * delta2 * (x.n * x.n * y.M2 + y.n * y.n * x.M2) / n2;
		result.M4 += (T) 4.0 * delta * (x.n * y.M3 - y.n * x.M3) / n;
		return result;
	}
};

//----------------------------------------------------------------------------------------------------------------------------------------
__global__
void	kernel_horizontal_comb (double*d_I,double*d_ImData,double*d_pix_plane,int ii,long height,long width,long m2,int*count)
{
	long index,rownum,colnum,i;
	int tid = threadIdx.x + blockDim.x * blockIdx.x;
	int threadcount = blockDim.x*gridDim.x;
	int temp;

	for (i=tid; i<(width*height); i=i+threadcount)
	{
		rownum = i/width;
		colnum = i%width;
		index = colnum + (width *rownum);
		if (fabs(d_I[index] + ii - m2) < 1)
		{
			temp=atomicAdd(count,1);
			d_pix_plane[temp]=d_ImData[ (width * ((long)d_I[index]-1)) +temp ];
		}
	}
}
//----------------------------------------------------------------------------------------------------------------------------------------
__global__
 void	kernel_vertical_comb(double*d_J,double*d_ImData,double*d_pix_plane,long ii,long height,long width,long n2,int*count)
{
	long index,rownum,colnum,i;
	long tid = threadIdx.x + blockDim.x * blockIdx.x;
	int threadcount = blockDim.x*gridDim.x;
	int temp;

	for (i=tid; i<(width*height); i=i+threadcount)
	{
		rownum = i/width;
		colnum = i%width;
		index = colnum + (width * rownum);
		if (fabs(d_J[index] + ii - n2) < 1)
		{
			temp=atomicAdd(count,1);
			d_pix_plane[temp]= d_ImData[((long)(d_J[index]-1)  + ((temp) * width) )];
		}
	}
}

__global__
void kernel_plus_45_degree(double *d_I,double *d_J1,double *d_ImData,double *d_pix_plane, int ii, long height, long width, int *count)
{
	long index,rownum,colnum,i;
	int tid = threadIdx.x + blockDim.x * blockIdx.x;
	int threadcount = blockDim.x*gridDim.x;
	int temp;

	for (i=tid; i<(width*height); i=i+threadcount)
	{
		rownum = i/width;
		colnum = i%width;
		index = colnum + (width *rownum);
		if (fabs(d_I[index] + ii - d_J1[index]) < 1)
		{
			temp=atomicAdd(count,1);
			d_pix_plane[temp]=d_ImData[ ( (long) d_I[index] * width) - (long)d_J1[index] ];
		}
	}
}

__global__
 void kernel_minus_45_degree(double *d_I, double *d_J, double *d_ImData, double *d_pix_plane, int ii, long height, long width, int *count)
{
	long index,rownum,colnum,i;
	int tid = threadIdx.x + blockDim.x * blockIdx.x;
	int threadcount = blockDim.x*gridDim.x;
	int temp;

	for (i=tid; i<(width*height); i=i+threadcount)
	{
		rownum = i/width;
		colnum = i%width;
		index = colnum + (width * rownum);
		if (fabs(d_I[index] + ii - d_J[index]) < 1)
		{
			temp=atomicAdd(count,1);
			d_pix_plane[temp] = d_ImData[ ((  (long)d_I[index]-1)*width)  + ((long)d_J[index]-1) ];
		}
	}
}

__global__
void kernel_copy_data (double *d_I,double *d_J,double *d_J1,long width,long height)
{
	long index,rownum,colnum,i;
	int tid = threadIdx.x + blockDim.x * blockIdx.x;
	int threadcount = blockDim.x*gridDim.x;

	for (i=tid; i<(width*height); i=i+threadcount)
	{
		rownum = i/width;
		colnum = i%width;
		index = colnum + (width *rownum);

		d_I[index]=(index%height)+1;
		d_J[index]=(index/height)+1;
		d_J1[((width * height) -1) - index]= d_J[index];
	}
}

//-----------------------------------------------------------------------------------------------------------------
int matr4moments_to_hist_gpu(double matr4moments[4][N_COMB_SAMPLES], std::vector<double> &vec, int vec_start){
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

/* This function is doning all allocations required and calling the CUDA kernels involved
   in statistical momentd
*/
extern "C"
int CombFirst4Moments2D_gpu(const ImageMatrix &Im, std::vector<double> &vec){
	double *d_I=NULL, *d_J=NULL, *d_ImData=NULL, *d_J1=NULL, *d_pix_plane=NULL;
	double *pix_data=NULL, *h_pix_plane=NULL;
	double z4[4]={0,0,0,0},z[4];
	double matr4moments[4][N_COMB_SAMPLES];
	long m,n,n2,m2;
	long a,ii;
	int matr4moments_index;
	int vec_count=0;
	readOnlyPixels pix_plane = Im.ReadablePixels();
	double step,m_pos;

	Moments4 tmpMoments;

#ifdef DISPLAY_TIME
	clock_t start_cpu,stop_cpu;
	double time_cpu;
	start_cpu=clock();
#endif

	#pragma unroll
	for (a = 0; a < 4; a++)    /* initialize */
		for (matr4moments_index = 0; matr4moments_index < N_COMB_SAMPLES; matr4moments_index++)
			matr4moments[a][matr4moments_index] = 0;

	m=Im.height;
	n=Im.width;
	long size=m*n;

	pix_data=(double *)pix_plane.data();
//----------------------------------------------------------- -----------------------------
	h_pix_plane= (double*)malloc(size * sizeof(double));
	if(h_pix_plane == NULL){
		printf("Error in malloc : at h_pix_plane\n");
		free_CombFirstMoments_Memory(h_pix_plane);
		return -1;
	}
//---------------------------- allocating memory on device ----------------------------------
	int count;
	int *d_count=NULL;
	if(checkCudaErrors(cudaMalloc ((void **)&d_count, sizeof(int))) != cudaSuccess){
		free_CombFirstMoments_Memory(h_pix_plane);
		cudaDeviceReset();
		return -1;
	}
//------------------------------------------------------------------------------------------
	if(checkCudaErrors(cudaMalloc ((void **)&d_pix_plane, size*sizeof(double))) !=cudaSuccess){
		free_CombFirstMoments_Memory(h_pix_plane);
		cudaDeviceReset();
		return -1;
	}
//-------------------------------------------------------------------------------------------
	if(checkCudaErrors(cudaMalloc ((void **)&d_I, size*sizeof(double))) != cudaSuccess){
		free_CombFirstMoments_Memory(h_pix_plane);
		cudaDeviceReset();
		return -1;
	}
	if(checkCudaErrors(cudaMalloc ((void **)&d_J, size*sizeof(double))) != cudaSuccess){
		free_CombFirstMoments_Memory(h_pix_plane);
		cudaDeviceReset();
		return -1;
	}
	if(checkCudaErrors(cudaMalloc ((void **)&d_J1, size*sizeof(double))) != cudaSuccess){
		free_CombFirstMoments_Memory(h_pix_plane);
		cudaDeviceReset();
		return -1;
	}
//-------------------------------------------------------------------------------------------------
	if(checkCudaErrors(cudaMalloc((void **)&d_ImData, size * sizeof(double))) != cudaSuccess){
		free_CombFirstMoments_Memory(h_pix_plane);
		cudaDeviceReset();
		return -1;
	}

	if(checkCudaErrors(cudaMemcpy(d_ImData, pix_data, size * sizeof(double), cudaMemcpyHostToDevice)) != cudaSuccess){
		free_CombFirstMoments_Memory(h_pix_plane);
		cudaDeviceReset();
		return -1;
	}
//--------------------------------------------------------------------------------------------------
	kernel_copy_data <<<NUMB,NUMT>>>(d_I,d_J,d_J1,n,m);
	if(checkCudaErrors(cudaDeviceSynchronize()) != cudaSuccess){
		free_CombFirstMoments_Memory(h_pix_plane);
		cudaDeviceReset();
		return -1;
	}

//--------------------------------------------------------------------------------------------------

	n2 = (int)(roundl(n/2));
	m2 = (int)(roundl(m/2));
//--------------------------------- operators required for the reduction ----------------------------
	summary_stats_unary_op<double> unary_op;
	summary_stats_binary_op<double> binary_op;
	summary_stats_data<double> init;
//---------------------------------------------------------------------------------------------------

	/* major diag -45 degrees */
	matr4moments_index = 0;
	step = (double)m / 10.0;
	if (step < 1.0) step = 1.0;
	for (m_pos = 1.0 - m; m_pos <= m; m_pos += step) {
		ii = round (m_pos);
		for (a = 0; a < 4; a++) matr4moments[a][matr4moments_index]=z4[a];

		count=0;/* clear count value*/
		if(checkCudaErrors(cudaMemcpy(d_count, &count, sizeof(int), cudaMemcpyHostToDevice)) != cudaSuccess){
			free_CombFirstMoments_Memory(h_pix_plane);
			cudaDeviceReset();
			return -1;
		}
		/*launching kernel*/
		kernel_minus_45_degree <<<NUMB,NUMT >>>(d_I, d_J, d_ImData, d_pix_plane, ii, m, n, d_count);
		if(checkCudaErrors(cudaDeviceSynchronize()) != cudaSuccess){
			free_CombFirstMoments_Memory(h_pix_plane);
			cudaDeviceReset();
			return -1;
	        }


		/*copy back results*/
		if(checkCudaErrors(cudaMemcpy( &count,d_count, sizeof(int), cudaMemcpyDeviceToHost)) != cudaSuccess ){
			free_CombFirstMoments_Memory(h_pix_plane);
			cudaDeviceReset();
			return -1;
		}

		if(checkCudaErrors(cudaMemcpy( h_pix_plane,d_pix_plane, size * sizeof(double), cudaMemcpyDeviceToHost)) !=cudaSuccess){
			free_CombFirstMoments_Memory(h_pix_plane);
			cudaDeviceReset();
			return -1;
		}

		/*use of thrust for reduction (calculating moments)*/
		thrust::device_vector<double> d_x(h_pix_plane, h_pix_plane+ count);
		init.initialize();
		summary_stats_data<double> result = thrust::transform_reduce(d_x.begin(), d_x.end(), unary_op, init, binary_op);

		z[0]=result.mean ;
		z[1]=result.variance();
		z[2]=result.skewness();
		z[3]=result.kurtosis();

		#pragma unroll
		for (a = 0; a < 4; a++) matr4moments[a][matr4moments_index] = z[a];
		matr4moments_index++;
	}
	vec_count=matr4moments_to_hist_gpu(matr4moments,vec,vec_count);

//------------------------------------------------------------------------------------------------------------------
	/* major diag +45 degrees */
	/* fliplr J */

	count=0;
	matr4moments_index = 0;
	step = (double)m / 10.0;
	if (step < 1.0) step = 1.0;
	for (m_pos = 1.0 - m; m_pos <= m; m_pos += step) {
		ii = round (m_pos);
		for (a = 0; a < 4; a++) matr4moments[a][matr4moments_index]=z4[a];

		count=0;/*clear count */
		if(checkCudaErrors(cudaMemcpy(d_count, &count, sizeof(int), cudaMemcpyHostToDevice)) != cudaSuccess){
			free_CombFirstMoments_Memory(h_pix_plane);
			cudaDeviceReset();
			return -1;
		}
		/*launching kernel*/
		kernel_plus_45_degree <<<NUMB,NUMT>>>(d_I, d_J1, d_ImData, d_pix_plane, ii, m, n, d_count);
		if(checkCudaErrors(cudaDeviceSynchronize()) != cudaSuccess){
			free_CombFirstMoments_Memory(h_pix_plane);
			cudaDeviceReset();
			return -1;
		}

		/*copy back results*/
		if(checkCudaErrors(cudaMemcpy( &count,d_count, sizeof(int), cudaMemcpyDeviceToHost)) != cudaSuccess){
			free_CombFirstMoments_Memory(h_pix_plane);
			cudaDeviceReset();
			return -1;
		}

		if(checkCudaErrors(cudaMemcpy( h_pix_plane,d_pix_plane, size * sizeof(double), cudaMemcpyDeviceToHost)) != cudaSuccess){
			free_CombFirstMoments_Memory(h_pix_plane);
			cudaDeviceReset();
			return -1;
		}

		/*use of thrust for calculating statistical moments*/
		thrust::device_vector<double> d_x(h_pix_plane, h_pix_plane+ count);
		init.initialize();
		summary_stats_data<double> result = thrust::transform_reduce(d_x.begin(), d_x.end(), unary_op, init, binary_op);

		/*store moments*/
		z[0]=result.mean ;
		z[1]=result.variance();
		z[2]=result.skewness();
		z[3]=result.kurtosis();

		for (a = 0; a < 4; a++) matr4moments[a][matr4moments_index] = z[a];
		matr4moments_index++;
	}
	vec_count=matr4moments_to_hist_gpu(matr4moments,vec,vec_count);
//--------------------------------------------------------------------------------------------------------------
	/* vertical comb*/
	count=0;
	matr4moments_index = 0;
	step = (double)m / 10.0;
	if (step < 1.0) step = 1.0;
	for (m_pos = 1.0 - m; m_pos <= m; m_pos += step) {
		ii = round (m_pos);
		for (a = 0; a < 4; a++) matr4moments[a][matr4moments_index]=z4[a];


		count=0;
		if(checkCudaErrors(cudaMemcpy(d_count, &count, sizeof(int), cudaMemcpyHostToDevice)) != cudaSuccess){
			free_CombFirstMoments_Memory(h_pix_plane);
			cudaDeviceReset();
			return -1;
		}

		kernel_vertical_comb <<<NUMB,NUMT>>>(d_J, d_ImData, d_pix_plane, ii, m, n, n2, d_count);
		if(checkCudaErrors(cudaDeviceSynchronize()) != cudaSuccess){
			free_CombFirstMoments_Memory(h_pix_plane);
			cudaDeviceReset();
			return -1;
	        }

		if(checkCudaErrors(cudaMemcpy( &count,d_count, sizeof(int), cudaMemcpyDeviceToHost)) != cudaSuccess){
			free_CombFirstMoments_Memory(h_pix_plane);
			cudaDeviceReset();
			return -1;
		}

		if(checkCudaErrors(cudaMemcpy( h_pix_plane,d_pix_plane, size * sizeof(double), cudaMemcpyDeviceToHost)) != cudaSuccess){
			free_CombFirstMoments_Memory(h_pix_plane);
			cudaDeviceReset();
			return -1;
		}

		thrust::device_vector<double> d_x(h_pix_plane, h_pix_plane+ count);
		init.initialize();
		summary_stats_data<double> result = thrust::transform_reduce(d_x.begin(), d_x.end(), unary_op, init, binary_op);

		/*store moments*/
		z[0]=result.mean ;
		z[1]=result.variance();
		z[2]=result.skewness();
		z[3]=result.kurtosis();

		#pragma unroll
		for (a = 0; a < 4; a++) matr4moments[a][matr4moments_index] = z[a];
		matr4moments_index++;
	}
	vec_count=matr4moments_to_hist_gpu(matr4moments,vec,vec_count);

//---------------------------------------------------------------------------------------------------------------
	/* horizontal comb */

	count=0;
	matr4moments_index = 0;
	step = (double)m / 10.0;
	if (step < 1.0) step = 1.0;
	for (m_pos = 1.0 - m; m_pos <= m; m_pos += step) {
		ii = round (m_pos);
		for (a = 0; a < 4; a++) matr4moments[a][matr4moments_index] = z4[a];

		count=0;/*clear count*/
		if(checkCudaErrors(cudaMemcpy(d_count, &count, sizeof(int), cudaMemcpyHostToDevice)) != cudaSuccess){
			free_CombFirstMoments_Memory(h_pix_plane);
			cudaDeviceReset();
			return -1;
		}

		/*launching kernel*/
		kernel_horizontal_comb <<<NUMB,NUMT>>>(d_I, d_ImData, d_pix_plane, ii, m, n, m2, d_count);
		if(checkCudaErrors(cudaDeviceSynchronize()) != cudaSuccess){
			free_CombFirstMoments_Memory(h_pix_plane);
			cudaDeviceReset();
			return -1;
		}

		/*copy back results*/
		if(checkCudaErrors(cudaMemcpy( &count,d_count, sizeof(int), cudaMemcpyDeviceToHost)) != cudaSuccess){
			free_CombFirstMoments_Memory(h_pix_plane);
			cudaDeviceReset();
			return -1;
		}

		if(checkCudaErrors(cudaMemcpy( h_pix_plane,d_pix_plane, size * sizeof(double), cudaMemcpyDeviceToHost)) != cudaSuccess){
			free_CombFirstMoments_Memory(h_pix_plane);
			cudaDeviceReset();
			return -1;
		}

		/*use of thrust for calculating moments*/
		thrust::device_vector<double> d_x(h_pix_plane, h_pix_plane+ count);
		init.initialize();
		summary_stats_data<double> result = thrust::transform_reduce(d_x.begin(), d_x.end(), unary_op, init, binary_op);

		/*store statistical moments*/
		z[0]=result.mean ;
		z[1]=result.variance();
		z[2]=result.skewness();
		z[3]=result.kurtosis();

		#pragma unroll
		for (a = 0; a < 4; a++) matr4moments[a][matr4moments_index] = z[a];
		matr4moments_index++;
	}
	vec_count=matr4moments_to_hist_gpu(matr4moments,vec,vec_count);

#ifdef DISPLAY_TIME
        stop_cpu=clock();
        time_cpu=(double)(stop_cpu-start_cpu)/CLOCKS_PER_SEC*1000;
        std::cout<<"Time to Moments (GPU):"<<time_cpu<<std::endl;
#endif

	/*Free memory*/
	free_CombFirstMoments_Memory(h_pix_plane);
	if(checkCudaErrors(cudaFree(d_I)) != cudaSuccess){
		cudaDeviceReset();
		return -1;
        }
	if(checkCudaErrors(cudaFree(d_J)) != cudaSuccess){
		cudaDeviceReset();
		return -1;
        }
	if(checkCudaErrors(cudaFree(d_J1)) != cudaSuccess){
		cudaDeviceReset();
		return -1;
        }
	if(checkCudaErrors(cudaFree(d_count)) != cudaSuccess){
		cudaDeviceReset();
		return -1;
        }
	if(checkCudaErrors(cudaFree(d_ImData)) != cudaSuccess){
		cudaDeviceReset();
		return -1;
        }
	if(checkCudaErrors(cudaFree(d_pix_plane)) != cudaSuccess){
		cudaDeviceReset();
		return -1;
        }

	return(vec_count);
}

void vd_Comb4Moments_gpu(std::vector<double> &in) {

	std::vector<double> temp = in; //deep copy
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

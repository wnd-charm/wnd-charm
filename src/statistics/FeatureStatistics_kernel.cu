/*FeatureStatics_kernel.cu                                                       */
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
/* Written by:  Lalit Chandivade           				  	 */
/* Email: lalit.chandivade@einfochips.com                                        */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <cuda.h>
#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/scan.h>
#include <thrust/fill.h>
#include <thrust/copy.h>

#include<common.h>
#include "cmatrix.h"

const dim3 numThreads(32,32);

void free_FeatureCentroid_gpu()
{
	cudaDeviceReset();
}

__global__
void set_mass_kernel(const double *d_ImData, unsigned long *d_x_mass, unsigned long *d_y_mass,
		     unsigned long *d_mass, double object_index, unsigned int height, unsigned int width)
{
	unsigned long row, col;

	row = threadIdx.y + (blockIdx.y * blockDim.y);
	col = threadIdx.x + (blockIdx.x * blockDim.x);

	if ( row < height && col < width) {
		if (d_ImData[col + row * width] == object_index) {
			atomicAdd((unsigned long long int *)d_x_mass, (unsigned long long int)col+1);
			atomicAdd((unsigned long long int *)d_y_mass, (unsigned long long int)row+1);
			atomicAdd((unsigned long long int *)d_mass, (unsigned long long int)1);
		}
	}
}

extern "C"
int FeatureCentroid_gpu(const ImageMatrix &Im, double *sum_areas, double *sum_dists, unsigned long *object_areas, double *centroid_dists,
				unsigned long *count, double *centroid_x, double *centroid_y) {
	unsigned int w = Im.width, h = Im.height;
	unsigned long x_mass=0,y_mass=0,mass=0;
	readOnlyPixels pix_plane = Im.ReadablePixels();
	const double *ImData = pix_plane.data();
	double *d_ImData;
	unsigned long *d_x_mass;
	unsigned long *d_y_mass;
	unsigned long *d_mass;
	unsigned long object_index;
	double x_centroid,y_centroid;
	dim3 numBlocks;

	if(checkCudaErrors(cudaMalloc((void **)&d_ImData, sizeof(double) * h * w)) != cudaSuccess)
	{
		free_FeatureCentroid_gpu();
		return -1;
	}

	if(checkCudaErrors(cudaMalloc((void **)&d_x_mass, sizeof(unsigned long))) != cudaSuccess)
	{
		free_FeatureCentroid_gpu();
		return -1;
	}

	if(checkCudaErrors(cudaMalloc((void **)&d_y_mass, sizeof(unsigned long))) != cudaSuccess)
	{
		free_FeatureCentroid_gpu();
		return -1;
	}

	if(checkCudaErrors(cudaMalloc((void **)&d_mass, sizeof(unsigned long))) != cudaSuccess)
	{
		free_FeatureCentroid_gpu();
		return -1;
	}

	if(checkCudaErrors(cudaMemcpy(d_ImData, ImData, sizeof(double) * h * w, cudaMemcpyHostToDevice)) != cudaSuccess)
	{
		free_FeatureCentroid_gpu();
		return -1;
	}
        numBlocks.x = (int) ceil ((float) w/numThreads.x);
        numBlocks.y =  (int) ceil ((float) h/numThreads.y);

	for (object_index = 0; object_index < *count; object_index++) {


		if(checkCudaErrors(cudaMemset(d_x_mass, 0, sizeof(unsigned long))) != cudaSuccess)
		{
			free_FeatureCentroid_gpu();
			return -1;
		}

		if(checkCudaErrors(cudaMemset(d_y_mass, 0, sizeof(unsigned long))) != cudaSuccess)
		{
			free_FeatureCentroid_gpu();
			return -1;
		}

		if(checkCudaErrors(cudaMemset(d_mass, 0, sizeof(unsigned long))) != cudaSuccess)
		{
			free_FeatureCentroid_gpu();
			return -1;
		}

		set_mass_kernel<<<numBlocks, numThreads>>>(d_ImData, d_x_mass, d_y_mass, d_mass, object_index+1, h, w);
		checkCudaErrors(cudaDeviceSynchronize());

		if(checkCudaErrors(cudaMemcpy(&x_mass, d_x_mass, sizeof(unsigned long), cudaMemcpyDeviceToHost)) != cudaSuccess)
		{
			free_FeatureCentroid_gpu();
			return -1;
		}

		if(checkCudaErrors(cudaMemcpy(&y_mass, d_y_mass, sizeof(unsigned long), cudaMemcpyDeviceToHost)) != cudaSuccess)
		{
			free_FeatureCentroid_gpu();
			return -1;
		}

		if(checkCudaErrors(cudaMemcpy(&mass, d_mass, sizeof(unsigned long), cudaMemcpyDeviceToHost)) != cudaSuccess)
		{
			free_FeatureCentroid_gpu();
			return -1;
		}

		if (x_centroid) x_centroid=(double)x_mass/(double)mass;
		if (y_centroid) y_centroid=(double)y_mass/(double)mass;

		object_areas[object_index] = mass;
		centroid_dists[object_index] = sqrt(pow(x_centroid-(*centroid_x),2) + pow(y_centroid - (*centroid_y),2));
		*sum_areas += object_areas[object_index];
		*sum_dists += centroid_dists[object_index];
	}

	//de-allocation of device memory
	if(checkCudaErrors(cudaFree(d_ImData))!= cudaSuccess)
	{
		free_FeatureCentroid_gpu();
		return -1;
	}

	if(checkCudaErrors(cudaFree(d_x_mass))!= cudaSuccess)
	{
		free_FeatureCentroid_gpu();
		return -1;
	}


	if(checkCudaErrors(cudaFree(d_y_mass))!= cudaSuccess)
	{
		free_FeatureCentroid_gpu();
		return -1;
	}

	if(checkCudaErrors(cudaFree(d_mass))!= cudaSuccess)
	{
		free_FeatureCentroid_gpu();
		return -1;
	}

	return 0;
}

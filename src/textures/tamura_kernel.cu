/* tamura_kernel.cu                                                         */
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
/*      Rutika Ubale                                                             */
/*      rutika.ubale@einfochips.com                                              */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include<cuda_runtime.h>
#include<stdio.h>
#include<common.h>
#include "math.h"
#include <sys/time.h>


__device__ double efficientLocalMean_dev (const long x,const long y,const long k, double * input_img, int rowsize, int colsize) {
        long k2 = k/2;

        long dimx = rowsize;
        long dimy = colsize;

        //wanting average over area: (y-k2,x-k2) ... (y+k2-1, x+k2-1)
        long starty = y-k2;
        long startx = x-k2;
        long stopy = y+k2-1;
        long stopx = x+k2-1;

        if (starty < 0) starty = 0;
        if (startx < 0) startx = 0;
        if (stopx > dimx-1) stopx = dimx-1;
        if (stopy > dimy-1) stopy = dimy-1;

        double unten, links, oben, obenlinks;

        if (startx-1 < 0) links = 0;
        else links = *(input_img+(stopy * dimx + startx-1));

        if (starty-1 < 0) oben = 0;
        else oben = *(input_img+((stopy-1) * dimx + startx));

        if ((starty-1 < 0) || (startx-1 <0)) obenlinks = 0;
        else obenlinks = *(input_img+((stopy-1) * dimx + startx-1));

        unten = *(input_img+(stopy * dimx + startx));

        long counter = (stopy-starty+1)*(stopx-startx+1);
        return (unten-links-oben+obenlinks)/counter;
}

__global__ void process_coarseness_ak_pix(double * output_ak,double * input_img,int colsize, int rowsize,long lenOf_ak)
{
    int index;
    int y  = threadIdx.x + blockIdx.x * blockDim.x;
    int x = threadIdx.y + blockIdx.y * blockDim.y;
    if(y < (colsize) && x < (rowsize))
    {
        index = y * rowsize + x ;
	output_ak[index] = efficientLocalMean_dev(x,y,lenOf_ak,input_img,rowsize,colsize);
    }
}

__global__ void process_coarseness_ek_pix(double * output_ak, double *output_ekh, double *output_ekv,int colsize, int rowsize,long lenOf_ek)
{
    int y  = threadIdx.x + blockIdx.x * blockDim.x;
    int x = threadIdx.y + blockIdx.y * blockDim.y;
    double input1,input2;
    int posx1 = x+lenOf_ek;
    int posx2 = x-lenOf_ek;
    int posy1 = y+lenOf_ek;
    int posy2 = y-lenOf_ek;
    if(y < (colsize) && x < (rowsize))
    {
	 if(posx1 < (int)rowsize && posx2 >= 0)
	 {
	     input1 = output_ak[y * rowsize + posx1];
	     input2 = output_ak[y * rowsize + posx2];
             output_ekh[y*rowsize+x] = fabs(input1 - input2);
	 }
	 else output_ekh[y*rowsize+x] = 0;

	 if(posy1 < (int)colsize && posy2 >= 0)
	 {
             input1 = output_ak[posy1 * rowsize + x];
	     input2 = output_ak[posy2 * rowsize + x];
	     output_ekv[y*rowsize+x] = fabs(input1 - input2);
	 }
	 else output_ekv[y*rowsize+x] = 0;
    }
}

void tamura_free()
{
 /** It Destroy all allocations and reset all state on the current device in the current process*/
    cudaDeviceReset();
}

int GPU_process_pix(double * laufendeSumme_hostPtr,double *Ak_pix_plane_host[],double *Ekh_pix_plane_host[],double *Ekv_pix_plane_host[],int yDim,int xDim,int sizeK)
{

    double *Ak_pix_plane_dev;
    double *Ekh_pix_plane_dev;
    double *Ekv_pix_plane_dev;
    double *laufende_devPtr;
    long lenOf_ak = 1,lenOf_ek = 1 ;
    cudaStream_t stream[sizeK];
    int i = 0;

    /** Threads and block configuration **/
    dim3 thread_dim (16 ,16 ,1);

    int y = ((yDim%thread_dim.x))? (yDim/thread_dim.x +1) : (yDim/thread_dim.x);
    int x = ((xDim%thread_dim.y))? (xDim/thread_dim.y +1) : (xDim/thread_dim.y);
    dim3 block_dim (y ,x ,1);

    /** Device memory alloaction **/
    if(checkCudaErrors(cudaMalloc( (void**)&laufende_devPtr ,yDim *  xDim *sizeof(double)))!= cudaSuccess)
    {
       tamura_free();
       return -1;
    }
    if(checkCudaErrors(cudaMalloc( (void**)&Ekh_pix_plane_dev ,yDim *  xDim * sizeK * sizeof(double)))!= cudaSuccess)
    {
       tamura_free();
       return -1;
    }
    if(checkCudaErrors(cudaMalloc( (void**)&Ekv_pix_plane_dev ,yDim *  xDim * sizeK * sizeof(double)))!= cudaSuccess)
    {
       tamura_free();
       return -1;
    }
    if(checkCudaErrors(cudaMalloc( (void**)&Ak_pix_plane_dev ,yDim *  xDim * sizeK * sizeof(double)))!= cudaSuccess)
    {
       tamura_free();
       return -1;
    }
    if(checkCudaErrors(cudaMemcpy(laufende_devPtr, laufendeSumme_hostPtr , yDim *  xDim * sizeof(double),cudaMemcpyHostToDevice))!= cudaSuccess)
    {
       tamura_free();
       return -1;
    }

    /** stream creation **/
    for (i = 0 ;i < sizeK; i++)
    {
        if(checkCudaErrors(cudaStreamCreate(&stream[i]))!= cudaSuccess)
        {
		tamura_free();
		return -1;
        }
    }

    /** Kernel launching for each streams  **/
    for (i = 0 ;i < sizeK; i++)
    {
           lenOf_ak *= 2;
	    process_coarseness_ak_pix<<<block_dim,thread_dim,0,stream[i]>>>(Ak_pix_plane_dev + (i * yDim* xDim),
									laufende_devPtr,yDim ,xDim,lenOf_ak);


	    process_coarseness_ek_pix<<<block_dim,thread_dim,0,stream[i]>>>(Ak_pix_plane_dev + (i * yDim* xDim),
                                                                        Ekh_pix_plane_dev + (i * yDim* xDim),
									Ekv_pix_plane_dev + (i * yDim* xDim),
									yDim ,xDim,lenOf_ek);
        lenOf_ek *= 2;
    }

    for (i = 0 ;i < sizeK; i++)
    {
	    if(checkCudaErrors(cudaMemcpyAsync(Ekh_pix_plane_host[i] , Ekh_pix_plane_dev + (i * yDim* xDim), yDim *  xDim * sizeof(double),cudaMemcpyDeviceToHost,stream[i]))!= cudaSuccess)
            {
		    tamura_free();
		    return -1;
            }
	    if(checkCudaErrors(cudaMemcpyAsync(Ekv_pix_plane_host[i] , Ekv_pix_plane_dev + (i * yDim* xDim), yDim *  xDim * sizeof(double),cudaMemcpyDeviceToHost,stream[i]))!= cudaSuccess)
            {
		    tamura_free();
		    return -1;
            }
    }
    /** wait for streams to complete their execution  **/
    for (i = 0 ;i < sizeK; i++)
    {
        cudaStreamSynchronize(stream[i]);
        if(checkCudaErrors(cudaStreamDestroy(stream[i]))!= cudaSuccess)
        {
		tamura_free();
		return -1;
        }

    }
    /** Free the devicde memory  **/
    if(checkCudaErrors(cudaFree(laufende_devPtr))!= cudaSuccess)
    {
	    tamura_free();
	    return -1;
    }
    if(checkCudaErrors(cudaFree(Ekh_pix_plane_dev))!= cudaSuccess)
    {
	    tamura_free();
	    return -1;
    }
    if(checkCudaErrors(cudaFree(Ekv_pix_plane_dev))!= cudaSuccess)
    {
	    tamura_free();
	    return -1;
    }
    if(checkCudaErrors(cudaFree(Ak_pix_plane_dev))!= cudaSuccess)
    {
	    tamura_free();
	    return -1;
    }
    return 0;
}

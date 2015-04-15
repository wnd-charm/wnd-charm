 /* common.cu                                                                     */
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

#ifndef _COMMON_CUDA_
#define _COMMON_CUDA_

#include<stdio.h>
#include <iostream>
#include<cuda.h>
#include <cublas_v2.h>
#include <cufft.h>

const char * _cudaGetErrorEnum(cudaError_t error);
const char * _cudaGetErrorEnum(cublasStatus_t error);
const char * _cudaGetErrorEnum(cufftResult error);
template< typename T >
int check(T result, char const *const func, const char *const file, int const line)
{
    if (result)
    {
        printf("CUDA error at %s:%d code=%d(%s) \"%s\" \n",
                file, line, static_cast<unsigned int>(result), _cudaGetErrorEnum(result), func);

        // Make sure we call CUDA Device Reset before exiting
        // exit(EXIT_FAILURE);
    }
          return result;
}
void init();

#define checkCudaErrors(val)           check ( (val), #val, __FILE__, __LINE__ )

#define cuda_initialize()              init()

#endif

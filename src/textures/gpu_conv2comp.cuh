/*     gpu_conv2comp.cuh                                                         */
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
#ifndef _CONVCOMP_CUDA_
#define _CONVCOMP_CUDA_

#include <iostream>
#include <stdio.h>
#include <cuda.h>

#define NUMT 1024
#define NUMB 1024



__device__ void compute_conv(int, int, double2 *, double *, double2 *, int *, int *, int, int, int, int, int, int);
__global__ void kernel_conv(double2 *, double *, double2 *, int *, int *, int, int, int, int, int, int);
void generate_vect(int *, int, int);
int min_dim(int, int);
void free_func_conv2comp(int *, int *);
extern "C" void gpu_conv2comp(double *, double *, double *, int, int, int, int);
#endif

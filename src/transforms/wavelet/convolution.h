/*
Copyright (c) 2006 Filip Wasilewski <filipwasilewski@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE


 * Originally from PyWavelets 
 * http://www.pybytes.com/pywavelets/
 * 
 * pretty much it just does the output size calculations
 * and the convultion step based on a specific kind of padding used
 * ZERO PADDING
 */

#ifndef _CONVOLUTION_H_
#define _CONVOLUTION_H_

#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>

enum MODE {
	   MODE_INVALID = -1,
       MODE_ZEROPAD = 0, // default
       MODE_SYMMETRIC,
	   MODE_ASYMMETRIC,
	   MODE_CONSTANT_EDGE,
       MODE_SMOOTH,
	   MODE_PERIODIC,
	   MODE_PERIODIZATION,
	   MODE_MAX
};

#define wtmalloc(size_t) malloc(size_t)
#define wtcalloc(len, size_t) calloc(len, size_t)
#define wtfree(ptr) free(ptr)


int dwt_buffer_length(int input_len, int filter_len, int mode);
int reconstruction_buffer_length(int coeffs_len, int filter_len);
int idwt_buffer_length(int coeffs_len, int filter_len, int mode);
int swt_buffer_length(int input_len);


int dwt_max_level(int input_len, int filter_len);


///////////////////////////////////////////////////////////////////////////////
//
// performs convolution of input with filter and downsamples by taking every
// step-th element from result.
//
// input	- input vector
// N		- input vector length
// filter	- filter vector
// F		- filter vector length
// output	- output vector
// step		- downsample step
// mode		- 

#define DTYPE double

int downsampling_convolution(const DTYPE* input, const int N, const double* filter, const int F, double* output, const int step, const int mode);
int allocating_downsampling_convolution(const DTYPE* input, const int N, const double* filter, const int F, double* output, const int step, const int mode);

#define convolution(data, data_len, filter, filter_len, output) downsampling_convolution(data, data_len, filter, filter_len, output, 1, MODE_ZEROPAD);


///////////////////////////////////////////////////////////////////////////////
//
// upsamples input signal by inserting zeros and convolves with filter.
// input: i0 i1 i2 i3 -> (upsampling) -> i0 0 i1 0 i2 0 i3 (0)
//
// input	- input vector
// N		- input vector length
// filter	- filter vector
// F		- filter vector length
// output	- output vector
// mode		- 


int upsampling_convolution_full(const DTYPE* input, const int N, const double* filter, const int F, double* output, const int O);
int upsampling_convolution_valid_sf(const DTYPE* input, const int N, const double* filter, const int F, double* output, const int O, const int mode);

//TODO
//int extended_filter_convolution(const DTYPE* input, const int N, const double* filter, const int F, double* output, int step, int mode);

#endif

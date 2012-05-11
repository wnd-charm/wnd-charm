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
 
#ifndef _WT_H_
#define _WT_H_

//#pragma inline_depth(2)

#include <memory.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

//#include "wavelets.h"
#include "convolution.h"
#include "Wavelet.h"

int d_dec_a(double input[], int input_len, Wavelet* wavelet, double output[], int output_len, int mode);
int d_dec_d(double input[], int input_len, Wavelet* wavelet, double output[], int output_len, int mode);

int d_rec_a(double coeffs_a[], int coeffs_len, Wavelet* wavelet, double output[], int output_len);
int d_rec_d(double coeffs_d[], int coeffs_len, Wavelet* wavelet, double output[],int  output_len);

int d_idwt(double coeffs_a[], int coeffs_a_len, double coeffs_d[], int coeffs_d_len, Wavelet* wavelet, double output[], int output_len, int mode, int correct_size);

#define max_swt_level(input_len) (int)floor(log(input_len)/log(2))

#endif

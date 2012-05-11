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


#include "Wavelet.h"
#include "wt.h"
//#include <iostream>
//using namespace std;

int d_dec_a(double input[], int input_len, Wavelet* wavelet, double output[], int output_len, int mode){
	int dwt_buf_len = dwt_buffer_length(input_len, wavelet->dec_len, mode);
	//printf("output len: %d, dwt buf len: %d\n",output_len,dwt_buf_len);
	if(output_len != dwt_buf_len){
//		cout << "There was a difference in what the data length should be. " <<
//		" This is unrecoverable and should not have happened" << endl;
		exit(0);
	}
	return downsampling_convolution(input, input_len, wavelet->analysisLow->coeff, wavelet->dec_len, output, 2, mode);
}

int d_dec_d(double input[], int input_len, Wavelet* wavelet, double output[], int output_len, int mode){
	if(output_len != dwt_buffer_length(input_len, wavelet->dec_len, mode))
		return -1;
	return downsampling_convolution(input, input_len, wavelet->analysisHigh->coeff, wavelet->dec_len, output, 2, mode);
}

int d_rec_a(double coeffs_a[], int coeffs_len, Wavelet* wavelet, double output[], int output_len){
	if(output_len != reconstruction_buffer_length(coeffs_len, wavelet->rec_len))
		return -1;
	
	return upsampling_convolution_full(coeffs_a, coeffs_len, wavelet->synthesisLow->coeff, wavelet->rec_len, output, output_len);
}

int d_rec_d(double coeffs_d[], int coeffs_len, Wavelet* wavelet, double output[],int  output_len){
	if(output_len != reconstruction_buffer_length(coeffs_len, wavelet->rec_len))
		return -1;
	
	return upsampling_convolution_full(coeffs_d, coeffs_len, wavelet->synthesisHigh->coeff, wavelet->rec_len, output, output_len);
}


int d_idwt(double coeffs_a[], int coeffs_a_len, double coeffs_d[], int coeffs_d_len, Wavelet* wavelet, double output[],
		 int output_len, int mode, int correct_size)
{
	int input_len;
	
	if(coeffs_a != NULL && coeffs_d != NULL){
		if(correct_size){
			if( (coeffs_a_len > coeffs_d_len ? coeffs_a_len - coeffs_d_len : coeffs_d_len-coeffs_a_len) > 1) 
			{ 	
//				goto errorIDWT;
                                exit(0);
			}
			input_len = coeffs_a_len>coeffs_d_len? coeffs_d_len : coeffs_a_len; // min
		} else {
			if(coeffs_a_len != coeffs_d_len) {
//				goto errorIDWT;
                                exit(0);
			}
			input_len = coeffs_a_len;
		}
	} else if(coeffs_a != NULL){
		input_len  = coeffs_a_len;
	} else if (coeffs_d != NULL){
		input_len = coeffs_d_len;
	} else {
//		goto errorIDWT;
                exit(0);
	}
	
	int supposed_len = idwt_buffer_length(input_len, wavelet->rec_len, mode);
		
	if(output_len != supposed_len) {
		goto errorIDWT;
	}
	memset(output, 0, output_len * sizeof(double));

	if(coeffs_a){
		if(upsampling_convolution_valid_sf(coeffs_a, input_len, wavelet->synthesisLow->coeff, wavelet->rec_len, output, output_len, mode) < 0) {
			goto errorIDWT;
		}
	}
	// +=
	if(coeffs_d){
		if(upsampling_convolution_valid_sf(coeffs_d, input_len, wavelet->synthesisHigh->coeff, wavelet->rec_len, output, output_len, mode) < 0) {
			goto errorIDWT;
		}
	}	
	
	return 0;

errorIDWT:
//	cout << "Invalid input size for output array. This is  " <<
//	    	    "unrecoverable and should not have happened " << endl;
    exit(0);

}

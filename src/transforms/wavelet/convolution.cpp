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
 
 
#include "convolution.h"
//#include <iostream>
//using namespace std;

int dwt_buffer_length(int input_len, int filter_len, int mode){
	if(input_len < 1 || filter_len < 1)
		return 0;
	switch(mode){
			case MODE_PERIODIZATION:
				return (int) ceil(input_len / 2.);
				break;
			default:
				return (int) floor((input_len + filter_len - 1) / 2.);
				break;
	}
}

int reconstruction_buffer_length(int coeffs_len, int filter_len){
	if(coeffs_len < 1 || filter_len < 1)
		return 0;
	return 2*coeffs_len+filter_len-2;
}

int idwt_buffer_length(int coeffs_len, int filter_len, int mode){
	if(coeffs_len < 0 || filter_len < 0) {
//		 cout << "Invalid coefficient arrays length for specified " <<
//	        " wavelet. Probably reconstructing not a result of dwt or " <<
//	        " the wavelet or mode is mistaken. " << endl;
	        exit(0);
	}

	switch(mode){
			case MODE_PERIODIZATION:
				return 2*coeffs_len;
				break;
			default:
				return 2*coeffs_len-filter_len+2;
	}
}

int dwt_max_level(int input_len, int filter_len){
	if(input_len < 1 || filter_len < 2)
		return 0;
	return (int)floor(log((double)input_len/(filter_len-1))/log(2.));
}


// -> swt - todo
int extended_filter_convolution(const DTYPE* input, const int N, const double* filter, const int F, double* output, int step, int mode)
{
	return -1;
}

int downsampling_convolution(const DTYPE* input, const int N, const double* filter, const int F, double* output, const int step, const int mode)
{
    int i, j, k, F_2, corr;
	int start;
    double sum, tmp;

    double* ptr_w = output;

	i = start = step-1; // first element taken from input is input[step-1]

    if(F <= N){

		///////////////////////////////////////////////////////////////////////
		// 0 - F-1	- sliding in filter

		// border distortion
		switch(mode)
		{
			case MODE_SYMMETRIC:
				for(i=start; i < F; i+=step){
					sum = 0;
					for(j = 0; j <= i; ++j)
						sum += filter[j]*input[i-j];
					k = i+1;
					for(j = i+1; j < F; ++j)
						sum += filter[j] * input[j-k];
					*(ptr_w++) = sum;
				}
				break;

			case MODE_ASYMMETRIC:
				for(i=start; i < F; i+=step){
					sum = 0;
					for(j = 0; j <= i; ++j)
						sum += filter[j]*input[i-j];
					k = i+1;
					for(j = i+1; j < F; ++j)
						sum += filter[j] * (input[0] - input[j-k]); // -=
					*(ptr_w++) = sum;
				}
				break;

			case MODE_CONSTANT_EDGE:
				for(i=start; i < F; i+=step){
					sum = 0;
					for(j = 0; j <= i; ++j)
						sum += filter[j]*input[i-j];
				
					k = i+1;
					for(j = i+1; j < F; ++j)
						sum += filter[j] * input[0];

					*(ptr_w++) = sum;
				}  
				break;

			case MODE_SMOOTH:
				for(i=start; i < F; i+=step){
					sum = 0;
					for(j = 0; j <= i; ++j)
						sum += filter[j]*input[i-j];
					tmp = input[0]-input[1];
					for(j = i+1; j < F; ++j){
						sum += filter[j] * (input[0] + tmp * (j-i)); 
					}
					*(ptr_w++) = sum;
				}  
				break;

			case MODE_PERIODIC:
				for(i=start; i < F; i+=step){
					sum = 0;
					for(j = 0; j <= i; ++j)
						sum += filter[j]*input[i-j];

					k = N+i;
					for(j = i+1; j < F; ++j)
						sum += filter[j] * input[k-j];

					*(ptr_w++) = sum;
				}  
				break;

			case MODE_PERIODIZATION:
				//extending by (F-2)/2 elements
				start = F/2;

				F_2 = F/2;
				corr = 0;
				for(i=start; i < F; i+=step){
					sum = 0;
					for(j = 0; j < i+1-corr; ++j)	// overlapping
						sum += filter[j]*input[i-j-corr];

					if(N%2){
						if(F-j){ // if something to extend
							sum += filter[j] * input[N-1];
							if(F-j){
								for(k = 2-corr; k <= F-j; ++k)
									sum += filter[j-1+k] * input[N-k+1];
							}
						}
					} else { // extra element from input   -> i0 i1 i2 [i2]
						for(k = 1; k <= F-j; ++k)
							sum += filter[j-1+k] * input[N-k];
					}
					*(ptr_w++) = sum;
				}
				break;
				
			case MODE_ZEROPAD:
			default:
				for(i=start; i < F; i+=step){
					sum = 0;
					for(j = 0; j <= i; ++j)
						sum += filter[j]*input[i-j];
					*(ptr_w++) = sum;
				}  
				
				break;
		}

		///////////////////////////////////////////////////////////////////////
        // F - N-1		- filter in input range
        for(; i < N; i+=step){					// input elements
            sum = 0;
            for(j = 0; j < F; ++j)				
                sum += input[i-j]*filter[j];
            *(ptr_w++) = sum;
        }  

		///////////////////////////////////////////////////////////////////////
		// N - N+F-1	- sliding out filter
		switch(mode)
		{
			case MODE_SYMMETRIC:
				for(; i < N+F-1; i += step){	// input elements
					sum = 0;
					k = i-N+1;					// 1, 2, 3 // overlapped elements
					for(j = k; j < F; ++j)					//TODO: j < F-_offset
						sum += filter[j]*input[i-j];

					for(j = 0; j < k; ++j)		// out of boundary		//TODO: j = _offset
						sum += filter[j]*input[N-k+j]; // j-i-1			0*(N-1), 0*(N-2) 1*(N-1)

					*(ptr_w++) = sum;
				}  
				break;

			case MODE_ASYMMETRIC:
				for(; i < N+F-1; i += step){	// input elements
					sum = 0;
					k = i-N+1;
					for(j = k; j < F; ++j)		// overlapped elements
						sum += filter[j]*input[i-j];

					for(j = 0; j < k; ++j)		// out of boundary
						sum += filter[j]*(input[N-1]-input[N-k-1+j]); // -=	j-i-1 
					*(ptr_w++) = sum;
				}  
				break;

			case MODE_CONSTANT_EDGE:
				for(; i < N+F-1; i += step){	// input elements
					sum = 0;
					k = i-N+1;
					for(j = k; j < F; ++j)		// overlapped elements
						sum += filter[j]*input[i-j];

					for(j = 0; j < k; ++j)		// out of boundary (filter elements [0, k-1])
						sum += filter[j]*input[N-1]; // input[N-1] = const

					*(ptr_w++) = sum;
				}  
				break;

			case MODE_SMOOTH:
				for(; i < N+F-1; i += step){	// input elements
					sum = 0;
					k = i-N+1; // 1, 2, 3, ...
					for(j = k; j < F; ++j)		// overlapped elements
						sum += filter[j]*input[i-j];

					tmp = input[N-1]-input[N-2];
					for(j = 0; j < k; ++j)		// out of boundary (filter elements [0, k-1])
						sum += filter[j] * (input[N-1] + tmp * (k-j));
					*(ptr_w++) = sum;
				}
				break;

			case MODE_PERIODIC:
				for(; i < N+F-1; i += step){	// input elements
					sum = 0;
					k = i-N+1;
					for(j = k; j < F; ++j)		// overlapped elements
						sum += filter[j]*input[i-j];
					for(j = 0; j < k; ++j)		// out of boundary (filter elements [0, k-1])
						sum += filter[j]*input[k-1-j];
					*(ptr_w++) = sum;
				}  
				break;

			case MODE_PERIODIZATION:
				
				for(; i < N-step + (F/2)+1 + N%2; i += step){	// input elements
					sum = 0;
					k = i-N+1;
					for(j = k; j < F; ++j)				// overlapped elements
						sum += filter[j]*input[i-j];
					
					if(N%2 == 0){
						for(j = 0; j < k; ++j){			// out of boundary (filter elements [0, k-1])
							sum += filter[j]*input[k-1-j];
						}
					} else {							// repeating extra element -> i0 i1 i2 [i2]
						for(j = 0; j < k-1; ++j)		// out of boundary (filter elements [0, k-1])
							sum += filter[j]*input[k-2-j];
						sum += filter[k-1] * input[N-1];
					}
					*(ptr_w++) = sum;
				}
				break;

			case MODE_ZEROPAD:
			default:
				for(; i < N+F-1; i += step){
					sum = 0;
					for(j = i-(N-1); j < F; ++j)
						sum += input[i-j]*filter[j];
					*(ptr_w++) = sum;
				}  
				break;
		}		


		return 0;

	} else {
		return allocating_downsampling_convolution(input, N, filter, F, output, step, mode);
    }
}

///////////////////////////////////////////////////////////////////////////////
//
// like downsampling_convolution, but with memory allocation
//

int allocating_downsampling_convolution(const DTYPE* input, const int N, const double* filter, const int F, double* output, const int step, const int mode)
{
    int i, j, F_minus_1, N_extended_len, N_extended_right_start;
	int start;
    double sum, tmp;
	DTYPE *buffer;
    double* ptr_w = output;
    
	F_minus_1 = F - 1;
	start = F_minus_1+step-1;

	if(mode != MODE_PERIODIZATION){
		N_extended_len = N + 2*F_minus_1;
		N_extended_right_start = N + F_minus_1;

		buffer = (double*) calloc(N_extended_len, sizeof(double));
		if(buffer == NULL)
			return -1;

		memcpy(buffer+F_minus_1, input, sizeof(double) * N);
	} else {
		N_extended_len = N + F-1;
		N_extended_right_start = N-1 + F/2;

		buffer = (double*) calloc(N_extended_len, sizeof(double));
		if(buffer == NULL)
			return -1;

		memcpy(buffer+F/2-1, input, sizeof(double) * N);
		start -= 1;
	}

	switch(mode){

		case MODE_PERIODIZATION:
			if(N%2){ // odd - repeat last element
				buffer[N_extended_right_start] = input[N-1];
				for(j = 1; j < F/2; ++j)
					buffer[N_extended_right_start+j] = buffer[F/2-2 + j]; // copy from beggining of `input` to right
				for(j = 0; j < F/2-1; ++j)								  // copy from 'buffer' to left
					buffer[F/2-2-j] =  buffer[N_extended_right_start-j];
			} else {
				for(j = 0; j < F/2; ++j)		
					buffer[N_extended_right_start+j] = input[j%N]; // copy from beggining of `input` to right
				for(j = 0; j < F/2-1; ++j)						   // copy from 'buffer' to left
					buffer[F/2-2-j] =  buffer[N_extended_right_start-1-j];
			}
			break;

		case MODE_SYMMETRIC:
			for(j = 0; j < N; ++j){
				buffer[F_minus_1-1-j] = input[j%N];
				buffer[N_extended_right_start+j] = input[N-1-(j%N)];
			}
			i=j;
			// use `buffer` as source
			for(; j < F_minus_1; ++j){
				buffer[F_minus_1-1-j] =  buffer[N_extended_right_start-1+i-j];
				buffer[N_extended_right_start+j] = buffer[F_minus_1+j-i];
			}
			break;

		case MODE_ASYMMETRIC:
			for(j = 0; j < N; ++j){
				buffer[F_minus_1-1-j] = input[0] - input[j%N];
				buffer[N_extended_right_start+j] = (input[N-1] - input[N-1-(j%N)]);
			}
			i=j;
			// use `buffer` as source
			for(; j < F_minus_1; ++j){
				buffer[F_minus_1-1-j] =  buffer[N_extended_right_start-1+i-j];
				buffer[N_extended_right_start+j] = buffer[F_minus_1+j-i];
			}
			break;

		case MODE_SMOOTH:
			if(N>1){
				tmp = input[0]-input[1];
				for(j = 0; j < F_minus_1; ++j)
					buffer[j] = input[0] +	(tmp * (F_minus_1-j));
				tmp = input[N-1]-input[N-2];
				for(j = 0; j < F_minus_1; ++j)
					buffer[N_extended_right_start+j] = input[N-1] + (tmp*j);
				break;
			}

		case MODE_CONSTANT_EDGE:
			for(j = 0; j < F_minus_1; ++j){
				buffer[j] = input[0];			
				buffer[N_extended_right_start+j] = input[N-1];
			}
			break;

		case MODE_PERIODIC:
			for(j = 0; j < F_minus_1; ++j)		
				buffer[N_extended_right_start+j] = input[j%N]; // copy from beggining of `input` to right
			
			for(j = 0; j < F_minus_1; ++j)				       // copy from 'buffer' to left
				buffer[F_minus_1-1-j] =  buffer[N_extended_right_start-1-j];			
			break;

		case MODE_ZEROPAD:
		default:
			//memset(buffer, 0, sizeof(double)*F_minus_1);
			//memset(buffer+N_extended_right_start, 0, sizeof(double)*F_minus_1);
			//memcpy(buffer+N_extended_right_start, buffer, sizeof(double)*F_minus_1);
			break;
	}


	///////////////////////////////////////////////////////////////////////
    // F - N-1		- filter in input range
   
	for(i=start; i < N_extended_len; i+=step){					// input elements
        sum = 0;
		for(j = 0; j < F; ++j){
            sum += buffer[i-j]*filter[j];
		}
        *(ptr_w++) = sum;
    }  
    
	free(buffer);
	return 0;
}

///////////////////////////////////////////////////////////////////////////////
// requires zero-filled output buffer
// output is larger than input
// performs "normal" convolution of "upsampled" input coeffs array with filter

int upsampling_convolution_full(const DTYPE* input, const int N, const double* filter, const int F, double* output, const int O){
	register int i;
	register int j;
	double *ptr_out;

	if(F<2) return -1;
	ptr_out = output + ((N-1) << 1);

	for(i = N-1; i >= 0; --i){
		// sliding in filter from the right (end of input)
		// i0 0  i1 0  i2 0
		//                f1 -> o1
		//             f1 f2 -> o2
		//          f1 f2 f3 -> o3
		for(j = 0; j < F; ++j)		
			ptr_out[j] += input[i] * filter[j]; // input[i] - const in loop
		ptr_out -= 2;
	}
	return 0;
}

///////////////////////////////////////////////////////////////////////////////
// performs IDWT for all modes (PERIODIZATION & others)
// 
// characteristic: - splits filter to even and odd elements
//                 - extends input for PERIODIZATION mode


int upsampling_convolution_valid_sf(const DTYPE* input, const int N, const double* filter, const int F, double* output, const int O, const int mode){
	double *ptr_out = output;
	double *filter_even, *filter_odd;
	DTYPE *periodization_buf = NULL, *periodization_buf_rear = NULL;
	DTYPE *ptr_base;
	double sum_even, sum_odd;
	int i, j, k, N_p = 0;
	int F_2 = F/2;	// F/2

	if(F%2) return -3;

	if(F_2 > N){
		if(mode == MODE_PERIODIZATION){
			N_p = F_2-1 +N;
			periodization_buf = (double*) wtcalloc(N_p, sizeof(DTYPE));
			k = (F_2-1)/2;

			for(i=k; i < k+N; ++i)
				periodization_buf[i] = input[(i-k)%N];
			//if(N%2)
			//	periodization_buf[i++] = input[N-1];
			periodization_buf_rear = periodization_buf+i-1;

			j = i-k;
			for(; i < N_p; ++i)
				periodization_buf[i] = periodization_buf[i-j];

			j = 0;
			for(i=k-1; i >= 0; --i){
				periodization_buf[i] = periodization_buf_rear[j];
				--j;
			}
			
			if(F_2%2==0){
				// cheap result fix
				ptr_out = (double*) wtcalloc(idwt_buffer_length(N, F, MODE_PERIODIZATION), sizeof(double));
				if(ptr_out == NULL)
					return -3;

				upsampling_convolution_valid_sf(periodization_buf, N_p, filter, F, ptr_out, O, MODE_ZEROPAD);

				for(i=2*N-1; i > 0; --i){
					output[i] += ptr_out[i-1];
				}
				output[0] += ptr_out[2*N-1];
				wtfree(ptr_out);
			} else {
				upsampling_convolution_valid_sf(periodization_buf, N_p, filter, F, output, O, MODE_ZEROPAD);
			}
			return 0;
		} else {
			return -2;	// invalid lengths
		}
	}

	// allocate memory for even and odd elements of filter
	filter_even = (double*) malloc(F_2 * sizeof(double));
	filter_odd = (double*) malloc(F_2 * sizeof(double));

	if(filter_odd == NULL || filter_odd == NULL){
		return -1;
	}

	// copy values
	for(i = 0; i < F_2; ++i){
		filter_even[i] = filter[i << 1];
		filter_odd[i] = filter[(i << 1) + 1];
	}

	///////////////////////////////////////////////////////////////////////////
	// MODE_PERIODIZATION 
	// this part is quite complicated

		if(mode == MODE_PERIODIZATION){

			k = F_2-1;

			N_p = F_2-1 + (int) ceil(k/2.); /*split filter len correct. +  extra samples*/
			
			if(N_p > 0){
				periodization_buf = (double*) calloc(N_p, sizeof(double));
				periodization_buf_rear = (double*) calloc(N_p, sizeof(double));
			
				if(k <= N){
					memcpy(periodization_buf + N_p - k, input, k * sizeof(double));		// copy from beginning of input to end of buffer
					for(i = 1; i <= (N_p - k); ++i)										// kopiowanie 'cykliczne' od koñca input 
						periodization_buf[(N_p - k) - i] = input[N - (i%N)];
					memcpy(periodization_buf_rear, input + N - k, k * sizeof(double));	// copy from end of input to begginning of buffer
					for(i = 0; i < (N_p - k); ++i)										// kopiowanie 'cykliczne' od pocz¹tku input
						periodization_buf_rear[k + i] = input[i%N];
				} else {
					//printf("see %d line in %s!!\n", __LINE__, __FILE__);
					// TODO: is this ever called? if yes check for errors
					for(i = 0; i < k; ++i)
						periodization_buf[(N_p - k) + i] = input[i % N];
					for(i = 1; i < (N_p - k); ++i)
						periodization_buf[(N_p - k) - i] = input[N - (i%N)];

					//for(i = 0; i < N_p; ++i)
					//	printf("%f ", periodization_buf[i]);
					//printf("--\n");
					//for(i = 0; i < N_p; ++i)
					//	printf("%f ", periodization_buf_rear[i]);
					//printf("\n");
				}

				ptr_base = periodization_buf + F_2 - 1;
				if(k%2 == 1){
					sum_odd = 0;
					for(j = 0; j < F_2; ++j)
						sum_odd += filter_odd[j] * ptr_base[-j];
					*(ptr_out++) += sum_odd;
					--k;
					if(k)
						upsampling_convolution_valid_sf(periodization_buf + 1, N_p-1, filter, F, ptr_out, O-1, MODE_ZEROPAD);
					ptr_out += k; // k0 - 1 // really move backward by 1
				} else if(k){
					upsampling_convolution_valid_sf(periodization_buf, N_p, filter, F, ptr_out, O, MODE_ZEROPAD);
					ptr_out += k;
				}
			}
		}
	// MODE_PERIODIZATION
	///////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////
	// perform _valid_ convolution (only when all filter_even and filter_odd elements are in range of input data)
	// this part is quite simple, no extra hacks
		ptr_base = (DTYPE*)input + F_2 - 1;
		for(i = 0; i < N-(F_2-1); ++i){	// sliding over signal from left to right
			
			sum_even = 0;
			sum_odd = 0;
			// iterate filters elements
			for(j = 0; j < F_2; ++j){
				sum_even += filter_even[j] * ptr_base[i-j];
				sum_odd += filter_odd[j] * ptr_base[i-j];
			}

			*(ptr_out++) += sum_even; 
			*(ptr_out++) += sum_odd;
		}
	//
	///////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////
	// MODE_PERIODIZATION
		if(mode == MODE_PERIODIZATION){
			if(N_p>0){
				k = F_2-1;
				if(k%2 == 1){
					if(F/2 <= N_p - 1) // k > 1 ?
						upsampling_convolution_valid_sf(periodization_buf_rear , N_p-1, filter, F, ptr_out, O-1, MODE_ZEROPAD);
					ptr_out += k; // move forward anyway -> see lower

					if(F_2%2 == 0){ // remaining one element
						ptr_base = periodization_buf_rear + N_p - 1;
						sum_even = 0;
						for(j = 0; j < F_2; ++j)
							sum_even += filter_even[j] * ptr_base[-j];
						*(--ptr_out) += sum_even; // move backward first
					}
				} else {
					if(k)
						upsampling_convolution_valid_sf(periodization_buf_rear, N_p, filter, F, ptr_out, O, MODE_ZEROPAD);
				}
			}
			if(periodization_buf != NULL) wtfree(periodization_buf);
			if(periodization_buf_rear != NULL) wtfree(periodization_buf_rear);
		}
	// MODE_PERIODIZATION
	///////////////////////////////////////////////////////////////////////////

	wtfree(filter_even);
	wtfree(filter_odd);
	return 0;
}

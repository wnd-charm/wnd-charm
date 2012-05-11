#include "WaveletLow.h"
#include "FilterSet.h"

// 5 level 
// bior5.5
// 1.0 x 10^-7
WaveletLow::WaveletLow(double lim)
{
	nsteps = 5;
	mode = WAVELET_LOW;
	
	double bior55LowDe [] = {
		0,  0,   0.03968708834741,   0.00794810863724,  -0.05446378846824,   0.34560528195603,
		0.73666018142821, 0.34560528195603,  -0.05446378846824,   0.00794810863724,
		0.03968708834741,  0
	};
	
	double bior55HighDe [] = {
		-0.01345670945912,  -0.00269496688011,   0.13670658466433,  -0.09350469740094,
		-0.47680326579848,   0.89950610974865,  -0.47680326579848,	-0.09350469740094,
		0.13670658466433,  -0.00269496688011,  -0.01345670945912, 0
	};
	
	double bior55LowRe [] = {
		 0.01345670945912,  -0.00269496688011,  -0.13670658466433,  -0.09350469740094,
		 0.47680326579848,   0.89950610974865,   0.47680326579848, -0.09350469740094,
		 -0.13670658466433,  -0.00269496688011,   0.01345670945912,                  0
	};
	
	double bior55HighRe [] = {
		0,  0,   0.03968708834741,  -0.00794810863724,  -0.05446378846824,
		-0.34560528195603,   0.73666018142821,  -0.34560528195603,
		-0.05446378846824,  -0.00794810863724,   0.03968708834741,   0
	};
	
	FilterSet * bior55 = new FilterSet(1, bior55LowDe, bior55HighDe,  12, -6, -5,
										bior55LowRe, bior55HighRe,  12, -5, -6 );

	filterset = bior55;
	
	if (lim > 0) {
		limit = lim;
	}
	else {
		limit = .0000001;
	}
	
	
	analysisLow = filterset->analysisLow;
	analysisHigh = filterset->analysisHigh;
	synthesisLow = filterset->synthesisLow;
	synthesisHigh = filterset->synthesisHigh;
	
	npad = 16;
	
	vanishing_moments_psi = 55/10;
	vanishing_moments_phi = -1;
	support_width = -1;
	orthogonal = 0;
	biorthogonal = 1;
	symmetry = SYMMETRIC;
	compact_support = 1;
	family_name = "Biorthogonal";
	short_name = "bior5.5";
	
	dec_len = rec_len = 12;
	dec_hi_offset = rec_lo_offset = 0;
	dec_lo_offset = rec_hi_offset = 12;

}

WaveletLow::~WaveletLow() {
	delete filterset;
}


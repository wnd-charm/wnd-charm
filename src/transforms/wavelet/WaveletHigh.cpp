#include "WaveletHigh.h"
#include "FilterSet.h"



// 5 level decomp
// bior3.1
// remove 1.0 x 10^-4.
WaveletHigh::WaveletHigh(double lim)
{
	nsteps = 5;
	mode = WAVELET_HIGH;

	double bior31LowDe [] = {
		 -0.35355339059327379, 1.0606601717798214, 1.0606601717798214, -0.35355339059327379
	};
	
	double bior31HighDe [] = {
		-0.17677669529663689, 0.53033008588991071, -0.53033008588991071, 0.17677669529663689
	};
	
	// change 1 & 3 
	double bior31LowRe [] = {
		 0.17677669529663689, 0.53033008588991071, 0.53033008588991071, 0.17677669529663689
	};
	
	// change 2 and 4
	double bior31HighRe [] = {
		-0.35355339059327379, -1.0606601717798214, 1.0606601717798214, 0.35355339059327379
	};
	
	FilterSet * bior31 = new FilterSet(1, bior31LowDe, bior31HighDe,  4, -2, -2,
										bior31LowRe, bior31HighRe,  4, -2, -2 );
	
	filterset = bior31;
	
	if (lim > 0) {
		limit = lim;
	}
	else {
		limit = .0001;
	}
	
	
	analysisLow = filterset->analysisLow;
	analysisHigh = filterset->analysisHigh;
	synthesisLow = filterset->synthesisLow;
	synthesisHigh = filterset->synthesisHigh;
	
	npad = max(analysisLow->size, analysisHigh->size);
	
	
	vanishing_moments_psi = 31/10;
	vanishing_moments_phi = -1;
	support_width = -1;
	orthogonal = 0;
	biorthogonal = 1;
	symmetry = SYMMETRIC;
	compact_support = 1;
	family_name = "Biorthogonal";
	short_name = "bior";
	
	dec_len = rec_len = 4;
	dec_hi_offset = rec_lo_offset = 0;
	
}

WaveletHigh::~WaveletHigh() {
	delete filterset;
}

#include "WaveletMedium.h"
#include "FilterSet.h"

// 5 level decomp
// bior3.1
// remove 1.0 x 10^-5.

#ifndef WIN32
#define max(a,b) (a) < (b) ? (b) : (a)
#endif

WaveletMedium::WaveletMedium(double lim)
{
	nsteps = 5;
	mode = WAVELET_MED;
	
	double bior31LowDe [] = {
		 -0.35355339059327,   1.06066017177982,   1.06066017177982,  -0.35355339059327
	};
	
	double bior31HighDe [] = {
		-0.17677669529664 ,  0.53033008588991,  -0.53033008588991,   0.17677669529664
	};
	
	double bior31LowRe [] = {
		 0.17677669529664 ,  0.53033008588991,   0.53033008588991,   0.17677669529664
	};
	
	double bior31HighRe [] = {
		-0.35355339059327 , -1.06066017177982 ,  1.06066017177982,   0.35355339059327
	};
	
	FilterSet * bior31 = new FilterSet(1, bior31LowDe, bior31HighDe,  4, -2, -2,
										bior31LowRe, bior31HighRe,  4, -2, -2 );
										
	filterset = bior31;
	
	if (lim > 0) {
		limit = lim;
	}
	else {
		limit = .00001;
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

WaveletMedium::~WaveletMedium() {
	delete filterset;
}


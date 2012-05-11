#include "Symlet5.h"
#include "FilterSet.h"


#define max(a, b)  (((a) > (b)) ? (a) : (b))

// 5 level decomp
// bior3.1
// remove 1.0 x 10^-4.
Symlet5::Symlet5(double lim, int steps)
{
	nsteps = 5;
        if (steps>0) nsteps=steps;

	double sym5LowDe [] = 	{0.027333068345077982, 0.029519490925774643, -0.039134249302383094, 0.1993975339773936, 0.72340769040242059, 0.63397896345821192, 0.016602105764522319, -0.17532808990845047, -0.021101834024758855, 0.019538882735286728};

	double sym5HighDe [] = {-0.019538882735286728, -0.021101834024758855, 0.17532808990845047, 0.016602105764522319, -0.63397896345821192, 0.72340769040242059, -0.1993975339773936, -0.039134249302383094, -0.029519490925774643, 0.027333068345077982};
	
	// change 1 & 3 
	double sym5LowRe [] = {0.019538882735286728, -0.021101834024758855, -0.17532808990845047, 0.016602105764522319, 0.63397896345821192, 0.72340769040242059, 0.1993975339773936, -0.039134249302383094, 0.029519490925774643, 0.027333068345077982};
	
	// change 2 and 4
	double sym5HighRe [] = {0.027333068345077982, -0.029519490925774643, -0.039134249302383094, -0.1993975339773936, 0.72340769040242059, -0.63397896345821192, 0.016602105764522319, 0.17532808990845047, -0.021101834024758855, -0.019538882735286728};
	
	FilterSet * sym5 = new FilterSet(1, sym5LowDe, sym5HighDe,  10, -5, -5,
										sym5LowRe, sym5HighRe,  10, -5, -5 );
	
	filterset = sym5;
	
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
	
	
	vanishing_moments_psi = 5;
	vanishing_moments_phi = 0;
	support_width = 2*5-1;
	orthogonal = 1;
	biorthogonal = 1;
	symmetry = NEAR_SYMMETRIC;
	compact_support = 1;
//	family_name = "Symlet";
//	short_name = "Symlet";
	
	dec_len = rec_len = 5 << 1;
	dec_hi_offset = rec_lo_offset = 0;
	
}

Symlet5::~Symlet5() {
	delete filterset;
}

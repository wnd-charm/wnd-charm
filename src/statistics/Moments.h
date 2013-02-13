#ifndef MomentsH
#define MomentsH

#include <stddef.h>
#include <cfloat> // DBL_MAX
#include <cmath>
class Moments {
  private:
  	double _min, _max, _mean, M2, M3, M4;
  	size_t _n;

  public:
  	Moments() { reset(); }
	void reset() { _mean = M2 = M3 = M4 = 0.0; _min = DBL_MAX; _max = -DBL_MAX; _n = 0;}
  	void add (const double x) {
  		size_t n1;
  		double delta, delta_n, delta_n2, term1;

		n1 = _n;
		_n = _n + 1;
		delta = x - _mean;
		delta_n = delta / _n;
		delta_n2 = delta_n * delta_n;
		term1 = delta * delta_n * n1;
		_mean = _mean + delta_n;
		M4 = M4 + term1 * delta_n2 * (_n*_n - 3*_n + 3) + 6 * delta_n2 * M2 - 4 * delta_n * M3;
		M3 = M3 + term1 * delta_n * (_n - 2) - 3 * delta_n * M2;
		M2 = M2 + term1;
		
		if (x > _max) _max = x;
		if (x < _min) _min = x;
  	}
  	
  	size_t n()         { return _n; }
  	double min()       { return _min; }
  	double max()       { return _max; }
  	double mean()      { return _mean; }
  	double std ()      { return (_n > 2 ? sqrt ( M2/(_n - 1) ) : 0.0); }
  	double var ()      { return (_n > 2 ? M2/(_n - 1) : 0.0); }
  	double skewness () { return (_n > 3 ? (sqrt (_n) * M3) / pow (M2, 1.5) : 0.0); }
  	double kurtosis () { return (_n > 4 ? (_n*M4) / (M2*M2) : 0.0); } // matlab does not subtract 3
  	void momentVector (double *z) { z[0] = mean(); z[1] = std(); z[2] = skewness(); z[3] = kurtosis(); }
};


#endif // MomentsH

#ifndef MomentsH
#define MomentsH

#include <stddef.h>
#include <cfloat> // DBL_MAX
#include <cmath>
#define MIN_VAL -FLT_MAX
#define MAX_VAL FLT_MAX

class Moments4 {
  private:
  	double _min, _max, _mean, M2, M3, M4;
  	size_t _n;

  public:
  	Moments4() { reset(); }
	void reset() { _mean = M2 = M3 = M4 = 0.0; _min = MAX_VAL; _max = MIN_VAL; _n = 0;}
  	inline double add (const double x) {
  		size_t n1;
  		double delta, delta_n, delta_n2, term1;
		if (std::isnan (x) || x > MAX_VAL || x < MIN_VAL) return (x);

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
		return (x);
  	}
 	
  	size_t n()         const { return _n; }
  	double min()       const { return _min; }
  	double max()       const { return _max; }
  	double mean()      const { return _mean; }
  	double std ()      const { return (_n > 2 ? sqrt ( M2/(_n - 1) ) : 0.0); }
  	double var ()      const { return (_n > 2 ? M2/(_n - 1) : 0.0); }
  	double skewness () const { return (_n > 3 ? (sqrt (_n) * M3) / pow (M2, 1.5) : 0.0); }
  	double kurtosis () const { return (_n > 4 ? (_n*M4) / (M2*M2) : 0.0); } // matlab does not subtract 3
  	void momentVector (double *z) { z[0] = mean(); z[1] = std(); z[2] = skewness(); z[3] = kurtosis(); }
};

// functor to call add on a reference using the () operator
// for example, using Eigen: ReadablePixels().unaryExpr (Moments4func(stats)).sum();
class Moments4func {
	Moments4 &moments;
	public:
		Moments4func (Moments4 &in_moments) : moments (in_moments) {in_moments.reset();}
	  	const double operator()(const double &x) const {
	  		return (moments.add (x));
	  	}
};

class Moments2 {
  private:
  	double _min, _max, _mean, M2;
  	size_t _n;

  public:
  	Moments2() { reset(); }
	void reset() { _mean = M2 = 0.0; _min = DBL_MAX; _max = -DBL_MAX; _n = 0;}
  	inline double add (const double x) {
  		size_t n1;
  		double delta, delta_n, term1;
		if (std::isnan (x) || x > MAX_VAL || x < MIN_VAL) return (x);

		n1 = _n;
		_n = _n + 1;
		delta = x - _mean;
		delta_n = delta / _n;
		term1 = delta * delta_n * n1;
		_mean = _mean + delta_n;
		M2 += term1;
		
		if (x > _max) _max = x;
		if (x < _min) _min = x;
		return (x);
  	}
  	
  	size_t n()    const { return _n; }
  	double min()  const { return _min; }
  	double max()  const { return _max; }
  	double mean() const { return _mean; }
  	double std () const { return (_n > 2 ? sqrt ( M2/(_n - 1) ) : 0.0); }
  	double var () const { return (_n > 2 ? (M2/(_n - 1)) : 0.0); }
  	void momentVector (double *z) const { z[0] = mean(); z[1] = std(); }
};

// functor to call add on a reference using the () operator
// for example, using Eigen: ReadablePixels().unaryExpr (Moments4func(stats)).sum();
// N.B.: The sum() in the expression above is to force evaluation of all of the coefficients.
// The return value of Eigen's unaryExpr is a unaryExpr, which doesn't actually do anything until its assigned to something.
class Moments2func {
	Moments2 &moments;
	public:
		Moments2func (Moments2 &in_moments) : moments (in_moments) {in_moments.reset();}
	  	const double operator()(const double &x) const {
	  		return (moments.add (x));
	  	}
};

#endif // MomentsH

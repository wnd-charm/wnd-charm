#ifndef __FEATURE_ALGORITHMS_H_
#define __FEATURE_ALGORITHMS_H_

#include "cmatrix.h"
#include <vector>
#include <string>


class ImageMatrix;

class FeatureAlgorithm {
	public:
		std::string name;
		int n_features;
		virtual std::vector<double> calculate( ImageMatrix * IN_matrix ) const { return std::vector<double>(); };
		void print_info() const;	
	protected:
		FeatureAlgorithm ();
		FeatureAlgorithm (const std::string &s,const int i);
		FeatureAlgorithm (const char *s,const int i);
};



class EmptyFeatureAlgorithm : public FeatureAlgorithm {
	public:
		virtual std::vector<double> calculate( ImageMatrix * IN_matrix ) const { return std::vector<double>(); };
		EmptyFeatureAlgorithm () : FeatureAlgorithm ("Empty", 0) {};
		EmptyFeatureAlgorithm (const std::string &s) : FeatureAlgorithm (s, 0) {};
		EmptyFeatureAlgorithm (const char *s) : FeatureAlgorithm (s, 0) {};
};

class ChebyshevFourierCoefficients : public FeatureAlgorithm {
	public:
		ChebyshevFourierCoefficients();
		virtual std::vector<double> calculate( ImageMatrix * IN_matrix ) const;
};

class ChebyshevCoefficients : public FeatureAlgorithm {
	public:
		ChebyshevCoefficients();
		virtual std::vector<double> calculate( ImageMatrix * IN_matrix ) const;
};

class ZernikeCoefficients : public FeatureAlgorithm {
	public:
		ZernikeCoefficients();
		virtual std::vector<double> calculate( ImageMatrix * IN_matrix ) const;
};

class HaralickTextures : public FeatureAlgorithm {
	public:
		HaralickTextures();
		virtual std::vector<double> calculate( ImageMatrix * IN_matrix ) const;
};

class MultiscaleHistograms : public FeatureAlgorithm {
	public:
		MultiscaleHistograms();
		virtual std::vector<double> calculate( ImageMatrix * IN_matrix ) const;
};

class TamuraTextures : public FeatureAlgorithm {
	public:
		TamuraTextures();
		virtual std::vector<double> calculate( ImageMatrix * IN_matrix ) const;
};

class CombFirstFourMoments : public FeatureAlgorithm {
	public:
		CombFirstFourMoments();
		virtual std::vector<double> calculate( ImageMatrix * IN_matrix ) const;
};

class RadonCoefficients : public FeatureAlgorithm {
	public:
		RadonCoefficients();
		virtual std::vector<double> calculate( ImageMatrix * IN_matrix ) const;
};

class FractalFeatures : public FeatureAlgorithm {
	public:
		FractalFeatures();
		virtual std::vector<double> calculate( ImageMatrix * IN_matrix ) const;
};

class PixelIntensityStatistics : public FeatureAlgorithm {
	public:
		PixelIntensityStatistics();
		virtual std::vector<double> calculate( ImageMatrix * IN_matrix ) const;
};

class EdgeFeatures : public FeatureAlgorithm {
	public:
		EdgeFeatures();
		virtual std::vector<double> calculate( ImageMatrix * IN_matrix ) const;
};

class ObjectFeatures : public FeatureAlgorithm {
	public:
		ObjectFeatures();
		virtual std::vector<double> calculate( ImageMatrix * IN_matrix ) const;
};

class InverseObjectFeatures : public FeatureAlgorithm {
	public:
		InverseObjectFeatures();
		virtual std::vector<double> calculate( ImageMatrix * IN_matrix ) const;
};

class GaborTextures : public FeatureAlgorithm {
	public:
		GaborTextures();
		virtual std::vector<double> calculate( ImageMatrix * IN_matrix ) const;
};

class GiniCoefficient : public FeatureAlgorithm {
	public:
		GiniCoefficient();
		virtual std::vector<double> calculate( ImageMatrix * IN_matrix ) const;
};

class ColorHistogram : public FeatureAlgorithm {
	public:
		ColorHistogram();
		virtual std::vector<double> calculate( ImageMatrix * IN_matrix ) const;
};

#endif //__FEATURE_ALGORITHMS_H_

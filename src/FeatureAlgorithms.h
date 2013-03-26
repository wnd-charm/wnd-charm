/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                               */
/* Copyright (C) 2007 Open Microscopy Environment                                */
/*       Massachusetts Institue of Technology,                                   */
/*       National Institutes of Health,                                          */
/*       University of Dundee                                                    */
/*                                                                               */
/*                                                                               */
/*                                                                               */
/*    This library is free software; you can redistribute it and/or              */
/*    modify it under the terms of the GNU Lesser General Public                 */
/*    License as published by the Free Software Foundation; either               */
/*    version 2.1 of the License, or (at your option) any later version.         */
/*                                                                               */
/*    This library is distributed in the hope that it will be useful,            */
/*    but WITHOUT ANY WARRANTY; without even the implied warranty of             */
/*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU          */
/*    Lesser General Public License for more details.                            */
/*                                                                               */
/*    You should have received a copy of the GNU Lesser General Public           */
/*    License along with this library; if not, write to the Free Software        */
/*    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA  */
/*                                                                               */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                               */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Written by:                                                                   */
/*      Christopher E. Coletta <colettace [at] mail [dot] nih [dot] gov>         */
/*      Ilya G. Goldberg <goldbergil [at] mail [dot] nih [dot] gov>              */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#ifndef __FEATURE_ALGORITHMS_H_
#define __FEATURE_ALGORITHMS_H_

#include <vector>
#include <string>
#include "Tasks.h"

// Forward declarations
class ImageMatrix;

class FeatureAlgorithm : public ComputationTask {
	public:
		int n_features;
		virtual std::vector<double> execute (const ImageMatrix &IN_matrix) const { return std::vector<double>(); };
		virtual void print_info() const;
		virtual bool register_task() const;
	protected:
		FeatureAlgorithm (const std::string &s,const int i) : ComputationTask (s, FeatureAlgorithmTask) {
			n_features = i;
		};
};

class EmptyFeatureAlgorithm : public FeatureAlgorithm {
	public:
		virtual std::vector<double> execute (const ImageMatrix &IN_matrix) const { return std::vector<double>(); };
		EmptyFeatureAlgorithm () : FeatureAlgorithm ("Empty", 0) {};
		EmptyFeatureAlgorithm (const std::string &s) : FeatureAlgorithm (s, 0) {};
};

class ChebyshevFourierCoefficients : public FeatureAlgorithm {
	public:
		ChebyshevFourierCoefficients();
		virtual std::vector<double> execute (const ImageMatrix &IN_matrix) const;
};

class ChebyshevCoefficients : public FeatureAlgorithm {
	public:
		ChebyshevCoefficients();
		virtual std::vector<double> execute (const ImageMatrix &IN_matrix) const;
};

class ZernikeCoefficients : public FeatureAlgorithm {
	public:
		ZernikeCoefficients();
		virtual std::vector<double> execute (const ImageMatrix &IN_matrix) const;
};

class HaralickTextures : public FeatureAlgorithm {
	public:
		HaralickTextures();
		virtual std::vector<double> execute (const ImageMatrix &IN_matrix) const;
};

class MultiscaleHistograms : public FeatureAlgorithm {
	public:
		MultiscaleHistograms();
		virtual std::vector<double> execute (const ImageMatrix &IN_matrix) const;
};

class TamuraTextures : public FeatureAlgorithm {
	public:
		TamuraTextures();
		virtual std::vector<double> execute (const ImageMatrix &IN_matrix) const;
};

class CombFirstFourMoments : public FeatureAlgorithm {
	public:
		CombFirstFourMoments();
		virtual std::vector<double> execute (const ImageMatrix &IN_matrix) const;
};

class RadonCoefficients : public FeatureAlgorithm {
	public:
		RadonCoefficients();
		virtual std::vector<double> execute (const ImageMatrix &IN_matrix) const;
};

class FractalFeatures : public FeatureAlgorithm {
	public:
		FractalFeatures();
		virtual std::vector<double> execute (const ImageMatrix &IN_matrix) const;
};

class PixelIntensityStatistics : public FeatureAlgorithm {
	public:
		PixelIntensityStatistics();
		virtual std::vector<double> execute (const ImageMatrix &IN_matrix) const;
};

class EdgeFeatures : public FeatureAlgorithm {
	public:
		EdgeFeatures();
		virtual std::vector<double> execute (const ImageMatrix &IN_matrix) const;
};

class ObjectFeatures : public FeatureAlgorithm {
	public:
		ObjectFeatures();
		virtual std::vector<double> execute (const ImageMatrix &IN_matrix) const;
};

class InverseObjectFeatures : public FeatureAlgorithm {
	public:
		InverseObjectFeatures();
		virtual std::vector<double> execute (const ImageMatrix &IN_matrix) const;
};

class GaborTextures : public FeatureAlgorithm {
	public:
		GaborTextures();
		virtual std::vector<double> execute (const ImageMatrix &IN_matrix) const;
};

class GiniCoefficient : public FeatureAlgorithm {
	public:
		GiniCoefficient();
		virtual std::vector<double> execute (const ImageMatrix &IN_matrix) const;
};

class ColorHistogram : public FeatureAlgorithm {
	public:
		ColorHistogram();
		virtual std::vector<double> execute (const ImageMatrix &IN_matrix) const;
};

#endif //__FEATURE_ALGORITHMS_H_

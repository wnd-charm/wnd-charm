/* FeatureAlgorithm.cpp */
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

#include "FeatureNames.h"
#include "FeatureAlgorithms.h"
#include "cmatrix.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
//start #including the functions directly once you start pulling them out of cmatrix
//#include "transforms/Chebyshev.h"

/* global variable */
extern int verbosity;

void FeatureAlgorithm::print_info() const {
	std::cout << typeLabel() << " '" << name << "' (" << n_features << " features) " << std::endl;
}

bool FeatureAlgorithm::register_task() const {
	return (FeatureNames::registerFeatureAlgorithm (this));
}

//===========================================================================
ChebyshevFourierCoefficients::ChebyshevFourierCoefficients() : FeatureAlgorithm ("Chebyshev-Fourier Coefficients", 32) {
//	cout << "Instantiating new " << name << " object." << endl;
}

std::vector<double> ChebyshevFourierCoefficients::execute (const ImageMatrix &IN_matrix) const {
	std::vector<double> coeffs;	
	if (verbosity > 3) std::cout << "calculating " << name << std::endl;
	coeffs.resize (n_features, 0);

	IN_matrix.ChebyshevFourierTransform2D(coeffs.data());

	return coeffs;
}

// Register a static instance of the class using a global bool
static bool ChebyshevFourierCoefficientsReg = ComputationTaskInstances::add (new ChebyshevFourierCoefficients);


//===========================================================================
ChebyshevCoefficients::ChebyshevCoefficients() : FeatureAlgorithm ("Chebyshev Coefficients", 32)  {
//	cout << "Instantiating new " << name << " object." << endl;
}

/**
 * Chebyshev Coefficients are calculated by performing a Chebyshev transform,
 * and generating a histogram of pixel intensities.
 *
 */
std::vector<double> ChebyshevCoefficients::execute (const ImageMatrix &IN_matrix) const {
	std::vector<double> coeffs;
	if (verbosity > 3) std::cout << "calculating " << name << std::endl;

	coeffs.resize (n_features, 0);
	IN_matrix.ChebyshevStatistics2D(coeffs.data(), 0, 32);

	return coeffs;
}

// Register a static instance of the class using a global bool
static bool ChebyshevCoefficientsReg = ComputationTaskInstances::add (new ChebyshevCoefficients);


//===========================================================================

ZernikeCoefficients::ZernikeCoefficients() : FeatureAlgorithm ("Zernike Coefficients", 72) {
	//cout << "Instantiating new " << name << " object." << endl;
}

std::vector<double> ZernikeCoefficients::execute (const ImageMatrix &IN_matrix) const {
	std::vector<double> coeffs;
	if (verbosity > 3) std::cout << "calculating " << name << std::endl;

	coeffs.resize (n_features, 0);

	long output_size;   // output size is normally 72

	IN_matrix.zernike2D(coeffs.data(), &output_size);

	return coeffs;
}

// Register a static instance of the class using a global bool
static bool ZernikeCoefficientsReg = ComputationTaskInstances::add (new ZernikeCoefficients);


//===========================================================================

HaralickTextures::HaralickTextures() : FeatureAlgorithm ("Haralick Textures", 28) {
	//cout << "Instantiating new " << name << " object." << endl;
}

std::vector<double> HaralickTextures::execute (const ImageMatrix &IN_matrix) const {
	std::vector<double> coeffs;
	if (verbosity > 3) std::cout << "calculating " << name << std::endl;

	coeffs.resize (n_features, 0);

	IN_matrix.HaralickTexture2D(0,coeffs.data());

	return coeffs;
}

// Register a static instance of the class using a global bool
static bool HaralickTexturesReg = ComputationTaskInstances::add (new HaralickTextures);


//===========================================================================

MultiscaleHistograms::MultiscaleHistograms() : FeatureAlgorithm ("Multiscale Histograms", 24) {
	//cout << "Instantiating new " << name << " object." << endl;
}

std::vector<double> MultiscaleHistograms::execute (const ImageMatrix &IN_matrix) const {
	std::vector<double> coeffs;
	if (verbosity > 3) std::cout << "calculating " << name << std::endl;

	coeffs.resize (n_features, 0);

	IN_matrix.MultiScaleHistogram(coeffs.data());

	return coeffs;
}

// Register a static instance of the class using a global bool
static bool MultiscaleHistogramsReg = ComputationTaskInstances::add (new MultiscaleHistograms);


//===========================================================================

TamuraTextures::TamuraTextures() : FeatureAlgorithm ("Tamura Textures", 6) {
	//cout << "Instantiating new " << name << " object." << endl;
}

std::vector<double> TamuraTextures::execute (const ImageMatrix &IN_matrix) const {
	std::vector<double> coeffs;
	if (verbosity > 3) std::cout << "calculating " << name << std::endl;

	coeffs.resize (n_features, 0);

	IN_matrix.TamuraTexture2D(coeffs.data());

	return coeffs;
}

// Register a static instance of the class using a global bool
static bool TamuraTexturesReg = ComputationTaskInstances::add (new TamuraTextures);


//===========================================================================

CombFirstFourMoments::CombFirstFourMoments() : FeatureAlgorithm ("Comb Moments", 48) {
	//cout << "Instantiating new " << name << " object." << endl;
}

std::vector<double> CombFirstFourMoments::execute (const ImageMatrix &IN_matrix) const {
	std::vector<double> coeffs;
	if (verbosity > 3) std::cout << "calculating " << name << std::endl;

	coeffs.resize (n_features, 0);

	IN_matrix.CombFirstFourMoments2D(coeffs.data());

	return coeffs;
}

// Register a static instance of the class using a global bool
static bool CombFirstFourMomentsReg = ComputationTaskInstances::add (new CombFirstFourMoments);


//===========================================================================

RadonCoefficients::RadonCoefficients() : FeatureAlgorithm ("Radon Coefficients", 12) {
	//cout << "Instantiating new " << name << " object." << endl;
}

std::vector<double> RadonCoefficients::execute (const ImageMatrix &IN_matrix) const {
	std::vector<double> coeffs;
	if (verbosity > 3) std::cout << "calculating " << name << std::endl;

	coeffs.resize (n_features, 0);

	IN_matrix.RadonTransform2D(coeffs.data());

	return coeffs;
}

// Register a static instance of the class using a global bool
static bool RadonCoefficientsReg = ComputationTaskInstances::add (new RadonCoefficients);


//===========================================================================
/* fractal 
   brownian fractal analysis 
   bins - the maximal order of the fractal
   output - array of the size k
   the code is based on: CM Wu, YC Chen and KS Hsieh, Texture features for classification of ultrasonic liver images, IEEE Trans Med Imag 11 (1992) (2), pp. 141Ð152.
   method of approaximation of CC Chen, JS Daponte and MD Fox, Fractal feature analysis and classification in medical imaging, IEEE Trans Med Imag 8 (1989) (2), pp. 133Ð142.
*/
FractalFeatures::FractalFeatures() : FeatureAlgorithm ("Fractal Features", 20) {
	//cout << "Instantiating new " << name << " object." << endl;
}

std::vector<double> FractalFeatures::execute (const ImageMatrix &IN_matrix) const {
	std::vector<double> coeffs;
	if (verbosity > 3) std::cout << "calculating " << name << std::endl;

	coeffs.resize (n_features, 0);

	int bins = n_features;
	int width = IN_matrix.width;
	int height = IN_matrix.height;
	readOnlyPixels IN_matrix_pix_plane = IN_matrix.ReadablePixels();
	int x, y, k, bin = 0;
	int K = ( ( width > height ) ? height : width) / 5; // MIN
	int step = (int) floor ( K / bins );
	if( step < 1 )
		step = 1;   // avoid an infinite loop if the image is small
	for( k = 1; k < K; k = k + step ) {  
		double sum = 0.0;
		for( x = 0; x < width; x++ )
			for( y = 0; y < height - k; y++ )
				sum += fabs( IN_matrix_pix_plane(y,x) - IN_matrix_pix_plane(y+k,x) );
		for( x = 0; x < width - k; x++ )
			for( y = 0; y < height; y++ )
				sum += fabs( IN_matrix_pix_plane(y,x) - IN_matrix_pix_plane(y,x + k) );
		if( bin < bins )
			coeffs[ bin++ ] = sum / ( width * ( width - k ) + height * ( height - k ) );    
	}

	return coeffs;
}

// Register a static instance of the class using a global bool
static bool FractalFeaturesReg = ComputationTaskInstances::add (new FractalFeatures);


//===========================================================================

PixelIntensityStatistics::PixelIntensityStatistics() : FeatureAlgorithm ("Pixel Intensity Statistics", 5) {
	//cout << "Instantiating new " << name << " object." << endl;
}

std::vector<double> PixelIntensityStatistics::execute (const ImageMatrix &IN_matrix) const {
	std::vector<double> coeffs;
	if (verbosity > 3) std::cout << "calculating " << name << std::endl;

	coeffs.resize (n_features, 0);
	
	Moments2 stats;
	IN_matrix.GetStats (stats);
	coeffs[0] = stats.mean();
	coeffs[1] = IN_matrix.get_median();
	coeffs[2] = stats.std();
	coeffs[3] = stats.min();
	coeffs[4] = stats.max();

	return coeffs;
}

// Register a static instance of the class using a global bool
static bool PixelIntensityStatisticsReg = ComputationTaskInstances::add (new PixelIntensityStatistics);

        
//===========================================================================

EdgeFeatures::EdgeFeatures() : FeatureAlgorithm ("Edge Features", 28) {
	//cout << "Instantiating new " << name << " object." << endl;
}

std::vector<double> EdgeFeatures::execute (const ImageMatrix &IN_matrix) const {
	std::vector<double> coeffs;
	if (verbosity > 3) std::cout << "calculating " << name << std::endl;

	coeffs.resize (n_features, 0);

	unsigned long EdgeArea = 0;
	double MagMean=0, MagMedian=0, MagVar=0, MagHist[8]={0,0,0,0,0,0,0,0}, DirecMean=0, DirecMedian=0, DirecVar=0, DirecHist[8]={0,0,0,0,0,0,0,0}, DirecHomogeneity=0, DiffDirecHist[4]={0,0,0,0};
	IN_matrix.EdgeStatistics(&EdgeArea, &MagMean, &MagMedian, &MagVar, MagHist, &DirecMean, &DirecMedian, &DirecVar, DirecHist, &DirecHomogeneity, DiffDirecHist, 8);


	int j, here = 0;
	coeffs[here++] = double( EdgeArea );

	for( j=0; j<4; j++ ){
		coeffs[here++] = DiffDirecHist[j];
	}
	for( j=0; j<8; j++ ){
		coeffs[here++] = DirecHist[j];
	}

	coeffs[here++] = DirecHomogeneity;
	coeffs[here++] = DirecMean;
	coeffs[here++] = DirecMedian;
	coeffs[here++] = DirecVar;

	for( j=0; j<8; j++ ){
		coeffs[here++] = MagHist[j];
	}

	coeffs[here++] = MagMean;
	coeffs[here++] = MagMedian;
	coeffs[here++] = MagVar;

	return coeffs;
}

// Register a static instance of the class using a global bool
static bool EdgeFeaturesReg = ComputationTaskInstances::add (new EdgeFeatures);


//===========================================================================

ObjectFeatures::ObjectFeatures() : FeatureAlgorithm ("Otsu Object Features", 34) {
	//cout << "Instantiating new " << name << " object." << endl;
}

std::vector<double> ObjectFeatures::execute (const ImageMatrix &IN_matrix) const {
	std::vector<double> coeffs;
	if (verbosity > 3) std::cout << "calculating " << name << std::endl;

	coeffs.resize (n_features, 0);

	unsigned long feature_count=0, AreaMin=0, AreaMax=0;
	long Euler=0;
	unsigned int AreaMedian=0,
			area_histogram[10]={0,0,0,0,0,0,0,0,0,0},
			dist_histogram[10]={0,0,0,0,0,0,0,0,0,0};

	double centroid_x=0, centroid_y=0, AreaMean=0, AreaVar=0, DistMin=0,
				 DistMax=0, DistMean=0, DistMedian=0, DistVar=0;

	IN_matrix.FeatureStatistics(&feature_count, &Euler, &centroid_x, &centroid_y,
			&AreaMin, &AreaMax, &AreaMean, &AreaMedian,
			&AreaVar, area_histogram, &DistMin, &DistMax,
			&DistMean, &DistMedian, &DistVar, dist_histogram, 10);

	int j, here = 0;

	for( j = 0; j < 10; j++ ){
		coeffs[here++] = area_histogram[j];
	}

	coeffs[here++] = AreaMax;
	coeffs[here++] = AreaMean;
	coeffs[here++] = AreaMedian;
	coeffs[here++] = AreaMin;
	coeffs[here++] = AreaVar;
	coeffs[here++] = centroid_x;
	coeffs[here++] = centroid_y;
	coeffs[here++] = feature_count;

	for( j = 0; j < 10; j++ ) {
		coeffs[here++] = dist_histogram[j];
	}

	coeffs[here++] = DistMax;
	coeffs[here++] = DistMean;
	coeffs[here++] = DistMedian;
	coeffs[here++] = DistMin;
	coeffs[here++] = DistVar;
	coeffs[here++] = Euler;

	return coeffs;
}

// Register a static instance of the class using a global bool
static bool ObjectFeaturesReg = ComputationTaskInstances::add (new ObjectFeatures);


//===========================================================================
InverseObjectFeatures::InverseObjectFeatures() : FeatureAlgorithm ("Inverse-Otsu Object Features", 34) {
	//cout << "Instantiating new " << name << " object." << endl;
}

std::vector<double> InverseObjectFeatures::execute (const ImageMatrix &IN_matrix) const {
	ImageMatrix InvMatrix;
	InvMatrix.copy (IN_matrix);
	InvMatrix.invert();
	static ObjectFeatures ObjFeaturesInst;
	return (ObjFeaturesInst.execute (InvMatrix));
}

// Register a static instance of the class using a global bool
static bool InverseObjectFeaturesReg = ComputationTaskInstances::add (new InverseObjectFeatures);


//===========================================================================

GaborTextures::GaborTextures() : FeatureAlgorithm ("Gabor Textures", 7) {
	//cout << "Instantiating new " << name << " object." << endl;
}

std::vector<double> GaborTextures::execute (const ImageMatrix &IN_matrix) const {
	std::vector<double> coeffs;
	if (verbosity > 3) std::cout << "calculating " << name << std::endl;

	coeffs.resize (n_features, 0);

	IN_matrix.GaborFilters2D(coeffs.data());
	return coeffs;
}

// Register a static instance of the class using a global bool
static bool GaborTexturesReg = ComputationTaskInstances::add (new GaborTextures);


//===========================================================================

/* gini
   compute the gini coefficient
   
   paper reference: Roberto G. Abraham, Sidney van den Bergh, Preethi Nair, A NEW APPROACH TO GALAXY MORPHOLOGY. I. ANALYSIS OF THE SLOAN DIGITAL SKY
        SURVEY EARLY DATA RELEASE, The Astrophysical Journal, vol. 588, p. 218-229, 2003.
*/
GiniCoefficient::GiniCoefficient() : FeatureAlgorithm ("Gini Coefficient", 1) {
	//cout << "Instantiating new " << name << " object." << endl;
}

std::vector<double> GiniCoefficient::execute (const ImageMatrix &IN_matrix) const {
	std::vector<double> coeffs;
	if (verbosity > 3) std::cout << "calculating " << name << std::endl;

	coeffs.resize(n_features, 0);

	long pixel_index, num_pixels;
	double *pixels, mean = 0.0, g = 0.0;
	long i, count = 0;
	double val;

	num_pixels = IN_matrix.height * IN_matrix.width;
	pixels = new double[ num_pixels ];

	readOnlyPixels IN_matrix_pix_plane = IN_matrix.ReadablePixels();
	for( pixel_index = 0; pixel_index < num_pixels; pixel_index++ ) {
		val = IN_matrix_pix_plane.array().coeff(pixel_index);
		if( val > 0 ) {
			pixels[ count ] = val;
			mean += val;
			count++;
		}
	}
	if( count > 0 )
		mean = mean / count;
	qsort( pixels, count, sizeof(double), compare_doubles );

	for( i = 1; i <= count; i++)
		g += (2. * i - count - 1.) * pixels[i-1];
	delete [] pixels;

	if( count <= 1 || mean <= 0.0 )
		coeffs[0] = 0.0;   // avoid division by zero
	else
		coeffs[0] = g / ( mean * count * ( count-1 ) );

	return coeffs;
}

// Register a static instance of the class using a global bool
static bool GiniCoefficientReg = ComputationTaskInstances::add (new GiniCoefficient);



//===========================================================================

/* Color Histogram
   compute the Color Histogram

*/
ColorHistogram::ColorHistogram() : FeatureAlgorithm ("Color Histogram", COLORS_NUM+1) {
	//cout << "Instantiating new " << name << " object." << endl;
}

std::vector<double> ColorHistogram::execute (const ImageMatrix &IN_matrix) const {
	std::vector<double> coeffs;
	if (verbosity > 3) std::cout << "calculating " << name << std::endl;

	coeffs.assign(n_features, 0);
	unsigned int x,y, width = IN_matrix.width, height = IN_matrix.height;
	HSVcolor hsv_pixel;
	unsigned long color_index=0;   
	double certainties[COLORS_NUM+1];

	readOnlyColors clr_plane = IN_matrix.ReadableColors();

	// find the colors
	for( y = 0; y < height; y++ ) {
		for( x = 0; x < width; x++ ) { 
			hsv_pixel = clr_plane (y, x);
			color_index = FindColor( hsv_pixel.h,  hsv_pixel.s, hsv_pixel.v, certainties );
			coeffs[ color_index ]++;
		}
	}
	/* normalize the color histogram */
	for (color_index = 0; color_index <= COLORS_NUM; color_index++)
		coeffs[color_index] /= (width*height);	 

	return coeffs;
}

// Register a static instance of the class using a global bool
static bool ColorHistogramReg = ComputationTaskInstances::add (new ColorHistogram);

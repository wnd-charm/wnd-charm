/* FeatureAlgorithm.cpp */

#include "FeatureAlgorithms.h"
#include "cmatrix.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
//start #including the functions directly once you start pulling them out of cmatrix
//#include "transforms/Chebyshev.h"

void FeatureAlgorithm::print_info() const {
        std::cout << "FeatureAlgorithm: " << name << " (" << n_features << " features) " << std::endl;
}

//===========================================================================
ChebyshevFourierCoefficients::ChebyshevFourierCoefficients() {
	name = "Chebyshev-Fourier Coefficients";
	n_features = 32;
	//cout << "Instantiating new " << name << " object." << endl;
}

std::vector<double> ChebyshevFourierCoefficients::calculate( ImageMatrix * IN_matrix )
{
	std::vector<double> coeffs;	
	std::cout << "calculating " << name << std::endl;
	if( IN_matrix == NULL ) {
		return coeffs;
	}
	coeffs.reserve(n_features-1);
	double temp_vec [32];
	int i;

	for( i = 0; i < n_features; i++ ) temp_vec[i] = 0;
	IN_matrix->ChebyshevFourierTransform2D(temp_vec);
	coeffs.assign( temp_vec, temp_vec + n_features);
	return coeffs;
}

//WNDCHARM_REGISTER_ALGORITHM(ChebyshevFourierCoefficients)

//===========================================================================
ChebyshevCoefficients::ChebyshevCoefficients() {
	name = "Chebyshev Coefficients";
	n_features = 32;
	//cout << "Instantiating new " << name << " object." << endl;
}

/**
 * Chebyshev Coefficients are calculated by performing a Chebyshev transform,
 * and generating a histogram of pixel intensities.
 *
 */
std::vector<double> ChebyshevCoefficients::calculate( ImageMatrix * IN_matrix )
{
	std::vector<double> coeffs;
	std::cout << "calculating " << name << std::endl;
	if( IN_matrix == NULL ) {
		return coeffs;
	}
	coeffs.reserve(n_features-1);
	double temp_vec [32];

	ImageMatrix * temp = IN_matrix->duplicate();
	for( int i = 0; i < n_features; i++ ) temp_vec[i] = 0;
	temp->ChebyshevStatistics2D(temp_vec, 0, 32);
	delete temp;
	coeffs.assign( temp_vec, temp_vec + n_features);

	return coeffs;
}

//WNDCHARM_REGISTER_ALGORITHM(ChebyshevCoefficients)

//===========================================================================

ZernikeCoefficients::ZernikeCoefficients() {
	name = "Zernike Coefficients";
	n_features = 72;
	//cout << "Instantiating new " << name << " object." << endl;
}

std::vector<double> ZernikeCoefficients::calculate( ImageMatrix * IN_matrix )
{
	std::vector<double> coeffs;
	std::cout << "calculating " << name << std::endl;
	if( IN_matrix == NULL ) {
		return coeffs;
	}
	coeffs.reserve(n_features-1);
	double temp_vec [72];
	int i;

	long output_size;   // output size is normally 72

	for( i = 0; i < n_features; i++ ) temp_vec[i] = 0;
	IN_matrix->zernike2D(temp_vec, &output_size);
	coeffs.assign( temp_vec, temp_vec + n_features);
	return coeffs;
}

//WNDCHARM_REGISTER_ALGORITHM(ZernikeCoefficients)

//===========================================================================

HaralickTextures::HaralickTextures() {
	name = "Haralick Textures";
	n_features = 28;
	//cout << "Instantiating new " << name << " object." << endl;
}

std::vector<double> HaralickTextures::calculate( ImageMatrix * IN_matrix )
{
	std::vector<double> coeffs;
	std::cout << "calculating " << name << std::endl;
	if( IN_matrix == NULL ) {
		return coeffs;
	}
	coeffs.reserve(n_features-1);
	double temp_vec [28];
	int i;

	for( i = 0; i < n_features; i++ ) temp_vec[i] = 0;
	IN_matrix->HaarlickTexture2D(0,temp_vec); // Note the misspelling
	coeffs.assign( temp_vec, temp_vec + n_features);
	return coeffs;
}

//WNDCHARM_REGISTER_ALGORITHM(HaralickTextures)

//===========================================================================

MultiscaleHistograms::MultiscaleHistograms() {
	name = "Multiscale Histograms";
	n_features = 24;
	//cout << "Instantiating new " << name << " object." << endl;
}

std::vector<double> MultiscaleHistograms::calculate( ImageMatrix * IN_matrix )
{
	std::vector<double> coeffs;
	std::cout << "calculating " << name << std::endl;
	if( IN_matrix == NULL ) {
		return coeffs;
	}
	coeffs.reserve(n_features-1);
	double temp_vec [24];
	int i;

	for( i = 0; i < n_features; i++ ) temp_vec[i] = 0;
	IN_matrix->MultiScaleHistogram(temp_vec);
	coeffs.assign( temp_vec, temp_vec + n_features);
	return coeffs;
}

//WNDCHARM_REGISTER_ALGORITHM(MultiscaleHistograms)

//===========================================================================

TamuraTextures::TamuraTextures() {
	name = "Tamura Textures";
	n_features = 6;
	//cout << "Instantiating new " << name << " object." << endl;
}

std::vector<double> TamuraTextures::calculate( ImageMatrix * IN_matrix )
{
	std::vector<double> coeffs;
	std::cout << "calculating " << name << std::endl;
	if( IN_matrix == NULL ) {
		return coeffs;
	}
	coeffs.reserve(n_features-1);
	double temp_vec [6];
	int i;

	for( i = 0; i < n_features; i++ ) temp_vec[i] = 0;
	IN_matrix->TamuraTexture2D(temp_vec);
	coeffs.assign( temp_vec, temp_vec + n_features);
	return coeffs;
}

//WNDCHARM_REGISTER_ALGORITHM(TamuraTextures)

//===========================================================================

CombFirstFourMoments::CombFirstFourMoments() {
	name = "Comb Moments";
	n_features = 48;
	//cout << "Instantiating new " << name << " object." << endl;
}

std::vector<double> CombFirstFourMoments::calculate( ImageMatrix * IN_matrix )
{
	std::vector<double> coeffs;
	std::cout << "calculating " << name << std::endl;
	if( IN_matrix == NULL ) {
		return coeffs;
	}
	coeffs.reserve(n_features-1);
	double temp_vec [48];
	int i;

	for( i = 0; i < n_features; i++ ) temp_vec[i] = 0;
	IN_matrix->CombFirstFourMoments2D(temp_vec);
	coeffs.assign( temp_vec, temp_vec + n_features);
	return coeffs;
}

//WNDCHARM_REGISTER_ALGORITHM(CombFirstFourMoments)

//===========================================================================

RadonCoefficients::RadonCoefficients() {
	name = "Radon Coefficients";
	n_features = 12;
	//cout << "Instantiating new " << name << " object." << endl;
}

std::vector<double> RadonCoefficients::calculate( ImageMatrix * IN_matrix )
{
	std::vector<double> coeffs;
	std::cout << "calculating " << name << std::endl;
	if( IN_matrix == NULL ) {
		return coeffs;
	}
	coeffs.reserve(n_features-1);
	double temp_vec [12];
	int i;

	for( i = 0; i < n_features; i++ ) temp_vec[i] = 0;
	IN_matrix->RadonTransform2D(temp_vec);
	coeffs.assign( temp_vec, temp_vec + n_features);
	return coeffs;
}

//WNDCHARM_REGISTER_ALGORITHM(RadonCoefficients)

//===========================================================================
/* fractal 
   brownian fractal analysis 
   bins - the maximal order of the fractal
   output - array of the size k
   the code is based on: CM Wu, YC Chen and KS Hsieh, Texture features for classification of ultrasonic liver images, IEEE Trans Med Imag 11 (1992) (2), pp. 141Ð152.
   method of approaximation of CC Chen, JS Daponte and MD Fox, Fractal feature analysis and classification in medical imaging, IEEE Trans Med Imag 8 (1989) (2), pp. 133Ð142.
*/
FractalFeatures::FractalFeatures() {
	name = "Fractal Features";
	n_features = 20;
	//cout << "Instantiating new " << name << " object." << endl;
}

std::vector<double> FractalFeatures::calculate( ImageMatrix * IN_matrix )
{
	std::vector<double> coeffs;
	std::cout << "calculating " << name << std::endl;
	if( IN_matrix == NULL ) {
		return coeffs;
	}
	coeffs.reserve(n_features-1);
	double temp_vec [20];
	int i;

	for( i = 0; i < n_features; i++ ) temp_vec[i] = 0;

	int bins = n_features;
	int width = IN_matrix->width;
	int height = IN_matrix->height;
	int x, y, k, bin = 0;
	int K = ( ( width > height ) ? height : width) / 5; // MIN
	int step = (long) floor ( K / bins );
	if( step < 1 )
		step = 1;   // avoid an infinite loop if the image is small
	for( k = 1; k < K; k = k + step )
	{  
		double sum = 0.0;
		for( x = 0; x < width; x++ )
			for( y = 0; y < height - k; y++ )
				sum += fabs( IN_matrix->pixel( x, y, 0 ).intensity - IN_matrix->pixel( x, y+k, 0 ).intensity );
		for( x = 0; x < width - k; x++ )
			for( y = 0; y < height; y++ )
				sum += fabs( IN_matrix->pixel( x, y, 0 ).intensity - IN_matrix->pixel( x + k, y, 0 ).intensity );
		if( bin < bins )
			temp_vec[ bin++ ] = sum / ( width * ( width - k ) + height * ( height - k ) );    
	}

	coeffs.assign( temp_vec, temp_vec + n_features);
	return coeffs;
}

//WNDCHARM_REGISTER_ALGORITHM(FractalFeatures)

//===========================================================================

PixelIntensityStatistics::PixelIntensityStatistics() {
	name = "Pixel Intensity Statistics";
	n_features = 5;
	//cout << "Instantiating new " << name << " object." << endl;
}

std::vector<double> PixelIntensityStatistics::calculate( ImageMatrix * IN_matrix )
{
	std::vector<double> coeffs;
	std::cout << "calculating " << name << std::endl;
	if( IN_matrix == NULL ) {
		return coeffs;
	}
	coeffs.reserve(n_features-1);

	double temp_vec[5];
	int j;
	for( j = 0; j < n_features; j++ ) temp_vec[j] = 0;

	IN_matrix->BasicStatistics(&temp_vec[0], &temp_vec[1], &temp_vec[2], &temp_vec[3], &temp_vec[4], NULL, 10);

	coeffs.assign( temp_vec, temp_vec + n_features);
	return coeffs;
}

//WNDCHARM_REGISTER_ALGORITHM(PixelIntensityStatistics)
        
//===========================================================================

EdgeFeatures::EdgeFeatures() {
	name = "Edge Features";
	n_features = 28;
	//cout << "Instantiating new " << name << " object." << endl;
}

std::vector<double> EdgeFeatures::calculate( ImageMatrix * IN_matrix )
{
	std::vector<double> coeffs;
	std::cout << "calculating " << name << std::endl;
	if( IN_matrix == NULL ) {
		return coeffs;
	}
	coeffs.reserve(n_features-1);

	long EdgeArea = 0;
	double MagMean=0, MagMedian=0, MagVar=0, MagHist[8]={0,0,0,0,0,0,0,0}, DirecMean=0, DirecMedian=0, DirecVar=0, DirecHist[8]={0,0,0,0,0,0,0,0}, DirecHomogeneity=0, DiffDirecHist[4]={0,0,0,0};

	IN_matrix->EdgeStatistics(&EdgeArea, &MagMean, &MagMedian, &MagVar, MagHist, &DirecMean, &DirecMedian, &DirecVar, DirecHist, &DirecHomogeneity, DiffDirecHist, 8);

	int j;

	double temp_vec[28];
	for( j = 0; j < n_features; j++ ) temp_vec[j] = 0;

	double * here = &temp_vec[0];
	*here = double( EdgeArea ); here++;

	for( j=0; j<4; j++ ){
		*here = DiffDirecHist[j]; here++;
	}
	for( j=0; j<8; j++ ){
		*here = DirecHist[j]; here++;
	}

	*here = DirecHomogeneity; here++;
	*here = DirecMean; here++;
	*here = DirecMedian; here++;
	*here = DirecVar; here++;

	for( j=0; j<8; j++ ){
		*here = MagHist[j]; here++;
	}

	*here = MagMean; here++;
	*here = MagMedian; here++;
	*here = MagVar; here++;

	coeffs.assign( temp_vec, temp_vec + n_features);
	return coeffs;
}

//WNDCHARM_REGISTER_ALGORITHM(EdgeFeatures)

//===========================================================================

ObjectFeatures::ObjectFeatures() {
	name = "Object Features";
	n_features = 34;
	//cout << "Instantiating new " << name << " object." << endl;
}

std::vector<double> ObjectFeatures::calculate( ImageMatrix * IN_matrix )
{
	std::vector<double> coeffs;
	std::cout << "calculating " << name << std::endl;
	if( IN_matrix == NULL ) {
		return coeffs;
	}
	coeffs.reserve(n_features-1);

	int feature_count=0, Euler=0, AreaMin=0, AreaMax=0, AreaMedian=0,
			area_histogram[10]={0,0,0,0,0,0,0,0,0,0},
			dist_histogram[10]={0,0,0,0,0,0,0,0,0,0};

	double centroid_x=0, centroid_y=0, AreaMean=0, AreaVar=0, DistMin=0,
				 DistMax=0, DistMean=0, DistMedian=0, DistVar=0;

	IN_matrix->FeatureStatistics(&feature_count, &Euler, &centroid_x, &centroid_y,
			NULL, &AreaMin, &AreaMax, &AreaMean, &AreaMedian,
			&AreaVar, area_histogram, &DistMin, &DistMax,
			&DistMean, &DistMedian, &DistVar, dist_histogram, 10);

	double temp_vec[34];
	int j = 0;

	for( j = 0; j < n_features; j++ ) temp_vec[j] = 0;

	double * here = &temp_vec[0];

	for( j = 0; j < 10; j++ ){
		*here = area_histogram[j]; here++;
	}

	*here = AreaMax; here++;
	*here = AreaMean; here++;
	*here = AreaMedian; here++;
	*here = AreaMin; here++;
	*here = AreaVar; here++;
	*here = centroid_x; here++;
	*here = centroid_y; here++;
	*here = feature_count; here++;

	for( j = 0; j < 10; j++ ) {
		*here = dist_histogram[j]; here++;
	}

	*here = DistMax; here++;
	*here = DistMean; here++;
	*here = DistMedian; here++;
	*here = DistMin; here++;
	*here = DistVar; here++;
	*here = Euler; here++;

	coeffs.assign( temp_vec, temp_vec + n_features);
	return coeffs;
}

//WNDCHARM_REGISTER_ALGORITHM(ObjectFeatures)

//===========================================================================

GaborTextures::GaborTextures() {
	name = "Gabor Textures";
	n_features = 7;
	//cout << "Instantiating new " << name << " object." << endl;
}

std::vector<double> GaborTextures::calculate( ImageMatrix * IN_matrix )
{
	std::vector<double> coeffs;
	std::cout << "calculating " << name << std::endl;
	if( IN_matrix == NULL ) {
		return coeffs;
	}
	coeffs.reserve(n_features-1);

	double temp_vec [7];
	int i;

	for( i = 0; i < n_features; i++ ) temp_vec[i] = 0;
	IN_matrix->GaborFilters2D(temp_vec);
	coeffs.assign( temp_vec, temp_vec + n_features);
	return coeffs;
}

//WNDCHARM_REGISTER_ALGORITHM(GaborTextures)

//===========================================================================

/* gini
   compute the gini coefficient
   
   paper reference: Roberto G. Abraham, Sidney van den Bergh, Preethi Nair, A NEW APPROACH TO GALAXY MORPHOLOGY. I. ANALYSIS OF THE SLOAN DIGITAL SKY
        SURVEY EARLY DATA RELEASE, The Astrophysical Journal, vol. 588, p. 218-229, 2003.
*/
GiniCoefficient::GiniCoefficient() {
	name = "Gini Coefficient";
	n_features = 1;
	//cout << "Instantiating new " << name << " object." << endl;
}

std::vector<double> GiniCoefficient::calculate( ImageMatrix * IN_matrix )
{
	std::vector<double> coeffs;
	std::cout << "calculating " << name << std::endl;
	if( IN_matrix == NULL ) {
		return coeffs;
	}
	coeffs.reserve(n_features-1);

	double temp_vec [1];
	int j;
	for( j = 0; j < n_features; j++ ) temp_vec[j] = 0;

	long pixel_index, num_pixels;
	double *pixels, mean = 0.0, g = 0.0;
	long i, count = 0;

	num_pixels = IN_matrix->height * IN_matrix->width * IN_matrix->depth;
	pixels = new double[ num_pixels ];

	for( pixel_index = 0; pixel_index < num_pixels; pixel_index++ )
		if( IN_matrix->data[ pixel_index ].intensity > 0 )
		{
			pixels[ count ] = IN_matrix->data[ pixel_index ].intensity;
			mean += IN_matrix->data[ pixel_index ].intensity;
			count++;
		}
	if( count > 0 )
		mean = mean / count;
	qsort( pixels, count, sizeof(double), compare_doubles );

	for( i = 1; i <= count; i++)
		g += (2. * i - count - 1.) * pixels[i-1];
	delete [] pixels;

	if( count <= 1 || mean <= 0.0 )
		temp_vec[0] = 0.0;   // avoid division by zero
	else
		temp_vec[0] = g / ( mean * count * ( count-1 ) );

	coeffs.assign( temp_vec, temp_vec + n_features);
	return coeffs;
}

//WNDCHARM_REGISTER_ALGORITHM(GiniCoefficient)

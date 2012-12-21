
#include <iostream> // used for debug output from instantiator methods
#include <cmath>
#include "cmatrix.h"
#include "FeatureTransforms.h"
#include "colors/FuzzyCalc.h" // for definition of compiler constant COLORS_NUM
#include "transforms/fft/bcb_fftw3/fftw3.h"

void Transform::print_info() {

}

EmptyTransform::EmptyTransform () {
	 Transform::name = "Empty";
};


ImageMatrix* EmptyTransform::transform( ImageMatrix * matrix_IN ) {
	std::cout << "Empty transform." << std::endl;
	return NULL;
}

//===========================================================================

FourierTransform::FourierTransform () {
	 Transform::name = "Fourier";
}

/* fft 2 dimensional transform */
// http://www.fftw.org/doc/
//TODO: The ImageMatrix::duplicate() function should really be const
ImageMatrix* FourierTransform::transform( ImageMatrix * matrix_IN ) {
	if( !matrix_IN )
		return NULL;
	
	std::cout << "Performing transform " << name << std::endl;
	ImageMatrix* matrix_OUT = new ImageMatrix(*matrix_IN);
	matrix_OUT->fft2();
	return matrix_OUT;
}



//WNDCHARM_REGISTER_TRANSFORM(FourierTransform)

//===========================================================================

ChebyshevTransform::ChebyshevTransform () {
	 Transform::name = "Chebyshev";
}


ImageMatrix* ChebyshevTransform::transform( ImageMatrix * matrix_IN ) {
	if( !matrix_IN )
		return NULL;	
	std::cout << "Performing transform " << name << std::endl;
	ImageMatrix* matrix_OUT = new ImageMatrix(*matrix_IN);

	matrix_OUT->ChebyshevTransform(0);
	return matrix_OUT;
}

//WNDCHARM_REGISTER_TRANSFORM(ChebyshevTransform)

//===========================================================================

WaveletTransform::WaveletTransform () {
	 Transform::name = "Wavelet";
};


ImageMatrix* WaveletTransform::transform( ImageMatrix * matrix_IN ) {
	if( !matrix_IN )
		return NULL;
	
	std::cout << "Performing transform " << name << std::endl;
	ImageMatrix* matrix_OUT = new ImageMatrix(*matrix_IN);
	matrix_OUT->Symlet5Transform();
	return matrix_OUT;
}

//WNDCHARM_REGISTER_TRANSFORM(WaveletTransform)

//===========================================================================

EdgeTransform::EdgeTransform () {
	 Transform::name = "Edge";
}


ImageMatrix* EdgeTransform::transform( ImageMatrix * matrix_IN ) {
	if( !matrix_IN )
		return NULL;
	
	std::cout << "Performing transform " << name << std::endl;
	ImageMatrix* matrix_OUT = new ImageMatrix(*matrix_IN);
	matrix_OUT->EdgeTransform();
	return matrix_OUT;
}

//WNDCHARM_REGISTER_TRANSFORM(EdgeTransform)

//===========================================================================

ColorTransform::ColorTransform () {
	 Transform::name = "Color";
}


ImageMatrix* ColorTransform::transform( ImageMatrix * matrix_IN ) {
	if( !matrix_IN )
		return NULL;
	
	std::cout << "Performing transform " << name << std::endl;
	double temp_vec [COLORS_NUM+1];

	ImageMatrix* matrix_OUT = new ImageMatrix(*matrix_IN);
	matrix_OUT->ColorTransform(temp_vec, 0);
	histogram_vals.assign(temp_vec, temp_vec+COLORS_NUM+1);
	return matrix_OUT;
}

//WNDCHARM_REGISTER_TRANSFORM(ColorTransform)

//===========================================================================

HueTransform::HueTransform () {
	 Transform::name = "Hue";
}


ImageMatrix* HueTransform::transform( ImageMatrix * matrix_IN ) {
	if( !matrix_IN )
		return NULL;
	
	std::cout << "Performing transform " << name << std::endl;
	ImageMatrix* matrix_OUT = new ImageMatrix(*matrix_IN);
	matrix_OUT->ColorTransform(NULL,1);
	return matrix_OUT;
}

//WNDCHARM_REGISTER_TRANSFORM(HueTransform)


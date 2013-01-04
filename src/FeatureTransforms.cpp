
#include <iostream> // used for debug output from instantiator methods
#include <cmath>
#include "cmatrix.h"
#include "FeatureTransforms.h"
#include "colors/FuzzyCalc.h" // for definition of compiler constant COLORS_NUM
#include "transforms/fft/bcb_fftw3/fftw3.h"

void Transform::print_info() {

}

SharedImageMatrix *Transform::getOutputIM ( const SharedImageMatrix * matrix_IN ) {
	SharedImageMatrix* matrix_OUT = new SharedImageMatrix;
	matrix_OUT->copy(*matrix_IN);
	matrix_OUT->source += "_" + name;
	return matrix_OUT;
};


EmptyTransform::EmptyTransform () {
	 Transform::name = "Empty";
};


SharedImageMatrix* EmptyTransform::transform( const SharedImageMatrix * matrix_IN ) {
	std::cout << "Empty transform." << std::endl;
	SharedImageMatrix* matrix_OUT = getOutputIM (matrix_IN);
	return matrix_OUT;
}

//===========================================================================

FourierTransform::FourierTransform () {
	 Transform::name = "Fourier";
}

/* fft 2 dimensional transform */
// http://www.fftw.org/doc/
SharedImageMatrix* FourierTransform::transform( const SharedImageMatrix * matrix_IN ) {
	if( !matrix_IN )
		return NULL;
	
	std::cout << "Performing transform " << name << std::endl;
	SharedImageMatrix* matrix_OUT = getOutputIM (matrix_IN);
	matrix_OUT->fft2();
	return matrix_OUT;
}



//WNDCHARM_REGISTER_TRANSFORM(FourierTransform)

//===========================================================================

ChebyshevTransform::ChebyshevTransform () {
	 Transform::name = "Chebyshev";
}


SharedImageMatrix* ChebyshevTransform::transform( const SharedImageMatrix * matrix_IN ) {
	if( !matrix_IN )
		return NULL;	
	std::cout << "Performing transform " << name << std::endl;
	SharedImageMatrix* matrix_OUT = getOutputIM (matrix_IN);

	matrix_OUT->ChebyshevTransform(0);
	return matrix_OUT;
}

//WNDCHARM_REGISTER_TRANSFORM(ChebyshevTransform)

//===========================================================================

WaveletTransform::WaveletTransform () {
	 Transform::name = "Wavelet";
};


SharedImageMatrix* WaveletTransform::transform( const SharedImageMatrix * matrix_IN ) {
	if( !matrix_IN )
		return NULL;
	
	std::cout << "Performing transform " << name << std::endl;
	SharedImageMatrix* matrix_OUT = getOutputIM (matrix_IN);
	matrix_OUT->Symlet5Transform();
	return matrix_OUT;
}

//WNDCHARM_REGISTER_TRANSFORM(WaveletTransform)

//===========================================================================

EdgeTransform::EdgeTransform () {
	 Transform::name = "Edge";
}


SharedImageMatrix* EdgeTransform::transform( const SharedImageMatrix * matrix_IN ) {
	if( !matrix_IN )
		return NULL;
	
	std::cout << "Performing transform " << name << std::endl;
	SharedImageMatrix* matrix_OUT = getOutputIM (matrix_IN);
	matrix_OUT->EdgeTransform();
	return matrix_OUT;
}

//WNDCHARM_REGISTER_TRANSFORM(EdgeTransform)

//===========================================================================

ColorTransform::ColorTransform () {
	 Transform::name = "Color";
}


SharedImageMatrix* ColorTransform::transform( const SharedImageMatrix * matrix_IN ) {
	if( !matrix_IN )
		return NULL;
	
	std::cout << "Performing transform " << name << std::endl;
	double temp_vec [COLORS_NUM+1];

	SharedImageMatrix* matrix_OUT = getOutputIM (matrix_IN);
	matrix_OUT->ColorTransform(temp_vec, 0);
	histogram_vals.assign(temp_vec, temp_vec+COLORS_NUM+1);
	return matrix_OUT;
}

//WNDCHARM_REGISTER_TRANSFORM(ColorTransform)

//===========================================================================

HueTransform::HueTransform () {
	 Transform::name = "Hue";
}


SharedImageMatrix* HueTransform::transform( const SharedImageMatrix * matrix_IN ) {
	if( !matrix_IN )
		return NULL;
	
	std::cout << "Performing transform " << name << std::endl;
	SharedImageMatrix* matrix_OUT = getOutputIM (matrix_IN);
	matrix_OUT->ColorTransform(NULL,1);
	return matrix_OUT;
}

//WNDCHARM_REGISTER_TRANSFORM(HueTransform)


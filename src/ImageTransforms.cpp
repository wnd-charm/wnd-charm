#include <iostream> // used for debug output from instantiator methods
#include <cmath>
#include <fcntl.h>

#include "cmatrix.h"
#include "ImageTransforms.h"
#include "colors/FuzzyCalc.h" // for definition of compiler constant COLORS_NUM
#include "transforms/fft/bcb_fftw3/fftw3.h"

/* global variable */
extern int verbosity;

void ImageTransform::print_info() {

}


//===========================================================================

void EmptyTransform::execute( const ImageMatrix * matrix_IN, ImageMatrix * matrix_OUT ) const {
	if( !( matrix_IN && matrix_OUT) )
		return;
	matrix_OUT->copy (*matrix_IN);
	if (verbosity > 3) std::cout << name << " transform." << std::endl;
}

//===========================================================================

FourierTransform::FourierTransform () : ImageTransform ("Fourier") {};

/* fft 2 dimensional transform */
// http://www.fftw.org/doc/
void FourierTransform::execute( const ImageMatrix * matrix_IN, ImageMatrix * matrix_OUT ) const {
	if( !( matrix_IN && matrix_OUT) )
		return;
	
	matrix_OUT->copy (*matrix_IN);
	if (verbosity > 3) std::cout << "Performing transform " << name << std::endl;
	matrix_OUT->fft2();
}



//WNDCHARM_REGISTER_TRANSFORM(FourierTransform)

//===========================================================================

ChebyshevTransform::ChebyshevTransform () : ImageTransform ("Chebyshev") {};

void ChebyshevTransform::execute( const ImageMatrix * matrix_IN, ImageMatrix * matrix_OUT ) const {
	if( !( matrix_IN && matrix_OUT) )
		return;	
	matrix_OUT->copy (*matrix_IN);
	if (verbosity > 3) std::cout << "Performing transform " << name << std::endl;

	matrix_OUT->ChebyshevTransform(0);
}

//WNDCHARM_REGISTER_TRANSFORM(ChebyshevTransform)

//===========================================================================

WaveletTransform::WaveletTransform () : ImageTransform ("Wavelet") {};

void WaveletTransform::execute( const ImageMatrix * matrix_IN, ImageMatrix * matrix_OUT ) const {
	if( !( matrix_IN && matrix_OUT) )
		return;
	
	matrix_OUT->copy (*matrix_IN);
	if (verbosity > 3) std::cout << "Performing transform " << name << std::endl;
	matrix_OUT->Symlet5Transform();
}

//WNDCHARM_REGISTER_TRANSFORM(WaveletTransform)

//===========================================================================

EdgeTransform::EdgeTransform () : ImageTransform ("Edge") {};

void EdgeTransform::execute( const ImageMatrix * matrix_IN, ImageMatrix * matrix_OUT ) const {
	if( !( matrix_IN && matrix_OUT) )
		return;
	
	matrix_OUT->copy (*matrix_IN);
	if (verbosity > 3) std::cout << "Performing transform " << name << std::endl;
	matrix_OUT->EdgeTransform();
}

//WNDCHARM_REGISTER_TRANSFORM(EdgeTransform)

//===========================================================================

ColorTransform::ColorTransform () : ImageTransform ("Color") {};

void ColorTransform::execute( const ImageMatrix * matrix_IN, ImageMatrix * matrix_OUT ) const {
	if( !( matrix_IN && matrix_OUT) )
		return;
	
	matrix_OUT->copy (*matrix_IN);
	if (verbosity > 3) std::cout << "Performing transform " << name << std::endl;
	matrix_OUT->ColorTransform();
}

//WNDCHARM_REGISTER_TRANSFORM(ColorTransform)

//===========================================================================

HueTransform::HueTransform () : ImageTransform ("Hue") {};

void HueTransform::execute( const ImageMatrix * matrix_IN, ImageMatrix * matrix_OUT ) const {
	if( !( matrix_IN && matrix_OUT) )
		return;
	
	matrix_OUT->copy (*matrix_IN);
	if (verbosity > 3) std::cout << "Performing transform " << name << std::endl;
	matrix_OUT->HueTransform();
}

//WNDCHARM_REGISTER_TRANSFORM(HueTransform)


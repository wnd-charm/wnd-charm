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
#include <iostream> // used for debug output from instantiator methods
#include <cmath>
#include <fcntl.h>

#include "cmatrix.h"
#include "FeatureNames.h"
#include "ImageTransforms.h"
#include "colors/FuzzyCalc.h" // for definition of compiler constant COLORS_NUM
#include "transforms/fft/bcb_fftw3/fftw3.h"

/* global variable */
extern int verbosity;


bool ImageTransform::register_task() const {
	return (FeatureNames::registerImageTransform (this));
}

//===========================================================================

void EmptyTransform::execute (const ImageMatrix &matrix_IN, ImageMatrix &matrix_OUT ) const {
	if (verbosity > 3) std::cout << name << " transform." << std::endl;
	matrix_OUT.copy (matrix_IN);
	matrix_OUT.finish();
}

//===========================================================================

FourierTransform::FourierTransform () : ImageTransform ("Fourier") {};

/* fft 2 dimensional transform */
// http://www.fftw.org/doc/
void FourierTransform::execute (const ImageMatrix &matrix_IN, ImageMatrix &matrix_OUT ) const {
	if (verbosity > 3) std::cout << "Performing transform " << name << std::endl;
	matrix_OUT.fft2(matrix_IN);
	matrix_OUT.finish();
}

// Register a static instance of the class using a global bool
static bool FourierTransformReg = ComputationTaskInstances::add (new FourierTransform);


//===========================================================================

ChebyshevTransform::ChebyshevTransform () : ImageTransform ("Chebyshev") {};

void ChebyshevTransform::execute (const ImageMatrix &matrix_IN, ImageMatrix &matrix_OUT ) const {
	if (verbosity > 3) std::cout << "Performing transform " << name << std::endl;
	matrix_OUT.ChebyshevTransform(matrix_IN, 0);
	matrix_OUT.finish();
}

// Register a static instance of the class using a global bool
static bool ChebyshevTransformReg = ComputationTaskInstances::add (new ChebyshevTransform);


//===========================================================================

WaveletTransform::WaveletTransform () : ImageTransform ("Wavelet") {};

void WaveletTransform::execute (const ImageMatrix &matrix_IN, ImageMatrix &matrix_OUT ) const {
	if (verbosity > 3) std::cout << "Performing transform " << name << std::endl;
	matrix_OUT.Symlet5Transform(matrix_IN);
	matrix_OUT.finish();
}

// Register a static instance of the class using a global bool
static bool WaveletTransformReg = ComputationTaskInstances::add (new WaveletTransform);


//===========================================================================

EdgeTransform::EdgeTransform () : ImageTransform ("Edge") {};

void EdgeTransform::execute (const ImageMatrix &matrix_IN, ImageMatrix &matrix_OUT ) const {
	if (verbosity > 3) std::cout << "Performing transform " << name << std::endl;
	matrix_OUT.EdgeTransform(matrix_IN);
	matrix_OUT.finish();
}

// Register a static instance of the class using a global bool
static bool EdgeTransformReg = ComputationTaskInstances::add (new EdgeTransform);


//===========================================================================

ColorTransform::ColorTransform () : ImageTransform ("Color") {};

void ColorTransform::execute (const ImageMatrix &matrix_IN, ImageMatrix &matrix_OUT ) const {
	if (verbosity > 3) std::cout << "Performing transform " << name << std::endl;
	matrix_OUT.ColorTransform(matrix_IN);
	matrix_OUT.finish();
}

// Register a static instance of the class using a global bool
static bool ColorTransformReg = ComputationTaskInstances::add (new ColorTransform);


//===========================================================================

HueTransform::HueTransform () : ImageTransform ("Hue") {};

void HueTransform::execute (const ImageMatrix &matrix_IN, ImageMatrix &matrix_OUT ) const {
	if (verbosity > 3) std::cout << "Performing transform " << name << std::endl;
	matrix_OUT.HueTransform(matrix_IN);
	matrix_OUT.finish();
}

// Register a static instance of the class using a global bool
static bool HueTransformReg = ComputationTaskInstances::add (new HueTransform);


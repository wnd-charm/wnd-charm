#include <iostream> // used for debug output from instantiator methods
#include <cmath>
#include <fcntl.h>
#include "sha1/sha1.h"
#include "b64/encode.h"
#include "cmatrix.h"
#include "FeatureTransforms.h"
#include "colors/FuzzyCalc.h" // for definition of compiler constant COLORS_NUM
#include "transforms/fft/bcb_fftw3/fftw3.h"


void SharedImageMatrix::SetShmemName ( const std::string &final_op ) {
	SHA1::SHA1 digest;
	uint8_t Message_Digest[SHA1::HashSize];

	for(size_t op = 0; op < operations.size(); op++) {
		digest.Input( (uint8_t *)(operations[op].data()), operations[op].length() );
	}
	if (final_op.length()) digest.Input( (uint8_t *)(final_op.data()), final_op.length() );
	digest.Result(Message_Digest);

	base64::encoder().encode ((char *)Message_Digest, SHA1::HashSize, shmem_name, false);
	std::cout <<         "         SharedImageMatrix::SetShmemName: " << shmem_name << std::endl;
	
}

void SharedImageMatrix::allocate (unsigned int w, unsigned int h) {
// 	std::ios::fmtflags stFlags = std::cout.flags();
// 	int stPrec = std::cout.precision();
// 	char stFill = std::cout.fill();
// 
// 	std::cout << "-------- called SharedImageMatrix::allocate (" << w << "," << h << ") on " << source << ", UID: ";
// 	for (unsigned int i = 0; i < sizeof (sourceUID); i++) cout << hex << std::setfill('0') << setw(2) << (int)sourceUID[i];
// 	std::cout << std::endl;
// 	std::cout.flags(stFlags);
// 	std::cout.precision(stPrec);
// 	std::cout.fill(stFill);

	ImageMatrix::allocate (w, h);
}

SharedImageMatrix *SharedImageMatrix::fromCache ( const std::string &final_op ) {
	std::ios::fmtflags stFlags = std::cout.flags();
	long stPrec = std::cout.precision();
	char stFill = std::cout.fill();

	std::cout <<         "-------- called SharedImageMatrix::fromCache (" << width << "," << height << ") on " << source << ", UID: ";
	for (unsigned int i = 0; i < sizeof (sourceUID); i++) cout << hex << std::setfill('0') << setw(2) << (int)sourceUID[i];
	std::cout << std::endl;
	for(size_t op = 0; op < operations.size(); op++) {
		if (! operations[op].compare(0,5,"Open_") ) {
			std::cout << "             Open_";
			for (size_t i = 0; i < sizeof (sourceUID); i++) cout << hex << std::setfill('0') << setw(2) << (int)((operations[op])[5+i]);
			std::cout << std::dec << std::endl;
		} else {
			std::cout << "             " << operations[op] << std::endl;
		}
	}
	if (final_op.length()) std::cout <<         "      param: " << final_op << std::endl;
	std::cout.flags(stFlags);
	std::cout.precision(stPrec);
	std::cout.fill(stFill);
	SetShmemName( final_op );

	return (NULL);
}

void SharedImageMatrix::Cache ( ) {
	std::ios::fmtflags stFlags = std::cout.flags();
	long stPrec = std::cout.precision();
	char stFill = std::cout.fill();

	std::cout <<         "-------- called SharedImageMatrix::Cache (" << width << "," << height << ") on " << source << ", UID: ";
	for (unsigned int i = 0; i < sizeof (sourceUID); i++) cout << hex << std::setfill('0') << setw(2) << (int)sourceUID[i];
	std::cout << std::endl;
	for(size_t op = 0; op < operations.size(); op++) {
		if (! operations[op].compare(0,5,"Open_") ) {
			std::cout << "             Open_";
			for (size_t i = 0; i < sizeof (sourceUID); i++) cout << hex << std::setfill('0') << setw(2) << (int)((operations[op])[5+i]);
			std::cout << std::dec << std::endl;
		} else {
			std::cout << "             " << operations[op] << std::endl;
		}
	}
	std::cout.flags(stFlags);
	std::cout.precision(stPrec);
	std::cout.fill(stFill);
	SetShmemName( "" );
}

int SharedImageMatrix::OpenImage(char *image_file_name,            // load an image of any supported format
		int downsample, rect *bounding_rect,
		double mean, double stddev) {

	int fildes = open (image_file_name, O_RDONLY);
	if (fildes > -1) {
		setSourceUID (fildes);
		close (fildes);
		operations.push_back ( std::string("Open_").append( (const char *)sourceUID,sizeof(sourceUID) ) );
	}

	// WARNING: Setting operations here seems somewhat brittle
	if (bounding_rect && bounding_rect->x >= 0)
		operations.push_back ( string_format ("Sub_%d_%d_%d_%d",
			bounding_rect->x, bounding_rect->y,
			bounding_rect->x+bounding_rect->w-1, bounding_rect->y+bounding_rect->h-1
		));
	if (downsample>0 && downsample<100)  /* downsample by a given factor */
		operations.push_back ( string_format ("DS_%lf_%lf",((double)downsample)/100.0,((double)downsample)/100.0) );
	if (mean>0)  /* normalize to a given mean and standard deviation */
		operations.push_back ( string_format ("Nstd_%lf_%lf",mean,stddev) );
	SharedImageMatrix *matrix_OUT = fromCache ( "" );
	
	if (!matrix_OUT) {
		operations.clear();
		int ret = ImageMatrix::OpenImage (image_file_name,downsample,bounding_rect,mean,stddev);
		Cache();
		return (ret);
	}
	return (1);

}


void Transform::print_info() {

}

SharedImageMatrix* Transform::transform( SharedImageMatrix * matrix_IN ) {
	SharedImageMatrix *matrix_OUT = matrix_IN->fromCache ( name );
	if (!matrix_OUT) {
		matrix_OUT = execute (matrix_IN);
		matrix_OUT->operations = matrix_IN->operations;
		matrix_OUT->operations.push_back (name);
		matrix_OUT->Cache();
	}

	return (matrix_OUT);
}


SharedImageMatrix *Transform::getOutputIM ( const SharedImageMatrix * matrix_IN ) {
	SharedImageMatrix* matrix_OUT = new SharedImageMatrix;
	matrix_OUT->copy(*matrix_IN);
	return matrix_OUT;
};

EmptyTransform::EmptyTransform () {
	 Transform::name = "Empty";
};


SharedImageMatrix* EmptyTransform::execute( const SharedImageMatrix * matrix_IN ) {
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
SharedImageMatrix* FourierTransform::execute( const SharedImageMatrix * matrix_IN ) {
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


SharedImageMatrix* ChebyshevTransform::execute( const SharedImageMatrix * matrix_IN ) {
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


SharedImageMatrix* WaveletTransform::execute( const SharedImageMatrix * matrix_IN ) {
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


SharedImageMatrix* EdgeTransform::execute( const SharedImageMatrix * matrix_IN ) {
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


SharedImageMatrix* ColorTransform::execute( const SharedImageMatrix * matrix_IN ) {
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


SharedImageMatrix* HueTransform::execute( const SharedImageMatrix * matrix_IN ) {
	if( !matrix_IN )
		return NULL;
	
	std::cout << "Performing transform " << name << std::endl;
	SharedImageMatrix* matrix_OUT = getOutputIM (matrix_IN);
	matrix_OUT->ColorTransform(NULL,1);
	return matrix_OUT;
}

//WNDCHARM_REGISTER_TRANSFORM(HueTransform)


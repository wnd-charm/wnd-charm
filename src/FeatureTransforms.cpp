
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


ImageMatrix* EmptyTransform::transform( ImageMatrix * matrix_IN )
{
	std::cout << "Empty transform." << std::endl;
	return NULL;
}

//===========================================================================

FourierTransform::FourierTransform () {
	 Transform::name = "Fourier";
};

/* fft 2 dimensional transform */
// http://www.fftw.org/doc/
//TODO: The ImageMatrix::duplicate() function should really be const
ImageMatrix* FourierTransform::transform( ImageMatrix * matrix_IN )
{
	if( !matrix_IN )
		return NULL;
	
	std::cout << "Performing transform " << name << std::endl;
	ImageMatrix* matrix_OUT = NULL;

	fftw_complex *out;
   double *in;
   fftw_plan p;
   long x,y,z;
	 matrix_OUT = matrix_IN->duplicate();
	 int width = matrix_OUT->width;
	 int height = matrix_OUT->height;
	 int depth = matrix_OUT->depth;
	 double val = 0;

   in = (double*) fftw_malloc(sizeof(double) * width*height*depth);
   out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * width*height*depth);

   if (depth==1) p=fftw_plan_dft_r2c_2d(width,height,in,out, FFTW_MEASURE /* FFTW_ESTIMATE */);
   else p=fftw_plan_dft_r2c_3d(width,height,depth,in,out, FFTW_MEASURE /* FFTW_ESTIMATE */);

   for (x=0;x<width;x++)
     for (y=0;y<height;y++)
       for (z=0;z<depth;z++)	   
         in[depth*height*x+y*depth+z]=matrix_OUT->pixel(x,y,z).intensity;

   fftw_execute(p); /* execute the transformation (repeat as needed) */

   if( depth == 1 ) {
		 // (to 56 including 56 )
		 long half_height = height / 2 + 1;

      // find the abs and angle
      for( x = 0 ; x < width; x++ ) {
        for( y = 0; y < half_height; y++ ) {
					// sqrt(real(X).^2 + imag(X).^2)
					val = sqrt( pow( out[half_height*x+y][0], 2 ) + pow( out[half_height*x+y][1], 2 ) );
          matrix_OUT->SetInt( x, y, 0, val );
				}
			}
   
     // complete the first column
     for( y = half_height; y < height; y++ )
       matrix_OUT->SetInt( 0, y, 0, matrix_OUT->pixel( 0, height - y, 0 ).intensity );

     // complete the rows
     for( y = half_height; y < height; y++ ) {
			 // 1 because the first column is already completed
       for( x = 1; x < width; x++ )
         matrix_OUT->SetInt( x, y, 0, matrix_OUT->pixel( width - x, height - y, 0 ).intensity );
		 }
   }
   else
   {
		 long half_depth=depth/2+1; 
      /* find the abs and angle */
      for (x=0;x<width;x++)
        for (y=0;y<height;y++)
          for (z=0;z<half_depth;z++)		
            matrix_OUT->SetInt(x,y,z,sqrt(pow(out[height*half_depth*x+half_depth*y+z][0],2)+pow(out[height*half_depth*x+half_depth*y+z][1],2)));    /* sqrt(real(X).^2 + imag(X).^2) */
   
      /* compute the first z */
      for (z=half_depth;z<depth;z++)
        for (y=0;y<height;y++)
          matrix_OUT->SetInt(0,y,z,matrix_OUT->pixel(0,y,depth-z).intensity);

      /* complete the rows */
	  for (z=half_depth;z<depth;z++)
        for (y=1;y<height;y++)   /* 1 because the first column is already completed */
          for (x=0;x<width;x++)   
            matrix_OUT->SetInt(x,y,z,matrix_OUT->pixel(width-x,height-y,depth-z).intensity);	  
   }
   
   fftw_destroy_plan(p);
   fftw_free(in);
   fftw_free(out);

   /* calculate the magnitude and angle */

   return matrix_OUT;
}


//int FourierTransform::transform(ImageMatrix * matrix_IN, ImageMatrix ** matrix_OUT_p)
//{
//	if( !matrix_IN ) {
//		return WC_INPUT_IMAGEMATRIX_NULL;
//
//	ImageMatrix* matrix_OUT = NULL;
//	
//	std::cout << "Performing transform " << name << std::endl;
//  *matrix_OUT_p = matrix_OUT = matrix_IN->duplicate();
//  matrix_OUT->fft2();
//	return WC_NO_ERROR;
//}


//struct FourierTransformRegistrar
//{ 
//  FourierTransformRegistrar()
//  {
//    FeatureNames *phonebook = FeatureNames::get_instance();
//		FourierTransform *tform_instance = new FourierTransform;
//    phonebook->register_transform( tform_instance->name, dynamic_cast<Transform*>( tform_instance ) );
//		std::cout << "Just registered the fourier transform." << std::endl;
//  }
//};
//static FourierTransformRegistrar FourierTransformRegistrar_instance;

//WNDCHARM_REGISTER_TRANSFORM(FourierTransform)

//===========================================================================

ChebyshevTransform::ChebyshevTransform () {
	 Transform::name = "Chebyshev";
};


ImageMatrix* ChebyshevTransform::transform( ImageMatrix * matrix_IN )
{
	if( !matrix_IN )
		return NULL;
	ImageMatrix* matrix_OUT = NULL;
	
	std::cout << "Performing transform " << name << std::endl;
  matrix_OUT = matrix_IN->duplicate();
  matrix_OUT->ChebyshevTransform(0);
	return matrix_OUT;
}

//WNDCHARM_REGISTER_TRANSFORM(ChebyshevTransform)

//===========================================================================

WaveletTransform::WaveletTransform () {
	 Transform::name = "Wavelet";
};


ImageMatrix* WaveletTransform::transform( ImageMatrix * matrix_IN )
{
	if( !matrix_IN )
		return NULL;
	ImageMatrix* matrix_OUT = NULL;
	
	std::cout << "Performing transform " << name << std::endl;
  matrix_OUT = matrix_IN->duplicate();
  matrix_OUT->Symlet5Transform();
	return matrix_OUT;
}

//WNDCHARM_REGISTER_TRANSFORM(WaveletTransform)

//===========================================================================

EdgeTransform::EdgeTransform () {
	 Transform::name = "Edge";
};


ImageMatrix* EdgeTransform::transform( ImageMatrix * matrix_IN )
{
	if( !matrix_IN )
		return NULL;
	ImageMatrix* matrix_OUT = NULL;
	
	std::cout << "Performing transform " << name << std::endl;
  matrix_OUT = matrix_IN->duplicate();
  matrix_OUT->EdgeTransform();
	return matrix_OUT;
}

//WNDCHARM_REGISTER_TRANSFORM(EdgeTransform)

//===========================================================================

ColorTransform::ColorTransform () {
	 Transform::name = "Color";
};


ImageMatrix* ColorTransform::transform( ImageMatrix * matrix_IN )
{
	if( !matrix_IN )
		return NULL;
	ImageMatrix* matrix_OUT = NULL;
	
	std::cout << "Performing transform " << name << std::endl;
	double temp_vec [COLORS_NUM+1];
  matrix_OUT = matrix_IN->duplicate();
  matrix_OUT->ColorTransform(temp_vec, 0);
	histogram_vals.assign(temp_vec, temp_vec+COLORS_NUM+1);
	return matrix_OUT;
}

//WNDCHARM_REGISTER_TRANSFORM(ColorTransform)

//===========================================================================

HueTransform::HueTransform () {
	 Transform::name = "Hue";
};


ImageMatrix* HueTransform::transform( ImageMatrix * matrix_IN )
{
	if( !matrix_IN )
		return NULL;
	ImageMatrix* matrix_OUT = NULL;
	
	std::cout << "Performing transform " << name << std::endl;
  matrix_OUT = matrix_IN->duplicate();
  matrix_OUT->ColorTransform(NULL,1);
	return matrix_OUT;
}

//WNDCHARM_REGISTER_TRANSFORM(HueTransform)


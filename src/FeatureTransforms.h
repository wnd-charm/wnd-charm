#ifndef __TRANSFORMS_H_
#define __TRANSFORMS_H_

#include <string>
#include <iomanip>
#include "wndchrm_error.h"
#include "cmatrix.h"

/*! SharedImageMatrix
* inherits from ImageMatrix to store data in a named mmap for memory sharing b/w processes
*/
class SharedImageMatrix: public ImageMatrix {
	public:
		void SetShmemName( const std::string &final_op );
		void allocate (unsigned int w, unsigned int h) ;
		SharedImageMatrix *fromCache ( const std::string &final_op );
		void Cache ();
		virtual int OpenImage(char *image_file_name,            // load an image of any supported format
			int downsample, rect *bounding_rect,
			double mean, double stddev);

		std::string shmem_name;
};

/*! Transform
 *  defines the interface for all inheriting transform classes
 *  Turns any class that inherits this interface into a singleton
 */
class Transform {
	public:
		SharedImageMatrix* transform( SharedImageMatrix * matrix_IN );
		virtual SharedImageMatrix* execute( const SharedImageMatrix * matrix_IN ) = 0;
		std::string name;
		void print_info();
		SharedImageMatrix* getOutputIM ( const SharedImageMatrix * matrix_IN );
	protected:
		Transform (const std::string &s) { name = s;}
		Transform (const char *s) { name = s;}
		Transform() {};
};

class EmptyTransform : public Transform {
	public:
		EmptyTransform ();
		virtual SharedImageMatrix* execute( const SharedImageMatrix * matrix_IN );
};


class FourierTransform : public Transform {
	public:
		FourierTransform();
		virtual SharedImageMatrix* execute( const SharedImageMatrix * matrix_IN );
};

class ChebyshevTransform: public Transform {
	public:
		ChebyshevTransform();
		virtual SharedImageMatrix* execute( const SharedImageMatrix * matrix_IN );
};

class WaveletTransform : public Transform {
	public:
		WaveletTransform();
		virtual SharedImageMatrix* execute( const SharedImageMatrix * matrix_IN );
};

class EdgeTransform : public Transform {
	public:
		EdgeTransform();
		virtual SharedImageMatrix* execute( const SharedImageMatrix * matrix_IN );
};

class ColorTransform : public Transform {
	public:
		ColorTransform();
		vector<double> histogram_vals;
		virtual SharedImageMatrix* execute( const SharedImageMatrix * matrix_IN );
};

class HueTransform : public Transform {
	public:
		HueTransform();
		virtual SharedImageMatrix* execute( const SharedImageMatrix * matrix_IN );
};

/*	
#define WNDCHARM_REGISTER_TRANSFORM(tform_name) \
struct tform_name##TransformRegistrar \
{ \
  tform_name##TransformRegistrar() \
  { \
    FeatureNames *phonebook = FeatureNames::get_instance(); \
		tform_name *tform_instance = new tform_name; \
		int retval = phonebook->register_transform( tform_instance->name, dynamic_cast<Transform*>( tform_instance ) ); \
  } \
}; \
static tform_name##TransformRegistrar tform_name##TransformRegistrar_instance;
//std::cout << "call to register_transform " << #tform_name << " returned " << retval << std::endl; \
*/
#endif


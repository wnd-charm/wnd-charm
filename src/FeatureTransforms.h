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
		void allocate (unsigned int w, unsigned int h) {
			std::cout << "-------- called SharedImageMatrix::allocate (" << w << "," << h << ") on " << source << ", UID: ";
			for (unsigned int i = 0; i < sizeof (sourceUID); i++) cout << hex << std::setfill('0') << setw(2) << (int)sourceUID[i];
			cout << endl;

			ImageMatrix::allocate (w, h);
		}
};

/*! Transform
 *  defines the interface for all inheriting transform classes
 *  Turns any class that inherits this interface into a singleton
 */
class Transform {
	public:
		virtual SharedImageMatrix* transform( const SharedImageMatrix * matrix_IN ) = 0;
		std::string name;
		void print_info();
		SharedImageMatrix* getOutputIM ( const SharedImageMatrix * matrix_IN );
	protected:
		Transform() {};
};

class EmptyTransform : public Transform {
	public:
		EmptyTransform (std::string &s) { name = s;}
		EmptyTransform (const char *s) { name = s;}
		EmptyTransform ();
		virtual SharedImageMatrix* transform( const SharedImageMatrix * matrix_IN );
};


class FourierTransform : public Transform {
	public:
		FourierTransform();
		virtual SharedImageMatrix* transform( const SharedImageMatrix * matrix_IN );
};

class ChebyshevTransform: public Transform {
	public:
		ChebyshevTransform();
		virtual SharedImageMatrix* transform( const SharedImageMatrix * matrix_IN );
};

class WaveletTransform : public Transform {
	public:
		WaveletTransform();
		virtual SharedImageMatrix* transform( const SharedImageMatrix * matrix_IN );
};

class EdgeTransform : public Transform {
	public:
		EdgeTransform();
		virtual SharedImageMatrix* transform( const SharedImageMatrix * matrix_IN );
};

class ColorTransform : public Transform {
	public:
		ColorTransform();
		vector<double> histogram_vals;
		virtual SharedImageMatrix* transform( const SharedImageMatrix * matrix_IN );
};

class HueTransform : public Transform {
	public:
		HueTransform();
		virtual SharedImageMatrix* transform( const SharedImageMatrix * matrix_IN );
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


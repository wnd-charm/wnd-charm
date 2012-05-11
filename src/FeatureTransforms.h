#ifndef __TRANSFORMS_H_
#define __TRANSFORMS_H_

#include <string>
#include "wndchrm_error.h"
#include "cmatrix.h"

/*! Transform
 *  defines the interface for all inheriting transform classes
 *  Turns any class that inherits this interface into a singleton
 */
class Transform {
	public:
		virtual ImageMatrix* transform( ImageMatrix * matrix_IN ) = 0;
		std::string name;
		void print_info();
	protected:
		Transform() {};
};

class EmptyTransform : public Transform {
	public:
		EmptyTransform (std::string &s) { name = s;}
		EmptyTransform (const char *s) { name = s;}
		EmptyTransform ();
		virtual ImageMatrix* transform( ImageMatrix * matrix_IN );
};


class FourierTransform : public Transform {
	public:
		FourierTransform();
		virtual ImageMatrix* transform( ImageMatrix * matrix_IN );
};

class ChebyshevTransform: public Transform {
	public:
		ChebyshevTransform();
		virtual ImageMatrix* transform( ImageMatrix * matrix_IN );
};

class WaveletTransform : public Transform {
	public:
		WaveletTransform();
		virtual ImageMatrix* transform( ImageMatrix * matrix_IN );
};

class EdgeTransform : public Transform {
	public:
		EdgeTransform();
		virtual ImageMatrix* transform( ImageMatrix * matrix_IN );
};

class ColorTransform : public Transform {
	public:
		ColorTransform();
		vector<double> histogram_vals;
		virtual ImageMatrix* transform( ImageMatrix * matrix_IN );
};

class HueTransform : public Transform {
	public:
		HueTransform();
		virtual ImageMatrix* transform( ImageMatrix * matrix_IN );
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


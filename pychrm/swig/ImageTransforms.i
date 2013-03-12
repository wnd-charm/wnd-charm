%module pychrm
%{
#include "ImageTransforms.h"
%}
%include "std_string.i"
%include "std_vector.i"
// Instantiate templates used by ImageTransforms
// The next two lines fix the "vector of const pointers" issue that appears in SWIG 1.3
// Taken from http://thread.gmane.org/gmane.comp.programming.swig/12098/focus=12100
%traits_swigtype(ImageTransform);
%fragment(SWIG_Traits_frag(ImageTransform));
namespace std {
   %template(ConstImageTransformPtrVector) vector<const ImageTransform *>;
}

%include "ImageTransforms.h"

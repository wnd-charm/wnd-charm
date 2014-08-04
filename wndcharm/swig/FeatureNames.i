%module pychrm

%{
/* Include in the generated wrapper file */
#include "FeatureNames.h"
%}
%include "std_string.i"
%include "std_vector.i"


// The %traits_swigtype and %fragment lines fix the "vector of const pointers" issue that appears in SWIG 1.3
// Taken from http://thread.gmane.org/gmane.comp.programming.swig/12098/focus=12100
%traits_swigtype(ImageTransform);
%fragment(SWIG_Traits_frag(ImageTransform));
namespace std {
   %template(ConstImageTransformPtrVector) vector<const ImageTransform *>;
}

%traits_swigtype(FeatureGroup);
%fragment(SWIG_Traits_frag(FeatureGroup));
namespace std {
   %template(ConstFeatureGroupPtrVector) vector<const FeatureGroup *>;
}

%include "FeatureNames.h"

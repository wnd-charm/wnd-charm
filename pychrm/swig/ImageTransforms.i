%module pychrm
%{
#include "ImageTransforms.h"
%}
%include "std_string.i"
%include "std_vector.i"
// Instantiate templates used by ImageTransforms
namespace std {
   %template(ConstImageTransformPtrVector) vector<const ImageTransform *>;
}

%include "ImageTransforms.h"

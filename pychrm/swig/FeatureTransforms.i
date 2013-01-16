%module pychrm

%{
#include "FeatureTransforms.h"
%}

%include "std_string.i"
%nodefaultctor Transform;
%newobject Transform::transform;
%include "FeatureTransforms.h"

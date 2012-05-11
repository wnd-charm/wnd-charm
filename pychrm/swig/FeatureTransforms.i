%module pychrm

%{
#include "FeatureTransforms.h"
%}

%include "std_string.i"
%nodefaultctor Transform;
%include "FeatureTransforms.h"

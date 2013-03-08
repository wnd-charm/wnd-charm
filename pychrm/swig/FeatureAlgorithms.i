%module pychrm

%{
#include "FeatureAlgorithms.h"
%}

%include "std_string.i"
%include "std_vector.i"
// Instantiate templates used by FeatureAlgorithms
namespace std {
   %template(DoubleVector) vector<double>;
}

%include "FeatureAlgorithms.h"

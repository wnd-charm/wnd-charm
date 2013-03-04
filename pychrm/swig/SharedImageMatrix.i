%module pychrm

%{
#include "SharedImageMatrix.h"
%}

%include "std_string.i"

%newobject SharedImageMatrix::transform;
%include "SharedImageMatrix.h"

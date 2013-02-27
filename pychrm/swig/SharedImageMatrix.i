%module pychrm

%{
#include "SharedImageMatrix.h"
%}
%newobject SharedImageMatrix::transform;
%include "SharedImageMatrix.h"

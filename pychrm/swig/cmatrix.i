%module pychrm

%{
#include "cmatrix.h"
%}
%newobject ImageMatrix::transform;
%include "cmatrix.h"

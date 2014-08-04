%module pychrm

%{
/* Include in the generated wrapper file */
/*typedef unsigned long size_t;*/
#include "Tasks.h"
%}
/* Tell SWIG about size_t */
/*typedef unsigned long size_t;*/

%include "std_string.i"
%include "std_vector.i"


namespace std {
   %template(DoubleVector) vector<double>;
}

// Instantiate templates used by Tasks
// The %traits_swigtype and %fragment lines fix the "vector of const pointers" issue that appears in SWIG 1.3
// Taken from http://thread.gmane.org/gmane.comp.programming.swig/12098/focus=12100

%traits_swigtype(ComputationTask);
%fragment(SWIG_Traits_frag(ComputationTask));
namespace std {
   %template(ConstComputationTaskPtrVector) vector<const ComputationTask *>;
}

%traits_swigtype(ComputationTaskNode);
%fragment(SWIG_Traits_frag(ComputationTaskNode));
namespace std {
   %template(ConstComputationTaskNodePtrVector) vector<const ComputationTaskNode *>;
}

%include "Tasks.h"

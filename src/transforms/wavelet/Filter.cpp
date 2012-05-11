#include <stddef.h>
#include "Filter.h"


Filter::Filter () {
	coeff = NULL; 
	size = 0; 
	firstIndex = 0;
}


Filter::Filter(int size, int filterFirst, double *coeff) { 
	init (size, filterFirst,  coeff); 
}

Filter::Filter (const Filter &filter) {
	coeff=0;
	copy(filter); 
}


Filter::~Filter ()
{
	if (coeff != NULL)
		delete [] coeff;
}

/*---------------------------------------------------------------------------*/

void Filter::init (int filterSize, int filterFirst, double *data)
{
   size = filterSize;
   firstIndex = filterFirst;
   center = -firstIndex;

   coeff = new double [size];
   if (data != NULL) {
     for (int i = 0; i < size; i++)
       coeff[i] = data[i];
   } else {
     for (int i = 0; i < size; i++)
       coeff[i] = 0;
   }
}

/*---------------------------------------------------------------------------*/

void Filter::copy (const Filter& filter)
{
  if (coeff != NULL)
    delete [] coeff;
  init (filter.size, filter.firstIndex, filter.coeff);
}

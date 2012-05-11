#ifndef FILTERSET_H_
#define FILTERSET_H_

//#include <iostream>
//using namespace std;
#include "Filter.h"

class FilterSet {
public:
  
		
  FilterSet (int symmetric, 
	     double *anLow, double *anHigh, int anLowSize, int anLowFirst, int anHighFirst,
	     double *synLow, double *synHigh, int synLowSize, int synLowFirst, int synHighFirst );
  
  ~FilterSet ();

   
  int symmetric;
  Filter *analysisLow, *analysisHigh, *synthesisLow,
    *synthesisHigh;

protected:
  void copy (const FilterSet& filterset);
};


#endif /*FILTERSET_H_*/

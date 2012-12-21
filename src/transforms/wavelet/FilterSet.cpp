#include "FilterSet.h"

FilterSet::FilterSet (int symmetric, 
	     double *anLow, double *anHigh, int anLowSize, int anLowFirst, int anHighFirst,
	     double *synLow, double *synHigh, int synLowSize, int synLowFirst, int synHighFirst  ) : 
             symmetric(symmetric) 
{

  analysisLow = new Filter (anLowSize, anLowFirst, anLow);
  analysisHigh = new Filter (anLowSize, anHighFirst, anHigh);

  synthesisLow = new Filter (synLowSize, synLowFirst, synLow);
  synthesisHigh = new Filter (synLowSize, synHighFirst, synHigh);
   
}

FilterSet::~FilterSet ()
{
  delete analysisLow;
  delete analysisHigh;
  delete synthesisLow;
  delete synthesisHigh;
}


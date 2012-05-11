#ifndef FILTER_H_
#define FILTER_H_

/*---------------------------------------------------------------------------*/
// Baseline Wavelet Transform Coder Construction Kit
//
// Geoff Davis
// gdavis@cs.dartmouth.edu
// http://www.cs.dartmouth.edu/~gdavis
//
// Copyright 1996 Geoff Davis 9/11/96
//
// Permission is granted to use this software for research purposes as
// long as this notice stays attached to this software.
//

//#include <iostream>
//using namespace std;

class Filter {
public:
	int size, firstIndex, center;
	double *coeff;

	Filter();
	Filter(int size, int filterFirst, double *coeff );
	Filter (const Filter &filter);
	~Filter ();

	void init (int size, int filterFirst, double *coeff);
	Filter& operator= (const Filter &filter) {copy(filter); return *this;};
	double& operator[] (int index) {return coeff[index-firstIndex];};

protected:
	void copy (const Filter &filter);
};

#endif /*FILTER_H_*/

/* Wavelet High - high compression/high error
 *  uses  5-level decomposition with the bior3.1
 *  wavelet and removing numbers smaller than 1.0 x 10^-4.
 *  In testing, this resulting in 128 time compression ration
 *  with a max error of 3.03E-05	and an average error of
 *  1.55E-06.
 *
*/

#ifndef WAVELETHIGH_H_
#define WAVELETHIGH_H_

#include "Wavelet.h"
//#include "DataGrid.h"

//#include <iostream>
//using namespace std;

class Symlet5:  public Wavelet
{
public:
	Symlet5(double lim, int steps);
	~Symlet5();

//	void transform2D(DataGrid * data);
//	void transform3D(DataGrid * data);

//	void inverseTransform2D(DataGrid * data);
//	void inverseTransform3D(DataGrid * data);

//	void implementInverse(double * out, int N);

};

#endif /*WAVELETHIGH_H_*/


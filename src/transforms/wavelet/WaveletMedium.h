/* Wavelet Medium - medium compression/mdeium error
 *  uses  5-level decomposition with the bior3.1
 *  wavelet and removing numbers smaller than 1.0 x 10^-5.
 *  In testing, this resulting in 60 time compression ratio
 *  with a max error of 3.54E-06  and an average error of 
 *  2.55E-07
 * 
 */


#ifndef WAVELETMED_H_
#define WAVELETMED_H_

#include "Wavelet.h"
#include "DataGrid.h"

class WaveletMedium :  public Wavelet
{
public:
	WaveletMedium(double lim);
	virtual ~WaveletMedium();
	
protected:
	void transform2D(DataGrid * data);
	void transform3D(DataGrid * data);
	
	void inverseTransform2D(DataGrid * data);
	void inverseTransform3D(DataGrid * data);
};

#endif /*WAVELETMED_H_*/

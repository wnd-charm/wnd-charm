/* Wavelet Medium - medium compression/mdeium error
 *  uses  5-level decomposition with the bior5.5
 *  wavelet and removing numbers smaller than 1.0 x 10^-7.
 *  In testing, this resulting in 25 time compression ratio
 *  with a max error of 7.35E-08  and an average error of 
 *  3.67E-09
 * 
 */


#ifndef WAVELETLOW_H_
#define WAVELETLOW_H_

#include "Wavelet.h"
#include "DataGrid.h"

class WaveletLow : public Wavelet
{
public:
	WaveletLow(double lim);
	virtual ~WaveletLow();

protected:
	void transform2D(DataGrid * data);
	void transform3D(DataGrid * data);
	
	void inverseTransform2D(DataGrid * data);
	void inverseTransform3D(DataGrid * data);
};

#endif /*WAVELETLOW_H_*/

#ifndef WAVELET_H_
#define WAVELET_H_

#include "DataGrid.h"
#include "Filter.h"
#include "FilterSet.h"
//#include "DataGrid.h"
#include "Common.h"
//#include "RemoteProcessor.h"

enum TRANSFORM_TYPE {
	   TRANSFORM_TYPE_INVALID = -1,
	   DWT,
	   IDWT,
	   SWT,
       RECONSTRUCTION,
	   TRANSFORM_TYPE_MAX
};

enum SYMMETRY {
	UNKNOWN = -1,
	ASYMMETRIC = 0,
	NEAR_SYMMETRIC = 1,
	SYMMETRIC = 2
};

class Wavelet
{
public:
	Wavelet();
	virtual ~Wavelet();
	void transform(DataGrid * data);
	void inverseTransform(DataGrid * data);

	
	int npad;
	int nsteps;
	double limit;
	int mode;
	
	FilterSet * filterset;
	Filter *analysisLow, *analysisHigh;    // H and G
	Filter *synthesisLow, *synthesisHigh;  // H~ and G~	
	
	int dec_len;	// length of decomposition filter
	int rec_len;	// length of reconstruction filter

	int dec_hi_offset;		// usually 0, but some filters are shifted in time
	int dec_lo_offset;
	int rec_hi_offset;		// - || -
	int rec_lo_offset;		// - || -

	int vanishing_moments_psi;
	int vanishing_moments_phi;
	int support_width;

	int symmetry:3;

	int orthogonal:1;
	int biorthogonal:1;
	int orthonormal:1;

	int compact_support:1;

	int _builtin:1;

	const char* family_name;
	const char* short_name;

//protected:
	void transform1D(DataGrid * data);
	void transform2D(DataGrid * data);
	void transform3D(DataGrid * data);
	
	void inverseTransform1D(DataGrid * data);
//	void inverseTransform2D(DataGrid * data);
//	void inverseTransform3D(DataGrid * data);

//	RemoteProcessor * remote;

};


#endif /*WAVELET_H_*/

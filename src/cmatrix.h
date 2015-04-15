/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                               */
/* Copyright (C) 2007 Open Microscopy Environment                                */
/*       Massachusetts Institue of Technology,                                   */
/*       National Institutes of Health,                                          */
/*       University of Dundee                                                    */
/*                                                                               */
/*                                                                               */
/*                                                                               */
/*    This library is free software; you can redistribute it and/or              */
/*    modify it under the terms of the GNU Lesser General Public                 */
/*    License as published by the Free Software Foundation; either               */
/*    version 2.1 of the License, or (at your option) any later version.         */
/*                                                                               */
/*    This library is distributed in the hope that it will be useful,            */
/*    but WITHOUT ANY WARRANTY; without even the implied warranty of             */
/*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU          */
/*    Lesser General Public License for more details.                            */
/*                                                                               */
/*    You should have received a copy of the GNU Lesser General Public           */
/*    License along with this library; if not, write to the Free Software        */
/*    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA  */
/*                                                                               */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                               */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Written by:  Lior Shamir <shamirl [at] mail [dot] nih [dot] gov>              */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


//---------------------------------------------------------------------------
#ifndef cmatrixH
#define cmatrixH
//---------------------------------------------------------------------------

// this turns off runtime assertions if defined
#undef NDEBUG
#include <assert.h>
#include <string> // for source field
#include <vector> // CEC, required for issue 2, see https://github.com/wnd-charm/wnd-charm/issues/2
#include "Eigen/Dense"
#include "colors/FuzzyCalc.h"
#include "statistics/Moments.h"
//#define min(a,b) (((a) < (b)) ? (a) : (b))
//#define max(a,b) (((a) < (b)) ? (b) : (a))


#define INF 10E200
#define EPSILON 10E-20

// Forward declarations
class ImageTransform;

typedef unsigned char byte;
typedef struct {
	byte r,g,b;
} RGBcolor;
typedef struct {
	byte h,s,v;
} HSVcolor;
typedef Eigen::Matrix< HSVcolor, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > MatrixXhsv;
typedef MatrixXhsv clrDataMat;

// the meaning of the color channels is specified by ColorMode, but it uses the HSVcolor structure for storage
// All color modes other than cmGRAY contain color planes as well as intensity planes
enum ColorModes { cmRGB, cmHSV, cmGRAY };


typedef Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > pixDataMat;
typedef Eigen::Map< pixDataMat, Eigen::Aligned > pixDataMap;
typedef Eigen::Map< clrDataMat, Eigen::Aligned > clrDataMap;
typedef pixDataMap pixData;
typedef clrDataMap clrData;

typedef const pixData &readOnlyPixels;
typedef const clrData &readOnlyColors;
typedef pixData &writeablePixels;
typedef clrData &writeableColors;
typedef struct {
	int x,y,w,h;
} rect;

//---------------------------------------------------------------------------
// global functions
#define MIN(a,b) (a<b?a:b)
#define MAX(a,b) (a>b?a:b)

static inline int compare_doubles (const void *a, const void *b) {
	double ret = *((double *)a) - *((double*)b);
	if (ret > 0) return 1;
	else if (ret < 0) return -1;
	else return 0;
}

static inline HSVcolor RGB2HSV(const RGBcolor rgb) {
	double r,g,b,h,max,min,delta;
	HSVcolor hsv;

	r = (double)(rgb.r) / 255;
	g = (double)(rgb.g) / 255;
	b = (double)(rgb.b) / 255;

	max = MAX (r, MAX (g, b)), min = MIN (r, MIN (g, b));
	delta = max - min;

	hsv.v = (byte)(max*240.0);
	if (max != 0.0)
		hsv.s = (byte)((delta / max)*240.0);
	else
		hsv.s = 0;
	if (hsv.s == 0) hsv.h = 0; //-1;
	else {
		h = 0;
		if (r == max)
		h = (g - b) / delta;
		else if (g == max)
		h = 2 + (b - r) / delta;
		else if (b == max)
		h = 4 + (r - g) / delta;
		h *= 60.0;
		if (h >= 360) h -= 360.0;
		if (h < 0.0) h += 360.0;
		hsv.h = (byte)(h *(240.0/360.0));
	}
	return(hsv);
}
static inline RGBcolor HSV2RGB(const HSVcolor hsv) {
	RGBcolor rgb;
	double R=0, G=0, B=0;
	double H, S, V;
	double i, f, p, q, t;

	H = hsv.h;
	S = (double)(hsv.s)/240;
	V = (double)(hsv.v)/240;
	if(S == 0 && H == 0) {R=G=B=V;}  /*if S=0 and H is undefined*/
	H = H*(360.0/240.0);
	if(H == 360) H=0;
	H = H/60;
	i = floor(H);
	f = H-i;
	p = V*(1-S);
	q = V*(1-(S*f));
	t = V*(1-(S*(1-f)));

	if(i==0) {R=V;  G=t;  B=p;}
	if(i==1) {R=q;  G=V;  B=p;}
	if(i==2) {R=p;  G=V;  B=t;}
	if(i==3) {R=p;  G=q;  B=V;}
	if(i==4) {R=t;  G=p;  B=V;}
	if(i==5) {R=V;  G=p;  B=q;}

	rgb.r = (byte)(R*255);
	rgb.g = (byte)(G*255);
	rgb.b = (byte)(B*255);
	return rgb;
}
static inline double RGB2GRAY(const RGBcolor rgb) {
	return((0.2989*rgb.r+0.5870*rgb.g+0.1140*rgb.b));
}

//---------------------------------------------------------------------------

class ImageMatrix {
private:
	pixData _pix_plane;                              // pixel plane data  
	clrData _clr_plane;                              // 3-channel color data
	bool _is_pix_writeable;
	bool _is_clr_writeable;
	double _median;
public:
	std::string source;                             // path of image source file
	enum ColorModes ColorMode;                       // can be cmRGB, cmHSV or cmGRAY
	unsigned short bits;                            // the number of intensity bits (8,16, etc)
	unsigned int width,height;                               // width and height of the picture
	Moments2 stats;        // min, max, mean, std computed in single pass, median in separate pass
	bool has_median;                     // if the median has been computed
	const double *data_ptr() const { return _pix_plane.data(); }
	double *writable_data_ptr() { return _pix_plane.data(); }
	inline writeablePixels WriteablePixels() {
		assert(_is_pix_writeable && "Attempt to write to read-only pixels");
		has_median = false;
		stats.reset();
		return _pix_plane;
	}
	inline writeableColors WriteableColors() {
		assert(_is_pix_writeable && "Attempt to write to read-only pixels");
		return _clr_plane;
	}
	inline readOnlyPixels ReadablePixels() const {
		return _pix_plane;
	}
	inline readOnlyColors ReadableColors() const {
		return _clr_plane;
	}
	inline readOnlyPixels ReadOnlyPixels() const {
		assert(!_is_pix_writeable && "Attempt to read from write-only pixels");
		return _pix_plane;
	}
	inline readOnlyColors ReadOnlyColors() const {
		assert(!_is_pix_writeable && "Attempt to read from write-only pixels");
		return _clr_plane;
	}
	void finish() {
		WriteablePixelsFinish();
		WriteableColorsFinish();
	}
	void WriteablePixelsFinish () {
		_is_pix_writeable = false;
	}
	void WriteableColorsFinish () {
		_is_clr_writeable = false;
	}
	int LoadTIFF(char *filename);                   // load from TIFF file
	int SaveTiff(char *filename);                   // save a matrix in TIF format
	virtual int OpenImage(char *image_file_name,            // load an image of any supported format
		int downsample, rect *bounding_rect,
		double mean, double stddev);
	// constructor helpers
	void init();
	void remap_pix_plane (double *ptr, const unsigned int w, const unsigned int h);
	void remap_clr_plane (HSVcolor *ptr, const unsigned int w, const unsigned int h);
	virtual void allocate (unsigned int w, unsigned int h);
	void copyFields(const ImageMatrix &copy);
	void copyData(const ImageMatrix &copy);
	void copy(const ImageMatrix &copy);
	void submatrix(const ImageMatrix &matrix,
		const unsigned int x1, const unsigned int y1, const unsigned int x2, const unsigned int y2);
	// N.B.: See note in implementation
	ImageMatrix () : _pix_plane (NULL,0,0), _clr_plane (NULL,0,0) {
		init();
	};
	virtual ~ImageMatrix();                                 // destructor

	virtual void transform (const ImageMatrix &matrix_IN, const ImageTransform *transform);

	void normalize(double min, double max, long range, double mean, double stddev); // normalized an image to either min/max or mean/stddev
	void to8bits (const ImageMatrix &matrix_IN);
	void flipV();                                   // flip an image around a vertical axis (left to right)
	void flipH();                                   // flip an image around a horizontal axis (upside down)
	void invert();                                  // invert the intensity of an image
	void Downsample (const ImageMatrix &matrix_IN, double x_ratio, double y_ratio);// down sample an image
	void Rotate (const ImageMatrix &matrix_IN, double angle);              // rotate an image by 90,180,270 degrees
	void convolve(const pixDataMat &filter);
	double update_median ();
	double get_median () const;
	void UpdateStats();
	void GetStats (Moments2 &moments2) const;
	inline double min() {
		if (! stats.n() > 0) UpdateStats();
		return (stats.min());
	}
	inline double max() {
		if (! stats.n() > 0) UpdateStats();
		return (stats.max());
	}
	inline double mean() {
		if (! stats.n() > 0) UpdateStats();
		return (stats.mean());
	}
	inline double std() {
		if (! stats.n() > 0) UpdateStats();
		return (stats.std());
	}
	inline double var() {
		if (! stats.n() > 0) UpdateStats();
		return (stats.var());
	}
	inline double median() {
		if (!has_median) update_median();
		return (_median);
	}
	void GetColorStatistics(double *hue_avg, double *hue_std, double *sat_avg, double *sat_std, double *val_avg, double *val_std, double *max_color, double *colors) const;
	void ColorTransform(const ImageMatrix &matrix_IN);
	void HueTransform(const ImageMatrix &matrix_IN);
	void histogram(double *bins,unsigned short nbins, bool imhist = false, const Moments2 &in_stats = Moments2()) const; // by default, based on computed min and max.
    double Otsu(bool dynamic_range=true) const;                                  /* Otsu grey threshold                  */
	void MultiScaleHistogram(double *out) const;
	//   double AverageEdge();
	void EdgeTransform(const ImageMatrix &matrix_IN);                           // gradient binarized using otsu threshold
	double fft2 (const ImageMatrix &matrix_IN, const ImageMatrix &matrix_OUT);
	void ChebyshevTransform (const ImageMatrix &matrix_IN, unsigned int N);
	void ChebyshevFourierTransform2D (double *coeff) const;
	void Symlet5Transform (const ImageMatrix &matrix_IN);
	void PrewittMagnitude2D (const ImageMatrix &matrix_IN);
	void PrewittDirection2D (const ImageMatrix &matrix_IN);
	void ChebyshevStatistics2D (double *coeff, unsigned int N, unsigned int nbins) const;
	int CombFirstFourMoments2D (std::vector<double> &vec) const;
	void EdgeStatistics (unsigned long *EdgeArea, double *MagMean, double *MagMedian, double *MagVar,
		double *MagHist, double *DirecMean, double *DirecMedian, double *DirecVar, double *DirecHist,
		double *DirecHomogeneity, double *DiffDirecHist, unsigned int nbins) const;
	void RadonTransform2D(double *vec) const;
	double OtsuBinaryMaskTransform (const ImageMatrix &matrix_IN);
	unsigned long BWlabel(int level);
	void centroid(double *x_centroid, double *y_centroid) const;
	void FeatureStatistics(unsigned long *count, long *Euler, double *centroid_x, double *centroid_y, unsigned long *AreaMin, unsigned long *AreaMax,
		double *AreaMean, unsigned int *AreaMedian, double *AreaVar, unsigned int *area_histogram,double *DistMin, double *DistMax,
		double *DistMean, double *DistMedian, double *DistVar, unsigned int *dist_histogram, unsigned int nbins
	) const;
	void GaborFilters2D(double *ratios) const;
	void HaralickTexture2D(double distance, double *out) const;
	void TamuraTexture2D(double *vec) const;
	void zernike2D(double *zvalues, long *output_size) const;

	// disable the copy constructor
private:
    ImageMatrix(const ImageMatrix &matrix) : _pix_plane (NULL,0,0), _clr_plane (NULL,0,0) {
		assert(false && "Attempt to use copy constructor");
	};
};

#ifdef GPU
extern "C" void Chebyshev2D_gpu(const ImageMatrix &Im, double *out, unsigned int N);
extern "C" void gpu_convolve(double *, double *, unsigned long, unsigned long, unsigned long, unsigned long, double *);
extern "C" int FeatureCentroid_gpu(const ImageMatrix &Im, double *sum_areas, double *sum_dists, unsigned long *object_area, double *centroid_dists,
                                unsigned long *count, double *centroid_x, double *centroid_y);
extern "C" void gpu_fft2(ImageMatrix *,const ImageMatrix&);
extern "C" int gpu_mb_zernike2D (const ImageMatrix &Im, double order, double rad, double *zvalues, long *output_size );
#endif
#endif

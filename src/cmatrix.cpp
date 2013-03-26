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

#include <vector>
#include <math.h>
#include <stdio.h>
#include "cmatrix.h"
#include "ImageTransforms.h"
#include "colors/FuzzyCalc.h"
#include "transforms/fft/bcb_fftw3/fftw3.h"
#include "transforms/chebyshev.h"
#include "transforms/ChebyshevFourier.h"
#include "transforms/wavelet/Symlet5.h"
#include "transforms/wavelet/DataGrid2D.h"
#include "transforms/wavelet/DataGrid3D.h"
#include "transforms/radon.h"
#include "statistics/CombFirst4Moments.h"
#include "statistics/FeatureStatistics.h"
#include "textures/gabor.h"
#include "textures/tamura.h"
#include "textures/haralick/haralick.h"
#include "textures/zernike/zernike.h"

#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/types.h> // for dev_t, ino_t
#include <fcntl.h>     // for O_RDONLY

#include <stdlib.h>
#include <string.h>
#include <tiffio.h>


using namespace std;
//-----------------------------------------------------------------------


//-----------------------------------------------------------------------



/* global variable */
extern int verbosity;




/* LoadTIFF
   filename -char *- full path to the image file
*/
int ImageMatrix::LoadTIFF(char *filename) {
	unsigned int h,w,x=0,y=0;
	unsigned short int spp=0,bps=0;
	TIFF *tif = NULL;
	unsigned char *buf8;
	unsigned short *buf16;
	RGBcolor rgb = {0,0,0};
	ImageMatrix R_matrix, G_matrix, B_matrix;
	Moments2 R_stats, G_stats, B_stats;

	TIFFSetWarningHandler(NULL);
	if( (tif = TIFFOpen(filename, "r")) ) {
		source = filename;

		TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
		width = w;
		TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
		height = h;
		TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bps);
		bits=bps;
		if ( ! (bits == 8 || bits == 16) ) return (0); // only 8 and 16-bit images supported.
		TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &spp);
		if (!spp) spp=1;  /* assume one sample per pixel if nothing is specified */
		// regardless of how the image comes in, the stored mode is HSV
		if (spp == 3) {
			ColorMode = cmHSV;
			// If the bits are > 8, we do the read into doubles so that later
			// we can scale the image to its actual signal range.
			if (bits > 8) {
				R_matrix.ColorMode = cmGRAY;
				R_matrix.allocate (width, height);
				G_matrix.ColorMode = cmGRAY;
				G_matrix.allocate (width, height);
				B_matrix.ColorMode = cmGRAY;
				B_matrix.allocate (width, height);
			}
		} else {
			ColorMode = cmGRAY;
		}
		if ( TIFFNumberOfDirectories(tif) > 1) return(0);   /* get the number of slices (Zs) */

		/* allocate the data */
		allocate (width, height);
		writeablePixels pix_plane = WriteablePixels();
		writeableColors clr_plane = WriteableColors();

		/* read TIFF header and determine image size */
		buf8 = (unsigned char *)_TIFFmalloc(TIFFScanlineSize(tif)*spp);
		buf16 = (unsigned short *)_TIFFmalloc( (tsize_t)sizeof(unsigned short)*TIFFScanlineSize(tif)*spp );
		for (y = 0; y < height; y++) {
			int col;
			if (bits==8) TIFFReadScanline(tif, buf8, y);
			else TIFFReadScanline(tif, buf16, y);
			x=0;col=0;
			while (x<width) {
				unsigned char byte_data;
				unsigned short short_data;
				double val=0;
				int sample_index;
				for (sample_index=0;sample_index<spp;sample_index++) {
					byte_data=buf8[col+sample_index];
					short_data=buf16[col+sample_index];
					if (bits==8) val=(double)byte_data;
					else val=(double)(short_data);
					if (spp==3 && bits > 8) {  /* RGB image */
						if (sample_index==0) R_matrix.WriteablePixels()(y,x) = R_stats.add (val);
						if (sample_index==1) G_matrix.WriteablePixels()(y,x) = G_stats.add (val);
						if (sample_index==2) B_matrix.WriteablePixels()(y,x) = B_stats.add (val);
					} else if (spp == 3) {
						if (sample_index==0) rgb.r = (unsigned char)(R_stats.add (val));
						if (sample_index==1) rgb.g = (unsigned char)(G_stats.add (val));
						if (sample_index==2) rgb.b = (unsigned char)(B_stats.add (val));
					}
				}
				if (spp == 3 && bits == 8) {
					clr_plane (y, x) = RGB2HSV(rgb);
				} else if (spp == 1) {
					pix_plane (y, x) = stats.add (val);
				}
				x++;
				col+=spp;
			}
		}
		// Do the conversion to unsigned chars based on the input signal range
		// i.e. scale global RGB min-max to 0-255
		if (spp == 3 && bits > 8) {
			size_t a, num = width*height;
			double RGB_min=0, RGB_max=0, RGB_scale=0;
			R_matrix.WriteablePixelsFinish();
			G_matrix.WriteablePixelsFinish();
			B_matrix.WriteablePixelsFinish();
			// Get the min and max for all 3 channels
			if (R_stats.min() <= G_stats.min() && R_stats.min() <= B_stats.min()) RGB_min = R_stats.min();
			else if (G_stats.min() <= R_stats.min() && G_stats.min() <= B_stats.min()) RGB_min = G_stats.min();
			else if (B_stats.min() <= R_stats.min() && B_stats.min() <= G_stats.min()) RGB_min = B_stats.min();
			if (R_stats.max() >= G_stats.max() && R_stats.max() >= B_stats.max()) RGB_max = R_stats.max();
			else if (G_stats.max() >= R_stats.max() && G_stats.max() >= B_stats.max()) RGB_max = G_stats.max();
			else if (B_stats.max() >= R_stats.max() && B_stats.max() >= G_stats.max()) RGB_max = B_stats.max();
			// Scale the clrData to the global min / max.
			RGB_scale = (255.0/(RGB_max-RGB_min));
			for (a = 0; a < num; a++) {
				rgb.r = (unsigned char)( (R_matrix.ReadablePixels().array().coeff(a) - RGB_min) * RGB_scale);
				rgb.g = (unsigned char)( (G_matrix.ReadablePixels().array().coeff(a) - RGB_min) * RGB_scale);
				rgb.b = (unsigned char)( (B_matrix.ReadablePixels().array().coeff(a) - RGB_min) * RGB_scale);
				clr_plane (y, x) = RGB2HSV(rgb);
			}
		}
		_TIFFfree(buf8);
		_TIFFfree(buf16);
		TIFFClose(tif);

// 		WriteableColorsFinish();
// 		WriteablePixelsFinish();
	} else return(0);

	return(1);
}

/*  SaveTiff
    Save a matrix in TIFF format (16 bits per pixel)
*/
int ImageMatrix::SaveTiff(char *filename) {
	readOnlyPixels pix_plane = ReadablePixels();

	unsigned int x,y;
	TIFF* tif = TIFFOpen(filename, "w");
	if (!tif) return(0);
	unsigned short *BufImage16 = new unsigned short[width*height];
	if (!BufImage16) {
		TIFFClose(tif);
		return (0);
	}
	unsigned char *BufImage8 = new unsigned char[width*height];
	if (!BufImage8) {
		TIFFClose(tif);
		delete BufImage16;
		return (0);
	}


	for (y = 0; y < height; y++)
		for (x = 0; x < width ; x++) {
			if (bits==16) BufImage16[x + (y * width)] = (unsigned short) ( pix_plane (y,x) );
			else BufImage8[x + (y * width)] = (unsigned char) ( pix_plane (y,x) );
		}

	TIFFSetField(tif,TIFFTAG_IMAGEWIDTH, width);
	TIFFSetField(tif,TIFFTAG_IMAGELENGTH, height);
	TIFFSetField(tif, TIFFTAG_PLANARCONFIG,1);
	TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, 1);
	TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, bits);
	TIFFSetField(tif, TIFFTAG_COMPRESSION, 1);
	TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);

	for (y = 0; y < height; y ++) {
		if (bits==16) TIFFWriteScanline (tif, &(BufImage16[y*width]), y,0 );
		else TIFFWriteScanline (tif, &(BufImage8[y*width]), y,0 );
	}

	TIFFClose(tif);
	delete BufImage16;
	delete BufImage8;
	return(1);
}

int ImageMatrix::OpenImage(char *image_file_name, int downsample, rect *bounding_rect, double mean, double stddev) {  
	int res=0;
	if (strstr(image_file_name,".tif") || strstr(image_file_name,".TIF")) {  
		res=LoadTIFF(image_file_name);
	}

	// add the image only if it was loaded properly
	if (res) {
		// compute features only from an area of the image
		if (bounding_rect && bounding_rect->x >= 0) {
			submatrix (*this, (unsigned int)bounding_rect->x, (unsigned int)bounding_rect->y,
				(unsigned int)bounding_rect->x+bounding_rect->w-1, (unsigned int)bounding_rect->y+bounding_rect->h-1
			);
		}
		if (downsample>0 && downsample<100)  /* downsample by a given factor */
			Downsample(*this, ((double)downsample)/100.0,((double)downsample)/100.0);   /* downsample the image */
		if (mean>0)  /* normalize to a given mean and standard deviation */
			normalize(-1,-1,-1,mean,stddev);
	}
	if (! source.length() ) source = image_file_name;
	WriteablePixelsFinish();
	WriteableColorsFinish();
	return(res);
}


// pseudo-constructor helpers (see below)

// This sets default values for the different constructors
void ImageMatrix::init() {
	width      = 0;
	height     = 0;
	_is_pix_writeable = _is_clr_writeable = false;

	stats.reset();
	has_median = false;
	_median    = 0;
	source     = "";
	ColorMode = cmGRAY;
	bits      = 8;
}

// This is using the fancy shmancy placement-new operator
// The object is created new at the same memory location it was before, so no new allocation happens
// This allows us to call Eigen's Map constructor with new parameters without actually re-allocating the object
// N.B.: THis does not do any memory allocation or deallocation, it simply assigns the passed in memory to an Eigen object.
void ImageMatrix::remap_pix_plane(double *ptr, const unsigned int w, const unsigned int h) {
	width  = 0;
	height = 0;
	// N.B. Eigen matrix parameter order is rows, cols, not X, Y
	new (&_pix_plane) pixData(ptr, h, w);
	width  = (unsigned int)_pix_plane.cols();
	height = (unsigned int)_pix_plane.rows();
	// FIXME: Should check here if the pointer is different than what it was.
	// May be possible to do a conservative re-mapping.
	_is_pix_writeable = true;
	has_median = false;
	stats.reset();
}
// Same as above for the color plane.
void ImageMatrix::remap_clr_plane(HSVcolor *ptr, const unsigned int w, const unsigned int h) {
	if (ColorMode == cmGRAY) return;
	width  = 0;
	height = 0;
	// N.B. Eigen matrix parameter order is rows, cols, not X, Y
	new (&_clr_plane) clrData(ptr, h, w);
	width  = (unsigned int)_clr_plane.cols();
	height = (unsigned int)_clr_plane.rows();
	// FIXME: Should check here if the pointer is different than what it was.
	// May be possible to do a conservative re-mapping.
	_is_clr_writeable = true;
}

// If the image are changed size, then reallocate.
// If the image changed color mode, reallocate.
// Ensure that anything that's reallocated is deallocated first.
void ImageMatrix::allocate (unsigned int w, unsigned int h) {

	if ((unsigned int) _pix_plane.cols() != w || (unsigned int)_pix_plane.rows() != h) {
		// These throw exceptions, which we don't catch (catch in main?)
		// FIXME: We could check for shrinkage and simply remap instead of allocating.
		if (verbosity > 7 && _pix_plane.data()) fprintf (stdout, "deallocating grayscale %p ",(void *)_pix_plane.data());
		if (_pix_plane.data()) Eigen::aligned_allocator<double>().deallocate (_pix_plane.data(), _pix_plane.size());
		remap_pix_plane (Eigen::aligned_allocator<double>().allocate (w * h), w, h);
		if (verbosity > 7 && _pix_plane.data()) fprintf (stdout, "allocated grayscale %p\n",(void *)_pix_plane.data());
	}

	// cleanup the color plane if it changed size, or if we have a gray image.
	if ( ColorMode == cmGRAY || (_pix_plane.data() && ((unsigned int)_clr_plane.cols() != w || (unsigned int)_clr_plane.rows() != h)) ) {
		if (verbosity > 7 && _clr_plane.data()) fprintf (stdout, "  deallocating color %p\n",(void *)_clr_plane.data());
		if (_clr_plane.data()) Eigen::aligned_allocator<HSVcolor>().deallocate (_clr_plane.data(), _clr_plane.size());
		remap_clr_plane (NULL, 0, 0);
	}

	// Allocate a new color plane if necessary.
	if (ColorMode != cmGRAY && ! (_clr_plane.data()) ) {
		// These throw exceptions, which we don't catch (catch in main?)
		// FIXME: We could check for shrinkage and simply remap instead of allocating.
		remap_clr_plane (Eigen::aligned_allocator<HSVcolor>().allocate (w * h), w, h);
		if (verbosity > 7 && _clr_plane.data()) fprintf (stdout, "  allocated color %p\n",(void *)_clr_plane.data());
	}
}

void ImageMatrix::copyFields(const ImageMatrix &copy) {
	width = copy.width;
	height = copy.height;
	stats  = copy.stats;
	has_median = copy.has_median;
	_median    = copy._median;
	source     = copy.source;
	ColorMode = copy.ColorMode;
	bits      = copy.bits;
}
void ImageMatrix::copyData(const ImageMatrix &copy) {
	// WriteablePixels() resets the stats
	Moments2 old_stats = copy.stats;
	bool old_has_median = has_median;

	allocate(copy.width, copy.height);
	WriteablePixels() = copy.ReadablePixels();
	if (ColorMode != cmGRAY) {
		WriteableColors() = copy.ReadableColors();
	}
	stats = old_stats;
	has_median = old_has_median;
}
void ImageMatrix::copy(const ImageMatrix &copy) {
	copyFields (copy);
	// Do this stuff last - potentially virtual method calls.
	copyData (copy);
}

void ImageMatrix::submatrix (const ImageMatrix &matrix, const unsigned int x1, const unsigned int y1, const unsigned int x2, const unsigned int y2) {
	unsigned int x0, y0;

	// verify that the image size is OK
	x0 = (x1 < 0 ? 0 : x1);
	y0 = (y1 < 0 ? 0 : y1);
	unsigned int new_width  = (x2 >= matrix.width  ? matrix.width  : x2 - x0 + 1);
	unsigned int new_height = (y2 >= matrix.height ? matrix.height : y2 - y0 + 1);

	copyFields (matrix);
	allocate (new_width, new_height);
	// Copy the Eigen matrixes
	// N.B. Eigen matrix parameter order is rows, cols, not X, Y
	WriteablePixels() = matrix.ReadablePixels().block(y0,x0,height,width);
	if (ColorMode != cmGRAY) {
		WriteableColors() = matrix.ReadableColors().block(y0,x0,height,width);
	}
}

/*
* There is only one simple constructor implemented (no copy constructors or other constructors).
* The reason for this is that a SharedImageMatrix subclass needs to override the allocate() method to make it shareable.
* any fancy constructor would have to use the allocate() method to do anything fancier than the most basic initialization.
* So, what's the problem with that?  Well, you can't use virtual methods (i.e. allocate()) in constructors, that's why.
* So, we are forcing the use pattern of first making a bare-bones ImageMatrix, then calling the pseudo-constructor helpers above
* to make copies, etc.  This way, we can safely sub-class ImageMatrix, and override the allocate() method.
* A side-benefit is that if there are unintended/implicit copies being made, they will now generate a run-time assertion.
* If there are prettier mechanisms around this limitation, feel free.
*/

// ImageMatrix::ImageMatrix() {
// 	std::cout << "-------- called ImageMatrix::ImageMatrix - empty" << std::endl;
// 	init();
// }


/* destructor
   Ensure that the pixels are finalized.
   Allocated pixel memory gets taken care of by the Eigen destructors.
*/
ImageMatrix::~ImageMatrix() {
	WriteablePixelsFinish();
	if (verbosity > 7 && _pix_plane.data()) fprintf (stdout, "deallocating grayscale %p\n",(void *)_pix_plane.data());
	if (_pix_plane.data()) Eigen::aligned_allocator<double>().deallocate (_pix_plane.data(), _pix_plane.size());
	remap_pix_plane (NULL, 0, 0);

	WriteableColorsFinish();
	if (verbosity > 7 && _clr_plane.data()) fprintf (stdout, "deallocating color %p\n",(void *)_clr_plane.data());
	if (_clr_plane.data()) Eigen::aligned_allocator<HSVcolor>().deallocate (_clr_plane.data(), _clr_plane.size());
	remap_clr_plane (NULL, 0, 0);
}

// This is a general transform method that applies the specified transform to the specified ImageMatrix,
// storing the result in the ImageMatrix it was called on
void ImageMatrix::transform (const ImageMatrix &matrix_IN, const ImageTransform *transform) {
	transform->execute (matrix_IN, *this);
}


/* to8bits
   convert an arbitrary-range matrix to an 8 bit range by scaling the signal range to 0.0 to 255.0
*/
void ImageMatrix::to8bits(const ImageMatrix &matrix_IN) {

	copyFields (matrix_IN);
	bits = 8;

	allocate (matrix_IN.width, matrix_IN.height);

	Moments2 local_stats;
	matrix_IN.GetStats (local_stats);
	double min_val = local_stats.min();
	double max_val = local_stats.max();

	double scale255 = (255.0/(max_val - min_val));

	WriteablePixels() = ( (matrix_IN.ReadablePixels().array() - min_val) * scale255 ).unaryExpr (Moments2func(stats));
}

/* flipV
   flip an image vertically
*/

void ImageMatrix::flipV() {
	bool old_has_median = has_median;
	Moments2 old_stats = stats;

	WriteablePixels() = ReadablePixels().rowwise().reverse();
	if (ColorMode != cmGRAY) {
		WriteableColors() = ReadableColors().rowwise().reverse();
	}
	// on its own, this operation doesn't affect the stats
	stats = old_stats;
	has_median = old_has_median;
}
/* flipH
   flip an image horizontally
*/
void ImageMatrix::flipH() {
	bool old_has_median = has_median;
	Moments2 old_stats = stats;

	WriteablePixels() = ReadablePixels().colwise().reverse();
	if (ColorMode != cmGRAY) {
		WriteableColors() = ReadableColors().colwise().reverse();
	}
	// on its own, this operation doesn't affect the stats
	stats = old_stats;
	has_median = old_has_median;
}

void ImageMatrix::invert() {
	double max_val = max(), min_val = min();
	WriteablePixels() = (max_val - ReadablePixels().array() + min_val).unaryExpr (Moments2func(stats));
}

/* Downsample
   down sample an image
   x_ratio, y_ratio -double- (0 to 1) the size of the new image comparing to the old one
   FIXME: Since this is done in-place, there is potential for aliasing (i.e. new pixel values interfering with old pixel values)
*/
void ImageMatrix::Downsample (const ImageMatrix &matrix_IN, double x_ratio, double y_ratio) {
	double x,y,dx,dy,frac;
	unsigned int new_x,new_y,a;
	HSVcolor hsv;


	if (x_ratio>1) x_ratio=1;
	if (y_ratio>1) y_ratio=1;
	dx=1/x_ratio;
	dy=1/y_ratio;

	if (dx == 1 && dy == 1) return;   /* nothing to scale */

	ImageMatrix copy_matrix;
	copy_matrix.copyFields (matrix_IN);
	copy_matrix.allocate (matrix_IN.width, matrix_IN.height);
	writeablePixels copy_pix_x = copy_matrix.WriteablePixels();
	writeableColors copy_clr_x = copy_matrix.WriteableColors();

	readOnlyPixels pix_plane_x = matrix_IN.ReadablePixels();
	readOnlyColors clr_plane_x = matrix_IN.ReadableColors();
 	unsigned int new_width = (unsigned int)(x_ratio*matrix_IN.width), new_height = (unsigned int)(y_ratio*matrix_IN.height),
 		old_width = matrix_IN.width, old_height = matrix_IN.height;

	// first downsample x
	for (new_y = 0; new_y < old_height; new_y++) {
		x = 0;
		new_x = 0;
		while (x < old_width) {
			double sum_i = 0;
			double sum_h = 0;
			double sum_s = 0;
			double sum_v = 0;

			/* the leftmost fraction of pixel */
			a = (unsigned int)(floor(x));
			frac = ceil(x)-x;
			if (frac > 0 && a < old_width) {
				sum_i += pix_plane_x(new_y,a) * frac;
				if (ColorMode != cmGRAY) {
					sum_h += clr_plane_x(new_y,a).h * frac;
					sum_s += clr_plane_x(new_y,a).s * frac;
					sum_v += clr_plane_x(new_y,a).v * frac;
				}
			} 

			/* the middle full pixels */
			for (a = (unsigned int)(ceil(x)); a < floor(x+dx); a++) {
				if (a < old_width) {
					sum_i += pix_plane_x(new_y,a);
					if (ColorMode != cmGRAY) {
						sum_h += clr_plane_x(new_y,a).h;
						sum_s += clr_plane_x(new_y,a).s;
						sum_v += clr_plane_x(new_y,a).v;
					}
				}
			}
			/* the right fraction of pixel */
			frac = x+dx - floor(x+dx);
			if (frac > 0 && a < old_width) {
				sum_i += pix_plane_x(new_y,a) * frac;
				if (ColorMode != cmGRAY) {
					sum_h += clr_plane_x(new_y,a).h * frac;
					sum_s += clr_plane_x(new_y,a).s * frac;
					sum_v += clr_plane_x(new_y,a).v * frac;
				}
			}

			copy_pix_x (new_y,new_x) = sum_i/(dx);
			if (ColorMode != cmGRAY) {
				hsv.h = (byte)(sum_h/(dx));
				hsv.s = (byte)(sum_s/(dx));
				hsv.v = (byte)(sum_v/(dx));
				copy_clr_x (new_y, new_x) = hsv;
			}

			x+=dx;
			new_x++;
		}
	}

	allocate (new_width, new_height);
	writeablePixels copy_pix_y = WriteablePixels();
	writeableColors copy_clr_y = WriteableColors();

	readOnlyPixels pix_plane_y = copy_matrix.ReadablePixels();
	readOnlyColors clr_plane_y = copy_matrix.ReadableColors();

	/* downsample y */
	for (new_x = 0; new_x < new_width; new_x++) {
		y = 0;
		new_y = 0;
		while (y < old_height) {
			double sum_i = 0;
			double sum_h = 0;
			double sum_s = 0;
			double sum_v = 0;

			a = (unsigned int)(floor(y));
			frac = ceil(y) - y;
			// take also the part of the leftmost pixel (if needed)
			if (frac > 0 && a < old_height) {
				sum_i += pix_plane_y(a,new_x) * frac;
				if (ColorMode != cmGRAY) {
					sum_h += clr_plane_y(a,new_x).h * frac;
					sum_s += clr_plane_y(a,new_x).s * frac;
					sum_v += clr_plane_y(a,new_x).v * frac;
				}
			}
			for (a = (unsigned int)(ceil(y)); a < floor(y+dy); a++) {
				if (a < old_height) {
					sum_i += pix_plane_y(a,new_x);
					if (ColorMode != cmGRAY) {
						sum_h += clr_plane_y(a,new_x).h;
						sum_s += clr_plane_y(a,new_x).s;
						sum_v += clr_plane_y(a,new_x).v;
					}
				}
			}
			frac=y+dy-floor(y+dy);
			if (frac > 0 && a < old_height) {
				sum_i += pix_plane_y(a,new_x) * frac;
				if (ColorMode != cmGRAY) {
					sum_h += clr_plane_y(a,new_x).h * frac;
					sum_s += clr_plane_y(a,new_x).s * frac;
					sum_v += clr_plane_y(a,new_x).v * frac;
				}
			}
			if (new_x < new_width && new_y < new_height) {
				copy_pix_y (new_y,new_x) = stats.add (sum_i/dy);
			}
			if (ColorMode != cmGRAY) {
				hsv.h = (byte)(sum_h/(dy));
				hsv.s = (byte)(sum_s/(dy));
				hsv.v = (byte)(sum_v/(dy));
				copy_clr_y (new_y, new_x) = hsv;
			}

			y+=dy;
			new_y++;
		}
	}
}


/* Rotate
   Rotate an image by 90, 120, or 270 degrees
   angle -double- (0 to 360) the degrees of rotation.  Only values of 90, 180, 270 are currently allowed
*/
void ImageMatrix::Rotate(const ImageMatrix &matrix_IN, double angle) {
	unsigned int new_width,new_height;
	Moments2 old_stats = stats;
	bool old_has_median = has_median;

	// Only deal with right angles
	if (! ( (angle == 90) || (angle == 180) || (angle == 270) ) ) return;

	// switch width/height if 90 or 270
	if ( (angle == 90) || (angle == 270) ) {
		new_width = matrix_IN.height;
		new_height = matrix_IN.width;
	} else {
		new_width = matrix_IN.width;
		new_height = matrix_IN.height;
	}

	// Copy fields from input
	copyFields (matrix_IN);

	allocate (new_width, new_height);

	// a 180 is simply a reverse of the matrix
	// a 90 is m.transpose().rowwise.reverse()
	// a 270 is m.transpose()
	switch ((int)angle) {
		case 90:
			WriteablePixels() = matrix_IN.ReadablePixels().transpose().rowwise().reverse();
		break;

		case 180:
			WriteablePixels() = matrix_IN.ReadablePixels().reverse();
		break;

		case 270:
			WriteablePixels() = matrix_IN.ReadablePixels().transpose();
		break;
	}
	
	// this operation does not change pixels statistics
	stats = old_stats;
	has_median = old_has_median;
}

// This pair of methods makes median-finding with and without caching for regular and const ImageMatrix objects
double ImageMatrix::update_median () {
	if (has_median) return _median;
	
	_median = get_median ();
	has_median = true;
	return _median;
}

double ImageMatrix::get_median () const {
	double median;
	size_t num = width * height;
	std::vector<double> v (num);

	readOnlyPixels pix_plane = ReadablePixels();
	for (size_t i = 0; i < num; i++) v[i] = pix_plane.array().coeff(i);
    size_t half = num / 2;
	if (num % 2 == 0) {
		nth_element(v.begin(), v.begin()+half, v.end());
		median = v[half];
		nth_element(v.begin(), v.begin()+half-1, v.end());
		median += (v[half-1]);
		median /= 2.0;
	} else {
		nth_element(v.begin(), v.begin()+half, v.end());
		median = v[half];
	}

	return median;
}

// This updates the sats and caches them in a non-const ImageMatrix.
void ImageMatrix::UpdateStats() {
	stats.reset();
	// note that sum() is there to force the evaluation of all of the coefficients.
	// hopefully the compiler won't optimize this away since we ignore the return.
	// Its not necessary to use sum() if using unaryExpr as part of an Eigen expression that gets assigned to something
	// e.g. invert and compute stats in one pass:
	// WriteablePixels() = (max_val - ReadablePixels().array() + min_val).unaryExpr (Moments2func(stats));
	ReadablePixels().unaryExpr (Moments2func(stats)).sum();
}

// This calculates sats and puts them in an externally supplied stats object, keeping the ImageMatrix const
// it uses the cached stats if they're valid.
void ImageMatrix::GetStats (Moments2 &moments2) const {
	if (stats.n() == width*height)
		moments2 = stats;
	else {
		// note that sum() is there to force the evaluation of all of the coefficients.
		// hopefully the compiler won't optimize this away since we ignore the return.
		// Its not necessary to use sum() if using unaryExpr as part of an Eigen expression that gets assigned to something
		// e.g. invert and compute stats in one pass:
		// WriteablePixels() = (max_val - ReadablePixels().array() + min_val).unaryExpr (Moments2func(stats));
		ReadablePixels().unaryExpr (Moments2func(moments2)).sum();
		// std::cout << "cache miss n=" << moments2.n() << std::endl;
	}
}

/* normalize the pixel values into a given range 
   min -double- the min pixel value (ignored if <0)
   max -double- the max pixel value (ignored if <0)
   range -long- nominal dynamic range (ignored if <0)
   n_mean -double- the mean of the normalized image (ignored if <0)
   n_std -double- the stddev of the normalized image (ignored if <0)
*/
void ImageMatrix::normalize(double n_min, double n_max, long n_range, double n_mean, double n_std) {
	unsigned int x,y;
	double val;
	writeablePixels pix_plane = WriteablePixels();

	/* normalized to n_min and n_max */
	if (n_min >= 0 && n_max > 0 && n_range > 0) {
		double norm_fact = n_range / (n_max-n_min);
		for (y = 0; y < height; y++) {
			for (x = 0; x < width; x++) {
				val = pix_plane (y,x);
				if (val < n_min) val = 0;
				else if (val > n_max) val = n_range;
				else val = norm_fact * (val - n_min);
				pix_plane (y,x) = stats.add(val);
			}
		}
	}

    /* normalize to n_mean and n_std */
	if (n_mean > 0) {
		double d_mean = mean() - n_mean, std_fact = (n_std > 0 ? n_std : 0)/std();
		double max_range = pow((double)2,bits)-1;
		for (y = 0; y < height; y++) {
			for (x = 0; x < width; x++) {
				val = pix_plane (y,x) - d_mean;
				if (n_std > 0)
					val = n_mean + (val-n_mean) * std_fact;
				if (val < 0) val = 0;
				else if (val > max_range) val = max_range;
				pix_plane (y,x) = stats.add (val);
			}
        }
	}	   
}

/* convolve
*/
void ImageMatrix::convolve(const pixDataMat &filter) {
	unsigned long x, y, xx, yy;
	long i, j;
	long height2=filter.rows()/2;
	long width2=filter.cols()/2;
	double tmp;

	std::string tmp_string;
	ImageMatrix temp;
	temp.copy (*this);
	temp.WriteablePixelsFinish();
	readOnlyPixels copy_pix_plane = temp.ReadablePixels();
	writeablePixels pix_plane = WriteablePixels();
	for (x = 0; x < width; ++x) {
		for (y = 0; y < height; ++y) {
			tmp=0.0;
			for (i = -width2; i <= width2; ++i) {
				xx=x+i;
				if (xx < width && xx >= 0) {
					for(j = -height2; j <= height2; ++j) {
						yy=y+j;
						if (yy >= 0 && yy < height) {
							tmp += filter (j+height2, i+width2) * copy_pix_plane(yy,xx);
						}
					}
				}
			}
			pix_plane (y,x) = stats.add(tmp);
		}
	}
}

/* find the basic color statistics
   hue_avg_p -double *- average hue
   hue_std_p -double *- standard deviation of the hue
   sat_avg_p -double *- average saturation
   sat_std_p -double *- standard deviation of the saturation
   val_avg_p -double *- average value
   val_std_p -double *- standard deviation of the value
   max_color_p -double *- the most popular color
   colors -double *- a histogram of colors
   if values are NULL - the value is not computed
*/

void ImageMatrix::GetColorStatistics(double *hue_avg_p, double *hue_std_p, double *sat_avg_p, double *sat_std_p, double *val_avg_p, double *val_std_p, double *max_color_p, double *colors) const {
	double hue_avg=0, hue_std=0, sat_avg=0, sat_std=0, val_avg=0, val_std=0;
	double delta, M2h=0, M2s=0, M2v=0;
	unsigned int a, x, y, pixel_index=0;
	unsigned long color_index=0;
	double max_val,pixel_num;
	double certainties[COLORS_NUM+1];
	byte h, s, v;
	readOnlyColors clr_plane = ReadableColors();

	pixel_num=height*width;

	/* calculate the average hue, saturation, value */
	if (hue_avg_p) *hue_avg_p=0;
	if (sat_avg_p) *sat_avg_p=0;
	if (val_avg_p) *val_avg_p=0;
	if (colors)
		for (a=0;a<=COLORS_NUM;a++)
			colors[a]=0;
	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			h = clr_plane(y,x).h;
			s = clr_plane(y,x).s;
			v = clr_plane(y,x).v;
			// This is Welford's cumulative mean+variance algorithm as reported by Knuth
			pixel_index++;
			// h
			delta = h - hue_avg;
			hue_avg += delta/pixel_index;
			M2h += delta * (h - hue_avg);
			// s
			delta = s - sat_avg;
			sat_avg += delta/pixel_index;
			M2s += delta * (s - sat_avg);
			// v
			delta = v - val_avg;
			val_avg += delta/pixel_index;
			M2v += delta * (v - val_avg);

			color_index=FindColor(h,s,v,certainties);
			colors[color_index]+=1;
		}
	}
	hue_std = sqrt ( M2h / (pixel_index - 1) );
	sat_std = sqrt ( M2s / (pixel_index - 1) );
	val_std = sqrt ( M2v / (pixel_index - 1) );

	if (hue_avg_p) *hue_avg_p = hue_avg;
	if (sat_avg_p) *sat_avg_p = sat_avg;
	if (val_avg_p) *val_avg_p = val_avg;
	if (hue_std_p) *hue_std_p = hue_std;
	if (sat_std_p) *sat_std_p = sat_std;
	if (val_std_p) *val_std_p = val_std;

	/* max color (the most common color in the image) */
	if (max_color_p) {
		*max_color_p=0;
		max_val=0.0;
		for (a = 0; a <= COLORS_NUM; a++) {
			if (colors[a] > max_val) {
				max_val=colors[a];
				*max_color_p=a;
			}
		}
	}
	/* colors */
	if (colors)
		for (a = 0; a <= COLORS_NUM; a++)
			colors[a]=colors[a]/pixel_num;
}

/* ColorTransform
   Transform a color image to a greyscale image such that each
   color_hist -double *- a histogram (of COLORS_NUM + 1 bins) of the colors. This parameter is ignored if NULL
   use_hue -int- 0 if classifying colors, 1 if using the hue component of the HSV vector
   grey level represents a different color
*/
void ImageMatrix::ColorTransform (const ImageMatrix &matrix_IN) {  
	unsigned int x,y; //,base_color;
	double cb_intensity;
	double max_range = pow((double)2,8)-1;
	HSVcolor hsv_pixel;
	unsigned long color_index=0;   
	double certainties[COLORS_NUM+1];

	copyFields (matrix_IN);
	// The result is an intensity image, so eliminate the color plane
	ColorMode = cmGRAY;
	allocate (matrix_IN.width, matrix_IN.height);
	writeablePixels out_plane = WriteablePixels();
	readOnlyColors clr_plane = matrix_IN.ReadableColors();

	// find the colors
	for( y = 0; y < height; y++ ) {
		for( x = 0; x < width; x++ ) { 
			hsv_pixel = clr_plane (y, x);
			color_index = FindColor( hsv_pixel.h,  hsv_pixel.s, hsv_pixel.v, certainties );
			// convert the color index to a greyscale value
			cb_intensity = int( ( max_range * color_index ) / COLORS_NUM );
			out_plane (y, x) = stats.add (cb_intensity);
		}
	}
}

void ImageMatrix::HueTransform (const ImageMatrix &matrix_IN) {  
	unsigned int x,y; //,base_color;

	copyFields (matrix_IN);
	// The result is an intensity image, so eliminate the color plane
	ColorMode = cmGRAY;
	allocate (matrix_IN.width, matrix_IN.height);
	writeablePixels out_plane = WriteablePixels();
	readOnlyColors clr_plane = matrix_IN.ReadableColors();

	// find the colors
	for( y = 0; y < height; y++ ) {
		for( x = 0; x < width; x++ ) { 
			out_plane (y, x) = stats.add (clr_plane (y, x).h);
		}
	}
}

/* get image histogram */
void ImageMatrix::histogram(double *bins,unsigned short nbins, bool imhist, const Moments2 &in_stats) const {
	unsigned long a, bin, num = width*height;
	double val,  h_min = INF, h_max = -INF, h_scale;
	readOnlyPixels pix_plane = ReadablePixels();

	/* find the minimum and maximum */
	if (imhist) {    /* similar to the Matlab imhist */
		h_min = 0;
		h_max = pow((double)2,bits)-1;
	} else if (in_stats.n() > 0) {
		h_min = in_stats.min();
		h_max = in_stats.max();
	} else {
		// to keep this const method from modifying the object, we use GetStats on a local Moments2 object
		Moments2 local_stats;
		GetStats (local_stats);
		h_min = local_stats.min();
		h_max = local_stats.max();
	}
	if (h_max-h_min > 0) h_scale = (double)nbins / double(h_max-h_min);
	else h_scale = 0;

	// initialize the bins
	memset(bins, 0, nbins * sizeof (double));

   // build the histogram
	for (a = 0; a < num; a++) {
		val = pix_plane.array().coeff(a);
		bin = (unsigned long)(( (val - h_min)*h_scale));
		if (bin >= nbins) bin = nbins-1;
		bins[bin] += 1.0;
	}

	return;
}

/* fft 2 dimensional transform */
// http://www.fftw.org/doc/
double ImageMatrix::fft2 (const ImageMatrix &matrix_IN) {
	fftw_plan p;
	unsigned int half_height = matrix_IN.height/2+1;

	copyFields (matrix_IN);
	allocate (matrix_IN.width, matrix_IN.height);
	writeablePixels out_plane = WriteablePixels();
	readOnlyPixels in_plane = matrix_IN.ReadablePixels();

	double *in = (double*) fftw_malloc(sizeof(double) * width*height);
 	fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * width*height);
	p = fftw_plan_dft_r2c_2d(width,height,in,out, FFTW_MEASURE); // FFTW_ESTIMATE: deterministic
	unsigned int x,y;
 	for (x=0;x<width;x++)
 		for (y=0;y<height;y++)
 			in[height*x+y]=in_plane.coeff(y,x);
 
 	fftw_execute(p);

	// The resultant image uses the modulus (sqrt(nrm)) of the complex numbers for pixel values
	unsigned long idx;
 	for (x=0;x<width;x++) {
 		for (y=0;y<half_height;y++) {
 			idx = half_height*x+y;
 			out_plane (y,x) = stats.add (sqrt( pow( out[idx][0],2)+pow(out[idx][1],2)));    // sqrt(real(X).^2 + imag(X).^2)
 		}
 	}

	// complete the first column
 	for (y=half_height;y<height;y++)
 		out_plane (y,0) = stats.add (out_plane (height - y, 0));

	// complete the rest of the columns
 	for (y=half_height;y<height;y++)
 		for (x=1;x<width;x++)   // 1 because the first column is already completed
 			out_plane (y,x) = stats.add (out_plane (height - y, width - x));

	// clean up
	fftw_destroy_plan(p);
	fftw_free(in);
	fftw_free(out);

// 
// 	// Doing this using the Eigen library
// 	// FIXME: This doesn't quite work - causes indirect segfault when reading the output into pix_plane
// 	using namespace Eigen; // just for this function
// 	double *in = (double*) fftw_malloc(sizeof(double) * width*height);
// 	fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * width*height);
// 	p = fftw_plan_dft_r2c_2d(width,height,in,out, FFTW_MEASURE); // FFTW_ESTIMATE: deterministic
// 	// Eigen takes care of converting the input from row major (pix_plane) to column major (Eigen and FFTW's default)
// 	Map<Array< double, Dynamic, Dynamic, ColMajor >, Aligned > in_m (&in[0],height,width);
// 	in_m = pix_plane;
// 
// 	fftw_execute(p);
// 
// 	// The resultant image uses the modulus (sqrt(norm)) of the complex numbers for pixel values
// 	// Map an Eigen matrix onto what was returned by FFTW.
// 	// Tricky bits specific to FFTW output:
// 	//  We are making two matrixes out the single output - one for the real component and one for complex.
// 	//  We are using strides to skip every other element of the output, starting at 0 for real, and 1 for complex.
// 	//  We are not using the bottom half, making it symmetrical around the horizontal midline instead, but we still stride over twice the height b/c its complex
// 	Map<const Array< double, Dynamic, Dynamic, ColMajor >, Aligned, Stride<Dynamic,2> > out_r (&out[0][0],half_height,width, Stride<Dynamic,2>(height*sizeof(fftw_complex),2));
// 	Map<const Array< double, Dynamic, Dynamic, ColMajor >, Aligned, Stride<Dynamic,2> > out_i (&out[0][1],half_height,width, Stride<Dynamic,2>(height*sizeof(fftw_complex),2));
// 	pix_plane.topRows(0,0,half_height,width) = (out_r.pow(2)+out_i.pow(2)).sqrt();
// 
// 	// complete the first column
// 	unsigned int  odd_pad = (height % 2 == 0 ? 0 : 1);
// 	pix_plane.block (half_height,0,half_height-odd_pad,1) = pix_plane.block (odd_pad,0,half_height-odd_pad,1).reverse();
// 
// 	// complete the rest of the columns
// 	pix_plane.block (half_height,1,half_height-odd_pad,width-1) = pix_plane.block (odd_pad,1,half_height-odd_pad,width-1).reverse();
// 	has_median = false;
// 	stats.reset();
// 
// 	// clean up
// 	fftw_destroy_plan(p);
// 	fftw_free(in);
// 	fftw_free(out);

	WriteablePixelsFinish();
	return(0);
}

/* chebyshev transform */
void ImageMatrix::ChebyshevTransform(const ImageMatrix &matrix_IN, unsigned int N) {
	double *out;
	unsigned int x,y;

	copyFields (matrix_IN);

	if (N<2)
		N = MIN( width, height );
	out=new double[height*N];
	Chebyshev2D(matrix_IN, out,N);
	width=N;
	height = MIN( height, N );   /* prevent error */

	allocate (width,height);
	writeablePixels pix_plane = WriteablePixels();
	for(y=0;y<height;y++)
		for(x=0;x<width;x++)
			pix_plane (y,x) = stats.add (out[y * width + x]);
	delete [] out;
	WriteablePixelsFinish();
}

/* chebyshev transform
   coeff -array of double- a pre-allocated array of 32 doubles
*/
void ImageMatrix::ChebyshevFourierTransform2D(double *coeff) const {
	ImageMatrix *matrix;
	const ImageMatrix *const_matrix;
	if( (width * height) > (300 * 300) ) {
		matrix = new ImageMatrix;
		matrix->copy (*this);
		matrix->Downsample(*this, MIN( 300.0/(double)width, 300.0/(double)height ), MIN( 300.0/(double)width, 300.0/(double)height ) );  /* downsample for avoiding memory problems */
		const_matrix = matrix;
	} else {
		const_matrix = this;
	}
	ChebyshevFourier2D(*const_matrix, 0, coeff,32);
	if (const_matrix != this) delete const_matrix;
}


/* Symlet5 transform */
void ImageMatrix::Symlet5Transform(const ImageMatrix &matrix_IN) {
	unsigned int x,y;
	DataGrid *grid=NULL;
	Symlet5 *Sym5;

	copyFields (matrix_IN);
	readOnlyPixels in_plane = matrix_IN.ReadablePixels();

	grid = new DataGrid2D(width,height,-1);
  
	for (y=0;y<height;y++)
		for(x=0;x<width;x++)
			grid->setData(x,y,-1,in_plane(y,x));
	Sym5=new Symlet5(0,1);
	Sym5->transform2D(grid);
	
	allocate (grid->getX(), grid->getY());
	writeablePixels out_plane = WriteablePixels();
	for (y=0;y<height;y++)
		for(x=0;x<width;x++)
			out_plane (y,x) = stats.add (grid->getData(x,y,-1));
		 
	delete Sym5;
	delete grid;
	WriteablePixelsFinish();
}

/* chebyshev statistics
   coeff -array of double- pre-allocated memory of 20 doubles
   nibs_num - (32 is normal)
*/
void ImageMatrix::ChebyshevStatistics2D(double *coeff, unsigned int N, unsigned int nbins) const {
	if (N<2) N=20;
	if (N>MIN(width,height)) N=MIN(width,height);   

	ImageMatrix tmp_transform;
	tmp_transform.ChebyshevTransform(*this, N);
	tmp_transform.histogram(coeff,nbins,0);
}

/* CombFirstFourMoments
   vec should be pre-alocated array of 48 doubles
*/
int ImageMatrix::CombFirstFourMoments2D(double *vec) const {
	int count;
	count = CombFirst4Moments2D (*this, vec);   
	vd_Comb4Moments (vec);   
	return (count);
}

/* Edge Transform */
void ImageMatrix::EdgeTransform (const ImageMatrix &matrix_IN) {
	unsigned int x,y;
	double max_x=0,max_y=0;

	copyFields (matrix_IN);
	allocate (matrix_IN.width, matrix_IN.height);

	readOnlyPixels in_plane = matrix_IN.ReadablePixels();
	writeablePixels out_plane = WriteablePixels();

	for (y = 0; y < height; y++)
		for (x = 0; x < width; x++) {
			max_x = max_y = 0;
			if (y > 0 && y < height-1) max_y=MAX(fabs(in_plane(y,x) - in_plane(y-1,x)), fabs(in_plane(y,x) - in_plane(y+1,x)));
			if (x > 0 && x < width-1)  max_x=MAX(fabs(in_plane(y,x) - in_plane(y,x-1)), fabs(in_plane(y,x) - in_plane(y,x+1)));
			out_plane(y,x) = stats.add (MAX(max_x,max_y));
		}
	WriteablePixelsFinish();
}

/* Prewitt gradient magnitude
   output - a pre-allocated matrix that will hold the output (the input matrix is not changed)
            output should be of the same size as the input matrix
*/
void ImageMatrix::PrewittMagnitude2D (const ImageMatrix &matrix_IN) {
	long x,y,i,j,w,h;
	double sumx,sumy;

	copyFields (matrix_IN);
	allocate (matrix_IN.width, matrix_IN.height);
	w = width;
	h = height;
	writeablePixels out_pix_plane = WriteablePixels();
	readOnlyPixels pix_plane = matrix_IN.ReadablePixels();

	for (x = 0; x < w; x++) {
		for (y = 0; y < h; y++) {
			sumx=0;
			sumy=0;		  
			for (j = y-1; j <= y+1; j++)
				if (j >= 0 && j < h && x-1 >= 0)
					sumx += pix_plane(j,x-1)*1;//0.3333;
			for (j = y-1; j <= y+1; j++)
				if (j >= 0 && j < h && x+1 < w)
					sumx += pix_plane(j,x+1)*-1;//-0.3333;
			for (i = x-1; i <= x+1; i++)
				if (i >= 0 && i < w && y-1 >= 0)
					sumy += pix_plane(y-1,i)*1;//-0.3333;
			for (i = x-1; i <= x+1; i++)
				if (i >= 0 && i < w && y+1 < h)
					sumy += pix_plane(y+1,i)*-1;//0.3333;
			out_pix_plane(y,x) = stats.add (sqrt(sumx*sumx+sumy*sumy));
		}
	}
	WriteablePixelsFinish();
}

/* Prewitt gradient direction
   output - a pre-allocated matrix that will hold the output (the input matrix is not changed)
            output should be of the same size as the input matrix
*/
void ImageMatrix::PrewittDirection2D(const ImageMatrix &matrix_IN) {
	long x,y,i,j,w,h;
	double sumx,sumy;

	copyFields (matrix_IN);
	allocate (matrix_IN.width, matrix_IN.height);
	w = width;
	h = height;

	writeablePixels out_pix_plane = WriteablePixels();
	readOnlyPixels pix_plane = matrix_IN.ReadablePixels();

	for (x = 0; x < w; x++)
		for (y = 0;y < h; y++) {
			sumx=0;
			sumy=0;
			for (j = y-1; j <= y+1; j++)
				if (j >= 0 && j < h && x-1 >= 0)
					sumx += pix_plane(j,x-1)*1;//0.3333;
			for (j = y-1; j <= y+1; j++)
				if (j >= 0 && j < h && x+1 < w)
					sumx += pix_plane(j,x+1)*-1;//-0.3333;
			for (i = x-1; i <= x+1; i++)
				if (i >= 0 && i < w && y-1 >= 0)
					sumy += pix_plane(y-1,i)*1;//-0.3333;
			for (i = x-1; i <= x+1; i++)
				if (i >= 0 && i < w && y+1 < h)
					sumy += pix_plane(y+1,i)*-1;//0.3333;
			if (sumy == 0 || fabs(sumy)<1/INF) out_pix_plane(y,x) = stats.add (3.1415926 * (sumx < 0 ? 1 : 0));
			else out_pix_plane(y,x) = stats.add (atan2(sumy,sumx));
		}
	WriteablePixelsFinish();
}

/* edge statistics */
//#define NUM_BINS 8
//#define NUM_BINS_HALF 4
/* EdgeArea -long- number of edge pixels
   MagMean -double- mean of the gradient magnitude
   MagMedian -double- median of the gradient magnitude
   MagVar -double- variance of the gradient magnitude
   MagHist -array of double- histogram of the gradient magnitude. array of size "nbins" should be allocated before calling the function
   DirecMean -double- mean of the gradient direction
   DirecMedian -double- median of the gradient direction
   DirecVar -double- variance of the gradient direction
   DirecHist -array of double- histogram of the gradient direction. array of size "nbins" should be allocated before calling the function
   DirecHomogeneity -double-
   DiffDirecHist -array of double- array of size nbins/2 should be allocated
*/

void ImageMatrix::EdgeStatistics(unsigned long *EdgeArea, double *MagMean, double *MagMedian, double *MagVar,
	double *MagHist, double *DirecMean, double *DirecMedian, double *DirecVar, double *DirecHist,
	double *DirecHomogeneity, double *DiffDirecHist, unsigned int nbins) const {

	unsigned int a,bin_index;
	double min_val,max_val,sum, level;

	ImageMatrix GradientMagnitude;
	GradientMagnitude.PrewittMagnitude2D (*this);

	readOnlyPixels GM_pix_plane = GradientMagnitude.ReadablePixels();
	ImageMatrix GradientDirection;
	GradientDirection.PrewittDirection2D (*this);

	/* find gradient statistics */
	*MagMedian = GradientMagnitude.update_median();
	*MagMean   = GradientMagnitude.mean();
	*MagVar    = GradientMagnitude.var();
	min_val   = GradientMagnitude.min();
	max_val   = GradientMagnitude.max();
	GradientMagnitude.histogram (MagHist, nbins, false);

	/* find the edge area (number of edge pixels) */
	*EdgeArea = 0;
	// level = min_val + ((max_val-min_val)/2.0);
	// level = Otsu();
	level = *MagMean;
	// level = min_val + ((max_val-min_val)/2.0);   // level=duplicate->OtsuBinaryMaskTransform()   // level=MagMean

	for (a = 0; a < GradientMagnitude.height*GradientMagnitude.width; a++)
		if (GM_pix_plane.array().coeff(a) > level) (*EdgeArea)+=1; /* find the edge area */
//   GradientMagnitude->OtsuBinaryMaskTransform();

	/* find direction statistics */
	*DirecMedian = GradientDirection.update_median();
	*DirecMean   = GradientDirection.mean();
	*DirecVar    = GradientDirection.var();
	min_val      = GradientDirection.min();
	max_val      = GradientDirection.max();
	GradientDirection.histogram (DirecHist, nbins, false);

	/* Calculate statistics about edge difference direction
	   Histogram created by computing differences amongst histogram bins at angle and angle+pi
	*/
	for (bin_index = 0; bin_index < (nbins/2); bin_index++)
		DiffDirecHist[bin_index] = fabs(DirecHist[bin_index]-DirecHist[bin_index+(int)(nbins/2)]);
	sum=0;
	for (bin_index = 0; bin_index < (nbins/2); bin_index++) {
		if (DirecHist[bin_index] + DirecHist[bin_index+(int)(nbins/2)] != 0)  /* protect from a numeric flaw */
			DiffDirecHist[bin_index] = DiffDirecHist[bin_index]/(DirecHist[bin_index]+DirecHist[bin_index+(int)(nbins/2)]);
		sum += (DirecHist[bin_index]+DirecHist[bin_index+(int)(nbins/2)]);
	}

	/* The fraction of edge pixels that are in the first two bins of the histogram measure edge homogeneity */
	if (sum > 0) *DirecHomogeneity = (DirecHist[0]+DirecHist[1])/sum;

}

/* radon transform
   vec -array of double- output column. a pre-allocated vector of the size 3*4=12
*/
void ImageMatrix::RadonTransform2D(double *vec) const {
	unsigned int x,y,val_index,output_size,vec_index;
	double *pixels,*ptr,bins[3];
	int angle,num_angles=4;
	double theta[4]={0,45,90,135};
	int rLast,rFirst;
	rLast = (int) ceil(sqrt(pow( (double)(width-1-(width-1)/2),2)+pow( (double)(height-1-(height-1)/2),2))) + 1;
	rFirst = -rLast;
	output_size=rLast-rFirst+1;
	readOnlyPixels pix_plane = ReadablePixels();

	ptr = new double[output_size*num_angles];
	for (val_index = 0; val_index < output_size*num_angles; val_index++)
		ptr[val_index] = 0;  /* initialize the output vector */

    pixels=new double[width*height];
    vec_index = 0;

	for (x = 0; x < width; x++)
		for (y = 0; y < height; y++)
			pixels[y+height*x] = pix_plane(y,x);

    radon(ptr,pixels, theta, height, width, (width-1)/2, (height-1)/2, num_angles, rFirst, output_size);

	for (angle = 0; angle < num_angles; angle++) {
		//radon(ptr,pixels, &theta, height, width, (width-1)/2, (height-1)/2, 1, rFirst, output_size);
		/* create histogram */
		double h_min = INF, h_max = -INF, h_scale;
		unsigned long bin, nbins = 3;
		/* find the minimum and maximum values */
		for (val_index = angle*output_size; val_index < (angle+1)*output_size; val_index++) {
			if (ptr[val_index] > h_max) h_max = ptr[val_index];
			if (ptr[val_index] < h_min) h_min = ptr[val_index];
		}

		if (h_max-h_min > 0) h_scale = (double)nbins / double(h_max - h_min);
		else h_scale = 0;

		memset(bins, 0, nbins * sizeof (double));
		for (val_index=angle*output_size;val_index<(angle+1)*output_size;val_index++) {
			bin = (unsigned long)(( (ptr[val_index] - h_min)*h_scale));
			if (bin >= nbins) bin = nbins-1;
			bins[bin] += 1.0;
		}

		for (bin = 0; bin < 3; bin++)
			vec[vec_index++] = bins[bin];
	}
	vd_RadonTextures(vec);
	delete [] pixels;
	delete [] ptr;
}

//-----------------------------------------------------------------------------------
/* Otsu
   Find otsu threshold
*/
double ImageMatrix::Otsu(bool dynamic_range) const {
     /* binarization by Otsu's method 
	based on maximization of inter-class variance */
#define OTSU_LEVELS 1024
	double hist[OTSU_LEVELS];
	double omega[OTSU_LEVELS];
	double myu[OTSU_LEVELS];
	double max_sigma, sigma[OTSU_LEVELS]; // inter-class variance
	int i;
	int threshold;
    double min_val,max_val; // pixel range

	if (!dynamic_range) {
		histogram(hist,OTSU_LEVELS,true);
		min_val = 0.0;
		max_val = pow(2.0,bits)-1;
	} else {
		// to keep this const method from modifying the object, we use GetStats on a local Moments2 object
		Moments2 local_stats;
		GetStats (local_stats);
		min_val = local_stats.min();
		max_val = local_stats.max();
		histogram(hist,OTSU_LEVELS,false);
	}
  
	// omega & myu generation
	omega[0] = hist[0] / (width * height);
	myu[0] = 0.0;
	for (i = 1; i < OTSU_LEVELS; i++) {
		omega[i] = omega[i-1] + (hist[i] / (width * height));
		myu[i] = myu[i-1] + i*(hist[i] / (width * height));
	}
  
	// maximization of inter-class variance
	threshold = 0;
	max_sigma = 0.0;
	for (i = 0; i < OTSU_LEVELS-1; i++) {
		if (omega[i] != 0.0 && omega[i] != 1.0)
			sigma[i] = pow(myu[OTSU_LEVELS-1]*omega[i] - myu[i], 2) / 
				(omega[i]*(1.0 - omega[i]));
		else
			sigma[i] = 0.0;
		if (sigma[i] > max_sigma) {
			max_sigma = sigma[i];
			threshold = i;
		}
	}

	// threshold is a histogram index - needs to be scaled to a pixel value.
	return ( (((double)threshold / (double)(OTSU_LEVELS-1)) * (max_val - min_val)) + min_val );
}

//-----------------------------------------------------------------------------------
/*
  OtsuBinaryMaskTransform
  Transforms an image to a binary image such that the threshold is otsu global threshold
*/
double ImageMatrix::OtsuBinaryMaskTransform (const ImageMatrix &matrix_IN) {
	double OtsuGlobalThreshold;

	copyFields (matrix_IN);
	allocate (matrix_IN.width, matrix_IN.height);

	writeablePixels out_plane = WriteablePixels();
	readOnlyPixels in_plane = matrix_IN.ReadablePixels();

	OtsuGlobalThreshold = matrix_IN.Otsu();

	/* classify the pixels by the threshold */
	for (unsigned int a = 0; a < width*height; a++)
		if (in_plane.array().coeff(a) > OtsuGlobalThreshold) (out_plane.array())(a) = stats.add (1);
		else (out_plane.array())(a) = stats.add (0);

	return(OtsuGlobalThreshold);
}

/*  BWlabel
    label groups of connected pixel (4 or 8 connected dependes on the value of the parameter "level").
    This is an implementation of the Matlab function bwlabel
    returned value -int- the number of objects found
*/
//--------------------------------------------------------
unsigned long ImageMatrix::BWlabel(int level) {
	return(bwlabel(*this,level));
}

//--------------------------------------------------------

void ImageMatrix::centroid(double *x_centroid, double *y_centroid) const {
	GlobalCentroid(*this,x_centroid,y_centroid);
}

//--------------------------------------------------------

/*
  FeatureStatistics
  Find feature statistics. Before calling this function the image should be transformed into a binary
  image using "OtsuBinaryMaskTransform".

  count -int *- the number of objects detected in the binary image
  Euler -int *- the euler number (number of objects - number of holes
  centroid_x -int *- the x coordinate of the centroid of the binary image
  centroid_y -int *- the y coordinate of the centroid of the binary image
  AreaMin -int *- the smallest area
  AreaMax -int *- the largest area
  AreaMean -int *- the mean of the areas
  AreaMedian -int *- the median of the areas
  AreaVar -int *- the variance of the areas
  DistMin -int *- the smallest distance
  DistMax -int *- the largest distance
  DistMean -int *- the mean of the distance
  DistMedian -int *- the median of the distances
  DistVar -int *- the variance of the distances

*/

int compare_ulongs (const void *a, const void *b) {
	if (*((unsigned long *)a) > *((unsigned long *)b)) return(1);
	if (*((unsigned long *)a) == *((unsigned long *)b)) return(0);
	return(-1);
}

void ImageMatrix::FeatureStatistics(unsigned long *count, long *Euler, double *centroid_x, double *centroid_y, unsigned long *AreaMin, unsigned long *AreaMax,
	double *AreaMean, unsigned int *AreaMedian, double *AreaVar, unsigned int *area_histogram,double *DistMin, double *DistMax,
	double *DistMean, double *DistMedian, double *DistVar, unsigned int *dist_histogram, unsigned int nbins
) const {
	unsigned long object_index, bin;
	double sum_areas,sum_dists;
	ImageMatrix BWImage;
	unsigned long *object_areas;
	double *centroid_dists, sum_dist, hist_scale;

	BWImage.OtsuBinaryMaskTransform(*this);
	BWImage.centroid(centroid_x,centroid_y);
	*count = BWImage.BWlabel(8);
	*Euler=EulerNumber(BWImage,8);

	// calculate the areas 
	sum_areas = 0;
	sum_dists = 0;
	object_areas = new unsigned long[*count];
	centroid_dists = new double[*count];
	for (object_index = 0; object_index < *count; object_index++) {
		double x_centroid,y_centroid;
		object_areas[object_index] = FeatureCentroid(BWImage, object_index+1, &x_centroid, &y_centroid);
		centroid_dists[object_index] = sqrt(pow(x_centroid-(*centroid_x),2) + pow(y_centroid - (*centroid_y),2));
		sum_areas += object_areas[object_index];
		sum_dists += centroid_dists[object_index];
	}

	// compute area statistics
	qsort(object_areas,*count,sizeof(unsigned long),compare_ulongs);
	*AreaMin=object_areas[0];
	*AreaMax=object_areas[*count-1];
	if (*count>0) *AreaMean=sum_areas/(*count);
	else *AreaMean=0;
	*AreaMedian=(unsigned int)(object_areas[(*count)/2]);
	memset (area_histogram, 0, nbins * sizeof (int));
	if ((*AreaMax - *AreaMin) > 0) hist_scale = (double)nbins / double(*AreaMax - *AreaMin);
	else hist_scale = 0;
	/* compute the variance and the histogram */
	sum_areas = 0;
	if (*AreaMax - *AreaMin > 0) {
		for (object_index = 0; object_index < *count; object_index++) {
			sum_areas += pow(object_areas[object_index] - *AreaMean,2);
			bin = (unsigned long)( ((double)object_areas[object_index] - *AreaMin) * hist_scale);
			if (bin >= nbins) bin = nbins - 1;
			area_histogram[bin] += 1;
		}
	}
	if (*count > 1) *AreaVar=sum_areas / ((*count)-1);
	else *AreaVar=sum_areas;

   /* compute distance statistics */
	qsort(centroid_dists,*count,sizeof(double),compare_doubles);
	*DistMin=centroid_dists[0];
	*DistMax=centroid_dists[*count-1];
	if (*count>0) *DistMean=sum_dists/(*count);
	else *DistMean=0;
	*DistMedian=centroid_dists[(*count)/2];
	memset (dist_histogram, 0, nbins * sizeof (int));
	if ((*DistMax - *DistMin) > 0) hist_scale = (double)nbins / double(*DistMax - *DistMin);
	else hist_scale = 0;

   /* compute the variance and the histogram */
	sum_dist=0;
	for (object_index = 0; object_index < *count; object_index++) {
		sum_dist += pow(centroid_dists[object_index] - *DistMean,2);
		bin = (unsigned long)((centroid_dists[object_index] - *DistMin) * hist_scale);
		if (bin >= nbins) bin = nbins - 1;
		dist_histogram[bin] += 1;
	}
	if (*count>1) *DistVar=sum_dist/((*count)-1);
	else *DistVar=sum_dist;

	delete [] object_areas;
	delete [] centroid_dists;
}

/* GaborFilters */
/* ratios -array of double- a pre-allocated array of double[7]
*/
void ImageMatrix::GaborFilters2D(double *ratios) const {
	GaborTextureFilters2D(*this, ratios);
}


/* haralick
   output -array of double- a pre-allocated array of 28 doubles
*/
void ImageMatrix::HaralickTexture2D(double distance, double *out) const {
	if (distance<=0) distance=1;
	haralick2D(*this,distance,out);
}

/* MultiScaleHistogram
   histograms into 3,5,7,9 bins
   Function computes signatures based on "multiscale histograms" idea.
   Idea of multiscale histogram came from the belief of a unique representativity of an
   image through infinite series of histograms with sequentially increasing number of bins.
   Here we used 4 histograms with number of bins being 3,5,7,9.
   out -array of double- a pre-allocated array of 24 bins
*/
void ImageMatrix::MultiScaleHistogram(double *out) const {
	int a;
	double max_val=0;
	Moments2 local_moments;
	GetStats (local_moments);
	histogram (out,        3, false, local_moments);
	histogram (&(out[3]),  5, false, local_moments);
	histogram (&(out[8]),  7, false, local_moments);
	histogram (&(out[15]), 9, false, local_moments);
	for (a = 0; a < 24; a++)
		if (out[a] > max_val) max_val = out[a];
	for (a = 0; a < 24; a++)
		out[a] /= max_val;
}

/* TamuraTexture
   Tamura texture signatures: coarseness, directionality, contrast
   vec -array of double- a pre-allocated array of 6 doubles
*/
void ImageMatrix::TamuraTexture2D(double *vec) const {
	Tamura3Sigs2D(*this,vec);
}

/* zernike
   zvalue -array of double- a pre-allocated array of double of a suficient size
                            (the actual size is returned by "output_size))
   output_size -* long- the number of enteries in the array "zvalues" (normally 72)
*/
void ImageMatrix::zernike2D(double *zvalues, long *output_size) const {
	mb_zernike2D(*this, 0, 0, zvalues, output_size);
}


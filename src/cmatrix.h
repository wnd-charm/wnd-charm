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

#include <vector>
#include <string> // for what_am_i definition
#ifdef WIN32
  #include <vcl.h>
#else  
  #include "colors/FuzzyCalc.h"
  //#define min(a,b) (((a) < (b)) ? (a) : (b))
  //#define max(a,b) (((a) < (b)) ? (b) : (a))
#endif

#define cmRGB 1
#define cmHSV 2

#define INF 10E200

using namespace std;

class FeatureGroup;

typedef unsigned char byte;

typedef struct RGBCOLOR
{  byte red,green,blue;
}RGBcolor;

typedef struct HSVCOLOR
{  byte hue,saturation,value;
}HSVcolor;

typedef union
{  RGBcolor RGB;
   HSVcolor HSV;
}color;


typedef struct PIX_DATA
{  color clr;
   double intensity;  /* normailized to (0,255) interval */
} pix_data;

typedef struct
{  int x,y,w,h;
}
rect;

int compare_doubles (const void *a, const void *b);

class ImageMatrix
{
  public:
   pix_data *data;                                 /* data of the colors                   */  
	 //std::string what_am_i;                          // informative label
	 int ColorMode;                                  /* can be cmRGB or cmHSV                */
   unsigned short bits;                            /* the number of intensity bits (8,16, etc) */
   int width,height,depth;                         /* width and height of the picture      */
#ifdef WIN32   
   int LoadImage(TPicture *picture,int ColorMode);
   int LoadBMP(char *filename,int ColorMode);      /* load from a bitmap file              */
#endif
   int LoadTIFF(char *filename);                   /* load from TIFF file                  */
   int SaveTiff(char *filename);                   /* save a matrix in TIF format          */
   int LoadPPM(char *filename, int ColorMode);     /* load from a PPM file                 */
   int OpenImage(char *image_file_name, int downsample, rect *bounding_rect, double mean, double stddev); /* load an image of any supported format */
    //void CmatrixMessage();
   ImageMatrix();                                  /* basic constructor                    */
   ImageMatrix(int width,int height,int depth);    /* construct a new empty matrix         */
   ImageMatrix(ImageMatrix *matrix,int x1, int y1, int x2, int y2, int z1, int z2);  /* create a new matrix which is part of the original one */
   ~ImageMatrix();                                 /* destructor */
   ImageMatrix *duplicate();                       /* create a new identical matrix        */
	 //void dump();                                    // dump the pixel intensities to a file for inspection
   pix_data pixel(int x,int y,int z);              /* get a pixel value                    */
   void set(int x,int y,int z, pix_data val);      /* assign a pixel value                 */       
   void SetInt(int x,int y,int z, double val);     /* set only the intensity of the pixel  */   
   void diff(ImageMatrix *matrix);                 /* compute the difference from another image */
   void normalize(double min, double max, long range, double mean, double stddev); /* normalized an image to either min/max or mean/stddev */
   void to8bits();
   void flip();                                    /* flip an image horizonatally          */
   void invert();                                  /* invert the intensity of an image     */
   void Downsample(double x_ratio, double y_ratio);/* down sample an image                 */
   ImageMatrix *Rotate(double angle);              /* rotate an image by 90,180,270 degrees*/
   void convolve(ImageMatrix *filter);
   void BasicStatistics(double *mean, double *median, double *std, double *min, double *max, double *histogram, int bins);
   void GetColorStatistics(double *hue_avg, double *hue_std, double *sat_avg, double *sat_std, double *val_avg, double *val_std, double *max_color, double *colors);
   void ColorTransform(double *color_hist, int use_hue);
   void histogram(double *bins,unsigned short bins_num, int imhist);
   double Otsu();                                  /* Otsu grey threshold                  */
   void MultiScaleHistogram(double *out);
//   double AverageEdge();
   void EdgeTransform();                           /* gradient binarized using otsu threshold */
   double fft2();
   void ChebyshevTransform(int N);
   void ChebyshevFourierTransform2D(double *coeff);
   void Symlet5Transform();
   void GradientMagnitude(int span);
   void GradientDirection2D(int span);
   void PerwittMagnitude2D(ImageMatrix *output);
   void PerwittDirection2D(ImageMatrix *output);
   void ChebyshevStatistics2D(double *coeff, int N, int bins_num);
   int CombFirstFourMoments2D(double *vec);
   void EdgeStatistics(long *EdgeArea, double *MagMean, double *MagMedian, double *MagVar, double *MagHist, double *DirecMean, double *DirecMedian, double *DirecVar, double *DirecHist, double *DirecHomogeneity, double *DiffDirecHist, int num_bins);
   void RadonTransform2D(double *vec);
   double OtsuBinaryMaskTransform();
   int BWlabel(int level);
   void centroid(double *x_centroid, double *y_centroid, double *z_centroid);
   void FeatureStatistics(int *count, int *Euler, double *centroid_x, double *centroid_y, double *centroid_z, int *AreaMin, int *AreaMax,
                                    double *AreaMean, int *AreaMedian, double *AreaVar, int *area_histogram,double *DistMin, double *DistMax,
                                    double *DistMean, double *DistMedian, double *DistVar, int *dist_histogram, int num_bins);
   void GaborFilters2D(double *ratios);
   void HaarlickTexture2D(double distance, double *out);
   void TamuraTexture2D(double *vec);
   void zernike2D(double *zvalues, long *output_size);
};


/* global functions */
HSVcolor RGB2HSV(RGBcolor rgb);
RGBcolor HSV2RGB(HSVcolor hsv);
TColor RGB2COLOR(RGBcolor rgb);
double COLOR2GRAY(TColor color);

#endif

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



#pragma hdrstop

#include <vector>
#include <math.h>
#include <stdio.h>
#include "cmatrix.h"
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
#include "textures/haarlick/haarlick.h"
#include "textures/zernike/zernike.h"

#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <time.h>
#include <sys/time.h>

#ifndef WIN32
#include <stdlib.h>
#include <string.h>
#include <tiffio.h>
#else
#include "libtiff32/tiffio.h"
#endif

#define MIN(a,b) (a<b?a:b)
#define MAX(a,b) (a>b?a:b)

using namespace std;

RGBcolor HSV2RGB(HSVcolor hsv)
{   RGBcolor rgb;
    float R=0, G=0, B=0;
    float H, S, V;
    float i, f, p, q, t;

    H=hsv.hue;
    S=(float)(hsv.saturation)/240;
    V=(float)(hsv.value)/240;
    if(S==0 && H==0) {R=G=B=V;}  /*if S=0 and H is undefined*/
    H=H*(360.0/240.0);
    if(H==360) {H=0;}
    H=H/60;
    i=floor(H);
    f=H-i;
    p=V*(1-S);
    q=V*(1-(S*f));
    t=V*(1-(S*(1-f)));

    if(i==0) {R=V;  G=t;  B=p;}
    if(i==1) {R=q;  G=V;  B=p;}
    if(i==2) {R=p;  G=V;  B=t;}
    if(i==3) {R=p;  G=q;  B=V;}
    if(i==4) {R=t;  G=p;  B=V;}
    if(i==5) {R=V;  G=p;  B=q;}

    rgb.red=(byte)(R*255);
    rgb.green=(byte)(G*255);
    rgb.blue=(byte)(B*255);
    return rgb;
}
//-----------------------------------------------------------------------
HSVcolor RGB2HSV(RGBcolor rgb)
{
  float r,g,b,h,max,min,delta;
  HSVcolor hsv;

  r=(float)(rgb.red) / 255;
  g=(float)(rgb.green) / 255;
  b=(float)(rgb.blue) / 255;

  max = MAX (r, MAX (g, b)), min = MIN (r, MIN (g, b));
  delta = max - min;

  hsv.value = (byte)(max*240.0);
  if (max != 0.0)
    hsv.saturation = (byte)((delta / max)*240.0);
  else
    hsv.saturation = 0;
  if (hsv.saturation == 0) hsv.hue = 0; //-1;
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
    hsv.hue = (byte)(h *(240.0/360.0));
  }
  return(hsv);
}


//--------------------------------------------------------------------------
TColor RGB2COLOR(RGBcolor rgb)
{  return((TColor)(rgb.blue*65536+rgb.green*256+rgb.red));
}

double COLOR2GRAY(TColor color1)
{  double r,g,b;

   r=(byte)(color1 & 0xFF);
   g=(byte)((color1 & 0xFF00)>>8);
   b=(byte)((color1 & 0xFF0000)>>16);

   return((0.3*r+0.59*g+0.11*b));
}


#ifdef WIN32

//--------------------------------------------------------------------------
int ImageMatrix::LoadImage(TPicture *picture,int ColorMode)
{  int a,b,x,y;
   pix_data pix;
   width=picture->Width;
   height=picture->Height;
   depth=1;   /* TPicture is a two-dimentional structure */
   bits=8;
   this->ColorMode=ColorMode;
   /* allocate memory for the image's pixels */
   data=new pix_data[width*height*depth];
   if (!data) return(0); /* memory allocation failed */
   /* load the picture */
   for (y=0;y<height;y++)
     for (x=0;x<width;x++)
     {  pix.clr.RGB.red=(byte)(picture->Bitmap->Canvas->Pixels[x][y] & 0xFF);               /* red value */
        pix.clr.RGB.green=(byte)((picture->Bitmap->Canvas->Pixels[x][y] & 0xFF00) >> 8);    /* green value */
        pix.clr.RGB.blue=(byte)((picture->Bitmap->Canvas->Pixels[x][y] & 0xFF0000) >> 16);  /* blue value */
        if (ColorMode==cmHSV) pix.clr.HSV=RGB2HSV(pix.clr.RGB);
        pix.intensity=COLOR2GRAY(picture->Bitmap->Canvas->Pixels[x][y]);
        set(x,y,0,pix);
     }
   return(1);
}
/*
void ImageMatrix::CmatrixMessage()
{
  printf( "you have succeeded.\n" ); 
}
*/
int ImageMatrix::LoadBMP(char *filename,int ColorMode)
{  TPicture *picture;
   int ret_val;
   picture = new TPicture;
   picture->LoadFromFile(filename);
   ret_val=LoadImage(picture,ColorMode);
   delete picture;
   return(ret_val);
}

#endif



/* LoadTIFF
   filename -char *- full path to the image file
*/
int ImageMatrix::LoadTIFF(char *filename)
{
//#ifndef WIN32
   unsigned long h,w,x,y,z;
   unsigned short int spp,bps;
   TIFF *tif = NULL;
   //tdata_t buf;
   unsigned char *buf8;
   unsigned short *buf16;
   double max_val;
   pix_data pix;
   TIFFSetWarningHandler(NULL);
   if( (tif = TIFFOpen(filename, "r")) )
   {
     TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
     width = w;
     TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
     height = h;
     TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bps);
     bits=bps;
     TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &spp);
     if (!spp) spp=1;  /* assume one sample per pixel if nothing is specified */
     if ((depth=TIFFNumberOfDirectories(tif))<0) return(0);   /* get the number of slices (Zs) */
     
     /* allocate the data */
     data=new pix_data[width*height*depth];
     if (!data) return(0); /* memory allocation failed */

     max_val=pow(2,bits)-1;
     /* read TIFF header and determine image size */
     buf8 = (unsigned char *)_TIFFmalloc(TIFFScanlineSize(tif)*spp);
     buf16 = (unsigned short *)_TIFFmalloc(TIFFScanlineSize(tif)*sizeof(unsigned short)*spp);
     for (z=0;z<depth;z++)
	 {  TIFFSetDirectory(tif,z);
        for (y = 0; y < height; y++)
        {   int col;
            if (bits==8) TIFFReadScanline(tif, buf8, y);
            else TIFFReadScanline(tif, buf16, y);
            x=0;col=0;
            while (x<width)
            { unsigned char byte_data;
              unsigned short short_data;
              double val=0;
              int sample_index;
              for (sample_index=0;sample_index<spp;sample_index++)
              {  byte_data=buf8[col+sample_index];
                 short_data=buf16[col+sample_index];
                 if (bits==8) val=(double)byte_data;
                 else val=(double)(short_data);
                 if (spp==3)  /* RGB image */
                 {  if (sample_index==0) pix.clr.RGB.red=(unsigned char)(255*(val/max_val));
                    if (sample_index==1) pix.clr.RGB.green=(unsigned char)(255*(val/max_val));
                    if (sample_index==2) pix.clr.RGB.blue=(unsigned char)(255*(val/max_val));
                    if ( ColorMode==cmHSV ) pix.clr.HSV=RGB2HSV(pix.clr.RGB);
                 }
              }
              if (spp==3) pix.intensity=COLOR2GRAY(RGB2COLOR(pix.clr.RGB));
              if (spp==1)
              {  pix.clr.RGB.red=(unsigned char)(255*(val/max_val));
                 pix.clr.RGB.green=(unsigned char)(255*(val/max_val));
                 pix.clr.RGB.blue=(unsigned char)(255*(val/max_val));
                 pix.intensity=val;
              }
              set(x,y,z,pix);		   
              x++;
              col+=spp;
            }
        }
     }
     _TIFFfree(buf8);
     _TIFFfree(buf16);
     TIFFClose(tif);
   }
   else return(0);
//#endif
   return(1);
}

/*  SaveTiff
    Save a matrix in TIFF format (16 bits per pixel)
*/
int ImageMatrix::SaveTiff(char *filename)
{
//#ifndef WIN32
   int x,y;
   TIFF* tif = TIFFOpen(filename, "w");
   unsigned short *BufImage16 = new unsigned short[width*height];
   unsigned char *BufImage8 = new unsigned char[width*height];

   if (!tif) return(0);

   for (y = 0; y < height; y++)
     for (x = 0; x < width ; x++)
     {  if (bits==16) BufImage16[x + (y * width)] = (unsigned short)(pixel(x,y,0).intensity);
        else BufImage8[x + (y * width)] = (unsigned char)(pixel(x,y,0).intensity);
     }

   TIFFSetField(tif,TIFFTAG_IMAGEWIDTH, width);
   TIFFSetField(tif,TIFFTAG_IMAGELENGTH, height);
   TIFFSetField(tif, TIFFTAG_PLANARCONFIG,1);
   TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, 1);
   TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, bits);
   TIFFSetField(tif, TIFFTAG_COMPRESSION, 1);
   TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);

   for (y = 0; y < height; y ++)
   {  if (bits==16) TIFFWriteScanline (tif, &(BufImage16[y*width]), y,0 );
      else TIFFWriteScanline (tif, &(BufImage8[y*width]), y,0 );
   }

   TIFFClose(tif);
   delete BufImage8;
   delete BufImage16;
//#endif
   return(1);
}

/* LoadPPM
   filename -char *- full path to the image file
*/

int ImageMatrix::LoadPPM(char *filename, int ColorMode)
{  FILE *fi;
   char ty[256],line[256];
   byte *buffer;
   int x, y, w, h, m=-1;
   pix_data pix;

   fi=fopen(filename,"r");  /* open the file */
   if (!fi) return(0);
   /* PPM header */
   fgets(line,sizeof(line),fi);
   while (line[0]=='#') fgets(line,sizeof(line),fi);
   sscanf(line, "%s %d %d %d", ty, &w, &h, &m);
   if (m>255 || m<=0)
   {  fgets(line,sizeof(line),fi);
      while (line[0]=='#') fgets(line,sizeof(line),fi);
      sscanf(line, "%d %d %d", &w, &h, &m);
      if (m>255 || m<=0)
      {  fgets(line,sizeof(line),fi);
         while (line[0]=='#') fgets(line,sizeof(line),fi);
         sscanf(line, "%d", &m);
      }
   }

   /* allocate memory */
   bits=(unsigned short)ceil(log(m+1)/log(2));
   width=w;
   height=h;
   depth=1;   /* PPM is a 2D format */
   data=new pix_data [width*height*depth];
   if (!data) return(0); /* memory allocation failed */

   /* read the pixel data */
   int res,length,index=0;
   length=width*height;
   if (strcmp(ty, "P6") == 0) length=3*length;
   buffer=new byte[length];
   res=fread(buffer,sizeof(byte),length,fi);
   if (res==0) return(0);
   /* color image */
   if (strcmp(ty, "P6") == 0)
   {  for (y = 0; y < height; y++)
        for (x = 0; x < width; x++)
        {  pix.clr.RGB.red=buffer[index++];
           pix.clr.RGB.green=buffer[index++];
           pix.clr.RGB.blue=buffer[index++];
           pix.intensity=COLOR2GRAY(RGB2COLOR(pix.clr.RGB));
           if (ColorMode==cmHSV) pix.clr.HSV=RGB2HSV(pix.clr.RGB);
		   set(x,y,0,pix);
        }
   }
   else /* greyscale image */
   if (strcmp(ty, "P5") == 0)
   {  for (y = 0; y < height; y++)
        for (x = 0; x < width; x++)
        {  pix.clr.RGB.red=buffer[index];
           pix.clr.RGB.green=buffer[index];
           pix.clr.RGB.blue=buffer[index++];
           pix.intensity=COLOR2GRAY(RGB2COLOR(pix.clr.RGB));
           if (ColorMode==cmHSV) pix.clr.HSV=RGB2HSV(pix.clr.RGB);
		   set(x,y,0,pix);		   
        }
   }
   delete [] buffer;
   fclose(fi);
   return(1);
}


int ImageMatrix::OpenImage(char *image_file_name, int downsample, rect *bounding_rect, double mean, double stddev)
{  
	int res=0;
	#ifdef WIN32
	if (strstr(image_file_name,".bmp") || strstr(image_file_name,".BMP"))
		res=LoadBMP(image_file_name,cmHSV);
	#endif
	if (strstr(image_file_name,".tif") || strstr(image_file_name,".TIF"))
	{  
		res=LoadTIFF(image_file_name);
		// if (res && matrix->bits==16) matrix->to8bits();
	}
	if (strstr(image_file_name,".ppm") || strstr(image_file_name,".PPM"))
		res=LoadPPM(image_file_name,cmHSV);
	if (strstr(image_file_name,".dcm") || strstr(image_file_name,".DCM"))
	{ 
		char buffer[512],temp_filename[64];
		sprintf(temp_filename,"tempimage%d.tif",rand() % 30000);  /* the getpid is to allow several processes to run from the same folder */
		sprintf(buffer,"convert %s %s",image_file_name,temp_filename);  
		system(buffer);
		res=LoadTIFF(temp_filename);
		if (res<=0) printf("Could not convert dcm to tiff\n");
		sprintf(buffer,"rm %s",temp_filename);
		system(buffer);	  
	}
	if (res)  /* add the image only if it was loaded properly */
	{  
		if (bounding_rect && bounding_rect->x>=0)    /* compute features only from an area of the image */
		{ 
			ImageMatrix *temp;
			temp=new ImageMatrix(this,bounding_rect->x,bounding_rect->y,bounding_rect->x+bounding_rect->w-1,bounding_rect->y+bounding_rect->h-1,0,depth-1);
			delete data;
			width=temp->width;height=temp->height;
			if (!(data=new pix_data[width*height*depth])) return(0);  /* allocate new memory */
			memcpy(data,temp->data,width*height*depth*sizeof(pix_data));		 
			//         for (int a=0;a<width*height*depth;a++)
			//		   data[a]=temp->data[a];
			delete temp;
		}
		if (downsample>0 && downsample<100)  /* downsample by a given factor */
			Downsample(((double)downsample)/100.0,((double)downsample)/100.0);   /* downsample the image */
		if (mean>0)  /* normalize to a given mean and standard deviation */
			normalize(-1,-1,-1,mean,stddev);
/*
		// keep track of what file this pixel plane came from
		what_am_i = image_file_name;
		size_t last_slash = what_am_i.find_last_of( '/' );
		if( last_slash != std::string::npos ) 
			what_am_i = what_am_i.substr( last_slash+1 );
	*/
	}
	return(res);
}

/* simple constructors */

ImageMatrix::ImageMatrix()
{
   data=NULL;
   width=0;
   height=0;
   depth=1;
   ColorMode=cmHSV;     /* set a diffult color mode */
}

ImageMatrix::ImageMatrix(int width, int height, int depth)
{  
   bits=8; /* set some default value */
   if (depth<1) depth=1;    /* make sure the image is at least two dimensional */
   ColorMode=cmHSV;
   this->width=width;
   this->height=height;
   this->depth=depth;
   data=new pix_data[width*height*depth];
   memset(data,0,width*height*depth*sizeof(pix_data));  /* initialize */
}

/* create an image which is part of the image
   (x1,y1) - top left
   (x2,y2) - bottom right
*/
ImageMatrix::ImageMatrix(ImageMatrix *matrix,int x1, int y1, int x2, int y2, int z1, int z2)
{  int x,y,z;
   bits=matrix->bits;
   ColorMode=matrix->ColorMode;
   /* verify that the image size is OK */
   if (x1<0) x1=0;
   if (y1<0) y1=0;
   if (z1<0) z1=0;   
   if (x2>=matrix->width) x2=matrix->width-1;
   if (y2>=matrix->height) y2=matrix->height-1;
   if (z2>=matrix->depth) z2=matrix->depth-1;

   width=x2-x1+1;
   height=y2-y1+1;
   depth=z2-z1+1;
   data=new pix_data[width*height*depth];
   
   for (z=z1;z<z1+depth;z++)
     for (y=y1;y<y1+height;y++)
       for (x=x1;x<x1+width;x++)
	     set(x-x1,y-y1,z-z1,matrix->pixel(x,y,z));

	 //what_am_i = matrix->what_am_i;
}

/* free the memory allocated in "ImageMatrix::LoadImage" */
ImageMatrix::~ImageMatrix()
{  if (data) delete [] data;
   data=NULL;
}

/* get a pixel value */
pix_data ImageMatrix::pixel(int x,int y,int z)
{  return(data[z*width*height+y*width+x]);
}

/* assigne a pixel value */
void ImageMatrix::set(int x,int y,int z, pix_data val)
{  data[z*width*height+y*width+x]=val;
}

/* assigne a pixel intensity only */
void ImageMatrix::SetInt(int x,int y,int z, double val)
{  data[z*width*height+y*width+x].intensity=val;
}

/* compute the difference from another image */
void ImageMatrix::diff(ImageMatrix *matrix)
{  int x,y,z;
   for (z=0;z<depth;z++)
     for (y=0;y<height;y++)
       for (x=0;x<width;x++)
       {  pix_data pix1,pix2;
          pix1=pixel(x,y,z);
          pix2=matrix->pixel(x,y,z);
          pix1.intensity=fabs(pix1.intensity-pix2.intensity);
          pix1.clr.RGB.red=(byte)fabs(pix1.clr.RGB.red-pix2.clr.RGB.red);		  
          pix1.clr.RGB.green=(byte)fabs(pix1.clr.RGB.green-pix2.clr.RGB.green);		  
          pix1.clr.RGB.blue=(byte)fabs(pix1.clr.RGB.blue-pix2.clr.RGB.blue);
          set(x,y,z,pix1);
       }
}


/* duplicate
   create another matrix the same as the first
*/
ImageMatrix* ImageMatrix::duplicate()
{ 
	ImageMatrix *new_matrix;
	new_matrix=new ImageMatrix;
	new_matrix->data=new pix_data[width*height*depth];
	if (!(new_matrix->data)) {
		fprintf (stderr,"Could not allocate memory for duplicate image\n");
		return(NULL); /* memory allocation failed */
	}
	new_matrix->width=width;
	new_matrix->height=height;
	new_matrix->depth=depth;
	new_matrix->bits=bits;
	new_matrix->ColorMode=ColorMode;
	memcpy(new_matrix->data,data,width*height*depth*sizeof(pix_data));

	//new_matrix->what_am_i = what_am_i;
	return(new_matrix);
}
/*
void ImageMatrix::dump( )
{
	
	time_t ltime;
	struct tm *Tm;
	struct timeval detail_time;

	ltime=time(NULL);
	Tm=localtime(&ltime);

	std::ostringstream filename;
	filename << what_am_i << "_" << (Tm->tm_year+1900) << '-' 
		<< setw(2) << setfill('0') << (Tm->tm_mon+1) << '-' 
		<< setw(2) << setfill('0') << Tm->tm_mday << "_" 
		<< setw(2) << setfill('0') << Tm->tm_hour <<  '-'
		<< setw(2) << setfill('0') << Tm->tm_min << '-'
		<< setw(2) << setfill('0') << Tm->tm_sec << ".";
	gettimeofday(&detail_time,NULL);
	filename << setw(2) <<  detail_time.tv_usec / 1000;
	filename << "pixel_dump.txt";
	std::cout << filename.str() << std::endl;
	ofstream pixel_dump_file ( filename.str().c_str(), ios::app );

	int count = 0;
	int numpix = width*height*depth;
	while (count < numpix) {
		pixel_dump_file << data[count].intensity << std::endl;
		++count;
	}
	pixel_dump_file.close();
}
*/

/* to8bits
   convert a 16 bit matrix to 8 bits
*/
void ImageMatrix::to8bits()
{  double max_val;
   if (bits==8) return;
   max_val=pow(2,bits)-1;
   bits=8;
   for (int a=0;a<width*height*depth;a++)
     data[a].intensity=255*(data[a].intensity/max_val);
}

/* flip
   flip an image horizontaly
*/
void ImageMatrix::flip()
{  int x,y,z;
   pix_data temp;
   for (z=0;z<depth;z++)   
     for (y=0;y<height;y++)
       for (x=0;x<width/2;x++)
       {  temp=pixel(x,y,z);
          set(x,y,z,pixel(width-x-1,y,z));
          set(width-x-1,y,z,temp);
       }
}

void ImageMatrix::invert()
{  	   
   for (int a=0;a<width*height*depth;a++)
     data[a].intensity=(pow(2,bits)-1)-data[a].intensity;
}

/* Downsample
   down sample an image
   x_ratio, y_ratio -double- (0 to 1) the size of the new image comparing to the old one
*/
void ImageMatrix::Downsample(double x_ratio, double y_ratio)
{  double x,y,dx,dy,frac;
   int new_x,new_y,a;
   pix_data pix;
   if (x_ratio>1) x_ratio=1;
   if (y_ratio>1) y_ratio=1;
   dx=1/x_ratio;
   dy=1/y_ratio;

   if (dx==1 && dy==1) return;   /* nothing to scale */
   
   for (int z=0;z<depth;z++)
   {  /* first downsample x */
      for (new_y=0;new_y<height;new_y++)
      { x=0;
        new_x=0;
        while (x<width)
        {  double sum_i=0;
           double sum_r=0;
           double sum_g=0;
           double sum_b=0;

           /* the leftmost fraction of pixel */
           a=(int)(floor(x));
           frac=ceil(x)-x;
           if (frac>0 && a<width)
           {  sum_i+=pixel(a,new_y,z).intensity*frac;
              sum_r+=pixel(a,new_y,z).clr.RGB.red*frac;
              sum_g+=pixel(a,new_y,z).clr.RGB.green*frac;
              sum_b+=pixel(a,new_y,z).clr.RGB.blue*frac;
           } 
           /* the middle full pixels */
           for (a=(int)(ceil(x));a<floor(x+dx);a=a+1)
           if (a<width)
           {  sum_i+=pixel(a,new_y,z).intensity;
              sum_r+=pixel(a,new_y,z).clr.RGB.red;
              sum_g+=pixel(a,new_y,z).clr.RGB.green;
              sum_b+=pixel(a,new_y,z).clr.RGB.blue;
           }
           /* the right fraction of pixel */
           frac=x+dx-floor(x+dx);
           if (frac>0 && a<width)
           {  sum_i+=pixel(a,new_y,z).intensity*frac;
              sum_r+=pixel(a,new_y,z).clr.RGB.red*frac;
              sum_g+=pixel(a,new_y,z).clr.RGB.green*frac;
              sum_b+=pixel(a,new_y,z).clr.RGB.blue*frac;
           }

           pix.intensity=sum_i/(dx);
           pix.clr.RGB.red=(byte)(sum_r/(dx));
           pix.clr.RGB.green=(byte)(sum_g/(dx));
           pix.clr.RGB.blue=(byte)(sum_b/(dx));
           set(new_x,new_y,z,pix);
		   
           x+=dx;
           new_x++;
        }
      }
 
      /* downsample y */
      for (new_x=0;new_x<x_ratio*width;new_x++)
      { y=0;
        new_y=0;
        while (y<height)
        {  double sum_i=0;
           double sum_r=0;
           double sum_g=0;
           double sum_b=0;

           a=(int)(floor(y));
           frac=ceil(y)-y;
           if (frac>0 && a<height)   /* take also the part of the leftmost pixel (if needed) */
           {  sum_i+=pixel(new_x,a,z).intensity*frac;
              sum_r+=pixel(new_x,a,z).clr.RGB.red*frac;
              sum_g+=pixel(new_x,a,z).clr.RGB.green*frac;
              sum_b+=pixel(new_x,a,z).clr.RGB.blue*frac;
           }
           for (a=(int)(ceil(y));a<floor(y+dy);a=a+1)
           if (a<height)
           {  sum_i+=pixel(new_x,a,z).intensity;
              sum_r+=pixel(new_x,a,z).clr.RGB.red;
              sum_g+=pixel(new_x,a,z).clr.RGB.green;
              sum_b+=pixel(new_x,a,z).clr.RGB.blue;
           }
           frac=y+dy-floor(y+dy);
           if (frac>0 && a<height)
           {  sum_i+=pixel(new_x,a,z).intensity*frac;
              sum_r+=pixel(new_x,a,z).clr.RGB.red*frac;
              sum_g+=pixel(new_x,a,z).clr.RGB.green*frac;
              sum_b+=pixel(new_x,a,z).clr.RGB.blue*frac;
           }

           pix.intensity=sum_i/(dy);
           pix.clr.RGB.red=(byte)(sum_r/(dy));
           pix.clr.RGB.green=(byte)(sum_g/(dy));
           pix.clr.RGB.blue=(byte)(sum_b/(dy));
           set(new_x,new_y,z,pix);

           y+=dy;
           new_y++;
        }
      }
   }

   width=(int)(x_ratio*width);
   height=(int)(y_ratio*height);
}


/* Rotate
   Rotate an image by 90, 120, or 270 degrees
   angle -double- (0 to 360) the degrees of rotation.  Only values of 90, 180, 270 are currently allowed
*/
ImageMatrix* ImageMatrix::Rotate(double angle) {
	ImageMatrix *new_matrix;
	int new_x,new_y,new_width,new_height;
	pix_data pix;

	// Only deal with right angles
	if (! ( (angle == 90) || (angle == 180) || (angle == 270) ) ) return (this);

	// Make a new image matrix
	new_matrix=new ImageMatrix;
	new_matrix->data=new pix_data[width*height*depth];
	if (!(new_matrix->data)) {
		fprintf (stderr,"Could not allocate memory for duplicate image\n");
		return(NULL); /* memory allocation failed */
	}
	new_matrix->depth=depth;
	new_matrix->bits=bits;
	new_matrix->ColorMode=ColorMode;
	// switch width/height if 90 or 270
	if ( (angle == 90) || (angle == 270) ) {
		new_matrix->width = new_width = height;
		new_matrix->height = new_height = width;
	} else {
		new_matrix->width = new_width = width;
		new_matrix->height = new_height = height;
	}

	// write the new pixels
	for (int z=0;z<depth;z++) {
		for (new_y=0;new_y<new_height;new_y++) {
			new_x=0;
			for (new_x=0;new_x<new_width;new_x++) {
				switch ((int)angle) {
					case 90:
						pix=pixel(new_y,height-new_x-1,z);
					break;

					case 180:
						pix=pixel(width-new_x-1,height-new_y-1,z);
					break;

					case 270:
						pix=pixel(width-new_y-1,new_x,z);
					break;
				}
	
				new_matrix->set(new_x,new_y,z,pix);
			} /* for x */
		} /* for y */
   } /* for z */

	return(new_matrix);
}


/* find basic intensity statistics */

int compare_doubles (const void *a, const void *b)
{
  if (*((double *)a) > *((double*)b)) return(1);
  if (*((double*)a) == *((double*)b)) return(0);
  return(-1);
}

/* BasicStatistics
   get basic statistical properties of the intensity of the image
   mean -double *- pre-allocated one double for the mean intensity of the image
   median -double *- pre-allocated one double for the median intensity of the image
   std -double *- pre-allocated one double for the standard deviation of the intensity of the image
   min -double *- pre-allocated one double for the minimum intensity of the image
   max -double *- pre-allocated one double for the maximal intensity of the image
   histogram -double *- a pre-allocated vector for the histogram. If NULL then histogram is not calculated
   nbins -int- the number of bins for the histogram
   
   if one of the pointers is NULL, the corresponding value is not computed.
*/
void ImageMatrix::BasicStatistics(double *mean, double *median, double *std, double *min, double *max, double *hist, int bins)
{  long pixel_index,num_pixels;
   double *pixels;
   double min1=INF,max1=-INF,mean_sum=0;
   
   num_pixels=height*width*depth;
   pixels=new double[num_pixels];

   /* compute the average, min and max */
   for (pixel_index=0;pixel_index<num_pixels;pixel_index++)
   {  mean_sum+=data[pixel_index].intensity;
      if (data[pixel_index].intensity>max1)
        max1=data[pixel_index].intensity;
      if (data[pixel_index].intensity<min1)
        min1=data[pixel_index].intensity;
      pixels[pixel_index]=data[pixel_index].intensity;		
   }
   if (max) *max=max1;	 
   if (min) *min=min1;
   if (mean || std) *mean=mean_sum/num_pixels;

   /* calculate the standard deviation */
   if (std)
   {  *std=0;
      for (pixel_index=0;pixel_index<num_pixels;pixel_index++)
        *std=*std+pow(data[pixel_index].intensity-*mean,2);
      *std=sqrt(*std/(num_pixels-1));
   }
   
   if (hist)  /* do the histogram only if needed */
     histogram(hist,bins,0);

   /* find the median */
   if (median)
   {  qsort(pixels,num_pixels,sizeof(double),compare_doubles);
      *median=pixels[num_pixels/2];
   }
   delete [] pixels;	     
}

/* normalize the pixel values into a given range 
   min -double- the min pixel value (ignored if <0)
   max -double- the max pixel value (ignored if <0)
   range -long- nominal dynamic range (ignored if <0)
   mean -double- the mean of the normalized image (ignored if <0)
   stddev -double- the stddev of the normalized image (ignored if <0)
*/
void ImageMatrix::normalize(double min, double max, long range, double mean, double stddev)
{  long x;
   /* normalized to min and max */
   if (min>=0 && max>0 && range>0)
     for (x=0;x<width*height*depth;x++)
     {  if (data[x].intensity<min) data[x].intensity=0;
        else if (data[x].intensity>max) data[x].intensity=range;
        else data[x].intensity=((data[x].intensity-min)/(max-min))*range;
	 }

    /* normalize to mean and stddev */
	if (mean>0)
	{   double original_mean,original_stddev;
        BasicStatistics(&original_mean, NULL, &original_stddev, NULL, NULL, NULL, 0);	
        for (x=0;x<width*height*depth;x++)		
        {  data[x].intensity-=(original_mean-mean);
           if (stddev>0)
             data[x].intensity=mean+(data[x].intensity-mean)*(stddev/original_stddev);
           if (data[x].intensity<0) data[x].intensity=0;
           if (data[x].intensity>pow(2,bits)-1) data[x].intensity=pow(2,bits)-1;  
        }
//BasicStatistics(&original_mean, NULL, &original_stddev, NULL, NULL, NULL, 0);		  
//printf("%f %f\n",original_mean,original_stddev);
	}	   
}

/* convolve
*/
void ImageMatrix::convolve(ImageMatrix *filter)
{ int x,y,z;
  ImageMatrix *copy;
  int height2=filter->height/2;
  int width2=filter->width/2;
  int depth2=filter->depth/2;
  double tmp;

  copy=duplicate();
  for (z=0;z<depth;z++)
    for(x=0;x<width;++x)
       for(y=0;y<height;++y)
       {
         tmp=0.0;
         for(int i=-width2;i<=width2;++i)
         {  int xx=x+i;
            if(xx<width && xx >= 0) {
            for(int j=-height2;j<=height2;++j) {
              int yy=y+j;
              if(int(yy)>=0 && yy < height) {
			  
for (int k=-depth2;k<=depth2;++k) {
  int zz=z+k;
  if (int(zz)>=0 && zz < depth) {
			  
                tmp+=filter->pixel(i+width2,j+height2,k+depth2).intensity*copy->pixel(xx,yy,zz).intensity;
} }				
              }
            }
          }
        }
        SetInt(x,y,z,tmp);
      }
  delete copy;
}

/* find the basic color statistics
   hue_avg -double *- average hue
   hue_std -double *- standard deviation of the hue
   sat_avg -double *- average saturation
   sat_std -double *- standard deviation of the saturation
   val_avg -double *- average value
   val_std -double *- standard deviation of the value
   max_color -double *- the most popular color
   colors -double *- a histogram of colors
   if values are NULL - the value is not computed
*/

void ImageMatrix::GetColorStatistics(double *hue_avg, double *hue_std, double *sat_avg, double *sat_std, double *val_avg, double *val_std, double *max_color, double *colors)
{  long a,color_index;
   color hsv;
   double max,pixel_num;
   float certainties[COLORS_NUM+1];

   pixel_num=height*width;

   /* calculate the average hue, saturation, value */
   if (hue_avg) *hue_avg=0;
   if (sat_avg) *sat_avg=0;
   if (val_avg) *val_avg=0;
   if (colors)
     for (a=0;a<=COLORS_NUM;a++)
       colors[a]=0;

   for (a=0;a<width*height*depth;a++)
   {  hsv=data[a].clr;
      if (hue_avg) *hue_avg+=hsv.HSV.hue;
      if (sat_avg) *sat_avg+=hsv.HSV.saturation;
      if (val_avg) *val_avg+=hsv.HSV.value;
       color_index=FindColor(hsv.HSV.hue,hsv.HSV.saturation,hsv.HSV.value,certainties);
       colors[color_index]+=1;
   }
   *hue_avg=*hue_avg/pixel_num;
   *sat_avg=*sat_avg/pixel_num;
   *val_avg=*val_avg/pixel_num;

   /* max color (the most common color in the image) */
   if (max_color)
   {  *max_color=0;
      max=0.0;
      for (a=0;a<=COLORS_NUM;a++)
        if (colors[a]>max)
        {  max=colors[a];
           *max_color=a;
        }
   }
   /* colors */
   if (colors)
     for (a=0;a<=COLORS_NUM;a++)
       colors[a]=colors[a]/pixel_num;

   /* standard deviation of hue, saturation and value */
   if (hue_std) *hue_std=0;
   if (sat_std) *sat_std=0;
   if (val_std) *val_std=0;
   for (a=0;a<width*height*depth;a++)
   {  hsv=data[a].clr;
      if (hue_std && hue_avg) *hue_std+=pow(hsv.HSV.hue-*hue_avg,2);
      if (sat_std && sat_avg) *sat_std+=pow(hsv.HSV.saturation-*sat_avg,2);
      if (val_std && val_avg) *val_std+=pow(hsv.HSV.value-*val_avg,2);
   }
   if (hue_std && hue_avg) *hue_std=sqrt(*hue_std/pixel_num);
   if (sat_std && sat_avg) *sat_std=sqrt(*sat_std/pixel_num);
   if (val_std && val_avg) *val_std=sqrt(*val_std/pixel_num);
}

/* ColorTransform
   Transform a color image to a greyscale image such that each
   color_hist -double *- a histogram (of COLOR_NUM + 1 bins) of the colors. This parameter is ignored if NULL
   use_hue -int- 0 if classifying colors, 1 if using the hue component of the HSV vector
   grey level represents a different color
*/
void ImageMatrix::ColorTransform(double *color_hist, int use_hue)
{  
	long x,y,z; //,base_color;
	double cb_intensity;
	double max_value;
	pix_data hsv_pixel;
	int color_index=0;   
	RGBcolor rgb;
	float certainties[COLORS_NUM+1];
	max_value=pow(2,bits)-1;
	// initialize the color histogram
	if( color_hist ) 
		for( color_index = 0; color_index <= COLORS_NUM; color_index++ )
			color_hist[color_index]=0;
	// find the colors
	for( z = 0; z < depth; z++ )
		for( y = 0; y < height; y++ )
			for( x = 0; x < width; x++ )
			{ 
				hsv_pixel = pixel( x, y, z );
				if( use_hue == 0 ) // not using hue
				{  
					color_index = FindColor( hsv_pixel.clr.HSV.hue,
					                         hsv_pixel.clr.HSV.saturation,
					                         hsv_pixel.clr.HSV.value,
					                         certainties );
					if( color_hist )
						color_hist[ color_index ] += 1;
					// convert the color index to a greyscale value
					cb_intensity = int( ( max_value * color_index ) / COLORS_NUM );
				} 
				else // using hue
					cb_intensity = hsv_pixel.clr.HSV.hue;
				rgb.red = rgb.green = rgb.blue = (byte)( 255 * ( cb_intensity / max_value ) );
				hsv_pixel.clr.HSV = RGB2HSV( rgb );
				hsv_pixel.intensity = cb_intensity;
				set( x, y, z, hsv_pixel );
			}
	/* normalize the color histogram */
	if (color_hist) 
		for (color_index=0;color_index<=COLORS_NUM;color_index++)
			color_hist[color_index]/=(width*height);	 
}

/* get image histogram */
void ImageMatrix::histogram(double *bins,unsigned short bins_num, int imhist)
{  long a;
   double min=INF,max=-INF;
   /* find the minimum and maximum */
   if (imhist==1)    /* similar to the Matlab imhist */
   {  min=0;
      max=pow(2,bits)-1;
   }
   else
   {  for (a=0;a<width*height*depth;a++)
      {  if (data[a].intensity>max)
           max=data[a].intensity;
         if (data[a].intensity<min)
           min=data[a].intensity;
      }
   }
   /* initialize the bins */
   for (a=0;a<bins_num;a++)
     bins[a]=0;

   /* build the histogram */
   for (a=0;a<width*height*depth;a++)   
     if (data[a].intensity==max) bins[bins_num-1]+=1;
     else bins[(int)(((data[a].intensity-min)/(max-min))*bins_num)]+=1;
	
   return;
}

/* fft 2 dimensional transform */
// http://www.fftw.org/doc/
double ImageMatrix::fft2()
{  fftw_complex *out;
   double *in;
   fftw_plan p;
   long x,y,z;

   in = (double*) fftw_malloc(sizeof(double) * width*height*depth);
   out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * width*height*depth);

   if (depth==1) p=fftw_plan_dft_r2c_2d(width,height,in,out, FFTW_MEASURE /* FFTW_ESTIMATE */);
   else p=fftw_plan_dft_r2c_3d(width,height,depth,in,out, FFTW_MEASURE /* FFTW_ESTIMATE */);

   for (x=0;x<width;x++)
     for (y=0;y<height;y++)
       for (z=0;z<depth;z++)	   
         in[depth*height*x+y*depth+z]=pixel(x,y,z).intensity;

   fftw_execute(p); /* execute the transformation (repeat as needed) */

   if (depth==1)
   {  long half_height=height/2+1;   /* (to 56 including 56 ) */
      /* find the abs and angle */
      for (x=0;x<width;x++)
        for (y=0;y<half_height;y++)
          SetInt(x,y,0,sqrt(pow(out[half_height*x+y][0],2)+pow(out[half_height*x+y][1],2)));    /* sqrt(real(X).^2 + imag(X).^2) */
   
     /* complete the first column */
     for (y=half_height;y<height;y++)
       SetInt(0,y,0,pixel(0,height-y,0).intensity);

     /* complete the rows */
     for (y=half_height;y<height;y++)
       for (x=1;x<width;x++)   /* 1 because the first column is already completed */
         SetInt(x,y,0,pixel(width-x,height-y,0).intensity);
   }
   else
   {  long half_depth=depth/2+1; 
      /* find the abs and angle */
      for (x=0;x<width;x++)
        for (y=0;y<height;y++)
          for (z=0;z<half_depth;z++)		
            SetInt(x,y,z,sqrt(pow(out[height*half_depth*x+half_depth*y+z][0],2)+pow(out[height*half_depth*x+half_depth*y+z][1],2)));    /* sqrt(real(X).^2 + imag(X).^2) */
   
      /* compute the first z */
      for (z=half_depth;z<depth;z++)
        for (y=0;y<height;y++)
          SetInt(0,y,z,pixel(0,y,depth-z).intensity);

      /* complete the rows */
	  for (z=half_depth;z<depth;z++)
        for (y=1;y<height;y++)   /* 1 because the first column is already completed */
          for (x=0;x<width;x++)   
            SetInt(x,y,z,pixel(width-x,height-y,depth-z).intensity);	  
   }
   
   fftw_destroy_plan(p);
   fftw_free(in);
   fftw_free(out);

   /* calculate the magnitude and angle */

   return(0);
}

/* chebyshev transform */
void ImageMatrix::ChebyshevTransform(int N)
{  double *out;
   int x,y,old_width;

   if (N<2)
		 N = MIN( width, height );
   out=new double[height*N];
   Chebyshev2D(this, out,N);

   old_width=width;  /* keep the old width to free the memory */
   width=N;
   height = MIN( height, N );   /* prevent error */

   for(y=0;y<height;y++)
     for(x=0;x<width;x++)
       SetInt(x,y,0,out[y*width+x]);
   delete [] out;
}

/* chebyshev transform
   coeff -array of double- a pre-allocated array of 32 doubles
*/
void ImageMatrix::ChebyshevFourierTransform2D(double *coeff)
{  ImageMatrix *matrix;
   matrix=duplicate();
   if( (width * height) > (300 * 300) )
		 matrix->Downsample( MIN( 300.0/(double)width, 300.0/(double)height ), MIN( 300.0/(double)width, 300.0/(double)height ) );  /* downsample for avoiding memory problems */
   ChebyshevFourier2D(matrix, 0, coeff,32);
   delete matrix;
}


/* Symlet5 transform */
void ImageMatrix::Symlet5Transform()
{  long x,y,z;
   DataGrid2D *grid2d=NULL;
   DataGrid3D *grid3d=NULL;
   DataGrid *grid;
   Symlet5 *Sym5;

   if (depth==1) grid2d = new DataGrid2D(width,height,-1);
   else grid3d=new DataGrid3D(width,height,depth);
   if (grid2d) grid=grid2d;
   else grid=grid3d;
   
   for (z=0;z<depth;z++)   
     for (y=0;y<height;y++)
       for(x=0;x<width;x++)
         grid->setData(x,y,z-(depth==1),pixel(x,y,z).intensity);
   Sym5=new Symlet5(0,1);
   if (depth==1) Sym5->transform2D(grid);
   else Sym5->transform3D(grid);

   delete [] data; /* free the old memory of the matrix */

   /* allocate new memory (new dimensions) and copy the values */
   width=grid->getX();
   height=grid->getY();
   if (depth>1) depth=grid->getZ();
   
   data=new pix_data[width*height*depth];
   for (z=0;z<depth;z++)   
     for (y=0;y<height;y++)
       for(x=0;x<width;x++)
         SetInt(x,y,z,grid->getData(x,y,z-(depth==1)));
		 
   delete Sym5;
   if (grid2d) delete grid2d;
   if (grid3d) delete grid3d;   
}

/* chebyshev statistics
   coeff -array of double- pre-allocated memory of 20 doubles
   nibs_num - (32 is normal)
*/
void ImageMatrix::ChebyshevStatistics2D(double *coeff, int N, int bins_num)
{
   if (N<2) N=20;
   if (N>MIN(width,height)) N=MIN(width,height);   
   ChebyshevTransform(N);
   histogram(coeff,bins_num,0);
}

/* CombFirstFourMoments
   vec should be pre-alocated array of 48 doubles
*/
int ImageMatrix::CombFirstFourMoments2D(double *vec)
{  int count;
   ImageMatrix *matrix;
   if (bits==16) 
   {  matrix=this->duplicate();
      matrix->to8bits();
   }
   else matrix=this;
   count=CombFirst4Moments2D(matrix, vec);   
   vd_Comb4Moments(vec);   
   if (bits==16) delete matrix;
   return(count);
}

/* Edge Transform */
void ImageMatrix::EdgeTransform()
{  long x,y,z;
   ImageMatrix *TempMatrix;
   TempMatrix=duplicate();
   for (y=0;y<TempMatrix->height;y++)
     for (x=0;x<TempMatrix->width;x++)
       for (z=0;z<TempMatrix->depth;z++)	 
       {  double max_x=0,max_y=0,max_z=0;
          if (y>0 && y<height-1) max_y=MAX(fabs(TempMatrix->pixel(x,y,z).intensity-TempMatrix->pixel(x,y-1,z).intensity),fabs(TempMatrix->pixel(x,y,z).intensity-TempMatrix->pixel(x,y+1,z).intensity));
          if (x>0 && x<width-1) max_x=MAX(fabs(TempMatrix->pixel(x,y,z).intensity-TempMatrix->pixel(x-1,y,z).intensity),fabs(TempMatrix->pixel(x,y,z).intensity-TempMatrix->pixel(x+1,y,z).intensity));
          if (z>0 && z<depth-1) max_z=MAX(fabs(TempMatrix->pixel(x,y,z).intensity-TempMatrix->pixel(x,y,z-1).intensity),fabs(TempMatrix->pixel(x,y,z).intensity-TempMatrix->pixel(x,y,z+1).intensity));
          SetInt(x,y,z,MAX(MAX(max_x,max_z),max_y));
       }

   /* use otsu global threshold to set edges to 0 or 1 */
/*   
   double OtsuGlobalThreshold,max_val;
   max_val=pow(2,bits)-1;
   OtsuGlobalThreshold=Otsu();
   for (y=0;y<height;y++)
     for (x=0;x<width;x++)
       if (data[x][y].intensity>OtsuGlobalThreshold*max_val) data[x][y].intensity=max_val;
       else data[x][y].intensity=0;
*/
   delete TempMatrix;
}

/* transform by gradient magnitude */
void ImageMatrix::GradientMagnitude(int span)
{  long x,y,z;
   //double sum;
   if (span==0) span=2;  /* make sure 0 is not a default */
   for (x=0;x<width-span;x++)
     for (y=0;y<height-span;y++)
       for (z=0;z<depth;z++)
       {  double sum=pow(pixel(x+span,y,z).intensity-pixel(x,y,z).intensity,2)+pow(pixel(x,y+span,z).intensity-pixel(x,y,z).intensity,2);
          if (z<depth-span) sum+=pow(pixel(x,y,z+span).intensity-pixel(x,y,z).intensity,2);
         SetInt(x,y,z,sqrt(sum));
       }
}

/* transform by gradient direction */
void ImageMatrix::GradientDirection2D(int span)
{  long x,y;
   if (span==0) span=2;  /* make sure 0 is not a default */
   for (x=0;x<width-span;x++)
     for (y=0;y<height-span;y++)
     {  if (pixel(x,y+span,0).intensity-pixel(x,y,0).intensity==0)
          SetInt(x,y,0,0);
          else SetInt(x,y,0,atan2(pixel(x+span,y,0).intensity-pixel(x,y,0).intensity,pixel(x,y+span,0).intensity-pixel(x,y,0).intensity));
     }
}

/* Perwitt gradient magnitude
   output - a pre-allocated matrix that will hold the output (the input matrix is not changed)
            output should be of the same size as the input matrix
*/
void ImageMatrix::PerwittMagnitude2D(ImageMatrix *output)
{  long x,y,z,i,j;
   double sumx,sumy;
   for (x=0;x<width;x++)
     for (y=0;y<height;y++)
       for (z=0;z<depth;z++)
       {  sumx=0;
          sumy=0;		  
          for (j=y-1;j<=y+1;j++)
            if (j>=0 && j<height && x-1>=0)
               sumx+=pixel(x-1,j,z).intensity*1;//0.3333;
          for (j=y-1;j<=y+1;j++)
            if (j>=0 && j<height && x+1<width)
              sumx+=pixel(x+1,j,z).intensity*-1;//-0.3333;
          for (i=x-1;i<=x+1;i++)
            if (i>=0 && i<width && y-1>=0)
              sumy+=pixel(i,y-1,z).intensity*1;//-0.3333;
          for (i=x-1;i<=x+1;i++)
            if (i>=0 && i<width && y+1<height)
              sumy+=pixel(i,y+1,z).intensity*-1;//0.3333;
          output->SetInt(x,y,z,sqrt(sumx*sumx+sumy*sumy));
       }
}

/* Perwitt gradient direction
   output - a pre-allocated matrix that will hold the output (the input matrix is not changed)
            output should be of the same size as the input matrix
*/
void ImageMatrix::PerwittDirection2D(ImageMatrix *output)
{  long x,y,z,i,j;
   double sumx,sumy;
   for (x=0;x<width;x++)
     for (y=0;y<height;y++)
       for (z=0;z<depth;z++)	 
       {  sumx=0;
          sumy=0;
          for (j=y-1;j<=y+1;j++)
            if (j>=0 && j<height && x-1>=0)
               sumx+=pixel(x-1,j,0).intensity*1;//0.3333;
          for (j=y-1;j<=y+1;j++)
            if (j>=0 && j<height && x+1<width)
              sumx+=pixel(x+1,j,0).intensity*-1;//-0.3333;
          for (i=x-1;i<=x+1;i++)
            if (i>=0 && i<width && y-1>=0)
              sumy+=pixel(i,y-1,0).intensity*1;//-0.3333;
          for (i=x-1;i<=x+1;i++)
            if (i>=0 && i<width && y+1<height)
              sumy+=pixel(i,y+1,0).intensity*-1;//0.3333;
          if (sumy==0 || fabs(sumy)<1/INF) output->SetInt(x,y,z,3.1415926*(sumx<0));
          else output->SetInt(x,y,z,atan2(sumy,sumx));
       }
}


/* edge statistics */
//#define NUM_BINS 8
//#define NUM_BINS_HALF 4
/* EdgeArea -long- number of edge pixels
   MagMean -double- mean of the gradient magnitude
   MagMedian -double- median of the gradient magnitude
   MagVar -double- variance of the gradient magnitude
   MagHist -array of double- histogram of the gradient magnitude. array of size "num_bins" should be allocated before calling the function
   DirecMean -double- mean of the gradient direction
   DirecMedian -double- median of the gradient direction
   DirecVar -double- variance of the gradient direction
   DirecHist -array of double- histogram of the gradient direction. array of size "num_bins" should be allocated before calling the function
   DirecHomogeneity -double-
   DiffDirecHist -array of double- array of size num_bins/2 should be allocated
*/

void ImageMatrix::EdgeStatistics(long *EdgeArea, double *MagMean, double *MagMedian, double *MagVar, double *MagHist, double *DirecMean, double *DirecMedian, double *DirecVar, double *DirecHist, double *DirecHomogeneity, double *DiffDirecHist, int num_bins)
{  ImageMatrix *GradientMagnitude,*GradientDirection;
   long a,bin_index;
   double min,max,sum,max_intensity;
   
   max_intensity=pow(bits,2)-1;
   
   GradientMagnitude=duplicate();
   PerwittMagnitude2D(GradientMagnitude);
   GradientDirection=duplicate();
   PerwittDirection2D(GradientDirection);

   /* find gradient statistics */
   GradientMagnitude->BasicStatistics(MagMean, MagMedian, MagVar, &min, &max, MagHist, num_bins);
   *MagVar=pow(*MagVar,2);

   /* find the edge area (number of edge pixels) */
   *EdgeArea=0;
//   level=min+(max-min)/2;   // level=duplicate->OtsuBinaryMaskTransform()   // level=MagMean

   for (a=0;a<GradientMagnitude->height*GradientMagnitude->width*GradientMagnitude->depth;a++)
     if (GradientMagnitude->data[a].intensity>max_intensity*0.5) (*EdgeArea)+=1; /* find the edge area */
//   GradientMagnitude->OtsuBinaryMaskTransform();

   /* find direction statistics */
   GradientDirection->BasicStatistics(DirecMean, DirecMedian, DirecVar, &min, &max, DirecHist, num_bins);
   *DirecVar=pow(*DirecVar,2);

   /* Calculate statistics about edge difference direction
      Histogram created by computing differences amongst histogram bins at angle and angle+pi */
   for (bin_index=0;bin_index<(int)(num_bins/2);bin_index++)
      DiffDirecHist[bin_index]=fabs(DirecHist[bin_index]-DirecHist[bin_index+(int)(num_bins/2)]);
   sum=0;
   for (bin_index=0;bin_index<(int)(num_bins/2);bin_index++)
   {  if (DirecHist[bin_index]+DirecHist[bin_index+(int)(num_bins/2)]!=0)  /* protect from a numeric flaw */
        DiffDirecHist[bin_index]=DiffDirecHist[bin_index]/(DirecHist[bin_index]+DirecHist[bin_index+(int)(num_bins/2)]);
      sum+=(DirecHist[bin_index]+DirecHist[bin_index+(int)(num_bins/2)]);
   }

   /* The fraction of edge pixels that are in the first two bins of the histogram measure edge homogeneity */
   if (sum>0) *DirecHomogeneity = (DirecHist[0]+DirecHist[1])/sum;

   delete GradientMagnitude;
   delete GradientDirection;
}

/* radon transform
   vec -array of double- output column. a pre-allocated vector of the size 3*4=12
*/
void ImageMatrix::RadonTransform2D(double *vec)
{   int x,y,val_index,output_size,vec_index,bin_index;
    double *pixels,*ptr,bins[3];
    int angle,num_angles=4;
    double theta[4]={0,45,90,135};
    //double min,max;
    int rLast,rFirst;
    rLast = (int) ceil(sqrt(pow(width-1-(width-1)/2,2)+pow(height-1-(height-1)/2,2))) + 1;
    rFirst = -rLast;
    output_size=rLast-rFirst+1;

    ptr=new double[output_size*num_angles];
    for (val_index=0;val_index<output_size*num_angles;val_index++)
      ptr[val_index]=0;  /* initialize the output vector */

    pixels=new double[width*height];
    vec_index=0;

    for (x=0;x<width;x++)
      for (y=0;y<height;y++)
        pixels[y+height*x]=pixel(x,y,0).intensity;

    radon(ptr,pixels, theta, height, width, (width-1)/2, (height-1)/2, num_angles, rFirst, output_size);

    for (angle=0;angle<num_angles;angle++)
    {  //radon(ptr,pixels, &theta, height, width, (width-1)/2, (height-1)/2, 1, rFirst, output_size);
       /* create histogram */
       double min=INF,max=-INF;
       /* find the minimum and maximum values */
       for (val_index=angle*output_size;val_index<(angle+1)*output_size;val_index++)
       {  if (ptr[val_index]>max) max=ptr[val_index];
          if (ptr[val_index]<min) min=ptr[val_index];
       }

       for (val_index=0;val_index<3;val_index++)   /* initialize the bins */
         bins[val_index]=0;
       for (val_index=angle*output_size;val_index<(angle+1)*output_size;val_index++)
         if (ptr[val_index]==max) bins[2]+=1;
         else bins[(int)(((ptr[val_index]-min)/(max-min))*3)]+=1;

       for (bin_index=0;bin_index<3;bin_index++)
         vec[vec_index++]=bins[bin_index];
    }
    vd_RadonTextures(vec);
    delete [] pixels;
    delete [] ptr;
}

//-----------------------------------------------------------------------------------
/* Otsu
   Find otsu threshold
*/
double ImageMatrix::Otsu()
{  long a; //,x,y;
   double hist[256],omega[256],mu[256],sigma_b2[256],maxval=-INF,sum,count;
   double max=pow(2,bits)-1;
   histogram(hist,256,1);
   omega[0]=hist[0]/(width*height);
   mu[0]=1*hist[0]/(width*height);
   for (a=1;a<256;a++)
   {  omega[a]=omega[a-1]+hist[a]/(width*height);
      mu[a]=mu[a-1]+(a+1)*hist[a]/(width*height);
   }
   for (a=0;a<256;a++)
   {  if (omega[a]==0 || 1-omega[a]==0) sigma_b2[a]=0;
      else sigma_b2[a]=pow(mu[255]*omega[a]-mu[a],2)/(omega[a]*(1-omega[a]));
      if (sigma_b2[a]>maxval) maxval=sigma_b2[a];
   }
   sum=0.0;
   count=0.0;
   for (a=0;a<256;a++)
     if (sigma_b2[a]==maxval)
     {  sum+=a;
        count++;
     }	 
   return((pow(2,bits)/256.0)*((sum/count)/max));
}

//-----------------------------------------------------------------------------------
/*
  OtsuBinaryMaskTransform
  Transforms an image to a binary image such that the threshold is otsu global threshold
*/
double ImageMatrix::OtsuBinaryMaskTransform()
{  double OtsuGlobalThreshold;
   double max=pow(2,bits)-1;

   OtsuGlobalThreshold=Otsu();

   /* classify the pixels by the threshold */
   for (long a=0;a<width*height*depth;a++)
     if (data[a].intensity>OtsuGlobalThreshold*max) data[a].intensity=1;
     else data[a].intensity=0;
	 
   return(OtsuGlobalThreshold);
}

/*  BWlabel
    label groups of connected pixel (4 or 8 connected dependes on the value of the parameter "level").
    This is an implementation of the Matlab function bwlabel
    returned value -int- the number of objects found
*/
//--------------------------------------------------------
int ImageMatrix::BWlabel(int level)
{
   return(bwlabel(this,level));
}

//--------------------------------------------------------

void ImageMatrix::centroid(double *x_centroid, double *y_centroid, double *z_centroid)
{
   GlobalCentroid(this,x_centroid,y_centroid,z_centroid);
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

int compare_ints (const void *a, const void *b)
{
  if (*((int *)a) > *((int *)b)) return(1);
  if (*((int *)a) == *((int *)b)) return(0);
  return(-1);
}

void ImageMatrix::FeatureStatistics(int *count, int *Euler, double *centroid_x, double *centroid_y, double *centroid_z, int *AreaMin, int *AreaMax,
                                    double *AreaMean, int *AreaMedian, double *AreaVar, int *area_histogram,double *DistMin, double *DistMax,
                                    double *DistMean, double *DistMedian, double *DistVar, int *dist_histogram, int num_bins)
{  int object_index,inv_count;
   double sum_areas,sum_dists;
   ImageMatrix *BWImage,*BWInvert,*temp;
   int *object_areas;
   double *centroid_dists,sum_dist;

   BWInvert=duplicate();   // check if the background is brighter or dimmer
   BWInvert->invert();
   BWInvert->OtsuBinaryMaskTransform();
   inv_count=BWInvert->BWlabel(8);
	 
   
   BWImage=duplicate();
   BWImage->OtsuBinaryMaskTransform();
   BWImage->centroid(centroid_x,centroid_y,centroid_z);
   *count=BWImage->BWlabel(8);
   if (inv_count>*count)
   {  temp=BWImage;
      BWImage=BWInvert;
      BWInvert=temp;
      *count=inv_count;
      BWImage->centroid(centroid_x,centroid_y,centroid_z);	  
   }
   delete BWInvert;
   *Euler=EulerNumber(BWImage,*count)+1;

   // calculate the areas 
   sum_areas=0;
   sum_dists=0;
   object_areas=new int[*count];
   centroid_dists=new double[*count];
   for (object_index=1;object_index<=*count;object_index++)
   {  double x_centroid,y_centroid,z_centroid;
      object_areas[object_index-1]=FeatureCentroid(BWImage, object_index, &x_centroid, &y_centroid,&z_centroid);
      centroid_dists[object_index-1]=sqrt(pow(x_centroid-(*centroid_x),2)+pow(y_centroid-(*centroid_y),2));
      sum_areas+=object_areas[object_index-1];
      sum_dists+=centroid_dists[object_index-1];
   }
   /* compute area statistics */
   qsort(object_areas,*count,sizeof(int),compare_ints);
   *AreaMin=object_areas[0];
   *AreaMax=object_areas[*count-1];
   if (*count>0) *AreaMean=sum_areas/(*count);
   else *AreaMean=0;
   *AreaMedian=object_areas[(*count)/2];
   for (object_index=0;object_index<num_bins;object_index++)
     area_histogram[object_index]=0;
   /* compute the variance and the histogram */
   sum_areas=0;
   if (*AreaMax-*AreaMin>0)
     for (object_index=1;object_index<=*count;object_index++)
     {  sum_areas+=pow(object_areas[object_index-1]-*AreaMean,2);
        if (object_areas[object_index-1]==*AreaMax) area_histogram[num_bins-1]+=1;
        else area_histogram[((object_areas[object_index-1]-*AreaMin)/(*AreaMax-*AreaMin))*num_bins]+=1;
     }
   if (*count>1) *AreaVar=sum_areas/((*count)-1);
   else *AreaVar=sum_areas;

   /* compute distance statistics */
   qsort(centroid_dists,*count,sizeof(double),compare_doubles);
   *DistMin=centroid_dists[0];
   *DistMax=centroid_dists[*count-1];
   if (*count>0) *DistMean=sum_dists/(*count);
   else *DistMean=0;
   *DistMedian=centroid_dists[(*count)/2];
   for (object_index=0;object_index<num_bins;object_index++)
     dist_histogram[object_index]=0;

   /* compute the variance and the histogram */
   sum_dist=0;
   for (object_index=1;object_index<=*count;object_index++)
   {  sum_dist+=pow(centroid_dists[object_index-1]-*DistMean,2);
      if (centroid_dists[object_index-1]==*DistMax) dist_histogram[num_bins-1]+=1;
      else dist_histogram[(int)(((centroid_dists[object_index-1]-*DistMin)/(*DistMax-*DistMin))*num_bins)]+=1;
   }
   if (*count>1) *DistVar=sum_dist/((*count)-1);
   else *DistVar=sum_dist;

   delete BWImage;
   delete [] object_areas;
   delete [] centroid_dists;
}

/* GaborFilters */
/* ratios -array of double- a pre-allocated array of double[7]
*/
void ImageMatrix::GaborFilters2D(double *ratios)
{  GaborTextureFilters2D(this, ratios);
}


/* haarlick
   output -array of double- a pre-allocated array of 28 doubles
*/
void ImageMatrix::HaarlickTexture2D(double distance, double *out)
{  if (distance<=0) distance=1;
   haarlick2D(this,distance,out);
}

/* MultiScaleHistogram
   histograms into 3,5,7,9 bins
   Function computes signatures based on "multiscale histograms" idea.
   Idea of multiscale histogram came from the belief of a unique representativity of an
   image through infinite series of histograms with sequentially increasing number of bins.
   Here we used 4 histograms with number of bins being 3,5,7,9.
   out -array of double- a pre-allocated array of 24 bins
*/
void ImageMatrix::MultiScaleHistogram(double *out)
{  int a;
   double max=0;
   histogram(out,3,0);
   histogram(&(out[3]),5,0);
   histogram(&(out[8]),7,0);
   histogram(&(out[15]),9,0);
   for (a=0;a<24;a++)
     if (out[a]>max) max=out[a];
   for (a=0;a<24;a++)
     out[a]=out[a]/max;
}

/* TamuraTexture
   Tamura texture signatures: coarseness, directionality, contrast
   vec -array of double- a pre-allocated array of 6 doubles
*/
void ImageMatrix::TamuraTexture2D(double *vec)
{
  Tamura3Sigs2D(this,vec);
}

/* zernike
   zvalue -array of double- a pre-allocated array of double of a suficient size
                            (the actual size is returned by "output_size))
   output_size -* long- the number of enteries in the array "zvalues" (normally 72)
*/
void ImageMatrix::zernike2D(double *zvalues, long *output_size)
{  mb_zernike2D(this, 0, 0, zvalues, output_size);
}

#pragma package(smart_init)





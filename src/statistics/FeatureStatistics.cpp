/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                               */
/*    Copyright (C) 2007 Open Microscopy Environment                             */
/*         Massachusetts Institue of Technology,                                 */
/*         National Institutes of Health,                                        */
/*         University of Dundee                                                  */
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


#include <string.h>

#include "FeatureStatistics.h"

typedef struct POINT1 {
	long x,y;
} point;

//---------------------------------------------------------------------------
/*  BWlabel
    label groups of 4-connected pixels.
    This is an implementation of the Matlab function bwlabel
*/
unsigned long bwlabel(ImageMatrix *Im, int level) {
	long x, y, base_x,base_y,stack_count, w = Im->width, h = Im->height;
	unsigned long group_counter = 1;
	point *stack=new point[w*h];
	pixData &pix_plane = Im->WriteablePixels();

	for (y = 0; y < h; y++) {
		for (x = 0; x < w; x++) {
			if ( pix_plane(y,x) == 1) {
				/* start a new group */
				group_counter++;
				pix_plane(y,x) = group_counter;
				stack[0].x=x;
				stack[0].y=y;
				stack_count=1;
				while (stack_count > 0) {
					base_x=stack[0].x;
					base_y=stack[0].y;

					if (base_x > 0 && pix_plane(base_y,base_x-1) == 1) {
						pix_plane(base_y,base_x-1) = group_counter;
						stack[stack_count].x=base_x-1;
						stack[stack_count].y=base_y;
						stack_count++;
					}

					if (base_x < w-1 && pix_plane(base_y,base_x+1) == 1) {
						pix_plane(base_y,base_x+1) = group_counter;
						stack[stack_count].x=base_x+1;
						stack[stack_count].y=base_y;
						stack_count++;
					}

					if (base_y > 0 && pix_plane(base_y-1,base_x) == 1) {
						pix_plane(base_y-1,base_x) = group_counter;
						stack[stack_count].x=base_x;
						stack[stack_count].y=base_y-1;
						stack_count++;
					}

					if (base_y < h-1 && pix_plane(base_y+1,base_x) == 1) {
						pix_plane(base_y+1,base_x) = group_counter;
						stack[stack_count].x=base_x;
						stack[stack_count].y=base_y+1;
						stack_count++;
					}

					/* look for 8 connected pixels */
					if (level==8) {
						if (base_x > 0 && base_y > 0 && pix_plane(base_y-1,base_x-1) == 1) {
							pix_plane(base_y-1,base_x-1) = group_counter;
							stack[stack_count].x=base_x-1;
							stack[stack_count].y=base_y-1;
							stack_count++;
						}

						if (base_x < w-1 && base_y > 0 && pix_plane(base_y-1,base_x+1) == 1) {
							pix_plane(base_y-1,base_x+1) = group_counter;
							stack[stack_count].x=base_x+1;
							stack[stack_count].y=base_y-1;
							stack_count++;
						}

						if (base_x > 0 && base_y < h-1 && pix_plane(base_y+1,base_x-1) == 1) {
							pix_plane(base_y+1,base_x-1) = group_counter;
							stack[stack_count].x=base_x-1;
							stack[stack_count].y=base_y+1;
							stack_count++;
						}

						if (base_x < w-1 && base_y < h-1 && pix_plane(base_y+1,base_x+1) == 1) {
							pix_plane(base_y+1,base_x+1) = group_counter;
							stack[stack_count].x=base_x+1;
							stack[stack_count].y=base_y+1;
							stack_count++;
						}
					}
			  
					stack_count-=1;
					memmove(stack,&(stack[1]),sizeof(point)*stack_count);
				}
			}
		}
	}

	/* now decrease every non-zero pixel by one because the first group was "2" */
	for (y=0;y<h;y++)
		for (x=0;x<w;x++)
			if (pix_plane(y,x) != 0)
				pix_plane(y,x) -= 1;

	delete [] stack;
	Im->WriteablePixelsFinish();
	return(group_counter-1);
}

/* the input should be a binary image */
void GlobalCentroid(ImageMatrix *Im, double *x_centroid, double *y_centroid) {
	unsigned int x,y,w = Im->width, h = Im->height;
	double x_mass=0,y_mass=0,mass=0;
	readOnlyPixels pix_plane = Im->ReadablePixels();

	for (y = 0; y < h; y++)
		for (x = 0; x < w; x++)
			if (pix_plane(y,x) > 0) {
				x_mass=x_mass+x+1;    /* the "+1" is only for compatability with matlab code (where index starts from 1) */
				y_mass=y_mass+y+1;    /* the "+1" is only for compatability with matlab code (where index starts from 1) */
				mass++;
			}
	if (mass) {
		*x_centroid=x_mass/mass;
		*y_centroid=y_mass/mass;
	} else *x_centroid=*y_centroid=0;
}

/* find the centroid of a certain feature
   the input image should be a bwlabel transform of a binary image
   the retruned value is the area of the feature
*/
unsigned long FeatureCentroid(ImageMatrix *Im, double object_index,double *x_centroid, double *y_centroid) {
	unsigned int x,y,w = Im->width, h = Im->height;
	unsigned long x_mass=0,y_mass=0,mass=0;
	readOnlyPixels pix_plane = Im->ReadablePixels();

	for (y = 0; y < h; y++)
		for (x = 0; x < w; x++)
			if (pix_plane(y,x) == object_index) {
				x_mass=x_mass+x+1;      /* the "+1" is only for compatability with matlab code (where index starts from 1) */
				y_mass=y_mass+y+1;      /* the "+1" is only for compatability with matlab code (where index starts from 1) */
				mass++;
			}
	if (x_centroid) *x_centroid=(double)x_mass/(double)mass;
	if (y_centroid) *y_centroid=(double)y_mass/(double)mass;
	return(mass);
}

/* the number of pixels that are above the threshold
   the input image is a binary image
*/
unsigned long area(ImageMatrix *Im) {
	unsigned int x,y,w = Im->width, h = Im->height;
	unsigned long sum=0;
	readOnlyPixels pix_plane = Im->ReadablePixels();
	for (y = 0; y < h; y++)
		for (x = 0; x < w; x++)
			sum += (pix_plane(y,x) > 0 ? 1 : 0);
	return(sum);
}

/* EulerNumber
   The input image should be a binary image
*/
long EulerNumber(ImageMatrix *Im, unsigned long FeatureNumber) {
	unsigned int x,y;
	unsigned long HolesNumber;
	ImageMatrix cp;
	cp.copy (*Im);
	pixData &cp_pix_plane = cp.WriteablePixels();

	/* inverse the image */
	for (y = 0; y < cp.height; y++)
		for (x = 0; x < cp.width; x++)
			if (cp_pix_plane(y,x) > 0.0)
				cp_pix_plane(y,x) = 0.0;
			else cp_pix_plane(y,x) = 1.0;

	HolesNumber=cp.BWlabel(8);
	return(FeatureNumber-HolesNumber-1);
}

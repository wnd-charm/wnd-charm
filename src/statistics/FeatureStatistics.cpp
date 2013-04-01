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
unsigned long bwlabel(ImageMatrix &Im, int level) {
	long x, y, base_x,base_y,stack_count, w = Im.width, h = Im.height;
	unsigned long group_counter = 1;
	point *stack=new point[w*h];
	pixData &pix_plane = Im.WriteablePixels();

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
	return(group_counter-1);
}

/* the input should be a binary image */
void GlobalCentroid(const ImageMatrix &Im, double *x_centroid, double *y_centroid) {
	unsigned int x,y,w = Im.width, h = Im.height;
	double x_mass=0,y_mass=0,mass=0;
	readOnlyPixels pix_plane = Im.ReadablePixels();

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
unsigned long FeatureCentroid(const ImageMatrix &Im, double object_index,double *x_centroid, double *y_centroid) {
	unsigned int x,y,w = Im.width, h = Im.height;
	unsigned long x_mass=0,y_mass=0,mass=0;
	readOnlyPixels pix_plane = Im.ReadablePixels();

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
unsigned long area(const ImageMatrix &Im) {
	unsigned int x,y,w = Im.width, h = Im.height;
	unsigned long sum=0;
	readOnlyPixels pix_plane = Im.ReadablePixels();
	for (y = 0; y < h; y++)
		for (x = 0; x < w; x++)
			sum += (pix_plane(y,x) > 0 ? 1 : 0);
	return(sum);
}

/* EulerNumber
   The input image should be a binary image
*/
long EulerNumber(const ImageMatrix &Im, int mode) {  
	unsigned long x, y;
	size_t i;
	// quad-pixel match patterns
	unsigned char Px[] = {
		// P1 - single pixel
		(1 << 3) | (0 << 2) |
		(0 << 1) | (0 << 0),
		(0 << 3) | (1 << 2) |
		(0 << 1) | (0 << 0),
		(0 << 3) | (0 << 2) |
		(1 << 1) | (0 << 0),
		(0 << 3) | (0 << 2) |
		(0 << 1) | (1 << 0),
		// P3 - 3-pixel
		(0 << 3) | (1 << 2) |
		(1 << 1) | (1 << 0),
		(1 << 3) | (0 << 2) |
		(1 << 1) | (1 << 0),
		(1 << 3) | (1 << 2) |
		(0 << 1) | (1 << 0),
		(1 << 3) | (1 << 2) |
		(1 << 1) | (0 << 0),
		// Pd - diagonals
		(1 << 3) | (0 << 2) |
		(0 << 1) | (1 << 0),
		(0 << 3) | (1 << 2) |
		(1 << 1) | (0 << 0)
	};
	unsigned char Imq;
	// Pattern match counters
	long C1 = 0, C3 = 0, Cd = 0;
	readOnlyPixels pix_plane = Im.ReadablePixels();
	
	assert ( (mode == 4 || mode == 8) && "Calling EulerNumber with mode other than 4 or 8");

	// update pattern counters by scanning the image.
	for (y = 1; y < Im.height; y++) {
		for (x = 1; x < Im.width; x++) {
			// Get the quad-pixel at this image location
			Imq = 0;
			if (pix_plane(y-1,x-1) > 0) Imq |=  (1 << 3);
			if (pix_plane(y-1,x  ) > 0) Imq |=  (1 << 2);
			if (pix_plane(y  ,x-1) > 0) Imq |=  (1 << 1);
			if (pix_plane(y  ,x  ) > 0) Imq |=  (1 << 0);
			// find the matching pattern
			for (i = 0; i < 10; i++) if (Imq == Px[i]) break;
			if      (i >= 0 && i <= 3) C1++;
			else if (i >= 4 && i <= 7) C3++;
			else if (i == 8 && i == 9) Cd++;
		}
	}
	
	if (mode == 4)
		return ( (C1 - C3 + (2*Cd)) / 4);
	else
		return ( (C1 - C3 - (2*Cd)) / 4);
}

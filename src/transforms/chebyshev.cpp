/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*                                                                               */
/*    Copyright (C) 2003 Open Microscopy Environment                             */
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


#pragma hdrstop

#include <stdlib.h>
#include <math.h>

#include "../cmatrix.h"
#include "chebyshev.h"

//---------------------------------------------------------------------------

void TNx(double *x, double *out, int N, int height)
{  int ix,iy;
   double *temp,*temp1;
//	if( max(abs(x(:))) > 1 )
//		error(':: Cheb. Polynomials Tn :: abs(arg)>1');
//	end;
   temp=new double[N*height];
   temp1=new double[N*height];

// 	T = cos((ones(size(x,2),1)*(0:(N-1))).*acos(x'*ones(1,N)));
//     	T(:,1)=ones(size(x'));

   // x'*ones(1,N)
   for (ix=0;ix<N;ix++)
     for (iy=0;iy<height;iy++)
       temp[iy*N+ix]=x[iy];
   // acos
   for (ix=0;ix<N;ix++)
     for (iy=0;iy<height;iy++)
       if (fabs(temp[iy*N+ix])>1) temp[iy*N+ix]=0;   /* protect from acos domain error */
       else temp[iy*N+ix]=acos(temp[iy*N+ix]);
   // ones(size(x,2),1)*(0:(N-1))
   for (ix=0;ix<N;ix++)
     for (iy=0;iy<height;iy++)
       temp1[iy*N+ix]=ix;
   //.*
   for (ix=0;ix<N;ix++)
     for (iy=0;iy<height;iy++)
       out[iy*N+ix]=temp[iy*N+ix]*temp1[iy*N+ix];
   //cos
   for (ix=0;ix<N;ix++)
     for (iy=0;iy<height;iy++)
       out[iy*N+ix]=cos(out[iy*N+ix]);

   for (iy=0;iy<height;iy++)
       out[iy*N+0]=1;

   delete [] temp1;
   delete [] temp;
}

void getChCoeff1D(double *f,double *out,double *x,int N,int width)
{  double *Tj,*tj;
   int jj,a;

   Tj=new double[width*N];
   tj=new double[width];
   TNx(x,Tj,N,width);
   for (jj=0;jj<N;jj++)
   {   int jx;
       jx=jj;
       for (a=0;a<width;a++)
         tj[a]=Tj[a*N+jj];
       if (!jx)
       {
         for (a=0;a<width;a++)
           tj[a]=tj[a]/(double)width;
       }
       else
       {
         for (a=0;a<width;a++)
           tj[a]=tj[a]*2/(double)width;
       }
       out[jj]=0;
       for (a=0;a<width;a++)
         out[jj]+=f[a]*tj[a]/2;
   }
   delete [] tj;
   delete [] Tj;
}

void getChCoeff(double *Im, double *out, double *x,int N,int width, int height)
{  int iy,ix;
   double *y,*y_out;
   y=new double[width];
   y_out=new double[N];

   for (iy=0;iy<height;iy++)
   {
     for (ix=0;ix<width;ix++)
       y[ix]=Im[iy*width+ix];
     getChCoeff1D(y,y_out,x,N,width);
     for (ix=0;ix<N;ix++)
       out[iy*N+ix]=y_out[ix];
   }

   delete [] y;
   delete [] y_out;
}

/* inputs:
IM - image
N - coefficient
width - width of the image
height - height of the image
*/
void Chebyshev2D(ImageMatrix *Im, double *out, int N)
{  double *x,*y;
   double *in;
   int a,i,j;

// Make a default value for coeficient order if it was not given as an input
//   if (N<=0)
//     N=min(Im->width,Im->height);

   x=new double[Im->width];
   y=new double[Im->height];

   for (a=0;a<Im->width;a++)
     x[a]=2*(double)(a+1) / (double)Im->width -1;
   for (a=0;a<Im->height;a++)
     y[a]=2*(double)(a+1) / (double)Im->height -1;

   in=new double[Im->width*Im->height];
   for (j=0;j<Im->height;j++)
     for (i=0;i<Im->width;i++)
       in[j*Im->width+i]=(double)Im->pixel(i,j,0).intensity;
   getChCoeff(in,out,x,N,Im->width,Im->height);
   /* transpose the matrix "out" into "in" */
   for (j=0;j<N;j++)
     for (i=0;i<Im->height/*Im->width*/;i++)
       in[j*Im->height+i]=out[i*N+j];

   getChCoeff(in,out,y,N,Im->height,N);

   delete [] in;
   delete [] x;
   delete [] y;
}

#pragma package(smart_init)




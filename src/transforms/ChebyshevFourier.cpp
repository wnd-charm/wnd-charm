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

//#include <vcl.h>

#include <math.h>
#ifndef WIN32
#include <stdio.h>
#endif
#include "ChebyshevFourier.h"

#define min(a, b)  (((a) < (b)) ? (a) : (b))
//---------------------------------------------------------------------------

double *ChebPol(double x,int N, double *out)
{  int a;
   if (x<-1) x=-1;      /* protect the acos function */
   if (x>1) x=1;
   for (a=0;a<N;a++)
     out[a]=a;
   for (a=0;a<N;a++)
     out[a]=cos(out[a]*acos(x));
   out[0]=0.5;
   return(out);
}



/*
ChebyshevFourier - Chebyshev Fourier transform
"coeff_packed" -array of doubles- a pre-allocated array of 32 doubles
*/
void ChebyshevFourier2D(ImageMatrix *Im, long N, double *coeff_packed, int packingOrder)
{  long a,b,m,n,x,y,nLast,NN,Nmax,ind;
   double *xx,*yy,*f,*r,*tmp,*img,*Tn,*c,*coeff;
   float **C;
   double min,max;
   long *kk;

   if (N==0) N=11;
   m=Im->height;
   n=Im->width;
   if (m*n>120*160)
   {  // image too large ==> downsample
   }

   yy=new double[m*n];
   xx=new double[m*n];
   img=new double[m*n];
   f=new double[m*n];
   r=new double[m*n];
   kk=new long[m*n];  /* the required size of kk is equal to nLast */

//   for (a=0;a<m;a++)
//     y[a]=-1+2/m);
//   for (a=0;a<n;a++)
//     x[a]=-1+2/n);
   for (x=0;x<n;x++)
     for (y=0;y<m;y++)
       xx[x*m+y]=-1+(double)x*(2/((double)n-1));

   for (x=0;x<n;x++)
     for (y=0;y<m;y++)
       yy[x*m+y]=-1+(double)y*(2/((double)m-1));

   for (y=0;y<m;y++)
     for (x=0;x<n;x++)
       img[x*m+y]=Im->pixel(x,y,0).intensity;
   /* convert cartesian to polar */
   nLast=0;
   for (ind=0;ind<m*n;ind++)
   {  r[ind]=sqrt(pow(xx[ind],2)+pow(yy[ind],2));
      if (yy[ind]==0) f[ind]=-3.14159265*(xx[ind]<0);     /* avoid atan2 domain error */
      else f[ind]=-1*atan2(yy[ind],xx[ind]);
      if (r[ind]<1)
         kk[nLast++]=ind;
   }  
   delete [] xx;
   delete [] yy;
   Nmax=(int)((min(m,n)-1)/2);
   if (N>Nmax) N=Nmax;
   NN = 2*N + 1;
   C=new float*[nLast];
   for (a=0;a<nLast;a++)
   {  C[a]=new float[NN*NN*2];                    /* *2 because of the complex numbers */
      for (b=0;b<NN*NN*2;b++)                   /* *2 because of the complex numbers */
        C[a][b]=0.0;   /* initialize */
   }

   c=new double[NN*NN*2];                         /* *2 because of the complex numbers */
   Tn=new double[NN];
   tmp=new double[NN*2];  /* *2 for the complex numbers */
   for (ind=0;ind<nLast;ind++)
   {  double ri,fi;
      int im,c_ind;
      c_ind=0;
      ri=r[kk[ind]];
      fi=f[kk[ind]];
      if (ri>1) continue;
      ChebPol(ri*2-1,NN,Tn);
      for (im=1;im<=NN;im++)
      {  double Ftrm;
         int mf;
         mf=im-1-N;
         if (!mf) Ftrm=0.5;
         else Ftrm=1.0;
         for (a=0;a<NN;a++)
         {  //tmp[a]=Ftrm*exp(/*-i**/mf*fi)*Tn[a];
           tmp[a*2]=Ftrm*cos(mf*fi)*Tn[a];                 /* *2 because of the complex numbers */
           tmp[a*2+1]=Ftrm*sin(-1*mf*fi)*Tn[a];
         }
         for (a=0;a<NN*2;a++)                                 /* *2 because of the complex numbers */
           c[c_ind+a]=tmp[a];
         c_ind=c_ind+NN*2;                                    /* *2 because of the complex numbers */
      }
      for (a=0;a<NN*NN*2;a++)                                 /* *2 because of the complex numbers */
       C[ind][a]=c[a];  /* can optimize this by moving it 4 lines upwards */
   }
   delete []tmp;
   delete [] c;

   min=INF;
   max=-INF;
   coeff = new double[NN*NN];
   for (a=0;a<NN*NN;a++)                                 /* *2 because of the complex numbers */
   { double sumr=0,sumi=0;
     for (b=0;b<nLast;b++)
     {  sumr+=C[b][a*2]*img[kk[b]];     //C[b][a]*img[kk[b]];
        sumi+=C[b][a*2+1]*img[kk[b]];                  /* *2 because of the complex numbers */
     }
     coeff[a]=sqrt(sumr*sumr+sumi*sumi);
     if (coeff[a]<min) min=coeff[a];
     if (coeff[a]>max) max=coeff[a];
   }

   for (a=0;a<packingOrder;a++)
     coeff_packed[a]=0;
   if (max!=min)
     for (a=0;a<NN*NN;a++)
     {  if (coeff[a]==max) coeff_packed[31]+=1;
        else coeff_packed[(int)(((coeff[a]-min)/(max-min))*32)]+=1;
     }

   delete [] coeff;
   for (a=0;a<nLast;a++)
     delete [] C[a];
   delete [] C;
   delete [] f;
   delete [] r;
   delete [] kk;
   delete [] img;
   delete [] Tn;
}



#pragma package(smart_init)

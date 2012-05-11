/****************************************************************************/
/*                                                                          */
/*                                                                          */
/*                              mb_Znl.c                                    */
/*                                                                          */
/*                                                                          */
/*                           Michael Boland                                 */
/*                            09 Dec 1998                                   */
/*                                                                          */     
/*  Revisions:                                                              */
/*  9-1-04 Tom Macura <tmacura@nih.gov> modified to make the code ANSI C    */
/*         and work with included complex arithmetic library from           */
/*         Numerical Recepies in C instead of using the system's C++ STL    */
/*         Libraries.                                                       */
/*                                                                          */
/*  1-29-06 Lior Shamir <shamirl (-at-) mail.nih.gov> modified "factorial"  */
/*  to a loop, replaced input structure with ImageMatrix class.             */
/****************************************************************************/


//---------------------------------------------------------------------------

#pragma hdrstop

#include <math.h>

#include "../../cmatrix.h"
#include "complex.h"

#include "zernike.h"

//---------------------------------------------------------------------------

double factorial(double n)
{  double ans=1;
   int a;
   if (n<0) return(0);

   for (a=1;a<=n;a++)
     ans=ans*a;
   return(ans);
}

/* mb_imgmoments
   calculates the moment MXY for IMAGE

*/
double mb_imgmoments(ImageMatrix *image, int x, int y)
{  double *xcoords,sum;
   xcoords=new double[image->width*image->height];
   int row,col;
   /* Generate a matrix with the x coordinates of each pixel. */
   for (row=0;row<image->height;row++)
     for (col=0;col<image->width;col++)
       xcoords[row*image->width+col]=pow(col+1,x);

   sum=0;
   /* Generate a matrix with the y coordinates of each pixel. */
   for (col=0;col<image->width;col++)
     for (row=0;row<image->height;row++)
     {
        if (y!=0)
        {  if (x==0) xcoords[row*image->width+col]=pow(row+1,y);
           else
           xcoords[row*image->width+col]=pow(col+1,y)*xcoords[row*image->width+col];
        }
        sum+=xcoords[row*image->width+col]*image->pixel(col,row,0).intensity;
     }

   delete [] xcoords;
   return(sum);
}


/* mb_Znl
  Zernike moment generating function.  The moment of degree n and
  angular dependence l for the pixels defined by coordinate vectors
  X and Y and intensity vector P.  X, Y, and P must have the same
  length
*/
void mb_Znl(long n, long l, double *X, double *Y, double *P, int size, double *out_r, double *out_i)
{ double x, y, p ;   /* individual values of X, Y, P */
  int i,m;
  fcomplex sum;              /* Accumulator for complex moments */
  fcomplex Vnl;              /* Inner sum in Zernike calculations */
  double* preal;             /* Real part of return value */
  double* pimag;             /* Imag part of return value */

  sum = Complex (0.0, 0.0);

  for(i = 0 ; i < size ; i++) {
    x = X[i] ;
    y = Y[i] ;
    p = P[i] ;

    Vnl = Complex (0.0, 0.0);
    for( m = 0; m <= (n-l)/2; m++) {
      double tmp = (pow((double)-1.0,(double)m)) * ( factorial(n-m) ) /
				( factorial(m) * (factorial((n - 2.0*m + l) / 2.0)) *
	  			(factorial((n - 2.0*m - l) / 2.0)) ) *
				( pow( sqrt(x*x + y*y), (double)(n - 2*m)) );

	  Vnl = Cadd (Vnl, RCmul(tmp, Rpolar(1.0, l*atan2(y,x))) );
      /*
       NOTE: This function did not work with the following:
        ...pow((x*x + y*y), (double)(n/2 -m))...
        perhaps pow does not work properly with a non-integer
        second argument.
       'not work' means that the output did not match the 'old'
        Zernike calculation routines.
      */
    }

    /* sum += p * conj(Vnl) ; */
	sum = Cadd(sum, RCmul(p, Conjg(Vnl)));
  }

  /* sum *= (n+1)/3.14159265 ; */
  sum = RCmul((n+1)/3.14159265, sum);


  /* Assign the returned value */
  *out_r = sum.r ;
  *out_i = sum.i ;
}


void mb_zernike2D(ImageMatrix *I, double D, double R, double *zvalues, long *output_size)
{  int rows,cols,n,l;
   double *Y,*X,*P,psum;
   double moment10,moment00,moment01;
   int x,y,size,size2,size3,a;

   if (D<=0) D=15;
   if (R<=0)
   {  rows=I->height;
      cols=I->width;
      R=rows/2;
   }
   Y=new double[rows*cols];
   X=new double[rows*cols];
   P=new double[rows*cols];

   /* Find all non-zero pixel coordinates and values */
   size=0;
   psum=0;
   for (y=0;y<rows;y++)
     for (x=0;x<cols;x++)
     if (I->pixel(x,y,0).intensity!=0)
     {  Y[size]=y+1;
        X[size]=x+1;
        P[size]=I->pixel(x,y,0).intensity;
        psum+=I->pixel(x,y,0).intensity;
        size++;
     }

   /* Normalize the coordinates to the center of mass and normalize
      pixel distances using the maximum radius argument (R) */
   moment10=mb_imgmoments(I,1,0);
   moment00=mb_imgmoments(I,0,0);
   moment01=mb_imgmoments(I,0,1);
   size2=0;
   for (a=0;a<size;a++)
   {  X[size2]=(X[a]-moment10/moment00)/R;
      Y[size2]=(Y[a]-moment01/moment00)/R;
      P[size2]=P[a]/psum;
      if (sqrt(X[size2]*X[size2]+Y[size2]*Y[size2])<=1) size2=size2++;
   }

   size3=0;
   for (n=0;n<=D;n++)
     for (l=0;l<=n;l++)
       if (((n-l) % 2) ==0)
       {  double preal,pimag;
          mb_Znl(n,l,X,Y,P,size2,&preal,&pimag);
          zvalues[size3++]=fabs(sqrt(preal*preal+pimag*pimag));
       }
   *output_size=size3;

   delete [] Y;
   delete [] X;
   delete [] P;

}



#pragma package(smart_init)



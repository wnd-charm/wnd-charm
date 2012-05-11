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

#include <math.h>
#include "CombFirst4Moments.h"


double round(double x)
{  double a;
   a=fabs(x-int(x));
   if (a<0.5) return(int(x));
   else return(int(x+1*(x>0)));
}

double skewness(double *vec, double avg, double std, int length)
{  int a;
   double s2,m3;

   if (length==0) return(0);

   s2=0.0;
   m3=0.0;
   for (a=0;a<length;a++)
   {  m3=m3+pow(vec[a]-avg,3);
      s2=s2+pow(vec[a]-avg,2);
   }
   m3=m3/length;
   s2=s2/length;

   if (s2==0) return(0);
   else return(m3/pow(s2,1.5));
}
//---------------------------------------------------------------------------

double kurtosis(double *vec, double avg, double std, int length)
{  int a;
   double s2,m4;
   m4=0.0;
   s2=0.0;

  if (length==0) return(0);

  for (a=0;a<length;a++)
  {  m4=m4+pow(vec[a]-avg,4);
     s2=s2+pow(vec[a]-avg,2);
  }
  m4=m4/length;
  s2=s2/length;

  if (s2==0) return(0);
  else return(m4/(s2*s2));

/*
   for (int a = 0; a < length; a++)
       sum += pow((vec[a] - avg) / std, 4.0);

   double coefficientOne = (length * (length + 1)) / ((length - 1) * (length - 2) * (length - 3));
   double termTwo = ((3 * pow(length - 1, 2.0)) / ((length - 2) * (length - 3)));
   // Calculate kurtosis
   kurtosis = (coefficientOne * accum) - termTwo;

   return kurtosis;
*/
}

//---------------------------------------------------------------------------

void get4scalMoments(double *vec, int vec_length, double *z)
{  int a;
   double sum,avg,s1;
   /* compute the std dev */
   sum=0.0;
   for (a=0;a<vec_length;a++)
     sum+=vec[a];
   avg=sum/(double)vec_length;
   sum=0.0;
   if (vec_length<=1) s1=0;
   else
   {  for (a=0;a<vec_length;a++)
       sum+=(vec[a]-avg)*(vec[a]-avg);
      s1=sqrt(sum/(double)(vec_length-1));
   }

   if (s1==0)
   {  z[0]=avg;z[1]=s1;z[2]=0;z[3]=0;
   }
   else
   {  z[0]=avg;z[1]=s1;z[2]=skewness(vec,avg,s1,vec_length);z[3]=kurtosis(vec,avg,s1,vec_length);
   }
}

//---------------------------------------------------------------------------

int matr4moments_to_hist(double matr4moments[4][21], double *vec, int vec_start)
{  int a,b,vec_index,bin_index;
   int nbins = 3;
   double bins[3];
   vec_index=vec_start;
   for (a=0;a<4;a++)
   {  double min=10E20,max=-10E20;
      for (b=0;b<nbins;b++) bins[b]=0;
      /* find min and max (for the bins) */
      for (b=0;b<21;b++)
      {  if (matr4moments[a][b]>max) max=matr4moments[a][b];
         if (matr4moments[a][b]<min) min=matr4moments[a][b];
      }
      /* find the bins */
      for (b=0;b<21;b++)
      {  if (matr4moments[a][b]==max) bin_index=nbins-1;
		 else bin_index=(int)(nbins*(matr4moments[a][b]-min)/(max-min));
		 if (bin_index>nbins-1) bin_index=nbins-1;  /* make sure to prevent an error */
		 bins[bin_index]++;
	  }	   
      /* add the bins to the vector */
      for (bin_index=0;bin_index<nbins;bin_index++)
        vec[vec_index++]=bins[bin_index];
   }
   return(vec_index);
}

//---------------------------------------------------------------------------
/*
void func(ImageMatrix *Im, double **I, double **J, double *tmp, int m, int n, double z[], double z4[],double matr4moments[][])
{
   int ii,x,y;
   for (ii=1-m;ii<m;ii=ii+round((double)m/10))
   {  int count=0;
      for (y=0;y<m;y++)
        for (x=0;x<n;x++)
        {  if (fabs(I[x][y]+ii-J[x][y])<1)
             tmp[count++]=Im->data[x][y].intensity;
        }
        if (count==0)
        {  for (a=0;a<4;a++)
             matr4moments[a][ii+m-1]=z4[a];
        }
        else
        {  double z[4];
           get4scalMoments(tmp,m*n,z);
           for (a=0;a<4;a++)
             matr4moments[a][ii+m-1]=z[a];
        }
   }
}
//---------------------------------------------------------------------------
*/

int CombFirst4Moments2D(ImageMatrix *Im, double *vec)
{  double **I,**J,**J1,*tmp,z[4],z4[4]={0,0,0,0};
   double matr4moments[4][21];
   int m,n,n2,m2;
   int a,x,y,ii;
   int matr4moments_index;
   int vec_count=0;

   for (a=0;a<4;a++)    /* initialize */
     for (matr4moments_index=0;matr4moments_index<21;matr4moments_index++)
       matr4moments[a][matr4moments_index]=0;

   m=Im->height;
   n=Im->width;
   I=new double*[n];
   J=new double*[n];
   J1=new double*[n];
   tmp=new double[m*n];
   for (a=0;a<n;a++)
   {  I[a]=new double[m];
      J[a]=new double[m];
      J1[a]=new double[m];
   }

   for (y=0;y<Im->height;y++)
     for (x=0;x<Im->width;x++)
     {  I[x][y]=y+1;
        J[x][y]=x+1;
     }

   n2=(int)(round(n/2));
   m2=(int)(round(m/2));

   /* major diag */
   matr4moments_index=0;
   for (ii=1-m;ii<=m;ii=ii+(int)(round((double)m/10)))
   {  int count=0;
      for (y=0;y<m;y++)
        for (x=0;x<n;x++)
        {  if (fabs(I[x][y]+ii-J[x][y])<1)
             tmp[count++]=Im->pixel(x,y,0).intensity;
        }
        if (count==0)
        {  for (a=0;a<4;a++)
             matr4moments[a][matr4moments_index]=z4[a];
        }
        else
        {  double z[4];
           get4scalMoments(tmp,count,z);
           for (a=0;a<4;a++)
             matr4moments[a][matr4moments_index]=z[a];
        }
      matr4moments_index++;
   }

//   func(Im, I, J,tmp, m, n, z[], z4[],matr4moments[][]);
   vec_count=matr4moments_to_hist(matr4moments,vec,vec_count);

   /* fliplr J */
   for (y=0;y<m;y++)
     for (x=0;x<n;x++)
       J1[x][y]=J[n-1-x][y];

//   matr4moments = [];

   matr4moments_index=0;
   for (ii=1-m;ii<=m;ii=ii+(int)(round((double)m/10)))
   {  int count=0;
      for (y=0;y<m;y++)
        for (x=0;x<n;x++)
        {  if (fabs(I[x][y]+ii-J1[x][y])<1)
             tmp[count++]=Im->pixel(x,y,0).intensity;
        }
        if (count==0)
        {  for (a=0;a<4;a++)
             matr4moments[a][matr4moments_index]=z4[a];
        }
        else
        {  double z[4];
           get4scalMoments(tmp,count,z);
           for (a=0;a<4;a++)
             matr4moments[a][matr4moments_index]=z[a];
        }
      matr4moments_index++;
   }

//   func(Im, I, J,tmp, m, n, z[], z4[],matr4moments[][]);
   vec_count=matr4moments_to_hist(matr4moments,vec,vec_count);

   /* vertical comb */
   matr4moments_index=0;
   for (ii=1-n;ii<=n;ii=ii+(int)(round((double)n/10)))
   {  int count=0;
      for (y=0;y<m;y++)
        for (x=0;x<n;x++)
        {  if (fabs(J[x][y]+ii-n2)<1)
             tmp[count++]=Im->pixel(x,y,0).intensity;
        }
        if (count==0)
        {  for (a=0;a<4;a++)
             matr4moments[a][matr4moments_index]=z4[a];
        }
        else
        {  double z[4];
           get4scalMoments(tmp,count,z);
           for (a=0;a<4;a++)
             matr4moments[a][matr4moments_index]=z[a];
        }
      matr4moments_index++;
   }
   vec_count=matr4moments_to_hist(matr4moments,vec,vec_count);

   /* horizontal comb */
   matr4moments_index=0;
   for (ii=1-m;ii<=m;ii=ii+(int)(round((double)m/10)))
   {  int count=0;
      for (y=0;y<m;y++)
        for (x=0;x<n;x++)
        {  if (fabs(I[x][y]+ii-m2)<1)
             tmp[count++]=Im->pixel(x,y,0).intensity;
        }
        if (count==0)
        {  for (a=0;a<4;a++)
             matr4moments[a][matr4moments_index]=z4[a];
        }
        else
        {  double z[4];
           get4scalMoments(tmp,count,z);
           for (a=0;a<4;a++)
             matr4moments[a][matr4moments_index]=z[a];
        }
      matr4moments_index++;
   }
   vec_count=matr4moments_to_hist(matr4moments,vec,vec_count);

   /* free the memory used by the function */
   for (a=0;a<n;a++)
   {  delete [] I[a];
      delete [] J[a];
      delete [] J1[a];
   }
   delete [] I;
   delete [] J;
   delete [] J1;
	 delete [] tmp;

   return(vec_count);
}



void vd_Comb4Moments(double *in)
{  double temp[48];
   int a;
   for (a=0;a<48;a++)
     temp[a]=in[a];

   in[0]=temp[45];
   in[1]=temp[46];
   in[2]=temp[47];
   in[3]=temp[36];
   in[4]=temp[37];
   in[5]=temp[38];
   in[6]=temp[42];
   in[7]=temp[43];
   in[8]=temp[44];
   in[9]=temp[39];
   in[10]=temp[40];
   in[11]=temp[41];
   in[12]=temp[33];
   in[13]=temp[34];
   in[14]=temp[35];
   in[15]=temp[24];
   in[16]=temp[25];
   in[17]=temp[26];
   in[18]=temp[30];
   in[19]=temp[31];
   in[20]=temp[32];
   in[21]=temp[27];
   in[22]=temp[28];
   in[23]=temp[29];
   in[24]=temp[9];
   in[25]=temp[10];
   in[26]=temp[11];
   in[27]=temp[0];
   in[28]=temp[1];
   in[29]=temp[2];
   in[30]=temp[6];
   in[31]=temp[7];
   in[32]=temp[8];
   in[33]=temp[3];
   in[34]=temp[4];
   in[35]=temp[5];
   in[36]=temp[21];
   in[37]=temp[22];
   in[38]=temp[23];
   in[39]=temp[12];
   in[40]=temp[13];
   in[41]=temp[14];
   in[42]=temp[18];
   in[43]=temp[19];
   in[44]=temp[20];
   in[45]=temp[15];
   in[46]=temp[16];
   in[47]=temp[17];
}

#pragma package(smart_init)

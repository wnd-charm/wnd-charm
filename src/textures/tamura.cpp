//---------------------------------------------------------------------------

#pragma hdrstop

#include <math.h>
#include "../statistics/CombFirst4Moments.h"
#include "tamura.h"



double contrast(ImageMatrix *image)
{ double *vec;
  int x,y;
  double avg,std,k;

  vec=new double[image->width*image->height];
  avg=0;
  for (x=0;x<image->width;x++)
    for (y=0;y<image->height;y++)
    {  vec[x*image->height+y]=image->pixel(x,y,0).intensity;
       avg+=image->pixel(x,y,0).intensity;
    }
  avg=avg/(image->width*image->height);
  std=0;
  for (x=0;x<image->width;x++)
    for (y=0;y<image->height;y++)
      std+=(image->pixel(x,y,0).intensity-avg)*(image->pixel(x,y,0).intensity-avg);
  std=sqrt(std/(image->width*image->height));
  k=kurtosis(vec, avg, std, image->width*image->height);
  delete [] vec;
  if (std<0.0000000001) return(0);
  else return(std/  pow(k/pow(std,4),0.25)  );
}


double directionality(ImageMatrix *image)
{ ImageMatrix *deltaH,*deltaV,*matrixH,*matrixV,*phi;
  double sum,sum_r;
  int a,x,y;
  unsigned int xdim=image->width;
  unsigned int ydim=image->height;
#define  NBINS 125
  double Hd[NBINS];

  deltaH=image->duplicate();
  deltaV=image->duplicate();

  matrixH=new ImageMatrix(3,3,1);
  matrixV=new ImageMatrix(3,3,1);


  //step1
  matrixH->SetInt(0,0,0,-1); matrixH->SetInt(0,1,0,-2); matrixH->SetInt(0,2,0,-1);
  matrixH->SetInt(2,0,0,1); matrixH->SetInt(2,1,0,2); matrixH->SetInt(2,2,0,-1);

  matrixV->SetInt(0,0,0,1); matrixH->SetInt(1,0,0,2); matrixH->SetInt(2,0,0,1);
  matrixV->SetInt(0,2,0,-1); matrixH->SetInt(1,2,0,-2); matrixH->SetInt(2,2,0,-1);

  deltaH->convolve(matrixH);
  deltaV->convolve(matrixV);

  //step2
  phi=new ImageMatrix(xdim,ydim,1);

  sum_r=0;
  for(y=0;y<(int)ydim;++y) {
    for(x=0;x<(int)xdim;++x) {

      if(deltaH->pixel(x,y,0).intensity>=0.0001)
      {  phi->SetInt(x,y,0,atan(deltaV->pixel(x,y,0).intensity/deltaH->pixel(x,y,0).intensity)+(M_PI/2.0+0.001)); //+0.001 because otherwise sometimes getting -6.12574e-17
         sum_r+=pow(deltaH->pixel(x,y,0).intensity,2)+pow(deltaV->pixel(x,y,0).intensity,2)+pow(phi->pixel(x,y,0).intensity,2);
      }
      else phi->SetInt(x,y,0,0.0);
    }
  }

  delete deltaH;
  delete deltaV;
  delete matrixH;
  delete matrixV;

  phi->histogram(Hd,NBINS,0);
  delete phi;

  double max=0.0;
  int fmx=0;
  for (a=0;a<NBINS;a++)
    if (Hd[a]>max)
    {  max=Hd[a];
       fmx=a;
    }

  sum=0;
  for (a=0;a<NBINS;a++)
    sum+=Hd[a]*pow(a+1-fmx,2);

  return(fabs(log(sum/sum_r+0.0000001)));

}


double efficientLocalMean(const int x,const int y,const int k, ImageMatrix *laufendeSumme)
{ int k2=k/2;

  int dimx=laufendeSumme->width;
  int dimy=laufendeSumme->height;

  //wanting average over area: (y-k2,x-k2) ... (y+k2-1, x+k2-1)
  int starty=y-k2;
  int startx=x-k2;
  int stopy=y+k2-1;
  int stopx=x+k2-1;

  if(starty<0) starty=0;
  if(startx<0) startx=0;
  if(stopx>dimx-1) stopx=dimx-1;
  if(stopy>dimy-1) stopy=dimy-1;

  double unten, links, oben, obenlinks;

  if(startx-1<0) links=0;
  else links=laufendeSumme->pixel(startx-1,stopy,0).intensity;

  if(starty-1<0) oben=0;
  else oben=laufendeSumme->pixel(startx,stopy-1,0).intensity;

  if((starty-1 < 0) || (startx-1 <0)) obenlinks=0;
  else obenlinks=laufendeSumme->pixel(startx-1,stopy-1,0).intensity;

  unten=laufendeSumme->pixel(startx,stopy,0).intensity;

//   cout << "obenlinks=" << obenlinks << " oben=" << oben << " links=" << links << " unten=" <<unten << endl;
  int counter=(stopy-starty+1)*(stopx-startx+1);
  return (unten-links-oben+obenlinks)/counter;
}


/* coarseness
   hist -array of double- a pre-allocated array of "nbins" enetries
*/
double coarseness(ImageMatrix *image, double *hist,int nbins)
{ int x,y,k,max;
#define K_VALUE 7
// K_VALUE can also be 5
  const unsigned int yDim=image->height;
  const unsigned int xDim=image->width;
  double sum=0.0;
  ImageMatrix *laufendeSumme,*Ak[K_VALUE],*Ekh[K_VALUE],*Ekv[K_VALUE],*Sbest;

//  laufendeSumme=new ImageMatrix;
  laufendeSumme=image->duplicate();

  // initialize for running sum calculation
  double links, oben, obenlinks;
  for(y=0;y<(int)yDim;++y)
  { for(x=0;x<(int)xDim;++x)
    { if(x<1) links=0;
      else links=laufendeSumme->pixel(x-1,y,0).intensity;

      if(y<1) oben=0;
      else oben=laufendeSumme->pixel(x,y-1,0).intensity;

      if(y<1 || x<1) obenlinks=0;
      else obenlinks=laufendeSumme->pixel(x-1,y-1,0).intensity;

      laufendeSumme->SetInt(x,y,0,image->pixel(x,y,0).intensity+links+oben-obenlinks);
    }
  }

  for (k=1;k<=K_VALUE;k++)
  {  Ak[k-1]=new ImageMatrix(xDim,yDim,1);
     Ekh[k-1]=new ImageMatrix(xDim,yDim,1);
     Ekv[k-1]=new ImageMatrix(xDim,yDim,1);
  }
  Sbest=new ImageMatrix(image->width,image->height,1);


  //step 1
  int lenOfk=1;
  for(k=1;k<=K_VALUE;++k)
  { lenOfk*=2;
    for(y=0;y<(int)yDim;++y)
      for(x=0;x<(int)xDim;++x)
        Ak[k-1]->SetInt(x,y,0,efficientLocalMean(x,y,lenOfk,laufendeSumme));
  }

  //step 2
  lenOfk=1;
  for(k=1;k<=K_VALUE;++k)
  { int k2=lenOfk;
    lenOfk*=2;
    for(y=0;y<(int)yDim;++y)
    { for(x=0;x<(int)xDim;++x)
      { int posx1=x+k2;
        int posx2=x-k2;

        int posy1=y+k2;
        int posy2=y-k2;
        if(posx1<(int)xDim && posx2>=0)
          Ekh[k-1]->SetInt(x,y,0,fabs(Ak[k-1]->pixel(posx1,y,0).intensity-Ak[k-1]->pixel(posx2,y,0).intensity));
          else Ekh[k-1]->SetInt(x,y,0,0);
        if(posy1<(int)yDim && posy2>=0)
          Ekv[k-1]->SetInt(x,y,0,fabs(Ak[k-1]->pixel(x,posy1,0).intensity-Ak[k-1]->pixel(x,posy2,0).intensity));
          else Ekv[k-1]->SetInt(x,y,0,0);
      }
    }
  }
  //step3
  for(y=0;y<(int)yDim;++y)
  { for(x=0;x<(int)xDim;++x)
    { double maxE=0;
      int maxk=0;
      for(int k=1;k<=K_VALUE;++k) {
        if(Ekh[k-1]->pixel(x,y,0).intensity>maxE)
        { maxE=Ekh[k-1]->pixel(x,y,0).intensity;
          maxk=k;
        }
        if(Ekv[k-1]->pixel(x,y,0).intensity>maxE)
        { maxE=Ekv[k-1]->pixel(x,y,0).intensity;
          maxk=k;
        }
      }
      Sbest->SetInt(x,y,0,maxk);
      sum+=maxk;
    }
  }

  /* calculate the average coarseness */
  if (yDim==32 || xDim==32) sum/=((xDim+1-32)*(yDim+1-32));     /* prevent division by zero */
  else sum/=((yDim-32)*(xDim-32));
  /* calculate the 3-bin histogram */
  Sbest->histogram(hist,nbins,0);
  /* normalize the 3-bin histogram */
  max=(int)-INF;
  for (k=0;k<nbins;k++)
    if (hist[k]>max) max=(int)(hist[k]);
  for (k=0;k<nbins;k++)
    hist[k]=hist[k]/max;
  /* free allocated memory */
  delete laufendeSumme;
  for (k=1;k<=K_VALUE;k++)
  {  delete Ak[k-1];
     delete Ekh[k-1];
     delete Ekv[k-1];
  }
  delete Sbest;
  return(sum);  /* return the mean coarseness */
}




/* Tamura3Sigs
   vec -array of double- a pre-allocated array of 6 doubles
*/
void Tamura3Sigs2D(ImageMatrix *Im, double *vec)
{  double temp[6];
   temp[0]=coarseness(Im,&(temp[1]),3);
   temp[4]=directionality(Im);
   temp[5]=contrast(Im);

   /* rearange the order of the value so it will fit OME */
   vec[0]=temp[1];
   vec[1]=temp[2];
   vec[2]=temp[3];
   vec[3]=temp[5];
   vec[4]=temp[4];
   vec[5]=temp[0];
}

//---------------------------------------------------------------------------

#pragma package(smart_init)

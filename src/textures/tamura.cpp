//---------------------------------------------------------------------------

#include <cmath>
#include <cfloat> // DBL_MAX
#include "../statistics/Moments.h"
#include "tamura.h"



double contrast(const ImageMatrix &image) {
	unsigned int x,y;
	double z[4];
	Moments4 statMoments;
	readOnlyPixels pix_plane = image.ReadablePixels();

	for (x = 0; x < image.width; x++) {
		for (y = 0; y < image.height; y++) {
			statMoments.add (pix_plane(y,x));
		}
	}

// alternative for above:
//	ReadablePixels().unaryExpr (Moments4func(statMoments)).sum();
	statMoments.momentVector(z);

	if (z[1] < DBL_EPSILON) return(0);
	else return(z[1] / pow (z[3] / pow (z[1], 4), 0.25));
}


#define  NBINS 125
double directionality(const ImageMatrix &image) {
	double sum,sum_r, val;
	long a;
	unsigned int x, xdim = image.width;
	unsigned int y, ydim = image.height;
	double Hd[NBINS];

	ImageMatrix deltaH;
	deltaH.copy(image);
	ImageMatrix deltaV;
	deltaV.copy(image);
	

	pixDataMat matrixH (3,3);
	matrixH.setZero();

	pixDataMat matrixV (3,3);
	matrixV.setZero();


	//step1
	matrixH(0,0) = -1; matrixH(1,0) = -2; matrixH(2,0) = -1;
	matrixH(0,2) =  1; matrixH(1,2) =  2; matrixH(2,2) = -1;

	matrixV(0,0) =  1; matrixH(0,1) =  2; matrixH(0,2) =  1;
	matrixV(2,0) = -1; matrixH(2,1) = -2; matrixH(2,2) = -1;

	deltaH.convolve(matrixH);
	deltaV.convolve(matrixV);
	deltaH.finish();
	deltaV.finish();
	readOnlyPixels deltaH_pix_plane = deltaH.ReadOnlyPixels();
	readOnlyPixels deltaV_pix_plane = deltaV.ReadOnlyPixels();
	

	//step2
	ImageMatrix phi;
	phi.allocate (xdim, ydim);
	writeablePixels phi_pix_plane = phi.WriteablePixels();
	Moments2 phi_stats;

	sum_r = 0;
	for (y = 0; y < ydim; ++y) {
		for (x = 0; x < xdim; ++x) {
			if (deltaH_pix_plane(y,x) >= 0.0001) {
				val = atan(deltaV_pix_plane(y,x) / deltaH_pix_plane(y,x))+(M_PI/2.0+0.001); //+0.001 because otherwise sometimes getting -6.12574e-17
				phi_pix_plane(y,x) = phi_stats.add (val);
				sum_r += pow(deltaH_pix_plane(y,x),2)+pow(deltaV_pix_plane(y,x),2)+pow(val,2);
			} else phi_pix_plane(y,x) = phi_stats.add (0.0);
		}
	}
	phi.stats = phi_stats;
	phi.finish();
	phi.histogram(Hd,NBINS,0);

	double max = 0.0;
	long fmx = 0;
	for (a = 0; a < NBINS; a++) {
		if (Hd[a] > max) {
			max = Hd[a];
			fmx = a;
		}
	}

	sum = 0;
	for (a = 0; a < NBINS; a++)
		sum += Hd[a]*pow((double)(a+1-fmx),2);

	return(fabs(log(sum/sum_r+0.0000001)));

}


double efficientLocalMean(const long x,const long y,const long k, const pixDataMat &laufendeSumme) {
	long k2 = k/2;

	long dimx = laufendeSumme.cols();
	long dimy = laufendeSumme.rows();

	//wanting average over area: (y-k2,x-k2) ... (y+k2-1, x+k2-1)
	long starty = y-k2;
	long startx = x-k2;
	long stopy = y+k2-1;
	long stopx = x+k2-1;

	if (starty < 0) starty = 0;
	if (startx < 0) startx = 0;
	if (stopx > dimx-1) stopx = dimx-1;
	if (stopy > dimy-1) stopy = dimy-1;

	double unten, links, oben, obenlinks;

	if (startx-1 < 0) links = 0;
	else links = laufendeSumme(stopy,startx-1);

	if (starty-1 < 0) oben = 0;
	else oben = laufendeSumme(stopy-1,startx);

	if ((starty-1 < 0) || (startx-1 <0)) obenlinks = 0;
	else obenlinks = laufendeSumme(stopy-1,startx-1);

	unten = laufendeSumme(stopy,startx);

	//   cout << "obenlinks = " << obenlinks << " oben = " << oben << " links = " << links << " unten = " <<unten << endl;
	long counter = (stopy-starty+1)*(stopx-startx+1);
	return (unten-links-oben+obenlinks)/counter;
}


/* coarseness
   hist -array of double- a pre-allocated array of "nbins" enetries
*/
// K_VALUE can also be 5
#define K_VALUE 7
double coarseness(const ImageMatrix &image, double *hist,unsigned int nbins) {
	unsigned int x,y,k;
	int max;
	const unsigned int yDim = image.height;
	const unsigned int xDim = image.width;
	double sum = 0.0;
	ImageMatrix *Sbest;
	pixDataMat laufendeSumme (yDim,xDim);
	pixDataMat *Ak[K_VALUE], *Ekh[K_VALUE], *Ekv[K_VALUE];

	readOnlyPixels image_pix_plane = image.ReadablePixels();

	// initialize for running sum calculation
	double links, oben, obenlinks;
	for(y = 0; y < yDim; ++y) {
		for(x = 0; x < xDim; ++x) {
			if(x < 1) links = 0;
			else links = laufendeSumme(y,x-1);

			if(y < 1) oben = 0;
			else oben = laufendeSumme(y-1,x);

			if(y < 1 || x < 1) obenlinks = 0;
			else obenlinks = laufendeSumme(y-1,x-1);

			laufendeSumme(y,x) = image_pix_plane(y,x) + links + oben - obenlinks;
		}
	}

	for (k = 1; k <= K_VALUE; k++) {
		Ak[k-1] = new pixDataMat(yDim,xDim);
		Ekh[k-1] = new pixDataMat(yDim,xDim);
		Ekv[k-1] = new pixDataMat(yDim,xDim);
	}
	Sbest = new ImageMatrix;
	Sbest->allocate (image.width,image.height);


	//step 1
	int lenOfk = 1;
	for(k = 1; k <= K_VALUE; ++k) {
		lenOfk *= 2;
		pixDataMat &Ak_pix_plane = *Ak[k-1];
		for(y = 0; y < yDim; ++y)
			for(x = 0; x < xDim; ++x)
				Ak_pix_plane(y,x) = efficientLocalMean(x,y,lenOfk,laufendeSumme);
	}

	//step 2
	lenOfk = 1;
	for(k = 1; k <= K_VALUE; ++k) {
		int k2 = lenOfk;
		lenOfk *= 2;
		pixDataMat &Ekh_pix_plane = *Ekh[k-1];
		pixDataMat &Ekv_pix_plane = *Ekv[k-1];
		pixDataMat &Ak_pix_plane = *Ak[k-1];
		for(y = 0; y < yDim; ++y) {
			for(x = 0; x < xDim; ++x) {
				int posx1 = x+k2;
				int posx2 = x-k2;

				int posy1 = y+k2;
				int posy2 = y-k2;
				if(posx1 < (int)xDim && posx2 >= 0)
					Ekh_pix_plane(y,x) = fabs(Ak_pix_plane(y,posx1) - Ak_pix_plane(y,posx2));
				else Ekh_pix_plane(y,x) = 0;
				if(posy1 < (int)yDim && posy2 >= 0)
					Ekv_pix_plane(y,x) = fabs(Ak_pix_plane(posy1,x) - Ak_pix_plane(posy2,x));
				else Ekv_pix_plane(y,x) = 0;
			}
		}
	}

	//step3
	writeablePixels Sbest_pix_plane = Sbest->WriteablePixels();
	Moments2 Sbest_stats;
	for(y = 0; y < yDim; ++y) {
		for(x = 0; x < xDim; ++x) {
			double maxE = 0;
			int maxk = 0;
			for(int k = 1; k <= K_VALUE; ++k) {
				double Ekh_val = (*(Ekh[k-1]))(y,x);
				double Ekv_val = (*(Ekv[k-1]))(y,x);
				if(Ekh_val > maxE) {
					maxE = Ekh_val;
					maxk = k;
				}
				if(Ekv_val > maxE) {
					maxE = Ekv_val;
					maxk = k;
				}
			}
			Sbest_pix_plane(y,x) = Sbest_stats.add(maxk);
			sum += maxk;
		}
	}
	Sbest->stats = Sbest_stats;
	Sbest->finish();

	/* calculate the average coarseness */
	if (yDim == 32 || xDim == 32) sum /= ((xDim+1-32)*(yDim+1-32));     /* prevent division by zero */
	else sum /= ((yDim-32)*(xDim-32));

	/* calculate the 3-bin histogram */
	Sbest->histogram(hist,nbins,0);

	/* normalize the 3-bin histogram */
	max = (int)-INF;
	for (k = 0; k < nbins; k++)
		if (hist[k] > max) max = (int)(hist[k]);
	for (k = 0; k < nbins; k++)
		hist[k] = hist[k]/max;

	/* free allocated memory */
	for (k = 1; k <= K_VALUE; k++) {
		delete Ak[k-1];
		delete Ekh[k-1];
		delete Ekv[k-1];
	}
	delete Sbest;
	return(sum);  /* return the mean coarseness */
}




/* Tamura3Sigs
   vec -array of double- a pre-allocated array of 6 doubles
*/
void Tamura3Sigs2D(const ImageMatrix &Im, double *vec) {
	double temp[6];
	// to keep this method from modifying the const Im object, we use GetStats on a local Moments2 object
	Moments2 local_stats;
	Im.GetStats (local_stats);
	double min_val = local_stats.min();
	double max_val = local_stats.max();
	ImageMatrix normImg;

	normImg.allocate (Im.width, Im.height);
	normImg.WriteablePixels() = ((Im.ReadablePixels().array() - min_val) / max_val).unaryExpr (Moments2func(normImg.stats));

	temp[0] = coarseness(normImg,&(temp[1]),3);
	temp[4] = directionality(normImg);
	temp[5] = contrast(normImg);

	/* rearange the order of the value so it will fit OME */
	// {"Coarseness_Hist_Bin_00","Coarseness_Hist_Bin_01","Coarseness_Hist_Bin_02","Contrast","Directionality","Total_Coarseness"}
	vec[0] = temp[1];
	vec[1] = temp[2];
	vec[2] = temp[3];
	vec[3] = temp[5];
	vec[4] = temp[4];
	vec[5] = temp[0];
}

//---------------------------------------------------------------------------

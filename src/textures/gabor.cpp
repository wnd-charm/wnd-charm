
#include <math.h>
#include "cmatrix.h"
#include "gabor.h"


//  conv2comp - conv2 when the smaller matrix is of complex numbers

//    DOUBLE *c;	/* Result matrix (ma+mb-1)-by-(na+nb-1) */
//    DOUBLE *a;	/* Larger matrix */
//    DOUBLE *b;	/* Smaller matrix */
//    INT ma;		/* Row size of a */
//    INT na;		/* Column size of a */
//    INT mb;		/* Row size of b */
//    INT nb;		/* Column size of b */
//    INT plusminus;	/* add or subtract from result */
//    int *flopcnt;	/* flop count */

void conv2comp(double *c, double *a, double *b, int na, int ma, int nb, int mb) {
	double *p,*q;	/* Pointer to elements in 'a' and 'c' matrices */
	double wr,wi;     	/* Imaginary and real weights from matrix b    */
	int mc,nc;
	int k,l,i,j;
	double *r;				/* Pointer to elements in 'b' matrix */

	mc = ma+mb-1;
	nc = (na+nb-1)*2;

	/* initalize the output matrix */
	for (j = 0; j < mc; ++j)     /* For each element in b */
		for (i = 0; i < nc; ++i)
			c[j*nc+i] = 0;

    /* Perform convolution */
	r = b;
	for (j = 0; j < mb; ++j) {    /* For each element in b */
		for (i = 0; i < nb; ++i) {
			wr = *(r++);			/* Get weight from b matrix */
			wi = *(r++);
			p = c + j*nc + i*2;                 /* Start at first row of a in c. */
			q = a;
			for (l = 0; l < ma; l++) {               /* For each row of a ... */
				for (k = 0; k < na; k++) {
					*(p++) += *(q) * wr;	        /* multiply by the real weight and add.      */
					*(p++) += *(q++) * wi;       /* multiply by the imaginary weight and add. */
				}
				p += (nb-1)*2;	                /* Jump to next row position of a in c */
				//		*flopcnt += 2*ma*na;
			}
		}
	}
}


//  conv2 - the conv2 matlab function

//    DOUBLE *c;	/* Result matrix (ma+mb-1)-by-(na+nb-1) */
//    DOUBLE *a;	/* Larger matrix */
//    DOUBLE *b;	/* Smaller matrix */
//    INT ma;		/* Row size of a */
//    INT na;		/* Column size of a */
//    INT mb;		/* Row size of b */
//    INT nb;		/* Column size of b */
//    INT plusminus;	/* add or subtract from result */
//    int *flopcnt;	/* flop count */

void conv2(double *c, double *a, double *b, int ma, int na, int mb, int nb, int plusminus) {
	double *p,*q;	/* Pointer to elements in 'a' and 'c' matrices */
	double w;		/* Weight (element of 'b' matrix) */
	int mc,nc;
	int k,l,i,j;
	double *r;				/* Pointer to elements in 'b' matrix */

	mc = ma+mb-1;
	nc = na+nb-1;

	/* Perform convolution */

	r = b;
	for (j = 0; j < nb; ++j) {			/* For each non-zero element in b */
		for (i = 0; i < mb; ++i) {
			w = *(r++);				/* Get weight from b matrix */
			if (w != 0.0) {
				p = c + i + j*mc;	/* Start at first column of a in c. */
				for (l = 0, q = a; l < na; l++) {		/* For each column of a ... */
					for (k = 0; k < ma; k++) {
						*(p++) += *(q++) * w * plusminus;	/* multiply by weight and add. */
					}
					p += mb - 1;	/* Jump to next column position of a in c */
				}
				//		*flopcnt += 2*ma*na;
			} /* end if */
		}
    }
}


/*
Creates a non-normalized Gabor filter
*/

//function Gex = Gabor(f0,sig2lam,gamma,theta,fi,n),
double *Gabor(double f0, double sig2lam, double gamma, double theta, double fi, int n) {
	double *tx,*ty;
	double lambda = 2*M_PI/f0;
	double cos_theta = cos(theta), sin_theta = sin(theta);
	double sig = sig2lam * lambda;
	double sum;
	double *Gex;
	int x,y;
	int nx = n;
	int ny = n;
	tx = new double[nx+1];
	ty = new double[ny+1];

	if (nx%2 > 0) {
		tx[0] = -((nx-1)/2);
		for (x = 1; x < nx; x++)
			tx[x] = tx[x-1]+1;
	} else {
		tx[0] = -(nx/2);
		for (x = 1; x <= nx; x++)
			tx[x] = tx[x-1]+1;
	}

	if (ny%2 > 0) {
		ty[0] = -((ny-1)/2);
		for (y = 1; y < ny; y++)
			ty[y] = ty[y-1]+1;
	} else {
		ty[0] = -(ny/2);
		for (y = 1; y <= ny; y++)
			ty[y] = ty[y-1]+1;
	}

	Gex = new double [n*n*2];

	sum = 0;
	for (y = 0; y < n; y++) {
		for (x = 0; x < n; x++) {
			double argm,xte,yte,rte,ge;
			xte = tx[x]*cos_theta+ty[y]*sin_theta;
			yte = ty[y]*cos_theta-tx[x]*sin_theta;
			rte = xte*xte+gamma*gamma*yte*yte;
			ge = exp(-1*rte/(2*sig*sig));
			argm = xte*f0+fi;
			Gex[y*n*2+x*2] = ge*cos(argm);             // ge .* exp(j.*argm);
			Gex[y*n*2+x*2+1] = ge*sin(argm);
			sum += sqrt(pow(Gex[y*n*2+x*2],2)+pow(Gex[y*n*2+x*2+1],2));
		}
	}

	/* normalize */
	for (y = 0; y < n; y++)
		for (x = 0; x < n*2; x += 1)
			Gex[y*n*2+x] = Gex[y*n*2+x]/sum;

	delete [] tx;
	delete [] ty;

	return(Gex);
}

/* Computes Gabor energy */
//Function [e2] = GaborEnergy(Im,f0,sig2lam,gamma,theta,n),
double *GaborEnergy(const ImageMatrix &Im, double f0, double sig2lam, double gamma, double theta, int n) {
	double *Gexp, *image, *c,*out;
	double fi = 0;
	unsigned int a,b,x,y;
	Gexp = Gabor(f0,sig2lam,gamma,theta,fi,n);
	readOnlyPixels pix_plane = Im.ReadablePixels();


	c = new double[(Im.width+n-1)*(Im.height+n-1)*2];
	image = new double[Im.width*Im.height];
	for (y = 0; y < Im.height; y++)
		for (x = 0; x < Im.width; x++)
			image[y*Im.width+x] = pix_plane(y,x);

	conv2comp(c, image, Gexp, Im.width, Im.height, n, n);

	out = new double[Im.height*Im.width];
	b = 0;
	for (y = (int)ceil((double)n/2); b < Im.height; y++) {
		a = 0;
		for (x = (int)ceil((double)n/2); a < Im.width; x++) {
			out[b*Im.width+a] = sqrt(pow(c[y*2*(Im.width+n-1)+x*2],2)+pow(c[y*2*(Im.width+n-1)+x*2+1],2));
			a++;
		}
		b++;
	}

	delete [] image;
	delete [] Gexp;
	delete [] c;
	return(out);
}

//---------------------------------------------------------------------------
/*
the output value is in "ratios" which is an array of 7 doubles
*/
void GaborTextureFilters2D(const ImageMatrix &Im, double *ratios) {
	double GRAYthr;
	/* parameters set up in complience with the paper */
	double gamma = 0.5, sig2lam = 0.56;
	int n = 38;
	double f0[7] = {1,2,3,4,5,6,7};       // frequencies for several HP Gabor filters
	double f0LP = 0.1;     // frequencies for one LP Gabor filter
	double theta = 3.14159265/2;
	double *e2LP;
	unsigned int ii;
	unsigned long originalScore = 0;

	ImageMatrix e2img;

	// compute the original score before Gabor
	e2LP = GaborEnergy(Im,f0LP,sig2lam,gamma,theta,n);
	e2img.remap_pix_plane(e2LP, Im.width, Im.height);
	readOnlyPixels pix_plane = e2img.ReadablePixels();
	// N.B.: for the base of the ratios, the threshold is 0.4 of max energy,
	// while the comparison thresholds are Otsu.
	originalScore = (pix_plane.array() > pix_plane.maxCoeff() * 0.4).count();

	// cleanup the mapping
	e2img.remap_pix_plane(NULL, 0, 0);
	delete [] e2LP;

	for (ii = 0; ii < 7; ii++) {
		double *e2;
		unsigned long afterGaborScore = 0;
		e2 = GaborEnergy(Im,f0[ii],sig2lam,gamma,theta,n);
		e2img.remap_pix_plane(e2, Im.width, Im.height);
		writeablePixels e2_pix_plane = e2img.WriteablePixels();
		e2_pix_plane.array() /= e2_pix_plane.maxCoeff();
		GRAYthr = e2img.Otsu();
		afterGaborScore = (e2_pix_plane.array() > GRAYthr).count();
		ratios[ii] = (double)afterGaborScore/(double)originalScore;

		e2img.remap_pix_plane(NULL, 0, 0);
		delete [] e2;
	}
}

/***************************************************************************
* ======================================================================
* Computer Vision/Image Processing Tool Project - Dr. Scott Umbaugh SIUE
* ======================================================================
*
*             File Name: CVIPtexture.h
*           Description: contains function prototypes, type names, constants,
*			 etc. related to libdataserv (Data Services Toolkit.)
*         Related Files: Imakefile, cvip_pgmtexture.c
*   Initial Coding Date: 6/19/96
*           Portability: Standard (ANSI) C
*             Credit(s): Steve Costello
*                        Southern Illinois University @ Edwardsville
*
** Copyright (C) 1993 SIUE - by Gregory Hance.
**
** Permission to use, copy, modify, and distribute this software and its
** documentation for any purpose and without fee is hereby granted, provided
** that the above copyright notice appear in all copies and that both that
** copyright notice and this permission notice appear in supporting
** documentation.  This software is provided "as is" without express or
** implied warranty.
**
****************************************************************************/
#ifndef _CVIP_texture
#define _CVIP_texture
/* [0] -> 0 degree, [1] -> 45 degree, [2] -> 90 degree, [3] -> 135 degree,
   [4] -> average, [5] -> range (max - min) */

typedef unsigned char u_int8_t;
   
typedef struct  {
	double ASM;          /*  (1) Angular Second Moment */
	double contrast;     /*  (2) Contrast */
	double correlation;  /*  (3) Correlation */
	double variance;     /*  (4) Variance */
	double IDM;		    /*  (5) Inverse Diffenence Moment */
	double sum_avg;	    /*  (6) Sum Average */
	double sum_var;	    /*  (7) Sum Variance */
	double sum_entropy;	/*  (8) Sum Entropy */
	double entropy;	    /*  (9) Entropy */
	double diff_var;	    /* (10) Difference Variance */
	double diff_entropy;	/* (11) Diffenence Entropy */
	double meas_corr1;	/* (12) Measure of Correlation 1 */
	double meas_corr2;	/* (13) Measure of Correlation 2 */
	double max_corr_coef; /* (14) Maximal Correlation Coefficient */
	} TEXTURE;

TEXTURE * Extract_Texture_Features(int distance, int angle,
		 		register u_int8_t **grays, unsigned int nrows, unsigned int ncols);

#endif

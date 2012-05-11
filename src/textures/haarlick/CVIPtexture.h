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
	float ASM;          /*  (1) Angular Second Moment */
	float contrast;     /*  (2) Contrast */
	float correlation;  /*  (3) Correlation */
	float variance;     /*  (4) Variance */
	float IDM;		    /*  (5) Inverse Diffenence Moment */
	float sum_avg;	    /*  (6) Sum Average */
	float sum_var;	    /*  (7) Sum Variance */
	float sum_entropy;	/*  (8) Sum Entropy */
	float entropy;	    /*  (9) Entropy */
	float diff_var;	    /* (10) Difference Variance */
	float diff_entropy;	/* (11) Diffenence Entropy */
	float meas_corr1;	/* (12) Measure of Correlation 1 */
	float meas_corr2;	/* (13) Measure of Correlation 2 */
	float max_corr_coef; /* (14) Maximal Correlation Coefficient */
	} TEXTURE;

TEXTURE * Extract_Texture_Features(int distance, int angle,
		 		register u_int8_t **grays, int rows, int cols, int max_val);

#endif

/*
 * Included in OME by Tom Macura so that Michael Boland's mb_Znl.cpp doesn't 
 * require system specific complex arithmetic library.
 *
 * http://www.csounds.com/developers/html/complex_8c-source.html
 * Code from Press, Teukolsky, Vettering and Flannery
 * Numerical Recipes in C, 2nd Edition, Cambridge 1992.
*/

#ifndef _NR_COMPLEX_H_
#define _NR_COMPLEX_H_

#ifndef _FCOMPLEX_DECLARE_T_
typedef struct FCOMPLEX {double r,i;} fcomplex;
#define _FCOMPLEX_DECLARE_T_
#endif /* _FCOMPLEX_DECLARE_T_ */

fcomplex Cadd(fcomplex a, fcomplex b);
fcomplex Csub(fcomplex a, fcomplex b);
fcomplex Cmul(fcomplex a, fcomplex b);
fcomplex Complex(double re, double im);
fcomplex Conjg(fcomplex z);
fcomplex Cdiv(fcomplex a, fcomplex b);
double Cabs(fcomplex z);
fcomplex Csqrt(fcomplex z);
fcomplex RCmul(double x, fcomplex a);
fcomplex Rpolar (double rho, double theta);
#endif /* _NR_COMPLEX_H_ */

//---------------------------------------------------------------------------

#ifndef zernikeH
#define zernikeH
//---------------------------------------------------------------------------

void mb_zernike2D(const ImageMatrix &I, double D, double R, double *zvalues, long *output_size);
void mb_zernike2D_OLD(const ImageMatrix &I, double D, double R, double *zvalues, long *output_size);

#endif

//---------------------------------------------------------------------------

#ifndef gaborH
#define gaborH
//---------------------------------------------------------------------------

void GaborTextureFilters2D(const ImageMatrix &Im, double *ratios);

void conv2(double *c, double *a, double *b, int ma, int na, int mb, int nb, int plusminus);

extern "C" void gpu_conv2comp(double *c, double *a, double *b, int na, int ma, int nb, int mb);
#endif

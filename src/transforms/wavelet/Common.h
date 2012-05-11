
#ifndef _COMMON_H_
#define _COMMON_H_

static const short WAVELET_HIGH = 1;
static const short WAVELET_MED = 2;
static const short WAVELET_LOW =  3;

static const bool COMPRESS  = true;
static const bool DECOMPRESS = false;

static const int PORT = 9876;
static const int slave1port = 9870;
static const int slave2port = 9871;

extern bool verbose;
extern bool parallel;
extern bool master;
extern char * slave1;
extern char * slave2;

// This covers bzero defines used in Wavelet.cpp, and DataGrid*D.cpp
#include <string.h>
#define bzero(b,len) (memset((b), '\0', (len)), (void) 0)


#endif


#ifndef __WNDCHRM_ERROR_H__
#define __WNDCHRM_ERROR_H__

#include <string>

enum WNDCHRM_ERROR {
	WC_UNINITIALIZED,
	WC_NO_ERROR,
	WC_IPP_NULL,
	WC_MM_FAIL_RECURSIVE_CALL,
	WC_TRANSFORM_FAIL,
	WC_EMPTY,
	WC_NOT_IMPLEMENTED,
	WC_INPUT_IMAGEMATRIX_NULL
};

extern int verbosity;
void catErrno ();
void catError (const char *fmt, ...);
void catError (const std::string &error);
int showError(int stop, const char *fmt, ...);
const std::string getErrorString ();
const char* translateError( WNDCHRM_ERROR return_val );

#endif // __WNDCHRM_ERROR_H__

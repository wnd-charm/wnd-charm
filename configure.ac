
dnl Process this file with autoconf to produce a configure script.
AC_INIT(pychrm, 0.1, christopher.coletta@nih.gov)
AC_CONFIG_MACRO_DIR(m4)
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_SRCDIR(src/cmatrix.cpp)

AM_INIT_AUTOMAKE
AC_PROG_CC


AC_PROG_CXX
AM_PROG_LEX
AC_PROG_YACC
AC_DISABLE_STATIC
AC_PROG_LIBTOOL
AX_PYTHON_DEVEL(>= '2.7')
AM_PATH_PYTHON
AX_PKG_SWIG(2.0, [], [ AC_MSG_ERROR([SWIG is required to build..]) ])
AC_PROG_SWIG
SWIG_ENABLE_CXX
SWIG_PYTHON

dnl ------------------------------------------------------------
dnl check for FFTW3
dnl ------------------------------------------------------------
AC_CHECK_LIB(fftw3,fftw_execute,[AC_DEFINE(HAVE_FFTW3,1,"Turn on FFTW3 support.")],
	[AC_MSG_ERROR([
    FFTW3 library not found.
    Please download+install FFTW 3.x from http://www.fftw.org/download.html
])]
)
AC_CHECK_HEADER(fftw3.h, [],
        AC_MSG_ERROR([
    the <fftw3.h> header file from FFTW3 could not be found.
    Please download+install FFTW 3.x from http://www.fftw.org/download.html
]))


dnl ------------------------------------------------------------
dnl check for libtiff
dnl ------------------------------------------------------------
AC_CHECK_LIB(tiff, TIFFReadScanline)
AC_CHECK_HEADER(tiffio.h,LIBTIFF='yes',AC_MSG_ERROR([
    *** WND_CHARM requires libtiff; (Required TIFF headers not found) ***
    Please download+install libtiff from http://www.libtiff.org/
]))
dnl DX_HTML_FEATURE(ON)
dnl DX_CHM_FEATURE(OFF)
dnl DX_CHI_FEATURE(OFF)
dnl DX_MAN_FEATURE(OFF)
dnl DX_RTF_FEATURE(OFF)
dnl DX_XML_FEATURE(OFF)
dnl DX_PDF_FEATURE(OFF)
dnl DX_PS_FEATURE(OFF)

AC_OUTPUT(Makefile src/Makefile pychrm/Makefile)

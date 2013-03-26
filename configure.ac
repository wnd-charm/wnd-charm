
dnl Process this file with autoconf to produce a configure script.
AC_INIT(pychrm, 0.1, christopher.coletta@nih.gov)
AC_CONFIG_MACRO_DIR(m4)
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_SRCDIR(src/cmatrix.cpp)

AM_INIT_AUTOMAKE
AC_PROG_CC

AM_CONFIG_HEADER(config.h)

AC_PROG_CXX
AM_PROG_LEX
AC_PROG_YACC
AC_DISABLE_STATIC
AC_PROG_LIBTOOL
AX_PYTHON_DEVEL(>= '2.6')
AM_PATH_PYTHON
AX_PKG_SWIG(1.3, [], [ AC_MSG_ERROR([SWIG is required to build..]) ])
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

dnl ------------------------------------------------------------
dnl check for TR1 stuff individually
dnl ------------------------------------------------------------
AC_DEFUN([AC_HEADER_TR1_UNORDERED_MAP], [
  AC_CACHE_CHECK(for tr1/unordered_map,
  ac_cv_cxx_tr1_unordered_map,
  [AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  AC_TRY_COMPILE([#include <tr1/unordered_map>], [using std::tr1::unordered_map;],
  ac_cv_cxx_tr1_unordered_map=yes, ac_cv_cxx_tr1_unordered_map=no)
  AC_LANG_RESTORE
  ])
  if test "$ac_cv_cxx_tr1_unordered_map" = yes; then
    AC_DEFINE(HAVE_TR1_UNORDERED_MAP,,[Define if tr1/unordered_map is present. ])
  fi
])
AC_HEADER_TR1_UNORDERED_MAP
AC_DEFUN([AC_HEADER_TR1_UNORDERED_SET], [
  AC_CACHE_CHECK(for tr1/unordered_set,
  ac_cv_cxx_tr1_unordered_set,
  [AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  AC_TRY_COMPILE([#include <tr1/unordered_set>], [using std::tr1::unordered_set;],
  ac_cv_cxx_tr1_unordered_set=yes, ac_cv_cxx_tr1_unordered_set=no)
  AC_LANG_RESTORE
  ])
  if test "$ac_cv_cxx_tr1_unordered_set" = yes; then
    AC_DEFINE(HAVE_TR1_UNORDERED_SET,,[Define if tr1/unordered_set is present. ])
  fi
])
AC_HEADER_TR1_UNORDERED_SET

dnl ------------------------------------------------------------
dnl check for native unordered_map and unordered_set individually
dnl ------------------------------------------------------------
AC_DEFUN([AC_HEADER_UNORDERED_MAP], [
  AC_CACHE_CHECK(for native unordered_map,
  ac_unordered_map,
  [AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  AC_TRY_COMPILE([#include <unordered_map>], [using std::unordered_map;],
  ac_unordered_map=yes, ac_unordered_map=no)
  AC_LANG_RESTORE
  ])
  if test "$ac_unordered_map" = yes; then
    AC_DEFINE(HAVE_UNORDERED_MAP,,[Define if native unordered_map is present. ])
  fi
])
AC_HEADER_UNORDERED_MAP
AC_DEFUN([AC_HEADER_UNORDERED_SET], [
  AC_CACHE_CHECK(for native unordered_set,
  ac_unordered_set,
  [AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  AC_TRY_COMPILE([#include <unordered_set>], [using std::unordered_set;],
  ac_unordered_set=yes, ac_unordered_set=no)
  AC_LANG_RESTORE
  ])
  if test "$ac_unordered_set" = yes; then
    AC_DEFINE(HAVE_UNORDERED_SET,,[Define if native unordered_set is present. ])
  fi
])
AC_HEADER_UNORDERED_SET

AC_ARG_ENABLE(class-prob-tsv,
 [  --enable-class-prob-tsv   When creating a report, output a tsv\
  containing the Average Class Probability Matrix data for the purpose \
  of computing Morphological Divergence Scores],
 [ AC_DEFINE(AVG_CLASS_PROB_TSV,,[Define optional output. ]) 
 ])


dnl DX_HTML_FEATURE(ON)
dnl DX_CHM_FEATURE(OFF)
dnl DX_CHI_FEATURE(OFF)
dnl DX_MAN_FEATURE(OFF)
dnl DX_RTF_FEATURE(OFF)
dnl DX_XML_FEATURE(OFF)
dnl DX_PDF_FEATURE(OFF)
dnl DX_PS_FEATURE(OFF)


AC_OUTPUT(Makefile src/Makefile pychrm/Makefile)

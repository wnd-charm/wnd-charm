#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([wndchrm], [1.60], [christopher.coletta@nih.gov])
AM_INIT_AUTOMAKE
AM_PATH_PYTHON
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

LT_INIT

# Checks for programs.
AC_PROG_CXX
AC_PROG_AWK
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lchrm':
AC_CHECK_LIB([chrm], [main])
# FIXME: Replace `main' with a function in `-lfftw3':
AC_CHECK_LIB([fftw3], [main])
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [main])
# FIXME: Replace `main' with a function in `-ltiff':
AC_CHECK_LIB([tiff], [main])

# Checks for header files.
AC_CHECK_HEADERS([fcntl.h limits.h memory.h stddef.h stdlib.h string.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT32_T
AC_TYPE_MODE_T
AC_C_RESTRICT
AC_TYPE_SIZE_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_HEADER_MAJOR
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRERROR_R
AC_FUNC_STRTOD
AC_CHECK_FUNCS([bzero floor ftruncate memmove memset mkdir pow select sqrt strchr strerror strrchr strstr strtol])

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

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT

AC_INIT(wndchrm,1.60)
AM_INIT_AUTOMAKE
AC_PROG_CC
AC_PROG_RANLIB

dnl set the default c++ flags
dnl Need to make the default CXXFLAGS be blank instead of -g -O2.  Use these in AM_CXXFLAGS instead
if test x"${CXXFLAGS+set}" = xset; then 
  # the user set CXXFLAGS; don't override it.
  # Much more clever would be to issue a warning if users' CXXFLAGS have optimization higher than 0
  # Cleverer still would be to figure out which versions of g++ barf when compiling FeatureNames.cpp with optimization.
	MY_CXXFLAGS="${CXXFLAGS}"
else 
	MY_CXXFLAGS=""
fi 
# This sets the default CXXFLAGS (user/environment/or built-in default)
AC_PROG_CXX
# And here we overwrite it based on our test above
CXXFLAGS=$MY_CXXFLAGS
# Put the usual defualt into AM_CXXFLAGS
AM_CXXFLAGS="-Wall -g -O2"

SVNREV=`( svnversion $srcdir | sed 's/:/_/' ) 2>/dev/null`
if { ! ( svnversion ) >/dev/null 2>&1 || test "$SVNREV" = "exported"; } ;
	then SVNREV=`cat $srcdir/SVNREV`
	else echo $SVNREV>$srcdir/SVNREV
fi
AC_SUBST(SVNREV)

AC_DEFINE_UNQUOTED([SVNREV],["$SVNREV"], [SVN Revision from svnversion command])
AC_DEFINE_UNQUOTED([PACKAGE_VERSION],["$VERSION-$SVNREV"], [Version string, including SVNREV])
AC_DEFINE_UNQUOTED([PACKAGE_STRING],["$PACKAGE_NAME $VERSION-$SVNREV"], [Package string with version+SVNREV])


AC_CONFIG_HEADERS(config.h)

AC_CHECK_HEADER(math.h, [],
        AC_MSG_ERROR([the <math.h> header file could not be found.]))

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
dnl check for C++0x
dnl ------------------------------------------------------------
# AC_COMPILE_STDCXX_OX
# AC_DEFUN([AC_COMPILE_STDCXX_0X], [
#   AC_CACHE_CHECK(if g++ supports C++0x features without additional flags,
#   ac_cv_cxx_compile_cxx0x_native,
#   [AC_LANG_SAVE
#   AC_LANG_CPLUSPLUS
#   AC_TRY_COMPILE([
#   template <typename T>
#     struct check
#     {
#       static_assert(sizeof(int) <= sizeof(T), "not big enough");
#     };
# 
#     typedef check<check<bool>> right_angle_brackets;
# 
#     int a;
#     decltype(a) b;
# 
#     typedef check<int> check_type;
#     check_type c;
#     check_type&& cr = c;],,
#   ac_cv_cxx_compile_cxx0x_native=yes, ac_cv_cxx_compile_cxx0x_native=no)
#   AC_LANG_RESTORE
#   ])
# 
#   AC_CACHE_CHECK(if g++ supports C++0x features with -std=c++0x,
#   ac_cv_cxx_compile_cxx0x_cxx,
#   [AC_LANG_SAVE
#   AC_LANG_CPLUSPLUS
#   ac_save_CXXFLAGS="$CXXFLAGS"
#   CXXFLAGS="$CXXFLAGS -std=c++0x"
#   AC_TRY_COMPILE([
#   template <typename T>
#     struct check
#     {
#       static_assert(sizeof(int) <= sizeof(T), "not big enough");
#     };
# 
#     typedef check<check<bool>> right_angle_brackets;
# 
#     int a;
#     decltype(a) b;
# 
#     typedef check<int> check_type;
#     check_type c;
#     check_type&& cr = c;],,
#   ac_cv_cxx_compile_cxx0x_cxx=yes, ac_cv_cxx_compile_cxx0x_cxx=no)
#   CXXFLAGS="$ac_save_CXXFLAGS"
#   AC_LANG_RESTORE
#   ])
# 
#   AC_CACHE_CHECK(if g++ supports C++0x features with -std=gnu++0x,
#   ac_cv_cxx_compile_cxx0x_gxx,
#   [AC_LANG_SAVE
#   AC_LANG_CPLUSPLUS
#   ac_save_CXXFLAGS="$CXXFLAGS"
#   CXXFLAGS="$CXXFLAGS -std=gnu++0x"
#   AC_TRY_COMPILE([
#   template <typename T>
#     struct check
#     {
#       static_assert(sizeof(int) <= sizeof(T), "not big enough");
#     };
# 
#     typedef check<check<bool>> right_angle_brackets;
# 
#     int a;
#     decltype(a) b;
# 
#     typedef check<int> check_type;
#     check_type c;
#     check_type&& cr = c;],,
#   ac_cv_cxx_compile_cxx0x_gxx=yes, ac_cv_cxx_compile_cxx0x_gxx=no)
#   CXXFLAGS="$ac_save_CXXFLAGS"
#   AC_LANG_RESTORE
#   ])
# 
#   if test "$ac_cv_cxx_compile_cxx0x_native" = yes ||
#      test "$ac_cv_cxx_compile_cxx0x_cxx" = yes ||
#      test "$ac_cv_cxx_compile_cxx0x_gxx" = yes; then
#     AC_DEFINE(HAVE_STDCXX_0X,,[Define if g++ supports C++0x features. ])
#   fi
#   
#   if test "$ac_cv_cxx_compile_cxx0x_cxx" = yes; then
#   	AM_CXXFLAGS="$AM_CXXFLAGS -std=c++0x"
#   elif test "$ac_cv_cxx_compile_cxx0x_gxx" = yes; then
#   	AM_CXXFLAGS="$AM_CXXFLAGS -std=gnu++0x"
#   fi
# ])
# AC_HEADER_UNORDERED_MAP
# AC_DEFUN([AC_HEADER_UNORDERED_MAP], [
#   AC_CACHE_CHECK(for unordered_map,
#   ac_cv_cxx_unordered_map,
#   [AC_REQUIRE([AC_COMPILE_STDCXX_0X])
#   AC_LANG_SAVE
#   AC_LANG_CPLUSPLUS
#   ac_save_CXXFLAGS="$CXXFLAGS"
#   CXXFLAGS="$CXXFLAGS -std=gnu++0x"
#   AC_TRY_COMPILE([#include <unordered_map>], [using std::unordered_map;],
#   ac_cv_cxx_unordered_map=yes, ac_cv_cxx_unordered_map=no)
#   CXXFLAGS="$ac_save_CXXFLAGS"
#   AC_LANG_RESTORE
#   ])
#   if test "$ac_cv_cxx_unordered_map" = yes; then
#     AC_DEFINE(HAVE_UNORDERED_MAP,,[Define if unordered_map is present. ])
#   fi
# ])
# AC_HEADER_UNORDERED_MAP
# # AC_HEADER_UNORDERED_SET
# AC_DEFUN([AC_HEADER_UNORDERED_SET], [
#   AC_CACHE_CHECK(for unordered_set,
#   ac_cv_cxx_unordered_set,
#   [AC_REQUIRE([AC_COMPILE_STDCXX_0X])
#   AC_LANG_SAVE
#   AC_LANG_CPLUSPLUS
#   ac_save_CXXFLAGS="$CXXFLAGS"
#   CXXFLAGS="$CXXFLAGS -std=gnu++0x"
#   AC_TRY_COMPILE([#include <unordered_set>], [using std::unordered_set;],
#   ac_cv_cxx_unordered_set=yes, ac_cv_cxx_unordered_set=no)
#   CXXFLAGS="$ac_save_CXXFLAGS"
#   AC_LANG_RESTORE
#   ])
#   if test "$ac_cv_cxx_unordered_set" = yes; then
#     AC_DEFINE(HAVE_UNORDERED_SET,,[Define if unordered_set is present. ])
#   fi
# ])
# AC_HEADER_UNORDERED_SET


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

# Write out our compile flags
AC_SUBST(CXXFLAGS)
AC_SUBST(AM_CXXFLAGS)

AC_OUTPUT(Makefile)

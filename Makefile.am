# echo "noinst_HEADERS = $(find Eigen  -not -path '*.svn/*'  -type f | perl -n -e 'chomp;s/^\.\///;push (@foo,$_); END{print join(" \\\n    ",@foo)}')" > Eigen/sources.am 
include src/Eigen/sources.am

EXTRA_DIST = SVNREV
distdir = $(PACKAGE)-$(VERSION).$(SVNREV)
am__tar = COPY_EXTENDED_ATTRIBUTES_DISABLE=1 COPYFILE_DISABLE=1 ${AMTAR} cof - "$$tardir"
AM_CPPFLAGS = -fPIC -Isrc

noinst_LIBRARIES = libchrm.a

libchrm_a_SOURCES = \
  src/unordered_map_dfn.h \
  src/FeatureAlgorithms.cpp \
  src/FeatureAlgorithms.h \
  src/FeatureNames.cpp \
  src/FeatureNames.h \
  src/ImageTransforms.cpp \
  src/ImageTransforms.h \
  src/OldFeatureNamesFileStr.h \
  src/Tasks.cpp \
  src/Tasks.h \
  src/cmatrix.cpp \
  src/cmatrix.h \
  src/colors/FuzzyCalc.cpp \
  src/colors/FuzzyCalc.h \
  src/gsl/specfunc.cpp \
  src/gsl/specfunc.h \
  src/statistics/CombFirst4Moments.cpp \
  src/statistics/CombFirst4Moments.h \
  src/statistics/FeatureStatistics.cpp \
  src/statistics/FeatureStatistics.h \
  src/statistics/Moments.h \
  src/textures/gabor.cpp \
  src/textures/gabor.h \
  src/textures/haralick/CVIPtexture.cpp \
  src/textures/haralick/CVIPtexture.h \
  src/textures/haralick/haralick.cpp \
  src/textures/haralick/haralick.h \
  src/textures/tamura.cpp \
  src/textures/tamura.h \
  src/textures/zernike/zernike.cpp \
  src/textures/zernike/zernike.h \
  src/transforms/ChebyshevFourier.cpp \
  src/transforms/ChebyshevFourier.h \
  src/transforms/chebyshev.cpp \
  src/transforms/chebyshev.h \
  src/transforms/fft/bcb_fftw3/fftw3.h \
  src/transforms/fft/fftw3.h \
  src/transforms/radon.cpp \
  src/transforms/radon.h \
  src/transforms/wavelet/Common.cpp \
  src/transforms/wavelet/Common.h \
  src/transforms/wavelet/DataGrid.h \
  src/transforms/wavelet/DataGrid2D.cpp \
  src/transforms/wavelet/DataGrid2D.h \
  src/transforms/wavelet/DataGrid3D.cpp \
  src/transforms/wavelet/DataGrid3D.h \
  src/transforms/wavelet/Filter.cpp \
  src/transforms/wavelet/Filter.h \
  src/transforms/wavelet/FilterSet.cpp \
  src/transforms/wavelet/FilterSet.h \
  src/transforms/wavelet/Symlet5.cpp \
  src/transforms/wavelet/Symlet5.h \
  src/transforms/wavelet/Wavelet.cpp \
  src/transforms/wavelet/Wavelet.h \
  src/transforms/wavelet/WaveletHigh.cpp \
  src/transforms/wavelet/WaveletHigh.h \
  src/transforms/wavelet/WaveletLow.cpp \
  src/transforms/wavelet/WaveletLow.h \
  src/transforms/wavelet/WaveletMedium.cpp \
  src/transforms/wavelet/WaveletMedium.h \
  src/transforms/wavelet/convolution.cpp \
  src/transforms/wavelet/convolution.h \
  src/transforms/wavelet/wt.cpp \
  src/transforms/wavelet/wt.h \
  src/wndchrm_src/TrainingSet.cpp \
  src/wndchrm_src/TrainingSet.h \
  src/wndchrm_src/WORMfile.cpp \
  src/wndchrm_src/WORMfile.h \
  src/wndchrm_src/signatures.cpp \
  src/wndchrm_src/signatures.h \
  src/wndchrm_error.cpp \
  src/wndchrm_error.h

libchrm_a_CXXFLAGS = -Wall -g -Os

bin_PROGRAMS=wndchrm

wndchrm_SOURCES = src/wndchrm_src/wndchrm.cpp

wndchrm_LDADD = libchrm.a -lm -ltiff -L. -lchrm -lfftw3

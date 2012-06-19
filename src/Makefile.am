
AM_YFLAGS = -d
AM_LFLAGS = -o$(LEX_OUTPUT_ROOT).c
AM_CPPFLAGS = -fPIC

lib_LTLIBRARIES = libchrm.la

libchrm_la_SOURCES = \
       colors/FuzzyCalc.cpp \
       colors/FuzzyCalc.h \
       statistics/CombFirst4Moments.cpp \
       statistics/CombFirst4Moments.h \
       statistics/FeatureStatistics.cpp \
       statistics/FeatureStatistics.h \
       textures/gabor.cpp \
       textures/gabor.h \
       textures/haarlick/CVIPtexture.cpp \
       textures/haarlick/CVIPtexture.h \
       textures/haarlick/haarlick.cpp \
       textures/haarlick/haarlick.h \
       textures/haarlick/mapkit.cpp \
       textures/haarlick/mapkit.h \
       textures/haarlick/mapkit_generic.cpp \
       textures/haarlick/mapkit_generic.h \
       textures/tamura.cpp \
       textures/tamura.h \
       textures/zernike/complex.cpp \
       textures/zernike/complex.h \
       textures/zernike/zernike.cpp \
       textures/zernike/zernike.h \
       transforms/ChebyshevFourier.cpp \
       transforms/ChebyshevFourier.h \
       transforms/chebyshev.cpp \
       transforms/chebyshev.h \
       transforms/fft/bcb_fftw3/fftw3.h \
       transforms/fft/fftw3.h \
       transforms/radon.cpp \
       transforms/radon.h \
       transforms/wavelet/Common.cpp \
       transforms/wavelet/Common.h \
       transforms/wavelet/convolution.cpp \
       transforms/wavelet/convolution.h \
       transforms/wavelet/DataGrid.h \
       transforms/wavelet/DataGrid2D.cpp \
       transforms/wavelet/DataGrid2D.h \
       transforms/wavelet/DataGrid3D.cpp \
       transforms/wavelet/DataGrid3D.h \
       transforms/wavelet/Filter.cpp \
       transforms/wavelet/Filter.h \
       transforms/wavelet/FilterSet.cpp \
       transforms/wavelet/FilterSet.h \
       transforms/wavelet/Symlet5.cpp \
       transforms/wavelet/Symlet5.h \
       transforms/wavelet/Wavelet.cpp \
       transforms/wavelet/Wavelet.h \
       transforms/wavelet/WaveletHigh.cpp \
       transforms/wavelet/WaveletHigh.h \
       transforms/wavelet/WaveletLow.cpp \
       transforms/wavelet/WaveletLow.h \
       transforms/wavelet/WaveletMedium.cpp \
       transforms/wavelet/WaveletMedium.h \
       transforms/wavelet/wt.cpp \
       transforms/wavelet/wt.h \
       cmatrix.cpp \
       cmatrix.h \
			 wndchrm_error.cpp \
			 wndchrm_error.h \
			 FeatureTransforms.cpp \
			 FeatureTransforms.h \
			 FeatureAlgorithms.cpp \
			 FeatureAlgorithms.h

libchrm_la_LIBADD = -ltiff -L. -lfftw3
libchrm_la_CXXFLAGS = -Wall -g -O2 -fpermissive

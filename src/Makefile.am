# echo "noinst_HEADERS = $(find Eigen  -not -path '*.svn/*'  -type f | perl -n -e 'chomp;s/^\.\///;push (@foo,$_); END{print join(" \\\n    ",@foo)}')" > Eigen/sources.am
include Eigen/sources.am

AM_CPPFLAGS = -fPIC

noinst_LIBRARIES = libchrm.a

libchrm_a_SOURCES = \
  unordered_map_dfn.h \
  FeatureAlgorithms.cpp \
  FeatureAlgorithms.h \
  FeatureNames.cpp \
  FeatureNames.h \
  ImageTransforms.cpp \
  ImageTransforms.h \
  OldFeatureNamesFileStr.h \
  Tasks.cpp \
  Tasks.h \
  cmatrix.cpp \
  cmatrix.h \
  colors/FuzzyCalc.cpp \
  colors/FuzzyCalc.h \
  gsl/specfunc.cpp \
  gsl/specfunc.h \
  statistics/CombFirst4Moments.cpp \
  statistics/CombFirst4Moments.h \
  statistics/FeatureStatistics.cpp \
  statistics/FeatureStatistics.h \
  statistics/Moments.h \
  textures/gabor.cpp \
  textures/gabor.h \
  textures/haralick/CVIPtexture.cpp \
  textures/haralick/CVIPtexture.h \
  textures/haralick/haralick.cpp \
  textures/haralick/haralick.h \
  textures/tamura.cpp \
  textures/tamura.h \
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
  transforms/wavelet/convolution.cpp \
  transforms/wavelet/convolution.h \
  transforms/wavelet/wt.cpp \
  transforms/wavelet/wt.h \
  wndchrm_src/TrainingSet.cpp \
  wndchrm_src/TrainingSet.h \
  wndchrm_src/WORMfile.cpp \
  wndchrm_src/WORMfile.h \
  wndchrm_src/signatures.cpp \
  wndchrm_src/signatures.h \
  wndchrm_error.cpp \
  wndchrm_error.h

GPU_source = \
  common.cu \
  transforms/chebyshev_kernel.cu \
  textures/gpu_conv2comp.cu \
  textures/gpu_conv2comp.cuh \
  gpu_convolve.cu \
  gpu_convolve.cuh \
  statistics/FeatureStatistics_kernel.cu \
  fft2.cu \
  fft2.cuh \
  textures/haralick/f14_maxcorr_gpu.cu \
  textures/tamura_kernel.cu \
  textures/tamura_cuda.h \
  textures/zernike/zernike_kernel.cu

libchrm_a_CXXFLAGS = -Wall

bin_PROGRAMS=wndchrm

wndchrm_SOURCES = wndchrm_src/wndchrm.cpp

wndchrm_LDADD = libchrm.a -lm -ltiff -L. -lchrm -lfftw3

# Check if CUDA enable
# If enable add GPU sources to libchrm_a_SOURCES
# Add Extra dependancies.
# Ad autotools dont support CUDA have to add LDADD Manually

if CUDA_ENABLE

libchrm_a_SOURCES += $(GPU_source)
CXXFLAGS += -DGPU

wndchrm_LDADD+=$(NVCCFLAGS)
EXTRA_wndchrm_DEPENDENCIES = common.o transforms/chebyshev_kernel.o textures/gpu_conv2comp.o
EXTRA_wndchrm_DEPENDENCIES += gpu_convolve.o statistics/FeatureStatistics_kernel.o fft2.o textures/haralick/f14_maxcorr_gpu.o
EXTRA_wndchrm_DEPENDENCIES += textures/tamura.o textures/zernike/zernike_kernel.o

.cu.o:
	$(NVCC) $(NVCC_INCLUDES) $(DEFAULT_INCLUDES) $(INCLUDES) $(NVCCFLAGS) -c -o $@ $<

wndchrm$(EXEEXT): $(wndchrm_OBJECTS) $(wndchrm_DEPENDENCIES) $(EXTRA_wndchrm_DEPENDENCIES)
	@rm -f wndchrm$(EXEEXT)
	$(NVCC) $(NVCC_FLAGS) $(EXTRA_wndchrm_DEPENDENCIES) $(wndchrm_OBJECTS) $(wndchrm_LDADD) $(LIBS) -o $@
endif

#END GPU

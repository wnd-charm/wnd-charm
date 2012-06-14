#!/usr/bin/env python

"""
setup.py file for SWIG-ified wndchrm
"""

from setuptools import setup, Extension


wndchrm_module = Extension('_pychrm',
	sources=['pychrm/swig/pychrm_wrap.cxx',
		'src/colors/FuzzyCalc.cpp',
		'src/statistics/CombFirst4Moments.cpp',
		'src/statistics/FeatureStatistics.cpp',
		'src/textures/gabor.cpp',
		'src/textures/haarlick/CVIPtexture.cpp',
		'src/textures/haarlick/haarlick.cpp',
		'src/textures/haarlick/mapkit.cpp',
		'src/textures/haarlick/mapkit_generic.cpp',
		'src/textures/tamura.cpp',
		'src/textures/zernike/complex.cpp',
		'src/textures/zernike/zernike.cpp',
		'src/transforms/ChebyshevFourier.cpp',
		'src/transforms/chebyshev.cpp',
		'src/transforms/radon.cpp',
		'src/transforms/wavelet/Common.cpp',
		'src/transforms/wavelet/convolution.cpp',
		'src/transforms/wavelet/DataGrid2D.cpp',
		'src/transforms/wavelet/DataGrid3D.cpp',
		'src/transforms/wavelet/Filter.cpp',
		'src/transforms/wavelet/FilterSet.cpp',
		'src/transforms/wavelet/Symlet5.cpp',
		'src/transforms/wavelet/Wavelet.cpp',
		'src/transforms/wavelet/WaveletHigh.cpp',
		'src/transforms/wavelet/WaveletLow.cpp',
		'src/transforms/wavelet/WaveletMedium.cpp',
		'src/transforms/wavelet/wt.cpp',
		'src/cmatrix.cpp',
		'src/wndchrm_error.cpp',
		'src/FeatureTransforms.cpp',
		'src/FeatureAlgorithms.cpp',
	],
	include_dirs=['src/'],
	libraries=['tiff','fftw3'],
)

setup (name = 'pychrm',
	version = '0.1',
	author      = "Chris Coletta, Ilya Goldberg",
	url = 'http://code.google.com/p/wnd-charm',
	description = """Python bindings for wnd-charm""",
	ext_modules = [wndchrm_module],
	#py_modules = ['pychrm', 'TrainingSet', 'FeatureNameMap', 'FeatureRegistration'],
	packages = ['pychrm'],
	install_requires=[
		'numpy',
		'stringformat >= 0.4',
	],
	dependency_links = [
        "https://github.com/igg/stringformat/tarball/0.4egg=stringformat-0.4"
    ],
)

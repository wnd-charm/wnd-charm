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
		'src/textures/haralick/CVIPtexture.cpp',
		'src/textures/haralick/haralick.cpp',
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
		'src/gsl/specfunc.cpp',
	],
	include_dirs=['src/'],
	libraries=['tiff','fftw3'],
)

import os
pkg_dir = os.path.join (os.path.dirname(os.path.realpath(__file__)),'pychrm')

# this sets the __version__ variable
# the pychrm/_version.py file contains a single line assigning the __version__ variable
# This file gets imported by __init__.py to set __version__ during package import
# Here we do NOT want to import the package we are building to set this variable!
# We DO want to set the version in setup() and have it all synchronized from a single place (pychrm/_version.py)
execfile(os.path.join (pkg_dir,'_version.py'))

# If we're building from svn, get the svn revision using the svnversion command
# The svn revision (if it exists) gets appended to the version string (e.g. 0.3-r123 or 0.2-r141-local)
# This mechanism allows for detection of locally modified (uncommitted) svn versions, and does not
# require any action other than normal svn updates/commits to register new version numbers
# The svn version is written to pychrm/_svn_version.py (which is not under svn control) for inclusion by __init__.py
# If we're not building from svn, then _svn_version.py is not written, and __version__ will be as in pychrm/_version.py
try:
	import subprocess
	import re
	svn_vers = subprocess.Popen(['svnversion', pkg_dir], stdout=subprocess.PIPE).communicate()[0].strip()
	rev_re = re.search(r'^(\d+)(.+)?$', svn_vers)
	if rev_re:
		if rev_re.group(2):
			svn_vers = 'r'+rev_re.group(1)+'-local'
		else:
			svn_vers = 'r'+rev_re.group(1)
		print "svn revision "+svn_vers
		# this construction matches what is done in __init__.py by importing both _version.py and _svn_version.py
		__version__ = __version__+'-'+svn_vers
		f = open(os.path.join(pkg_dir,'_svn_version.py'), 'w+')
		f.write ("__svn_version__ = '"+svn_vers+"'\n")
		f.close()
except:
	pass

setup (name = 'pychrm',
	version = __version__,
	author      = "Chris Coletta, Ilya Goldberg",
	url = 'http://code.google.com/p/wnd-charm',
	description = """Python bindings for wnd-charm""",
	ext_modules = [wndchrm_module],
	#py_modules = ['pychrm', 'FeatureSet', 'FeatureNameMap', 'FeatureRegistration'],
	packages = ['pychrm'],
	install_requires=[
		'numpy',
		'scipy',
		'matplotlib',
		'argparse',
	],
)

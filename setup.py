#!/usr/bin/env python

"""
setup.py file for SWIG-ified wndchrm
"""

from setuptools import setup, Extension

# Need to subclass setuptools/distutils objects to change build order
# such that dynamically-generated swig files will be present when MANIFEST is created.
# See http://stackoverflow.com/a/21236111/1017549

from distutils.command.build import build
from setuptools.command.install import install

class CustomBuild( build):
    def run(self):
        self.run_command('build_ext')
        build.run(self)

class CustomInstall(install):
    def run(self):
        self.run_command('build_ext')
        self.do_egg_install()

import os
pkg_dir = os.path.join( os.path.dirname( os.path.realpath(__file__)), 'wndcharm' )

# this sets the __version__ variable
# the wndcharm/_version.py file contains a single line assigning the __version__ variable
# This file gets imported by __init__.py to set __version__ during package import
# Here we do NOT want to import the package we are building to set this variable!
# We DO want to set the version in setup() and have it all synchronized from a single place (wndcharm/_version.py)
execfile(os.path.join (pkg_dir,'_version.py'))

try:
    from subprocess import check_output

    # Get the git hash for this commit
    # If this returns non-zero, a.k.a. got a git repo, throws CalledProcessError
    git_hash = check_output(['git', 'rev-parse', '--short', 'HEAD']).strip()

    # Check for local modifications
    if check_output(['git', 'diff-index', '--name-only', 'HEAD']).strip():
        git_hash += 'localmod'

    print "git hash " + git_hash
    # this construction matches what is done in __init__.py by importing
    # both _version.py and _git_hash.py use "normalized" semantic version string (a.k.a., dots)
    __version__ = __version__+ '+' + git_hash
    with open( os.path.join( pkg_dir, '_git_hash.py' ), 'w+' ) as f:
	f.write( "__git_hash__ = '{0}'\n".format( git_hash) )
except:
    pass

# Since there's a large C++ underpinning for the wndcharm Python API, run autotools to test build environment
import os
if not os.path.exists( 'config.h' ):
    cmd = os.getcwd() + os.sep + 'configure'
    import subprocess
    p = subprocess.call( [cmd] )

    if p != 0:
	print "Error running configure script"
	import sys
	sys.exit(p)

wndchrm_module = Extension('_wndcharm',
	sources=[
		'wndcharm/swig/wndcharm.i',
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
		'src/ImageTransforms.cpp',
		'src/FeatureAlgorithms.cpp',
		'src/Tasks.cpp',
		'src/FeatureNames.cpp',
		'src/gsl/specfunc.cpp',
	],
	include_dirs=['./','src/', '/usr/local/include'],
	swig_opts=['-threads', '-c++', '-I./', '-I./src', '-outdir', 'wndcharm'],
	libraries=['tiff','fftw3'],
)

setup (
	cmdclass={'build': CustomBuild, 'install': CustomInstall},
	name = 'wndcharm',
	version = __version__,
	author      = "Chris Coletta, Ilya Goldberg",
	url = 'https://github.com/wnd-charm/wnd-charm',
	description = """Python bindings for wnd-charm""",
	license = 'LGPLv2',
	ext_modules = [wndchrm_module],
	packages = ['wndcharm'],
	install_requires=['numpy', 'scipy'],
)

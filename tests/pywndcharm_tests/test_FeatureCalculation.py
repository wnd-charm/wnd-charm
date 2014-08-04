#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Copyright (C) 2014 National Institutes of Health
# All rights reserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#

import sys
if sys.version_info < (2, 7):
    import unittest2 as unittest
else:
    import unittest

from wndcharm.FeatureSet import FeatureSet_Discrete, FisherFeatureWeights,\
        DiscreteImageClassificationResult, Signatures

from os.path import dirname, realpath, join

pychrm_test_dir = dirname( realpath( __file__ ) ) #WNDCHARM_HOME/tests/pywndchrm_tests
wndchrm_test_dir = join( dirname( pychrm_test_dir ), 'wndchrm_tests' )

test_dir = wndchrm_test_dir

class TestFeatureCalculation( unittest.TestCase ):
	"""Feature Calculation"""

	epsilon = 0.00001

	sig_file_path = join( test_dir,'010067_301x300-l_precalculated.sig' )
	test_tif_path = join( test_dir,'010067_301x300.tif' )

	# An alternate target for sig calculation:
	# sig_file = os.path.join (test_dir,'t1_s01_c05_ij-l_precalculated.sig')
	# test_tif = os.path.join (test_dir,'t1_s01_c05_ij.tif')

	# --------------------------------------------------------------------------
	@unittest.skip( "Not doing anything with this right now" )
	def test_ProfileLargeFeatureSet( self ):
		"""Profiling for calculating sigs"""

		import cProfile
		import tempfile
		import pstats
		prof = tempfile.NamedTemporaryFile()
		cmd = 'Signatures.LargeFeatureSet( "{0}" )'.format( self.test_tif_path )
		cProfile.run( cmd, prof.name, 'time')
		p = pstats.Stats(prof.name)
		p.sort_stats('time').print_stats(5)
		prof.close()

		# FIXME: Actually do some checking of the profile results

	# --------------------------------------------------------------------------
	def test_LargeFeatureSetGrayscale( self ):
		"""Large feature set, grayscale image"""
		reference_sample = Signatures.NewFromSigFile( self.sig_file_path,
			image_path=self.test_tif_path )
		target_sample = Signatures.LargeFeatureSet( self.test_tif_path )

		for target_val, res_val in zip( reference_sample.values, target_sample.values ):
			self.assertAlmostEqual( target_val, res_val, delta=self.epsilon )


if __name__ == '__main__':
	unittest.main()

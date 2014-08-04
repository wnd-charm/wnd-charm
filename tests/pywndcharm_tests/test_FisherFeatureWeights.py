#!/usr/bin/env python
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

from os.path import dirname, realpath, join

pychrm_test_dir = dirname( realpath( __file__ ) ) #WNDCHARM_HOME/tests/pywndchrm_tests
wndchrm_test_dir = join( dirname( pychrm_test_dir ), 'wndchrm_tests' )
test_dir = wndchrm_test_dir

from wndcharm.FeatureSet import FeatureSet_Discrete, FisherFeatureWeights

class TestFisherFeatureWeights( unittest.TestCase ):
	"""Fisher score calculation"""

	epsilon = 0.00001

	test_fit_path = join( test_dir,'test-l.fit' )
	test_normalized_fit_path = join( test_dir, 'test_fit-l-normalized.fit' )
	test_feat_weight_path = join( test_dir,'test_fit-l.weights' )
	
	# --------------------------------------------------------------------------
	def test_NewFromFeatureSet( self ):
		"""Fisher score calculation"""

		feature_set = FeatureSet_Discrete.NewFromFitFile( self.test_fit_path )
		feature_set.Normalize()
		result_weights = FisherFeatureWeights.NewFromFeatureSet( feature_set )

		# test weights generated from test-l.fit:
		# wndchrm classify -l -f1.0 -vtest_fit-l.weights test-l.fit test-l.fit 
		target_weights = FisherFeatureWeights.NewFromFile( self.test_feat_weight_path )

	 	for target_val, res_val in zip( target_weights.values, result_weights.values ):
			self.assertAlmostEqual( target_val, res_val, delta=self.epsilon )

	# --------------------------------------------------------------------------
	def test_Normalize( self ):
		"""FIXME: THIS TEST BELONGS IN TEST_FEATURESET.PY"""

		from numpy.testing import assert_allclose
		result_fs = FeatureSet_Discrete.NewFromFitFile( self.test_fit_path )
		result_fs.Normalize()
		target_fs = FeatureSet_Discrete.NewFromFitFile( self.test_normalized_fit_path )

		assert_allclose( result_fs.data_matrix, target_fs.data_matrix, rtol=self.epsilon )

if __name__ == '__main__':
	unittest.main()

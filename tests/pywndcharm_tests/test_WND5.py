#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Copyright (C) 2013 University of Dundee & Open Microscopy Environment.
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
#

import sys
if sys.version_info < (2, 7):
    import unittest2 as unittest
else:
    import unittest

import re
import numpy as np

from os.path import sep, dirname, realpath
from tempfile import mkdtemp
from shutil import rmtree

pychrm_test_dir = dirname( realpath( __file__ ) ) #WNDCHARM_HOME/tests/pychrm_tests
wndchrm_test_dir = dirname( pychrm_test_dir ) + sep + 'wndchrm_tests'

test_dir = wndchrm_test_dir

from wndcharm.FeatureSet import FeatureSet_Discrete, FisherFeatureWeights,\
        DiscreteBatchClassificationResult

epsilon = 0.00001

# Define paths to original files

test_sig_path = os.path.join( test_dir,'t1_s01_c05_ij-l_precalculated.sig' )
test_fit_path = os.path.join( test_dir,'test-l.fit' )
test_feat_wght_path = os.path.join( test_dir,'test_fit-l.weights' )
test_tif_path = os.path.join( test_dir,'t1_s01_c05_ij.tif' )

# Here are the correct values that Python API needs to return

# wndchrm classify -l -f1.0 test-l.fit test/2cell/t45_s06_c06_ij.tif
# test/2cell/t45_s06_c06_ij.tif	6.87e-28	0.614	0.386	*	2cell	2.771
# wndchrm classify -l -f0.14765 test-l.fit test/2cell/t45_s06_c06_ij.tif
# test/2cell/t45_s06_c06_ij.tif	1.34e-27	0.637	0.363	*	2cell	2.727
# wndchrm classify -l -f0.0685 test-l.fit test/2cell/t45_s06_c06_ij.tif
# test/2cell/t45_s06_c06_ij.tif	3.4e-27	0.649	0.351	*	2cell	2.701

correct_marg_probs_2919_feats = [0.047, 0.953]
correct_marg_probs_431_feats = [0.039, 0.961]
correct_marg_probs_200_feats = [0.032, 0.968]


class TestWND5Classification( unittest.TestCase ):
	"""
	WND5 Classification
	"""

	def setUp(self):
		test_sample = Signatures.NewFromSigFile( test_sig_path, test_tif_path )
		feature_set = FeatureSet_Discrete.NewFromFitFile( test_fit_path )
		feature_set.Normalize()

		self.all_weights = FisherFeatureWeights.NewFromFile( test_feat_wght_path )
		self.feature_set = feature_set.FeatureReduce( all_weights.names )
		self.test_sample.Normalize( feature_set )

	def test_WND5_all_features( self ):
		"""
		Checks WND5 Classification with large (2919) feature set
		"""
		result = DiscreteImageClassificationResult.NewWND5(
				self.feature_set, self.all_weights, self.test_sample )

		result_marg_probs = [ round( val, 3 ) \
				for val in result.marginal_probabilities ]

		for target_val, res_val in zip( correct_marg_probs_2919_feats, result_marg_probs ):  
			self.AssertAlmostEqual( target_val, res_val, delta=epsilon )

	def test_WND5_15percent_threshold( self ):
		"""
		Checks WND5 Classification with standard 15% threshold (431) feats
		of large feature set.
		"""

		weights = self.all_weights.Threshold( 431 )
		feat_set = self.feature_set.FeatureReduce( weights.names )
		sample = self.test_sample.FeatureReduce( weights.names )

		result = DiscreteImageClassificationResult.NewWND5( feat_set, weights, sample )

		result_marg_probs = [ round( val, 3 ) \
				for val in result.marginal_probabilities ]

		for target_val, res_val in zip( correct_marg_probs_431_feats, result_marg_probs ):  
			self.AssertAlmostEqual( target_val, res_val, delta=epsilon )

	def test_WND5_200_feats_threshold( self ):
		"""
		Checks WND5 Classification with 200 feature threshold.
		"""
		weights = self.all_weights.Threshold( 200 )
		feat_set = self.feature_set.FeatureReduce( weights.names )
		sample = self.test_sample.FeatureReduce( weights.names )

		result = DiscreteImageClassificationResult.NewWND5( feat_set, weights, sample )

		result_marg_probs = [ round( val, 3 ) \
				for val in result.marginal_probabilities ]

		for target_val, res_val in zip( correct_marg_probs_200_feats, result_marg_probs ):  
			self.AssertAlmostEqual( target_val, res_val, delta=epsilon )


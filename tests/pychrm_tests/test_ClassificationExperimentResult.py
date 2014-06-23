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


# =====================================================================
# Continuous
# =====================================================================
from pychrm.ArtificialFeatureSets import CreateArtificialFeatureSet_Continuous
from pychrm.FeatureSet import ContinuousClassificationExperimentResult,\
				ContinuousFeatureWeights

class TESTINGContinuousClassificationExperimentResult( unittest.TestCase ):
	"""Test various functions from the DiscreteClassificationExperimentResult class."""

	# ---------------------------------------------------------------------
	def test_NewShuffleSplitLeastSquares(self):
		"""CONTINUOUS SHUFFLE SPLIT LEAST SQUARES"""

		kwargs = { 'n_iter' : 5, 'name' : "Continuous Shuffle Split Least Squares POSITIVE CONTROL",
               'quiet' : True }

		fs = CreateArtificialFeatureSet_Continuous()
		exp = ContinuousClassificationExperimentResult.NewShuffleSplit( fs, **kwargs )
		exp.GenerateStats()

		self.assertIs( len(exp), kwargs['n_iter'] )

		# Positive control - Artificial data with defaults should corellate almost perfectly
		self.assertAlmostEqual( exp.pearson_coeff, 1.0, delta=0.01 )

		# Negative control - take the bottom quintile of the artificial features
		# which ARE functions of ground truth but should score low on linear correlation,
		# e.g., sin, x^2, etc.
		max_allowable_pearson_coeff = 0.15
		fs = CreateArtificialFeatureSet_Continuous( num_features_per_signal_type = 5 )
		all_features = ContinuousFeatureWeights.NewFromFeatureSet(fs)
		quintile = int( len(all_features) / 5 )
		crappy_features = all_features.Slice( quintile*4, len( all_features ) )
		crap_featureset = fs.FeatureReduce( crappy_features.names )
		kwargs = { 'n_iter' : 5, 'name' : "Continuous Shuffle Split Least Squares NEGATIVE CONTROL",
						'quiet' : True  }
		exp = ContinuousClassificationExperimentResult.NewShuffleSplit( crap_featureset, **kwargs )
		exp.GenerateStats()
		self.assertAlmostEqual( exp.pearson_coeff, 0.0, delta=max_allowable_pearson_coeff )

	# --------------------------------------------------------------------
	def test_NewShuffleSplitVoting(self):
		"""CONTINUOUS SHUFFLE SPLIT VOTING METHOD"""

		kwargs = { 'n_iter' : 5, 'name' : "Continuous Shuffle Split Voting-Regression POSITIVE CONTROL",
               'classifier' : 'voting', 'quiet' : True }

		fs = CreateArtificialFeatureSet_Continuous()
		exp = ContinuousClassificationExperimentResult.NewShuffleSplit( fs, **kwargs )
		exp.GenerateStats()

		self.assertIs( len(exp), kwargs['n_iter'] )

		# Positive control - Artificial data with defaults should corellate almost perfectly
		self.assertAlmostEqual( exp.pearson_coeff, 1.0, delta=0.01 )

		# Negative control - take the bottom quintile of the artificial features
		# which ARE functions of ground truth but should score low on linear correlation,
		# e.g., sin, x^2, etc.
		max_allowable_pearson_coeff = 0.3
		fs = CreateArtificialFeatureSet_Continuous( num_features_per_signal_type = 5 )
		all_features = ContinuousFeatureWeights.NewFromFeatureSet(fs)
		quintile = int( len(all_features) / 5 )
		crappy_features = all_features.Slice( quintile*4, len( all_features ) )
		crap_featureset = fs.FeatureReduce( crappy_features.names )
		kwargs = { 'n_iter' : 5, 'name' : "Continuous Shuffle Split Voting-Regression NEGATIVE CONTROL",
						'classifier' : 'voting', 'quiet' : True  }
		exp = ContinuousClassificationExperimentResult.NewShuffleSplit( crap_featureset, **kwargs )
		exp.GenerateStats()
		self.assertAlmostEqual( exp.pearson_coeff, 0.0, delta=max_allowable_pearson_coeff )

	# -------------------------------------------------------------------
	def test_PerSampleStatistics(self):
		"""CONTINUOUS PER-SAMPLE STATISTICS"""

		train_size = 16
		test_size = 4
		n_samples = train_size + test_size
		n_iter = 50 # Lots of iterations to make sure every sample is tested

		kwargs = { 'name' : "Continuous PerSample Statistics", 'quiet' : True,
						'n_iter' : n_iter, 'train_size' : train_size, 'test_size' : test_size }

		fs = CreateArtificialFeatureSet_Continuous( n_samples=n_samples,
						num_features_per_signal_type=2)
		exp = ContinuousClassificationExperimentResult.NewShuffleSplit( fs, **kwargs )
		exp.GenerateStats()

		# Capture output from PerSampleStatistics
		from StringIO import StringIO
		out = StringIO()
		try:
				exp.PerSampleStatistics( output_stream=out )
		except Exception as e:
				m = 'Error in experiment.PredictedValueAnalysis: %s' % e
				message += m + '\n'
				self.fail( m )

		# Count the number of lines
		# 3 header lines + 2*num_samples + n_iter*test_size
		per_sample_output = out.getvalue().splitlines()
		self.assertEqual( len(per_sample_output), 3 + 2*n_samples + n_iter*test_size )

# =====================================================================
# Discrete
# =====================================================================
from pychrm.ArtificialFeatureSets import CreateArtificialFeatureSet_Discrete
from pychrm.FeatureSet import FeatureSet_Discrete, DiscreteClassificationExperimentResult

class TESTINGDiscreteClassificationExperimentResult( unittest.TestCase ):
	"""Test various functions from the DiscreteClassificationExperimentResult class."""
	
	# -------------------------------------------------------------------
	def test_PerSampleStatisticsWITHOUTPredictedValue(self):
		"""DISCRETE ShuffleSplit/PerSampleStatistics w/ mini binucleate test set (no predicted value)"""

		fs = FeatureSet_Discrete.NewFromFitFile( '../wndchrm_tests/test-l.fit' )
		exp = DiscreteClassificationExperimentResult.NewShuffleSplit( fs, quiet=True )
		exp.PerSampleStatistics()
		self.assertTrue(True)

	# -------------------------------------------------------------------
	def test_PerSampleStatisticsWITHPredictedValue(self):
		"""DISCRETE PerSampleStatistics with numeric predicted value"""

		fs = CreateArtificialFeatureSet_Discrete()
		exp = DiscreteClassificationExperimentResult.NewShuffleSplit( fs, quiet=True )
		exp.PerSampleStatistics()
		self.assertTrue(True)


if __name__ == '__main__':
		unittest.main()

"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Copyright (C) 2015 National Institutes of Health

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    See the GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 Written by:    Christopher Coletta (github.com/colettace)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"""

import sys
if sys.version_info < (2, 7):
    import unittest2 as unittest
else:
    import unittest


# =====================================================================
# Continuous
# =====================================================================
from wndcharm.ArtificialFeatureSpace import CreateArtificialFeatureSpace_Continuous
from wndcharm.FeatureSpacePredictionExperiment import FeatureSpaceRegressionExperiment
from wndcharm.FeatureWeights import PearsonFeatureWeights

class TESTINGFeatureSpaceRegressionExperiment( unittest.TestCase ):
    """Test various functions from the FeatureSpaceRegressionExperiment class."""

    # ---------------------------------------------------------------------
    #@unittest.skip('')
    def test_NewShuffleSplitLeastSquares(self):
        """CONTINUOUS SHUFFLE SPLIT LEAST SQUARES"""

        fs_kwargs = {}
        fs_kwargs['name'] = "CONTINUOUS PerSampleStatistics_TESTFS"
        fs_kwargs['n_samples'] = 100
        fs_kwargs['num_features_per_signal_type'] = 5
        fs_kwargs['initial_noise_sigma'] = 5
        fs_kwargs['noise_gradient'] = 5
        fs_kwargs['n_samples_per_group'] = 1
        fs_kwargs['random_state'] = 43
        fs_kwargs['singularity'] = True
        fs_kwargs['clip'] = True

        fs = CreateArtificialFeatureSpace_Continuous( **fs_kwargs )

        ss_kwargs = {}
        ss_kwargs['n_iter'] = 5
        ss_kwargs['name'] = "Continuous Shuffle Split Least Squares POSITIVE CONTROL"
        ss_kwargs['quiet'] = True
        ss_kwargs['random_state'] = 43
        exp = FeatureSpaceRegressionExperiment.NewShuffleSplit( fs, **ss_kwargs )

        exp.GenerateStats()
        #exp.Print()

        # len( exp ) is supposed to be the number of batch results (split results)
        self.assertIs( len(exp), ss_kwargs['n_iter'] )

        # Positive control - Artificial data with defaults should corellate almost perfectly
        self.assertAlmostEqual( exp.pearson_coeff, 1.0, delta=0.02 )

        # Negative control - take the bottom quintile of the artificial features
        # which ARE functions of ground truth but should score low on linear correlation,
        # e.g., sin, x^2, etc.

        # With LSTSQ regression of noise features, pearson coeffs tend to be around -0.34 +/- .045
        max_allowable_pearson_coeff = 0.4

        temp_normalized_fs = fs.Normalize( inplace=False )
        ranked_nonzero_features = \
            PearsonFeatureWeights.NewFromFeatureSpace( temp_normalized_fs ).Threshold(_all='nonzero')

        quintile = int( len( ranked_nonzero_features ) / 5 )
        crappy_features = ranked_nonzero_features[ quintile*4 : len( ranked_nonzero_features ) ]
        #crappy_features.Print()
        crap_featureset = fs.FeatureReduce( crappy_features, inplace=False )

        ss_kwargs['name'] = "Continuous Shuffle Split Least Squares NEGATIVE CONTROL"
        exp = FeatureSpaceRegressionExperiment.NewShuffleSplit( crap_featureset, **ss_kwargs )
        exp.GenerateStats()
        exp.PerSampleStatistics()
        #exp.Print()
        self.assertAlmostEqual( exp.pearson_coeff, 0.0, delta=max_allowable_pearson_coeff )

    # --------------------------------------------------------------------
    #@unittest.skip('')
    def test_NewShuffleSplitLinearMultivariateRegression(self):
        """CONTINUOUS SHUFFLE SPLIT LINEAR MULTIVARIATE METHOD"""

        fs_kwargs = {}
        fs_kwargs['name'] = "CONTINUOUS PerSampleStatistics_TESTFS"
        fs_kwargs['n_samples'] = 100
        fs_kwargs['num_features_per_signal_type'] = 5
        fs_kwargs['initial_noise_sigma'] = 5
        fs_kwargs['noise_gradient'] = 5
        fs_kwargs['n_samples_per_group'] = 1
        fs_kwargs['random_state'] = 43
        fs_kwargs['singularity'] = True
        fs_kwargs['clip'] = False

        fs = CreateArtificialFeatureSpace_Continuous( **fs_kwargs )

        ss_kwargs = {}
        ss_kwargs['n_iter'] = 5
        ss_kwargs['name'] = "Continuous Shuffle Split Multivariate-Regression POSITIVE CONTROL"
        ss_kwargs['quiet'] = True
        ss_kwargs['random_state'] = 43
        ss_kwargs['classifier'] = 'linear'
        exp = FeatureSpaceRegressionExperiment.NewShuffleSplit( fs, **ss_kwargs )

        exp.GenerateStats()
        #exp.Print()

        self.assertIs( len(exp), ss_kwargs['n_iter'] )

        # Positive control - Artificial data with defaults should corellate almost perfectly
        self.assertAlmostEqual( exp.pearson_coeff, 1.0, delta=0.03 )

        # Negative control - take the bottom quintile of the artificial features
        # which ARE functions of ground truth but should score low on linear correlation,
        # e.g., sin, x^2, etc.

        # Voting method with crap features tends to be around 0.14 +/- 0.04
        max_allowable_pearson_coeff = 0.2

        temp_normalized_fs = fs.Normalize( inplace=False )
        ranked_nonzero_features = \
            PearsonFeatureWeights.NewFromFeatureSpace( temp_normalized_fs ).Threshold(_all='nonzero')

        quintile = int( len( ranked_nonzero_features ) / 5 )
        crappy_features = ranked_nonzero_features[ quintile*4 : len( ranked_nonzero_features ) ]
        #crappy_features.Print()
        crap_featureset = fs.FeatureReduce( crappy_features )

        ss_kwargs['name'] = "Continuous Shuffle Split Linear Multivariate-Regression NEGATIVE CONTROL",
        exp = FeatureSpaceRegressionExperiment.NewShuffleSplit( crap_featureset, **ss_kwargs )
        exp.GenerateStats()
        #exp.Print()
        self.assertAlmostEqual( exp.pearson_coeff, 0.0, delta=max_allowable_pearson_coeff )

    # -------------------------------------------------------------------
    #@unittest.skip('')
    def test_PerSampleStatistics(self):
        """Testing ContinuousClassificationExperimentResult.PerSampleStatistics()

        Goal is to check the aggregating functionality of PerSampleStatistics"""

        # create a small FeatureSet with not a lot of samples and not a lot of features
        # to enable quick classification

        fs_kwargs = {}
        fs_kwargs['name'] = "CONTINUOUS PerSampleStatistics_TESTFS"
        fs_kwargs['n_samples'] = n_samples = 20
        fs_kwargs['num_features_per_signal_type'] = 2 # small on purpose, to make test fast
        fs_kwargs['initial_noise_sigma'] = 75
        fs_kwargs['noise_gradient'] = 25
        fs_kwargs['n_samples_per_group'] = 1
        fs_kwargs['random_state'] = 42
        fs_kwargs['singularity'] = False
        fs_kwargs['clip'] = False

        fs = CreateArtificialFeatureSpace_Continuous( **fs_kwargs )

        ss_kwargs = {}
        ss_kwargs['name'] = "Continuous PerSampleStatistics ShuffleSplit"
        ss_kwargs['quiet'] = True
        ss_kwargs['n_iter'] = n_iter = 50 # do a lot of iterations so that all samples will be hit
        ss_kwargs['train_size'] = train_size = 16
        ss_kwargs['test_size' ] = test_size = 4
        ss_kwargs['random_state'] = 42
        exp = FeatureSpaceRegressionExperiment.NewShuffleSplit( fs, **ss_kwargs )

        exp.GenerateStats()

        # Capture output from PerSampleStatistics
        from StringIO import StringIO
        out = StringIO()
        try:
            exp.PerSampleStatistics( output_stream=out )
        except Exception as e:
            m = 'Error in experiment.PredictedValueAnalysis: %s' % e
            self.fail( m )

        #print out
        # Count the number of lines
        # 3 header lines + 2*num_samples + n_iter*test_size
        #per_sample_output = out.getvalue().splitlines()
        #self.assertEqual( len(per_sample_output), 3 + 2*n_samples + n_iter*test_size )
        self.assertTrue(True)


# =====================================================================
# Discrete
# =====================================================================
from wndcharm.ArtificialFeatureSpace import CreateArtificialFeatureSpace_Discrete
from wndcharm.FeatureSpacePredictionExperiment import FeatureSpaceClassificationExperiment

class TESTINGFeatureSpaceClassificationExperiment( unittest.TestCase ):
    """Test various functions from the DiscreteClassificationExperimentResult class."""

    # -------------------------------------------------------------------
    #@unittest.skip('')
    def test_PerSampleStatisticsWITHOUTPredictedValue(self):
        """DISCRETE ShuffleSplit/PerSampleStatistics w/ no predicted value"""

        # CAN'T USE THIS, SINCE THE CLASS NAMES ARE INTERPOLATABLE
        # 2-class, 10 samples per class
        #fs = FeatureSet_Discrete.NewFromFitFile( '../wndchrm_tests/test-l.fit' )

        fs_kwargs = {}
        fs_kwargs['name'] = "DISCRETE PerSampleStatistics No Pred Values"
        fs_kwargs['n_samples'] = n_samples = 20
        fs_kwargs['n_classes'] = 2
        fs_kwargs['num_features_per_signal_type'] = 10 # small on purpose, to make test fast
        fs_kwargs['noise_gradient'] = 50
        fs_kwargs['initial_noise_sigma'] = 75
        fs_kwargs['n_samples_per_group'] = 1
        fs_kwargs['random_state'] = 42
        fs_kwargs['interpolatable'] = False
        fs_kwargs['singularity'] = False
        fs_kwargs['clip'] = False
        fs = CreateArtificialFeatureSpace_Discrete( **fs_kwargs )

        ss_kwargs = {}
        ss_kwargs['name'] = "Discrete PerSampleStatistics ShuffleSplit No Pred Values"
        ss_kwargs['quiet'] = True
        ss_kwargs['n_iter'] = n_iter = 1
        ss_kwargs['train_size'] = train_size = 8 # per-class
        ss_kwargs['test_size' ] = test_size = 2 # per-class
        ss_kwargs['random_state'] = 42
        exp = FeatureSpaceClassificationExperiment.NewShuffleSplit( fs, **ss_kwargs )

        ss_kwargs['lda'] = True
        exp = FeatureSpaceClassificationExperiment.NewShuffleSplit( fs, **ss_kwargs )
        #Print calls self.GenereateStats()
        #from os import devnull
        exp.Print( )#output_stream=devnull )
        exp.PerSampleStatistics( )#output_stream=devnull )
        self.assertTrue(True)

    # -------------------------------------------------------------------
    #@unittest.skip('')
    def test_PerSampleStatisticsWITHPredictedValue(self):
        """DISCRETE PerSampleStatistics with numeric predicted value"""

        fs_kwargs = {}
        fs_kwargs['name'] = "DISCRETE PerSampleStatistics WITH Pred Values"
        fs_kwargs['n_samples'] = n_samples = 40
        fs_kwargs['n_classes'] = 2
        fs_kwargs['num_features_per_signal_type'] = 10 # small on purpose, to make test fast
        fs_kwargs['noise_gradient'] = 50
        fs_kwargs['initial_noise_sigma'] = 75
        fs_kwargs['n_samples_per_group'] = 1
        fs_kwargs['random_state'] = 42
        fs_kwargs['interpolatable'] = True
        fs_kwargs['singularity'] = False
        fs_kwargs['clip'] = False
        fs = CreateArtificialFeatureSpace_Discrete( **fs_kwargs )

        ss_kwargs = {}
        ss_kwargs['name'] = "Discrete PerSampleStatistics ShuffleSplit WITH Pred Values"
        ss_kwargs['quiet'] = True
        ss_kwargs['n_iter'] = n_iter = 10
        ss_kwargs['train_size'] = train_size = 8 # per-class
        ss_kwargs['test_size' ] = test_size = 2 # per-class
        ss_kwargs['random_state'] = 42
        exp = FeatureSpaceClassificationExperiment.NewShuffleSplit( fs, **ss_kwargs )

        ss_kwargs['lda'] = True
        exp = FeatureSpaceClassificationExperiment.NewShuffleSplit( fs, **ss_kwargs )
        #Print calls self.GenereateStats()
        #from os import devnull
        exp.Print( )#output_stream=devnull )
        exp.PerSampleStatistics( )#output_stream=devnull )
        self.assertTrue(True)

    # -------------------------------------------------------------------
    def test_GridSearches(self):

        fs_kwargs = {}
        fs_kwargs['name'] = "DISCRETE PerSampleStatistics WITH Pred Values"
        fs_kwargs['n_samples'] = n_samples = 500
        fs_kwargs['n_classes'] = 10
        fs_kwargs['num_features_per_signal_type'] = 10 # small on purpose, to make test fast
        fs_kwargs['noise_gradient'] = 5
        fs_kwargs['initial_noise_sigma'] = 75
        fs_kwargs['n_samples_per_group'] = 1
        fs_kwargs['random_state'] = 42
        fs_kwargs['interpolatable'] = True
        fs_kwargs['singularity'] = False
        fs_kwargs['clip'] = False
        fs = CreateArtificialFeatureSpace_Discrete( **fs_kwargs )

        ss_kwargs = {}
        ss_kwargs['feature_space'] = fs
        ss_kwargs['quiet'] = False
        ss_kwargs['n_iter'] = n_iter = 50
        ss_kwargs['random_state'] = 42
        exp = FeatureSpaceClassificationExperiment.NumSamplesGridSearch( **ss_kwargs )


if __name__ == '__main__':
    unittest.main()

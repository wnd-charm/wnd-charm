"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                               
 Copyright (C) 2015 National Institutes of Health 

    This library is free software; you can redistribute it and/or              
    modify it under the terms of the GNU Lesser General Public                 
    License as published by the Free Software Foundation; either               
    version 2.1 of the License, or (at your option) any later version.         
                                                                               
    This library is distributed in the hope that it will be useful,            
    but WITHOUT ANY WARRANTY; without even the implied warranty of             
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU          
    Lesser General Public License for more details.                            
                                                                               
    You should have received a copy of the GNU Lesser General Public           
    License along with this library; if not, write to the Free Software        
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA  
                                                                               
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                               
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 Written by:  Christopher Coletta (github.com/colettace)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"""

import sys
if sys.version_info < (2, 7):
    import unittest2 as unittest
else:
    import unittest

from os.path import dirname, realpath, join

pychrm_test_dir = dirname( realpath( __file__ ) ) #WNDCHARM_HOME/tests/pywndchrm_tests
wndchrm_test_dir = join( dirname( pychrm_test_dir ), 'wndchrm_tests' )
test_dir = wndchrm_test_dir

from wndcharm.FeatureSpace import FeatureSpace
from wndcharm.FeatureWeights import FisherFeatureWeights
from wndcharm.FeatureVector import FeatureVector
from wndcharm.SingleSamplePrediction import SingleSampleClassification

class TestWND5Classification( unittest.TestCase ):
    """WND5 Classification"""

    # --------------------------------------------------------------------------
    def test_WND5_all_features( self ):
        epsilon = 0.00001

        # Define paths to original files
        test_sig_path = join( test_dir,'t1_s01_c05_ij-l_precalculated.sig' )
        test_fit_path = join( test_dir,'test-l.fit' )
        test_feat_wght_path = join( test_dir,'test_fit-l.weights' )
        test_tif_path = join( test_dir,'t1_s01_c05_ij.tif' )

        # Here are the correct values that Python API needs to return:
        # wndchrm classify -l -f0.75 test-l.fit t1_s01_c05_ij.tif 
        # t1_s01_c05_ij.tif    1.6e-27    0.083    0.917    *    4cell    3.835
        # wndchrm classify -l test-l.fit t1_s01_c05_ij.tif
        # t1_s01_c05_ij.tif    3.19e-27    0.076    0.924    *    4cell    3.848
        # wndchrm classify -l -f0.05 test-l.fit t1_s01_c05_ij.tif
        # t1_s01_c05_ij.tif    1.06e-26    0.066    0.934    *    4cell    3.869

        correct_marg_probs = {}
        correct_marg_probs[2189] = [0.083, 0.917]
        correct_marg_probs[438] = [0.076, 0.924]
        correct_marg_probs[146] = [0.066, 0.934]

        # Load the original files once and only once for all this class's tests
        feature_set = FeatureSpace.NewFromFitFile( test_fit_path )
        fs1 = feature_set.feature_names
        feature_set.Normalize()
        fs2 = feature_set.feature_names
        self.assertSequenceEqual( fs1, fs2 )

        test_sample = FeatureVector( source_filepath=test_tif_path, long=True )
        test_sample.LoadSigFile( test_sig_path )
        self.assertSequenceEqual( feature_set.feature_names, test_sample.feature_names )
        test_sample.Normalize( feature_set )

        all_weights = FisherFeatureWeights.NewFromFile( test_feat_wght_path )

        def Check( num_feats ):
            weights = all_weights.Threshold( num_feats )
            feat_set = feature_set.FeatureReduce( weights )
            sample = test_sample.FeatureReduce( weights )
            result = SingleSampleClassification.NewWND5( feat_set, weights, sample )
            result_marg_probs = [ round( val, 3 ) \
                    for val in result.marginal_probabilities ]
            for target_val, res_val in zip( correct_marg_probs[ num_feats ], result_marg_probs ):
                self.assertAlmostEqual( target_val, res_val, delta=epsilon )

        for num_feats in correct_marg_probs:
            Check( num_feats )

if __name__ == '__main__':
    unittest.main()

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
        # wndchrm classify -l -f1.0 test-l.fit t1_s01_c05_ij.tif 
        # t1_s01_c05_ij.tif    2.09e-27    0.047    0.953    *    4cell    3.906
        # wndchrm classify -l -f0.14765 test-l.fit t1_s01_c05_ij.tif
        # t1_s01_c05_ij.tif    4.29e-27    0.039    0.961    *    4cell    3.923
        # wndchrm classify -l -f0.0685 test-l.fit t1_s01_c05_ij.tif
        # t1_s01_c05_ij.tif    9.05e-27    0.032    0.968    *    4cell    3.936

        correct_marg_probs = {}
        correct_marg_probs[2321] = [0.047, 0.953]
        correct_marg_probs[431] = [0.039, 0.961]
        correct_marg_probs[200] = [0.032, 0.968]

        # Load the original files once and only once for all this class's tests
        feature_set = FeatureSpace.NewFromFitFile( test_fit_path )
        fs1 = feature_set.featurenames_list
        feature_set.Normalize()
        fs2 = feature_set.featurenames_list
        self.assertSequenceEqual( fs1, fs2 )

        test_sample = FeatureVector.NewFromSigFile( test_sig_path, test_tif_path )
        self.assertSequenceEqual( feature_set.featurenames_list, test_sample.featurenames_list )
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

        Check( 2321 )
        Check( 431 )
        Check( 200 )

if __name__ == '__main__':
    unittest.main()

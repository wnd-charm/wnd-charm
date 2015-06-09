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

import numpy as np

from wndcharm.FeatureSpace import FeatureSpace
from wndcharm.FeatureWeights import PearsonFeatureWeights
from wndcharm.FeatureSpacePrediction import FeatureSpaceRegression
from wndcharm.ArtificialFeatureSpace import CreateArtificialFeatureSpace_Continuous

class TestCreateArtificialFeatureSpace( unittest.TestCase ):
    """
    Test CreateArtificialFeatureSpace function
    """

    def test_MultivariateLinearFitOnFitNoTiling( self ):

        fake_continuous = CreateArtificialFeatureSpace_Continuous( n_samples=100,
                num_features_per_signal_type=5, noise_gradient=5, initial_noise_sigma=10,
                n_samples_per_group=1 )
 
        fake_continuous.Normalize( quiet=True )
        reduced_fw = PearsonFeatureWeights.NewFromFeatureSpace( fake_continuous ).Threshold()
        reduced_fs = fake_continuous.FeatureReduce( reduced_fw )
        batch_result = FeatureSpaceRegression.NewMultivariateLinear(
                test_set=reduced_fs, feature_weights=reduced_fw, quiet=True )

    def test_LeastSquaresFitOnFitLeaveOneOutNoTiling( self ):

        fake_continuous = CreateArtificialFeatureSpace_Continuous( n_samples=100,
                num_features_per_signal_type=5, noise_gradient=5, initial_noise_sigma=10,
                n_samples_per_group=1 )
 
        normalized_fs = fake_continuous.Normalize( inplace=False, quiet=True )
        reduced_fw = PearsonFeatureWeights.NewFromFeatureSpace( normalized_fs ).Threshold()
        reduced_fs = fake_continuous.FeatureReduce( reduced_fw )

        batch_result = FeatureSpaceRegression.NewLeastSquares(
            training_set=reduced_fs, test_set=None, feature_weights=reduced_fw, quiet=True )

if __name__ == '__main__':
    unittest.main()

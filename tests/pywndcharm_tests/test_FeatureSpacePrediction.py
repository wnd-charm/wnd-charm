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

import re
import numpy as np

from os.path import sep, dirname, realpath
from tempfile import mkdtemp
from shutil import rmtree

pychrm_test_dir = dirname( realpath( __file__ ) ) #WNDCHARM_HOME/tests/pychrm_tests
wndchrm_test_dir = dirname( pychrm_test_dir ) + sep + 'wndchrm_tests'

from wndcharm.FeatureSpace import FeatureSpace
from wndcharm.FeatureWeights import FisherFeatureWeights, PearsonFeatureWeights
from wndcharm.FeatureSpacePrediction import FeatureSpaceClassification, FeatureSpaceRegression
from wndcharm.FeatureSpacePredictionExperiment import FeatureSpaceClassificationExperiment

class TestFeatureSpaceClassification( unittest.TestCase ):
    """
    Test the classification functionality
    """
    def test_FitOnFit( self ):
        """Uses a curated subset of the IICBU 2008 Lymphoma dataset, preprocessed as follows:
        auto-deconvolved, eosin channel only, tiled 5x6, 3 classes, 10 imgs per class,
        300 samples per class.
        """

        # Inflate the zipped test fit into a temp file
        import zipfile
        zipped_file_path = pychrm_test_dir + sep + 'lymphoma_iicbu2008_subset_EOSIN_ONLY_t5x6_v3.2features.fit.zip'
        zf = zipfile.ZipFile( zipped_file_path, mode='r' )
        tempdir = mkdtemp()
        zf.extractall( tempdir )

        try:
            fitfilepath = tempdir + sep + zf.namelist()[0]

            # Do fit on fit WITHOUT tiling and compare with fit on fit results
            # generated with wndchrm 1.60
            fs = FeatureSpace.NewFromFitFile( fitfilepath ).Normalize( inplace=True, quiet=True )
            #fs = FeatureSpace.NewFromFitFile( wndchrm_test_dir + sep + 'test-l.fit' )
            #fs.ToFitFile( 'temp.fit' )
            fw = FisherFeatureWeights.NewFromFeatureSpace( fs ).Threshold()
            fs.FeatureReduce( fw, inplace=True )
            #fw.Print()
            #fs.Print(verbose=True)
            pychrm_res = FeatureSpaceClassification.NewWND5( fs, fs, fw )
#
#            import cProfile as pr
#            #import profile as pr
#            import tempfile
#            import pstats
#            prof = tempfile.NamedTemporaryFile()
#            cmd = 'no_tile_pychrm_result = DiscreteBatchClassificationResult.New( reduced_fs, reduced_fs, fw )'
#            pr.runctx( cmd, globals(), locals(), prof.name)
#            p = pstats.Stats(prof.name)
#            p.sort_stats('time').print_stats(20)
#            prof.close()

            self.maxDiff = None

            html_path = pychrm_test_dir + sep + 'lymphoma_iicbu2008_subset_eosin_t5x6_v3.2feats_REFERENCE_RESULTS_900_samples_TRAINING_ERROR.html' 
            wres = FeatureSpaceClassificationExperiment.NewFromHTMLReport( html_path )
            wc_batch_result = wres.individual_results[0] # only 1 split in fit-on-fit

            # This takes WAY too long:
            #self.assertSequenceEqual( wc_batch_result.individual_results, pychrm_res.individual_results )
            wc_result = np.empty( (3* len( wc_batch_result.individual_results ) ) )
            for i, single_result in enumerate( wc_batch_result.individual_results ):
                wc_result[ i*3 : (i+1)*3 ] = single_result.marginal_probabilities

            pc_result = np.empty( (3* len( pychrm_res.individual_results ) ) )
            for i, single_result in enumerate( pychrm_res.individual_results ):
                # HTML report only has 3 decimal places
                pc_result[ i*3 : (i+1)*3 ] = \
                    [ float( "{0:0.3f}".format( val ) ) for val in single_result.marginal_probabilities ]

            from numpy.testing import assert_allclose
            assert_allclose( actual=pc_result, desired=wc_result, atol=0.003 )

            #wc_batch_result.Print()
            #pres.Print()

            # ==========================================================
            # Now do the same with tiling, reusing fs from before:
            fs.Update( tile_rows=5, tile_cols=6, num_samples_per_group=30 )
            with_tile_pychrm_result = FeatureSpaceClassification.NewWND5( fs, fs, fw )
            html_path = pychrm_test_dir + sep + 'lymphoma_iicbu2008_subset_eosin_t5x6_v3.2feats_REFERENCE_RESULTS_30_samples_tiled_TRAINING_ERROR.html' 
            with_tile_wndchrm_result = \
              FeatureSpaceClassificationExperiment.NewFromHTMLReport( html_path ).individual_results[0]

            #self.assertSequenceEqual( with_tile_pychrm_result.tiled_results, with_tile_wndchrm_result.individual_results )
            wc_result = np.empty( (3* len( with_tile_wndchrm_result.individual_results ) ) )
            for i, single_result in enumerate( with_tile_wndchrm_result.individual_results ):
                wc_result[ i*3 : (i+1)*3 ] = single_result.marginal_probabilities

            pc_result = np.empty( (3* len( with_tile_pychrm_result.tiled_results ) ) )
            for i, single_result in enumerate( with_tile_pychrm_result.tiled_results ):
                # HTML report only has 3 decimal places
                pc_result[ i*3 : (i+1)*3 ] = \
                    [ float( "{0:0.3f}".format( val ) ) for val in single_result.marginal_probabilities ]

            assert_allclose( actual=pc_result, desired=wc_result, atol=0.003 )

        finally:
            rmtree( tempdir )

    def test_TiledTrainTestSplit( self ):
        """Uses a fake FeatureSpace"""

        from wndcharm.ArtificialFeatureSpace import CreateArtificialFeatureSpace_Discrete
        fs_kwargs = {}
        fs_kwargs['name'] = "DiscreteArtificialFS 10-class"
        fs_kwargs['n_samples'] = 1000
        fs_kwargs['n_classes'] = 10 # 100 samples per class
        fs_kwargs['num_features_per_signal_type'] = 25
        fs_kwargs['initial_noise_sigma'] = 40
        fs_kwargs['noise_gradient'] = 20
        fs_kwargs['n_samples_per_group'] = 4 # 25 images, 2x2 tiling scheme
        fs_kwargs['interpolatable'] = True
        fs_kwargs['random_state'] = 43
        fs_kwargs['singularity'] = False
        fs_kwargs['clip'] = False

        fs = CreateArtificialFeatureSpace_Discrete( **fs_kwargs )

        train, test = fs.Split( random_state=False, quiet=True )
        train.Normalize( inplace=True, quiet=True )
        fw = FisherFeatureWeights.NewFromFeatureSpace( train ).Threshold()

        train.FeatureReduce( fw, inplace=True )
        test.FeatureReduce( fw, inplace=True, quiet=True ).Normalize( train, inplace=True, quiet=True )

        batch_result = FeatureSpaceClassification.NewWND5( train, test, fw  )
        batch_result.GenerateStats()
            
if __name__ == '__main__':
    unittest.main()

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

from wndcharm.FeatureVector import FeatureVector, GenerateFeatureComputationPlan, \
        IncompleteFeatureSetError
from wndcharm.utils import compare

from os.path import dirname, sep, realpath, join, abspath, splitext, basename
from tempfile import mkdtemp
from shutil import rmtree

pychrm_test_dir = dirname( realpath( __file__ ) ) #WNDCHARM_HOME/tests/pywndchrm_tests
wndchrm_test_dir = join( dirname( pychrm_test_dir ), 'wndchrm_tests' )

test_dir = wndchrm_test_dir

class TestFeatureCalculation( unittest.TestCase ):
    """Feature Calculation"""

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
        cmd = 'FeatureVector( source_filepath="{0}", long=True ).GenerateFeatures()'.format( self.test_tif_path )
        cProfile.run( cmd, prof.name, 'time')
        p = pstats.Stats(prof.name)
        p.sort_stats('time').print_stats(5)
        prof.close()

#Loaded features from file /Users/chris/src/wnd-charm/tests/wndchrm_tests/010067_301x300-l_precalculated.sig
#.Fri Jan 23 15:15:35 2015    /var/folders/cr/vsd9_15x6xbc3np6rvx12mqm0000gp/T/tmpf7Oof0
#
#         17595 function calls in 14.956 seconds
#
#   Ordered by: internal time
#   List reduced from 55 to 5 due to restriction <5>
#
#   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
#        1   14.926   14.926   14.926   14.926 {_wndcharm.FeatureComputationPlanExecutor_run}
#        1    0.005    0.005   14.956   14.956 <string>:1(<module>)
#     2920    0.004    0.000    0.004    0.000 {_wndcharm.SwigPyIterator_next}
#     2922    0.004    0.000    0.004    0.000 {method 'format' of 'str' objects}
#        1    0.003    0.003   14.951   14.951 /Users/chris/src/wnd-charm/build/lib.macosx-10.9-x86_64-2.7/wndcharm/FeatureSet.py:1126(GenerateFeatures)
        # FIXME: Actually do some checking of the profile results

    # --------------------------------------------------------------------------
    def test_LargeFeatureSetGrayscale( self ):
        """Large feature set, grayscale image"""
        reference_sample = FeatureVector.NewFromSigFile( self.sig_file_path,
            image_path=self.test_tif_path )

        target_sample = FeatureVector( source_filepath=self.test_tif_path,
            long=True).GenerateFeatures( write_to_disk=False )

#        This doesn't work since the ranges of features are so wide
#        Tried using relative tolerance, but no dice:
#        from numpy.testing import assert_allclose
#        assert_allclose( reference_sample.values, target_sample.values, rtol=1e-3 )

        # Remember we're reading these values in from strings. and the ranges are so wide
        # you only have 6 sig figs. Better apples to apples comparison is to 
        # compare strings.

        from wndcharm.utils import compare
        
        self.assertTrue( compare( target_sample.values, reference_sample.values ) )

    # --------------------------------------------------------------------------
    def test_LoadSubsetFromFile( self ):
        """Calculate one feature family, store to sig, load sig, and use to create larger fs"""

        img_filename = "lymphoma_eosin_channel_MCL_test_img_sj-05-3362-R2_001_E.tif"
        orig_img_filepath = pychrm_test_dir + sep + img_filename

        full_list = list( ('Pixel Intensity Statistics () [3]',
                'Pixel Intensity Statistics (Fourier ()) [3]',) )

        tempdir = mkdtemp()

        from shutil import copy
        
        try:
            # copy the tiff to the tempdir so the .sig files end up there too
            copy( orig_img_filepath, tempdir )
            input_img_path = tempdir + sep + img_filename

            kwargs = {}
            kwargs[ 'source_filepath' ] = input_img_path
            kwargs[ 'tile_num_cols' ] = 6
            kwargs[ 'tile_num_rows' ] = 5
            kwargs[ 'tiling_scheme' ] = '5x6'
            kwargs[ 'tile_col_index' ] = 0
            kwargs[ 'tile_row_index' ] = 0
            kwargs[ 'feature_names' ] = full_list[1:]

            fv1 = FeatureVector( **kwargs ).GenerateFeatures(quiet=False)

            # modify the sig value and write to sig file to make sure subsequent loading
            # used the value from disk and not recalculated it:
            fv1.values[0] = -9999
            fv1.ToSigFile(quiet=False)

            # Now, ask for more features:
            kwargs[ 'feature_names' ] = full_list
            fv2 = FeatureVector( **kwargs )
            with self.assertRaises( IncompleteFeatureSetError ):
                fv2.LoadSigFile()

            #import pdb; pdb.set_trace()
            fv2.GenerateFeatures()
            #self.assertEqual( fv1.values[0], fv2.values[0] )

        finally:
            rmtree( tempdir )



from wndcharm.FeatureSpace import FeatureSpace
from wndcharm.FeatureWeights import FisherFeatureWeights
from wndcharm.SingleSamplePrediction import SingleSampleClassification
from wndcharm.FeatureSpacePrediction import FeatureSpaceClassification
from wndcharm.utils import SampleImageTiles

from sys import exit
#from skimage import io
import numpy as np

class TestSampleImageTiles( unittest.TestCase ):

    def test_HeatMap_w_FeatureComputationPlan( self ):
        """Classification results using SampleImageTiles method and FOF should be the same.
        """

        # chris@NIA-LG-01778617 ~/src/wnd-charm/tests/pywndcharm_tests
        # $ tiffinfo lymphoma_eosin_channel_MCL_test_img_sj-05-3362-R2_001_E.tif
        # TIFF Directory at offset 0x18ea9c (1632924)
        #   Image Width: 1388 Image Length: 1040
        #   Bits/Sample: 8
        #   Compression Scheme: LZW
        #   Photometric Interpretation: min-is-black
        #   Samples/Pixel: 1
        #   Rows/Strip: 5
        #   Planar Configuration: single image plane

        # 5x6 tiling scheme => tile dims 208 x 231.33 each
        scan_x = 231
        scan_y = 208

        #num_features = 200

        # Inflate the zipped test fit into a temp file
        import zipfile
        zipped_file_path = pychrm_test_dir + sep + 'lymphoma_iicbu2008_subset_EOSIN_ONLY_t5x6_v3.2features.fit.zip'
        zf = zipfile.ZipFile( zipped_file_path, mode='r' )
        tempdir = mkdtemp()
        zf.extractall( tempdir )

        fitfilepath = tempdir + sep + zf.namelist()[0]

        img_filename = "lymphoma_eosin_channel_MCL_test_img_sj-05-3362-R2_001_E.tif"
        orig_img_filepath = pychrm_test_dir + sep + img_filename

        from shutil import copy
        
        try:
            # copy the tiff to the tempdir so the .sig files end up there too
            copy( orig_img_filepath, tempdir )
            input_image_path = tempdir + sep + img_filename

            fs = FeatureSpace.NewFromFitFile( fitfilepath, tile_num_rows=5, tile_num_cols=6 )
            fs.Normalize( inplace=True, quiet=True )

            #fs = FeatureSpace.NewFromFitFile( fitfilepath ).Normalize( inplace=True, quiet=True )
            fw = FisherFeatureWeights.NewFromFeatureSpace( fs ).Threshold( 0.05 )
            fs.FeatureReduce( fw, inplace=True )

            # Remember computation plan includes all features from all required
            # feature algorithm families
            comp_plan = GenerateFeatureComputationPlan( fw.feature_names )

            # create the tile image iterator
            image_iter = SampleImageTiles( input_image_path, scan_x, scan_y, True)
            print "Number of samples = " + str( image_iter.samples )

            base, ext = splitext( input_image_path )

            feature_vector_list = []
            # iterate over the image, classifying each tile
            for i, sample in enumerate( image_iter.sample() ):
                kwargs = {}
                kwargs[ 'name' ] = input_image_path
                kwargs[ 'source_filepath' ] = sample
                kwargs[ 'feature_names' ] = fw.feature_names
                kwargs[ 'feature_computation_plan' ] = comp_plan
                kwargs[ 'long' ] = True
                kwargs[ 'tile_num_cols' ] = image_iter.tiles_x
                kwargs[ 'tile_num_rows' ] = image_iter.tiles_y
                kwargs[ 'tiling_scheme' ] = '{0}x{1}'.format( image_iter.tiles_x, image_iter.tiles_y )
                kwargs[ 'tile_col_index' ] = image_iter.current_col
                kwargs[ 'tile_row_index' ] = image_iter.current_row
                kwargs[ 'samplegroupid' ] = 0

                # Setting feature_names initiates the feature reduce from
                # the larger set of features that comes back from computation
                #kwargs[ 'feature_names' ] = fw.feature_names
                # if these are set, then the code will try to take a ROI of a ROI:
                #kwargs[ 'x' ] = image_iter.current_x
                #kwargs[ 'y' ] = image_iter.current_y
                #kwargs[ 'w' ] = image_iter.tile_width
                #kwargs[ 'h' ] = image_iter.tile_height

                samp_feats = FeatureVector( **kwargs ).GenerateFeatures( quiet=False, write_to_disk=False )
                feature_vector_list.append( samp_feats )

            fs_kwargs = {}
            fs_kwargs[ 'feature_vectors_list' ] = feature_vector_list
            fs_kwargs[ 'num_samples' ] = len( feature_vector_list )
            fs_kwargs[ 'num_features' ] = len( fw )
            fs_kwargs[ 'name' ] = input_image_path
            fs_kwargs[ 'num_samples_per_group' ] = fs.num_samples_per_group

            test_set = FeatureSpace.NewFromListOfFeatureVectors( **fs_kwargs )
            test_set.Normalize( reference_features=fs, inplace=True )
            result = FeatureSpaceClassification.NewWND5( fs, test_set, fw )
            result.Print()

        finally:
            rmtree( tempdir )


if __name__ == '__main__':
    unittest.main()

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

from os.path import sep, dirname, realpath
import filecmp
from tempfile import mkdtemp
from shutil import rmtree

import numpy as np

pychrm_test_dir = dirname( realpath( __file__ ) ) #WNDCHARM_HOME/tests/pychrm_tests
wndchrm_test_dir = dirname( pychrm_test_dir ) + sep + 'wndchrm_tests'

from wndcharm.FeatureSpace import FeatureSpace
from wndcharm.FeatureWeights import FisherFeatureWeights
from wndcharm.FeatureSpacePrediction import FeatureSpaceClassification
from wndcharm.ArtificialFeatureSpace import CreateArtificialFeatureSpace_Discrete
from wndcharm.FeatureSpacePredictionExperiment import FeatureSpaceClassificationExperiment
from wndcharm.visualization import PredictedValuesGraph
try:
    import matplotlib
    HasMatplotlib = True
except ImportError:
    HasMatplotlib = False

class TestGraphs( unittest.TestCase ):
    """Test WND-CHARM's graph-making functionality."""
    
    fs_kwargs = {}
    fs_kwargs['name'] = "DiscreteArtificialFS 10-class"
    fs_kwargs['n_samples'] = 1000
    fs_kwargs['n_classes'] = 10
    fs_kwargs['num_features_per_signal_type'] = 25
    fs_kwargs['initial_noise_sigma'] = 40
    fs_kwargs['noise_gradient'] = 20
    fs_kwargs['n_samples_per_group'] = 1
    fs_kwargs['interpolatable'] = True
    fs_kwargs['random_state'] = 43
    fs_kwargs['singularity'] = False
    fs_kwargs['clip'] = False

    fs = CreateArtificialFeatureSpace_Discrete( **fs_kwargs )

    train_set, test_set = fs.Split( random_state=False, quiet=True )
    train_set.Normalize( quiet=True )
    fw = FisherFeatureWeights.NewFromFeatureSpace( train_set ).Threshold()

    reduced_train_set = train_set.FeatureReduce( fw )
    reduced_test_set = test_set.FeatureReduce( fw )
    reduced_test_set.Normalize( reduced_train_set, quiet=True )

    batch_result = FeatureSpaceClassification.NewWND5(
            reduced_train_set, reduced_test_set, fw, quiet=True )

    def setUp( self ):
        self.tempdir = mkdtemp()

    def tearDown( self ):
        rmtree( self.tempdir )

    def CompareGraphs( self, graph, testfilename ):
        """Helper function to check output graphs"""

        # Uncoment to see what graph looks like!
        #graph.SaveToFile( testfilename + 'GRAPH.png' )

        # We used to output the graphs to a png file and do a binary diff on a reference png
        # but there are superficial differences between matplotlib versions that result in
        # the points still being in the right place, but the font is slightly larger,
        # or the text is subtlely offset. So now, we interrogate the matplotlib.figure
        # object and retrieve its coordinates and check them against blessed numpy arrays
        # saved to a npy file.

        axessubplot = graph.figure.gca()

        if len( axessubplot.lines ) > 0:
            # line plot
            try:
                all_coords = np.dstack( tuple( [group._path._vertices for group in axessubplot.lines] ) )
            except AttributeError:
                # older version of matplotlib didn't include leading underscore in attribute
                # "_vertices"
                all_coords = np.dstack( tuple( [group._path.vertices for group in axessubplot.lines] ) )
        elif len( axessubplot.collections ) > 0:
            # scatter plot
            all_coords = np.dstack( tuple( [group._offsets for group in axessubplot.collections] ) )
        else:
            self.fail("Graph doesn't have any lines nor points")

        # uncomment to replace old coords
        #np.save( testfilename, all_coords )
        #from os.path import splitext
        #testfilename_base, ext = splitext( testfilename )
        #np.save( testfilename_base + 'NEW.npy', all_coords )
        reference_array = np.load( testfilename )

        if not np.array_equal( all_coords, reference_array ):
            if not np.allclose( all_coords, reference_array ):
                errmsg = 'Reference graph "{0}" coordinates '.format(testfilename) + \
                    'do not concur with coordinates generated by this test.'
                self.fail( errmsg )

    @unittest.skipIf( HasMatplotlib, "Skipped if matplotlib IS installed" )
    def test_ErrMsgIfMatplotibNotInstalled( self ):
        """Fail gracefully with informative message if matplotlib"""

        graph = PredictedValuesGraph( self.batch_result )
        with self.assertRaises( ImportError ):
            graph.RankOrderedPredictedValuesGraph()
        with self.assertRaises( ImportError ):
            graph.KernelSmoothedDensityGraph()

    @unittest.skipUnless( HasMatplotlib, "Skipped if matplotlib IS NOT installed" )
    def test_RankOrderedFromBatchClassificationResult( self ):
        """Rank Ordered Predicted values graph from a single split"""

        testfilename = 'test_graph_rank_ordered_interpolated_discrete.npy'
        graph = PredictedValuesGraph( self.batch_result )
        graph.RankOrderedPredictedValuesGraph()
        self.CompareGraphs( graph, testfilename )

    @unittest.skipUnless( HasMatplotlib, "Skipped if matplotlib IS NOT installed" )
    def test_KernelSmoothedFromBatchClassificationResult( self ):
        """Kernel Smoothed Probability density graph from a single split"""

        testfilename = 'test_graph_kernel_smoothed.npy'
        graph = PredictedValuesGraph( self.batch_result )
        graph.KernelSmoothedDensityGraph()
        self.CompareGraphs( graph, testfilename )

    @unittest.skipUnless( HasMatplotlib, "Skipped if matplotlib IS NOT installed" )
    def test_FromDiscreteClassificationExperimentResults( self ):
        """Rank Ordered Predicted values graph from an experiment result (multiple splits)"""

        testfilename = 'test_graph_rank_ordered_experiment.npy'

        # Make a smaller featureset to do multiple splits
        fs_kwargs = {}
        fs_kwargs['name'] = "DiscreteArtificialFS RANK ORDERED SHUFFLE SPLIT"
        fs_kwargs['n_samples'] = 100 # smaller
        fs_kwargs['n_classes'] = 5 # smaller, 20 samples per class
        fs_kwargs['num_features_per_signal_type'] = 10 # smaller
        fs_kwargs['initial_noise_sigma'] = 50
        fs_kwargs['noise_gradient'] = 20
        fs_kwargs['n_samples_per_group'] = 1
        fs_kwargs['interpolatable'] = True
        fs_kwargs['random_state'] = 42
        fs_kwargs['singularity'] = False
        fs_kwargs['clip'] = False

        small_fs = CreateArtificialFeatureSpace_Discrete( **fs_kwargs )

        ss_kwargs = {}
        ss_kwargs['quiet'] = True
        ss_kwargs['n_iter'] = n_iter = 10
        ss_kwargs['train_size'] = train_size = 18 # per-class
        ss_kwargs['test_size' ] = test_size = 2 # per-class
        ss_kwargs['random_state'] = 42
        exp = FeatureSpaceClassificationExperiment.NewShuffleSplit( small_fs, **ss_kwargs )
        graph = PredictedValuesGraph( exp, use_averaged_results=False )
        graph.RankOrderedPredictedValuesGraph()
        self.CompareGraphs( graph, testfilename )

    @unittest.skipUnless( HasMatplotlib, "Skipped if matplotlib IS NOT installed" )
    def test_FromHTML( self ):
        """Rank Ordered Predicted values graph from an experiment result (multiple splits)"""

        testfilename = 'test_graph_fromHTML.npy'
         # Inflate the zipped html file into a temp file
        import zipfile

        #zipped_file_path = pychrm_test_dir + sep + 'c_elegans_terminal_bulb.html'
        #import zlib
        #zf = zipfile.ZipFile( zipped_file_path + '.zip', mode='w' )
        #zf.write( zipped_file_path, compress_type=zipfile.ZIP_DEFLATED )
        #zf.close()
 
        zipped_file_path = pychrm_test_dir + sep + 'c_elegans_terminal_bulb.html.zip'
        zf = zipfile.ZipFile( zipped_file_path, mode='r' )
        zf.extractall( self.tempdir )
        htmlfilepath = self.tempdir + sep + zf.namelist()[0]
        graph = PredictedValuesGraph.NewFromHTMLReport( htmlfilepath, use_averaged_results=False ) 
        graph.RankOrderedPredictedValuesGraph()

        self.CompareGraphs( graph, testfilename )

    @unittest.skipUnless( HasMatplotlib, "Skipped if matplotlib IS NOTinstalled" )
    def test_IfNotInterpolatable( self ):
        """You can't graph predicted values if the classes aren't interpolatable."""

        testfilename = 'ShouldntBeGraphable.png'
        small_fs = CreateArtificialFeatureSpace_Discrete( 
                        n_samples=20, n_classes=2, random_state=42, interpolatable=False )
        train_set, test_set = small_fs.Split( random_state=False, quiet=True )
        train_set.Normalize()

        fw = FisherFeatureWeights.NewFromFeatureSpace( train_set ).Threshold()
        reduced_train_set = train_set.FeatureReduce( fw )
        reduced_test_set = test_set.FeatureReduce( fw )
        test_set.Normalize( train_set, quiet=True )

        batch_result = FeatureSpaceClassification.NewWND5(
                                    reduced_train_set, reduced_test_set, fw, quiet=True )
        with self.assertRaises( ValueError ):
            graph = PredictedValuesGraph( batch_result )

if __name__ == '__main__':
    unittest.main()
